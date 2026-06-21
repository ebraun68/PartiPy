#!/usr/bin/env python3
"""
label_tree.py - Root a Newick tree on the first taxon in --tax-order
                 and label all non-trivial internal nodes with their
                 Jakobsen hex bipartition codes.

The tree is re-rooted so that the first taxon in --tax-order is the sole
member of one side of the root split (i.e. it becomes the outgroup).  All
non-trivial internal nodes (both sides of the split contain >= 2 taxa) are
then labeled with the Jakobsen hex code of their clade bipartition, computed
using the bit ordering in --tax-order.

The labeled, rooted Newick is written to stdout.  If -o/--outfile is given,
it is also written to PREFIX.nwk.

Usage:
  python3 label_tree.py TREE.nwk --tax-order FILE [-o PREFIX]
"""

import argparse
import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

try:
    from partcompat import (
        read_tax_order,
        _canonicalize,
        _coded_to_bitmask,
        bitmask_to_jakobsen_hex,
    )
except ImportError as e:
    sys.exit(f"ERROR: Cannot import partcompat.py from {_HERE}: {e}")

# ---------------------------------------------------------------------------
# Newick parsing
# ---------------------------------------------------------------------------

def strip_comments(nwk):
    return re.sub(r'\[.*?\]', '', nwk, flags=re.DOTALL)


def tokenize(nwk):
    buf = []
    for ch in nwk:
        if ch in '(),;':
            if buf:
                yield ''.join(buf).strip()
                buf = []
            if ch != ';':
                yield ch
        else:
            buf.append(ch)
    if buf:
        s = ''.join(buf).strip()
        if s:
            yield s


def parse_label_token(token):
    """Split 'name:length' token into (name, length_or_None)."""
    token = token.strip()
    if ':' in token:
        left, right = token.rsplit(':', 1)
        try:
            length = float(right)
        except ValueError:
            length = None
        try:
            float(left)
            name = ''
        except ValueError:
            name = left.strip()
        return name, length
    else:
        try:
            float(token)
            return '', None
        except ValueError:
            return token, None


def parse_newick(text):
    """
    Parse a Newick string into a tree of dicts.

    Each node dict has:
      'name'          : str (leaf name or internal label, may be empty)
      'branch_length' : float or None
      'children'      : list of child node dicts
    """
    tokens = list(tokenize(strip_comments(text)))
    pos = [0]

    def peek():
        return tokens[pos[0]] if pos[0] < len(tokens) else None

    def consume():
        t = tokens[pos[0]]
        pos[0] += 1
        return t

    def parse_node():
        if peek() == '(':
            consume()   # '('
            children = [parse_node()]
            while peek() == ',':
                consume()
                children.append(parse_node())
            if peek() != ')':
                raise ValueError(
                    f"Expected ')' at position {pos[0]}, got {peek()!r}"
                )
            consume()   # ')'
            name, bl = '', None
            if peek() not in (')', ',', None):
                name, bl = parse_label_token(consume())
            return {'name': name, 'branch_length': bl, 'children': children}
        else:
            tok = consume() if peek() not in (')', ',', None) else ''
            name, bl = parse_label_token(tok)
            return {'name': name, 'branch_length': bl, 'children': []}

    return parse_node()


def get_leaves(node):
    if not node['children']:
        return [node['name']]
    leaves = []
    for child in node['children']:
        leaves.extend(get_leaves(child))
    return leaves

# ---------------------------------------------------------------------------
# Re-rooting
# ---------------------------------------------------------------------------

def find_leaf_path(node, target, path):
    """DFS for leaf `target`; builds path from root to target (inclusive)."""
    path.append(node)
    if not node['children'] and node['name'] == target:
        return True
    for child in node['children']:
        if find_leaf_path(child, target, path):
            return True
    path.pop()
    return False


def reroot(root, outgroup_name):
    """
    Re-root the tree so that `outgroup_name` is the sole child on one side
    of the new root.

    Algorithm:
      1. Find the path from the current root to the outgroup leaf.
      2. Detach the outgroup leaf from its parent.
      3. Reverse all edges from the outgroup parent back up to the old root,
         preserving branch lengths on each reversed edge.
      4. If the old root is now degree-1 (was degree-2 originally), absorb
         it by summing its branch length with its sole remaining child's.
      5. Build a new root with two children: the outgroup leaf and the
         modified subtree.
    """
    path = []
    if not find_leaf_path(root, outgroup_name, path):
        sys.exit(f"ERROR: outgroup taxon '{outgroup_name}' not found in tree")
    # path = [old_root, ..., outgroup_parent, outgroup_leaf]
    n = len(path)

    outgroup_leaf = path[-1]
    og_bl = outgroup_leaf['branch_length']

    # Step 1: detach outgroup leaf from its parent
    outgroup_parent = path[-2]
    outgroup_parent['children'] = [
        c for c in outgroup_parent['children'] if c is not outgroup_leaf
    ]

    # Step 2: reverse edges from outgroup_parent up toward old root
    for i in range(n - 2, 0, -1):
        node   = path[i]       # was a child; becomes a parent
        parent = path[i - 1]   # was a parent; becomes a child
        edge_len = node['branch_length']   # length of edge node -> parent
        parent['children'] = [c for c in parent['children'] if c is not node]
        parent['branch_length'] = edge_len
        node['branch_length'] = None
        node['children'].append(parent)

    # Step 3: absorb old root if it is now degree-1
    old_root = path[0]
    if len(old_root['children']) == 1:
        sole = old_root['children'][0]
        if old_root['branch_length'] is not None:
            sole['branch_length'] = (
                (sole['branch_length'] or 0.0) + old_root['branch_length']
            )
        if n >= 3:
            path[1]['children'] = [
                sole if c is old_root else c
                for c in path[1]['children']
            ]

    # Step 4: determine the subtree that becomes the second child of new root
    if n == 2:
        # Outgroup was a direct child of old root
        remaining = [c for c in old_root['children'] if c is not outgroup_leaf]
        if len(remaining) == 1:
            subtree = remaining[0]   # absorbed
        else:
            subtree = old_root       # still has multiple children
        subtree['branch_length'] = None
    else:
        subtree = outgroup_parent
        subtree['branch_length'] = None

    outgroup_leaf['branch_length'] = og_bl
    return {
        'name': '', 'branch_length': None,
        'children': [outgroup_leaf, subtree],
    }

# ---------------------------------------------------------------------------
# Bipartition labeling
# ---------------------------------------------------------------------------

def label_nodes(node, order, n_tax, prefix=''):
    """
    Post-order traversal: label each non-trivial internal node with its
    Jakobsen hex bipartition code (optionally with a prefix string).
    Returns the frozenset of leaf names in this node's clade.
    """
    if not node['children']:
        return frozenset([node['name']])

    clade = frozenset()
    for child in node['children']:
        clade |= label_nodes(child, order, n_tax, prefix=prefix)

    n_clade = len(clade)
    if 2 <= n_clade <= n_tax - 2:
        raw   = [1 if t in clade else 0 for t in order]
        canon = _canonicalize(raw)
        bm    = _coded_to_bitmask(canon)
        node['name'] = prefix + bitmask_to_jakobsen_hex(bm, n_tax)
    else:
        node['name'] = ''

    return clade


def count_labels(node):
    """Count non-empty labels on internal nodes."""
    if not node['children']:
        return 0
    return (1 if node.get('name') else 0) + sum(
        count_labels(c) for c in node['children']
    )

# ---------------------------------------------------------------------------
# Newick serialization
# ---------------------------------------------------------------------------

def newick_str(node, is_root=False):
    """Serialize a node dict to a Newick string fragment."""
    if not node['children']:
        s = node['name'] or ''
        if node['branch_length'] is not None:
            s += f":{node['branch_length']}"
        return s

    children_str = ','.join(newick_str(c) for c in node['children'])
    label = node.get('name', '')
    s = f"({children_str}){label}"
    if not is_root and node['branch_length'] is not None:
        s += f":{node['branch_length']}"
    if is_root:
        s += ';'
    return s

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            'Root a Newick tree on the first taxon in --tax-order and '
            'label all non-trivial internal nodes with their Jakobsen hex '
            'bipartition codes.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
The first taxon in --tax-order is used as the outgroup; the tree is re-rooted
so that this taxon is the sole member of one side of the root split.  All
non-trivial internal nodes (both sides of the split contain >=2 taxa) are then
labeled with the Jakobsen hex bipartition code of their clade.

Hex codes are computed using the bit ordering in --tax-order, so they match
the codes produced by partcompat.py, newick_to_biparts.py, and all other
PartiPy programs given the same --tax-order file.

The labeled Newick is written to stdout and, if -o/--outfile is given,
also to PREFIX.nwk.

Example:
  python3 label_tree.py tree.nwk --tax-order taxa.txt -o labeled
        ''',
    )
    p.add_argument('tree',
                   help='Input Newick tree file')
    p.add_argument('--tax-order', required=True,
                   help='Taxon order file (one name per line).  The first '
                        'taxon is used as the outgroup for rooting.  Bit '
                        'positions in the hex codes follow this order.')
    p.add_argument('-o', '--outfile', default=None,
                   help='Optional output prefix; writes PREFIX.nwk')
    p.add_argument('--label-prefix', default='',
                   metavar='STR',
                   help='String prepended to every hex label (default: none). '
                        'Use --label-prefix b_ or --label-prefix _ to work '
                        'around FigTree 1.4.4 parsing hex labels as numbers.')
    return p.parse_args()

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if not os.path.isfile(args.tree):
        sys.exit(f"ERROR: Tree file not found: {args.tree}")

    order = read_tax_order(args.tax_order)
    n_tax = len(order)
    if n_tax == 0:
        sys.exit(f"ERROR: No taxa found in {args.tax_order}")
    print(f"Taxon order: {n_tax} taxa", file=sys.stderr)

    with open(args.tree) as f:
        nwk_text = f.read().strip()

    try:
        root = parse_newick(nwk_text)
    except Exception as e:
        sys.exit(f"ERROR: Could not parse Newick tree: {e}")

    # Validate leaf set vs --tax-order
    tree_taxa  = set(get_leaves(root))
    order_set  = set(order)

    in_tree_not_order = tree_taxa  - order_set
    in_order_not_tree = order_set  - tree_taxa

    if in_tree_not_order or in_order_not_tree:
        lines = ["ERROR: taxon sets in tree and --tax-order do not match."]
        if in_tree_not_order:
            lines.append(
                f"  In tree but not in --tax-order ({len(in_tree_not_order)}):"
            )
            for t in sorted(in_tree_not_order):
                lines.append(f"    {t}")
        if in_order_not_tree:
            lines.append(
                f"  In --tax-order but not in tree ({len(in_order_not_tree)}):"
            )
            for t in sorted(in_order_not_tree):
                lines.append(f"    {t}")
        sys.exit('\n'.join(lines))

    outgroup = order[0]
    print(f"Outgroup (first in --tax-order): {outgroup}", file=sys.stderr)

    new_root = reroot(root, outgroup)

    label_nodes(new_root, order, n_tax, prefix=args.label_prefix)
    n_labeled = count_labels(new_root)
    print(f"Labeled {n_labeled} non-trivial internal node(s)",
          file=sys.stderr)

    result = newick_str(new_root, is_root=True)
    print(result)

    if args.outfile:
        out_path = args.outfile + '.nwk'
        with open(out_path, 'w') as f:
            f.write(result + '\n')
        print(f"Wrote {out_path}", file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
