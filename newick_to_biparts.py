#!/usr/bin/env python3

import os
import sys
import argparse
import re
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
try:
    from partcompat import bitmask_to_jakobsen_hex
    _HAVE_PARTCOMPAT = True
except ImportError:
    _HAVE_PARTCOMPAT = False


def strip_comments(nwk):
    """Remove bracketed comments e.g. [&comment]."""
    return re.sub(r'\[.*?\]', '', nwk)


def tokenize(nwk):
    """Yield tokens: '(', ')', ',', or a label string."""
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


def parse_newick(nwk):
    """
    Parse Newick into a nested list structure.
    Leaves are strings (taxon names); internal nodes are lists of children.
    Branch lengths and support values are stripped.
    """
    tokens = list(tokenize(strip_comments(nwk)))
    pos = [0]

    def strip_label(s):
        s = s.strip()
        if ':' in s:
            s = s[:s.index(':')]
        return s.strip()

    def parse_node():
        if tokens[pos[0]] == '(':
            pos[0] += 1
            children = [parse_node()]
            while tokens[pos[0]] == ',':
                pos[0] += 1
                children.append(parse_node())
            assert tokens[pos[0]] == ')', \
                f"Expected ')' at token {pos[0]}, got '{tokens[pos[0]]}'"
            pos[0] += 1
            # consume optional internal label / support value
            if pos[0] < len(tokens) and tokens[pos[0]] not in '(),':
                pos[0] += 1
            return children
        else:
            label = strip_label(tokens[pos[0]])
            pos[0] += 1
            return label

    return parse_node()


def get_leaves(node):
    if isinstance(node, str):
        return [node]
    leaves = []
    for child in node:
        leaves.extend(get_leaves(child))
    return leaves


def get_bipartitions(node, n_total):
    """
    Recursively collect non-trivial bipartitions (both sides >= 2 taxa).
    Returns (frozenset_of_clade_leaves, list_of_bipartition_frozensets).
    """
    if isinstance(node, str):
        return frozenset([node]), []

    child_clades = []
    biparts = []
    for child in node:
        clade, child_biparts = get_bipartitions(child, n_total)
        child_clades.append(clade)
        biparts.extend(child_biparts)

    this_clade = frozenset().union(*child_clades)
    clade_size = len(this_clade)

    # Non-trivial: both sides of the split have >= 2 taxa
    if 2 <= clade_size <= n_total - 2:
        biparts.append(this_clade)

    return this_clade, biparts


def clade_to_bitstring(clade, order):
    """Convert a frozenset of taxon names to a canonical binary string."""
    bits = [1 if t in clade else 0 for t in order]
    if bits[0] == 0:
        bits = [1 - b for b in bits]
    return ''.join(str(b) for b in bits)


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            'Extract non-trivial bipartitions from a Newick tree. '
            'Output is written to stdout and is suitable for direct use with '
            'partcompat.py -t/--tree-biparts or -p/--plot-biparts. '
            'Multiple outputs can be concatenated; duplicates and # comment '
            'lines are handled automatically by partcompat.py.'
        ),
    )
    p.add_argument('tree',
                   help='Newick tree file')
    p.add_argument('tax_order', nargs='?', default=None,
                   help='Optional taxon order file (one name per line). '
                        'Bit positions in the output follow this order. '
                        'If omitted, taxa are ordered by first appearance '
                        'in the Newick. Should match the --tax-order file '
                        'used with partcompat.py.')
    p.add_argument('--hex', action='store_true',
                   help='Include a commented section listing each bipartition '
                        'binary string alongside its Jakobsen hex index. '
                        'Useful for cross-referencing with partcompat.py TSV '
                        'output. Can be placed anywhere on the command line.')
    return p.parse_args()


def main():
    args = parse_args()

    if not os.path.isfile(args.tree):
        sys.exit(f"ERROR: Tree file not found: {args.tree}")
    with open(args.tree) as f:
        nwk = f.read().strip()

    tree     = parse_newick(nwk)
    all_taxa = get_leaves(tree)
    n_total  = len(set(all_taxa))

    # Taxon order: from file if supplied, else order of first appearance
    if args.tax_order:
        with open(args.tax_order) as f:
            order = [l.strip() for l in f if l.strip()]
        missing = set(all_taxa) - set(order)
        if missing:
            print(f"WARNING: taxa in tree but not in order file: {missing}",
                  file=sys.stderr)
        order = [t for t in order if t in set(all_taxa)]
    else:
        seen  = set()
        order = []
        for t in all_taxa:
            if t not in seen:
                order.append(t)
                seen.add(t)

    _, biparts = get_bipartitions(tree, n_total)

    # Deduplicate: a split and its complement are the same bipartition
    seen_bits      = set()
    unique_biparts = []
    for clade in biparts:
        bs = clade_to_bitstring(clade, order)
        if bs not in seen_bits:
            seen_bits.add(bs)
            unique_biparts.append((bs, clade))

    # Output
    print("# Bipartitions extracted from Newick tree")
    print(f"# {len(order)} taxa in bit-string order:")
    for i, t in enumerate(order):
        print(f"#   bit {i+1:2d}: {t}")
    print(f"# {len(unique_biparts)} non-trivial internal bipartitions")

    if args.hex:
        if not _HAVE_PARTCOMPAT:
            print("# WARNING: partcompat.py not found; hex indices unavailable",
                  file=sys.stderr)
        else:
            n_tax = len(order)
            print("#")
            print("# Hexadecimal bipartition indices:")
            for bs, _ in unique_biparts:
                bm      = int(bs, 2)
                hex_str = bitmask_to_jakobsen_hex(bm, n_tax)
                print(f"#   {bs}	{hex_str}")
            print("#")

    print()
    for bs, _ in unique_biparts:
        print(bs)


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
