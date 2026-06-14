#!/usr/bin/env python3
"""
compare_trees.py - Compare internal branch lengths between two Newick trees,
                   or compare a tree against Lento-plot support values.

Two modes (mutually exclusive):
  Tree-vs-tree:  -1 TREE1.nwk  -2 TREE2.nwk  --tax-order FILE
  Tree-vs-Lento: -1 TREE1.nwk  -l LENTO.tsv  --tax-order FILE

Produces:
  PREFIX.tsv        branch-length comparison table
  PREFIX.svg        paired bar plot (both bars upward, normalized by default)
  PREFIX.trunc.svg  (with --truncate N) first N pairs only

Usage:
  python3 compare_trees.py -1 TREE1.nwk -2 TREE2.nwk --tax-order FILE -o PREFIX
  python3 compare_trees.py -1 TREE1.nwk -l LENTO.tsv  --tax-order FILE -o PREFIX
"""

import argparse
import csv
import re
import sys
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

try:
    from partcompat import (
        read_tax_order,
        read_tree_biparts,
        _canonicalize,
        _coded_to_bitmask,
        _bitmask_to_coded,
        bitmask_to_binary_str,
        bitmask_to_jakobsen_hex,
    )
except ImportError as e:
    sys.exit(f"ERROR: Cannot import partcompat.py from {_HERE}: {e}")

# ---------------------------------------------------------------------------
# Newick parsing
# ---------------------------------------------------------------------------

def strip_comments(nwk):
    return re.sub(r'\[.*?\]', '', nwk)


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


def parse_label(token):
    token = token.strip()
    if ':' in token:
        left, right = token.rsplit(':', 1)
        try:
            length = float(right)
        except ValueError:
            length = None
        left = left.strip()
        try:
            float(left)
            name = ''
        except ValueError:
            name = left
        return name, length
    else:
        try:
            float(token)
            return '', None
        except ValueError:
            return token, None


def parse_newick(nwk):
    tokens = list(tokenize(strip_comments(nwk)))
    pos = [0]

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
            branch_length = None
            if pos[0] < len(tokens) and tokens[pos[0]] not in '(),':
                _, branch_length = parse_label(tokens[pos[0]])
                pos[0] += 1
            return {'children': children, 'name': '', 'branch_length': branch_length}
        else:
            name, branch_length = parse_label(tokens[pos[0]])
            pos[0] += 1
            return {'children': [], 'name': name, 'branch_length': branch_length}

    return parse_node()


def get_leaves(node):
    if not node['children']:
        return [node['name']]
    leaves = []
    for child in node['children']:
        leaves.extend(get_leaves(child))
    return leaves


def extract_internal_branches(node, n_total, order):
    if not node['children']:
        return frozenset([node['name']]), []

    child_clades, branches = [], []
    for child in node['children']:
        clade, child_branches = extract_internal_branches(child, n_total, order)
        child_clades.append(clade)
        branches.extend(child_branches)

    this_clade = frozenset().union(*child_clades)
    clade_size = len(this_clade)

    if 2 <= clade_size <= n_total - 2:
        bits  = [1 if t in this_clade else 0 for t in order]
        canon = _canonicalize(bits)
        bm    = _coded_to_bitmask(canon)
        branches.append((bm, node['branch_length']))

    return this_clade, branches


def load_tree(path, order, n_tax, label):
    """Parse a Newick file and return {bitmask: branch_length}."""
    if not os.path.isfile(path):
        sys.exit(f"ERROR: Tree file not found: {path}")
    with open(path) as f:
        nwk = f.read().strip()
    tree     = parse_newick(nwk)
    all_taxa = get_leaves(tree)
    order_set = set(order)

    missing = order_set - set(all_taxa)
    extra   = set(all_taxa) - order_set
    if missing:
        print(f"WARNING: {label}: {len(missing)} taxa in --tax-order absent "
              f"from tree: {sorted(missing)}", file=sys.stderr)
    if extra:
        print(f"WARNING: {label}: {len(extra)} taxa in tree absent from "
              f"--tax-order (ignored): {sorted(extra)}", file=sys.stderr)

    _, branches = extract_internal_branches(tree, n_tax, order)
    seen = {}
    for bm, length in branches:
        if bm not in seen:
            seen[bm] = length
    total = sum(v for v in seen.values() if v is not None)
    print(f"{label} ({path}): {len(seen)} non-trivial internal branches, "
          f"total internal length = {total:.6f}", file=sys.stderr)
    return seen

# ---------------------------------------------------------------------------
# Lento TSV reader
# ---------------------------------------------------------------------------

def load_lento(path, n_tax):
    """
    Read support_norm values from a partcompat.py / sample_loci.py TSV.
    Returns {bitmask: support_norm}.
    """
    result = {}
    try:
        with open(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                bp = row.get('bipartition_binary', '').strip().lstrip("'")
                if not bp or not all(c in '01' for c in bp):
                    continue
                try:
                    sup = float(row.get('support_norm', 0))
                except (ValueError, KeyError):
                    continue
                if len(bp) != n_tax:
                    continue
                bm = int(bp, 2)
                result[bm] = sup
    except Exception as e:
        sys.exit(f"ERROR: Cannot read Lento TSV '{path}': {e}")

    total = sum(result.values())
    print(f"Lento TSV ({path}): {len(result)} bipartitions, "
          f"total support_norm = {total:.4f}", file=sys.stderr)
    return result

# ---------------------------------------------------------------------------
# SVG paired bar plot
# ---------------------------------------------------------------------------

def _xml_escape(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def write_comparison_svg(path, sorted_bm, vals1, vals2, n_tax, order,
                         label1, label2, normalize, show_hex,
                         max_cols=None):
    """
    Paired upward bar plot.

    sorted_bm : ordered list of bitmask ints
    vals1/2   : dicts {bitmask -> float}  (raw values, before normalization)
    label1/2  : legend labels
    normalize : if True, scale vals2 so sum(vals2) == sum(vals1)
    show_hex  : draw Jakobsen hex labels below dot schematic
    max_cols  : if not None, truncate to first max_cols pairs
    """
    pairs = sorted_bm[:max_cols] if max_cols is not None else sorted_bm
    n_pairs = len(pairs)
    if n_pairs == 0:
        print("WARNING: no bipartitions to plot.", file=sys.stderr)
        return

    # Normalization
    total1 = sum((vals1.get(bm) or 0.0) for bm in sorted_bm)
    total2 = sum((vals2.get(bm) or 0.0) for bm in sorted_bm)
    if normalize and total1 > 0 and total2 > 0:
        scale2 = total1 / total2
    else:
        scale2 = 1.0

    # Layout constants
    BAR_W           = 10    # width of each bar
    BAR_GAP         = 2     # gap between the two bars in a pair
    PAIR_GAP        = 6     # gap between consecutive pairs
    PAIR_W          = 2 * BAR_W + BAR_GAP + PAIR_GAP
    MAX_BAR_H       = 120
    SCHEMATIC_ROW_H = 10
    SCHEMATIC_GAP   = 8
    HEX_FONT        = 7
    HEX_GAP         = 6
    MARGIN_BASE     = 10
    MARGIN_LEFT     = 120
    MARGIN_RIGHT    = 20
    AXIS_H          = 2
    LABEL_FONT      = 10
    TITLE_FONT      = 11

    # Colors
    COL1 = '#1565c0'   # blue for tree1
    COL2 = '#e65100'   # orange for tree2/Lento

    if show_hex:
        max_hex_len = max((len(bitmask_to_jakobsen_hex(bm, n_tax))
                           for bm in pairs), default=1)
        hex_label_h = max_hex_len * HEX_FONT
    else:
        hex_label_h = 0

    MARGIN_TOP  = MARGIN_BASE + hex_label_h + (HEX_GAP if show_hex else 0)
    schematic_h = n_tax * SCHEMATIC_ROW_H
    total_w     = MARGIN_LEFT + n_pairs * PAIR_W + MARGIN_RIGHT
    total_h     = (MARGIN_TOP
                   + MAX_BAR_H + AXIS_H
                   + SCHEMATIC_GAP + schematic_h
                   + MARGIN_BASE)

    # Find max value for scaling (coerce None branch lengths to 0.0)
    def _v1(bm): return vals1.get(bm) or 0.0
    def _v2(bm): return (vals2.get(bm) or 0.0) * scale2

    max_val = max(
        max((_v1(bm) for bm in sorted_bm), default=0.0),
        max((_v2(bm) for bm in sorted_bm), default=0.0),
        0.001,
    )

    def scale(v):
        return (v / max_val) * MAX_BAR_H

    axis_y        = MARGIN_TOP + MAX_BAR_H
    hex_anchor_y  = MARGIN_TOP
    schematic_y0  = axis_y + AXIS_H + SCHEMATIC_GAP
    dot_r         = SCHEMATIC_ROW_H * 0.35

    svg = []
    svg.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{total_w}" height="{total_h}" '
        f'font-family="Helvetica,Arial,sans-serif">'
    )
    svg.append(
        f'<rect width="{total_w}" height="{total_h}" fill="white" stroke="none"/>'
    )

    # Axis line
    svg.append(
        f'<line x1="{MARGIN_LEFT}" y1="{axis_y}" '
        f'x2="{MARGIN_LEFT + n_pairs * PAIR_W}" y2="{axis_y}" '
        f'stroke="#333" stroke-width="1.5"/>'
    )

    # Y-axis scale ticks (50% and 100%)
    for frac in [0.5, 1.0]:
        ty = axis_y - scale(max_val * frac)
        svg.append(
            f'<line x1="{MARGIN_LEFT - 4}" y1="{ty:.1f}" '
            f'x2="{MARGIN_LEFT}" y2="{ty:.1f}" '
            f'stroke="#666" stroke-width="1"/>'
        )
        svg.append(
            f'<text x="{MARGIN_LEFT - 6}" y="{ty + 3:.1f}" '
            f'text-anchor="end" font-size="8" fill="#555">'
            f'{max_val * frac:.4f}</text>'
        )

    # Y-axis label
    mid_y = axis_y - MAX_BAR_H // 2
    svg.append(
        f'<text transform="rotate(-90,{MARGIN_LEFT - 60},{mid_y})" '
        f'x="{MARGIN_LEFT - 60}" y="{mid_y + 4}" '
        f'text-anchor="middle" font-size="{TITLE_FONT}" fill="#333">'
        f'Branch length</text>'
    )

    # Legend
    lx = MARGIN_LEFT
    ly = MARGIN_TOP - hex_label_h - (HEX_GAP if show_hex else 0) - 4
    if ly < 4:
        ly = 4
    svg.append(
        f'<rect x="{lx}" y="{ly - 7}" width="{BAR_W}" height="8" '
        f'fill="{COL1}" stroke="none"/>'
        f'<text x="{lx + BAR_W + 4}" y="{ly}" '
        f'font-size="{LABEL_FONT}" fill="#333">'
        f'{_xml_escape(label1)}</text>'
    )
    lx2 = lx + BAR_W + 4 + max(60, len(label1) * 7)
    svg.append(
        f'<rect x="{lx2}" y="{ly - 7}" width="{BAR_W}" height="8" '
        f'fill="{COL2}" stroke="none"/>'
        f'<text x="{lx2 + BAR_W + 4}" y="{ly}" '
        f'font-size="{LABEL_FONT}" fill="#333">'
        f'{_xml_escape(label2)}</text>'
    )

    # Bars
    for i, bm in enumerate(pairs):
        pair_x = MARGIN_LEFT + i * PAIR_W
        x1     = pair_x
        x2     = pair_x + BAR_W + BAR_GAP

        v1 = vals1.get(bm) or 0.0
        v2 = (vals2.get(bm) or 0.0) * scale2
        h1 = scale(v1)
        h2 = scale(v2)

        if h1 > 0:
            svg.append(
                f'<rect x="{x1}" y="{axis_y - h1:.1f}" '
                f'width="{BAR_W}" height="{h1:.1f}" '
                f'fill="{COL1}" stroke="none"/>'
            )
        if h2 > 0:
            svg.append(
                f'<rect x="{x2}" y="{axis_y - h2:.1f}" '
                f'width="{BAR_W}" height="{h2:.1f}" '
                f'fill="{COL2}" stroke="none"/>'
            )

    # Taxon labels
    for t_idx, name in enumerate(order):
        cy = schematic_y0 + t_idx * SCHEMATIC_ROW_H + SCHEMATIC_ROW_H // 2
        svg.append(
            f'<text x="{MARGIN_LEFT - 4}" y="{cy + 3}" '
            f'text-anchor="end" font-size="{LABEL_FONT}" fill="#222">'
            f'{_xml_escape(name)}</text>'
        )

    # Dot schematic (centered under each pair)
    for i, bm in enumerate(pairs):
        bin_str = bitmask_to_binary_str(bm, n_tax)
        cx = MARGIN_LEFT + i * PAIR_W + BAR_W + BAR_GAP // 2
        for t_idx, bit in enumerate(bin_str):
            cy   = schematic_y0 + t_idx * SCHEMATIC_ROW_H + SCHEMATIC_ROW_H // 2
            fill = '#222' if bit == '1' else 'white'
            svg.append(
                f'<circle cx="{cx}" cy="{cy}" r="{dot_r:.1f}" '
                f'fill="{fill}" stroke="#555" stroke-width="0.8"/>'
            )

    # Hex labels (centered under pair, reading upward)
    if show_hex:
        for i, bm in enumerate(pairs):
            hex_str = bitmask_to_jakobsen_hex(bm, n_tax)
            cx = MARGIN_LEFT + i * PAIR_W + BAR_W + BAR_GAP // 2
            svg.append(
                f'<text transform="rotate(-90,{cx},{hex_anchor_y})" '
                f'x="{cx}" y="{hex_anchor_y}" '
                f'text-anchor="end" font-size="{HEX_FONT}" fill="#444">'
                f'{_xml_escape(hex_str)}</text>'
            )

    svg.append('</svg>')
    with open(path, 'w') as f:
        f.write('\n'.join(svg))

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            'Compare internal branch lengths between two Newick trees, or '
            'compare a tree against Lento-plot support values from '
            'partcompat.py / sample_loci.py.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Modes (mutually exclusive):
  Tree-vs-tree:  -1 TREE1.nwk  -2 TREE2.nwk  --tax-order FILE
  Tree-vs-Lento: -1 TREE1.nwk  -l LENTO.tsv  --tax-order FILE

In tree-vs-Lento mode, support_norm values from the Lento TSV replace tree2
branch lengths.  The same normalization (sum-to-equal) is applied.

Examples:
  python3 compare_trees.py -1 tree1.nwk -2 tree2.nwk --tax-order taxa.txt -o cmp
  python3 compare_trees.py -1 tree1.nwk -l results.tsv --tax-order taxa.txt -o cmp
  python3 compare_trees.py -1 tree1.nwk -2 tree2.nwk --tax-order taxa.txt \\
      -o cmp --truncate 24 --no-normalize
        ''',
    )
    p.add_argument('-1', '--tree1',     required=True,
                   help='First Newick tree file')
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument('-2', '--tree2',
                     help='Second Newick tree file (tree-vs-tree mode)')
    src.add_argument('-l', '--lento',
                     help='Lento TSV from partcompat.py or sample_loci.py '
                          '(tree-vs-Lento mode; support_norm used as "lengths")')
    p.add_argument('--tax-order',       required=True,
                   help='Taxon order file (one name per line)')
    p.add_argument('-o', '--outfile',   required=True,
                   help='Output prefix; produces PREFIX.tsv and PREFIX.svg')
    p.add_argument('--no-normalize',    action='store_true',
                   help='Do not scale values to equal sums before plotting '
                        '(default: normalize so both totals are equal)')
    p.add_argument('--no-hex',          action='store_true',
                   help='Suppress Jakobsen hex labels above the plot')
    p.add_argument('--omit-biparts',
                   help='File of bipartitions to omit from SVG output '
                        '(same binary/hex format as --tree-biparts in '
                        'partcompat.py). Omitted bipartitions still appear '
                        'in the TSV. Bipartitions in this file absent from '
                        'the data are silently ignored.')
    p.add_argument('--no-svg',          action='store_true',
                   help='Suppress SVG output (TSV is always written)')
    p.add_argument('--truncate',        type=int, default=None,
                   metavar='N',
                   help='Also write PREFIX.trunc.svg showing only the first '
                        'N pairs (same ordering as the full plot; requires '
                        'SVG output)')
    return p.parse_args()

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args  = parse_args()
    order = read_tax_order(args.tax_order)
    n_tax = len(order)
    if n_tax == 0:
        sys.exit(f"ERROR: No taxa found in {args.tax_order}")
    print(f"Taxon order: {n_tax} taxa", file=sys.stderr)

    # ── Load tree 1 ───────────────────────────────────────────────────────────
    branches1 = load_tree(args.tree1, order, n_tax, 'Tree 1')
    total1    = sum(v for v in branches1.values() if v is not None)
    if total1 == 0:
        print("WARNING: Tree 1 total internal branch length is 0.",
              file=sys.stderr)

    # ── Load tree 2 or Lento TSV ──────────────────────────────────────────────
    lento_mode = args.lento is not None
    if lento_mode:
        branches2 = load_lento(args.lento, n_tax)
        label2    = 'Lento support'
        source2   = args.lento
    else:
        branches2 = load_tree(args.tree2, order, n_tax, 'Tree 2')
        label2    = 'Tree 2'
        source2   = args.tree2
    label1 = 'Tree 1'

    total2 = sum((v or 0.0) for v in branches2.values())
    if total2 == 0:
        print(f"WARNING: {label2} total is 0.", file=sys.stderr)

    # ── Build union and sort ──────────────────────────────────────────────────
    all_bm = set(branches1.keys()) | set(branches2.keys())

    def sort_key(bm):
        in1 = bm in branches1
        in2 = bm in branches2
        cat = 0 if (in1 and in2) else (1 if in1 else 2)
        l1  = (branches1.get(bm) or 0.0)
        return (cat, -l1)

    sorted_bm = sorted(all_bm, key=sort_key)

    # ── Write TSV ─────────────────────────────────────────────────────────────
    tsv_path = args.outfile + '.tsv'
    col2_name = 'support_norm' if lento_mode else 'length_tree2'
    prop2_name = 'prop_support' if lento_mode else 'prop_tree2'

    with open(tsv_path, 'w') as f:
        f.write('\t'.join([
            'bipartition_binary', 'jakobsen_hex',
            'length_tree1', col2_name,
            'prop_tree1', prop2_name,
            'present_in', 'grouping',
        ]) + '\n')

        for bm in sorted_bm:
            bin_str = bitmask_to_binary_str(bm, n_tax)
            hex_str = bitmask_to_jakobsen_hex(bm, n_tax)
            in1     = bm in branches1
            in2     = bm in branches2
            l1      = (branches1.get(bm) or 0.0)
            l2      = (branches2.get(bm) or 0.0)
            if l1 is None: l1 = 0.0
            if l2 is None: l2 = 0.0
            p1      = (l1 / total1) if total1 > 0 else 0.0
            p2      = (l2 / total2) if total2 > 0 else 0.0
            if in1 and in2:
                present = 'both'
            elif in1:
                present = 'tree1_only'
            else:
                present = f'{("lento" if lento_mode else "tree2")}_only'
            group1   = [order[i] for i, c in enumerate(bin_str) if c == '1']
            group0   = [order[i] for i, c in enumerate(bin_str) if c == '0']
            grouping = f"({','.join(group1)})|({','.join(group0)})"
            f.write('\t'.join([
                "'" + bin_str, hex_str,
                f'{l1:.8f}', f'{l2:.8f}',
                f'{p1:.8f}', f'{p2:.8f}',
                present, grouping,
            ]) + '\n')
    print(f"Wrote {tsv_path}", file=sys.stderr)

    # ── Write SVG(s) ──────────────────────────────────────────────────────────
    normalize  = not args.no_normalize
    show_hex   = not args.no_hex

    # Load omit set (SVG only; TSV is unaffected)
    omit_bitmasks = set()
    if args.omit_biparts:
        omit_bitmasks = read_tree_biparts(args.omit_biparts, n_tax)
        # Bipartitions in omit file absent from data are silently ignored
        n_omit_present = sum(1 for bm in omit_bitmasks if bm in all_bm)
        print(f"Omit filter: {len(omit_bitmasks)} supplied, "
              f"{n_omit_present} present in data (will be excluded from SVG)",
              file=sys.stderr)

    # Apply omit filter to the ordered list used for SVG
    svg_bm = [bm for bm in sorted_bm if bm not in omit_bitmasks]

    if not args.no_svg:
        svg_path = args.outfile + '.svg'
        write_comparison_svg(
            svg_path, svg_bm, branches1, branches2,
            n_tax, order, label1, label2,
            normalize=normalize, show_hex=show_hex,
        )
        print(f"Wrote {svg_path}", file=sys.stderr)

        if args.truncate is not None:
            trunc_path = args.outfile + '.trunc.svg'
            write_comparison_svg(
                trunc_path, svg_bm, branches1, branches2,
                n_tax, order, label1, label2,
                normalize=normalize, show_hex=show_hex,
                max_cols=args.truncate,
            )
            print(f"Wrote {trunc_path}", file=sys.stderr)

    # ── Summary ───────────────────────────────────────────────────────────────
    n_both  = sum(1 for bm in all_bm if bm in branches1 and bm in branches2)
    n_only1 = sum(1 for bm in all_bm if bm in branches1 and bm not in branches2)
    n_only2 = sum(1 for bm in all_bm if bm not in branches1 and bm in branches2)

    print(f"\nTree 1: {args.tree1}")
    print(f"  Non-trivial internal branches: {len(branches1)}")
    print(f"  Total internal branch length:  {total1:.6f}")
    print(f"\n{label2}: {source2}")
    print(f"  Bipartitions / branches:       {len(branches2)}")
    print(f"  Total:                         {total2:.6f}")
    print(f"\nShared bipartitions:             {n_both}")
    print(f"Tree 1 only:                     {n_only1}")
    print(f"{label2} only:                   {n_only2}")
    print(f"Total unique bipartitions:       {len(all_bm)}")
    if normalize and total1 > 0 and total2 > 0:
        print(f"\nNormalization scale factor ({label2}): {total1/total2:.6f}")


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
