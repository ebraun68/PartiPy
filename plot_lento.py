#!/usr/bin/env python3
"""
plot_lento.py - Generate a LentoPlot SVG from a partcompat.py / sample_loci.py TSV.

Reads pre-computed support and conflict values from a TSV file and produces a
LentoPlot SVG without rerunning any analysis.  Accepts TSV files from
partcompat.py, sample_loci.py (summary.tsv), or analyze_loci.py.

Usage:
  python3 plot_lento.py -i PREFIX.tsv --tax-order FILE -o PREFIX [options]
"""

import argparse
import csv
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

try:
    from partcompat import (
        read_tax_order,
        read_tree_biparts,
        read_bipart_order,
        build_column_order,
        write_svg,
        bitmask_to_binary_str,
        bitmask_to_jakobsen_hex,
        _write_tree_bipart_tsv,
        sort_key,
    )
except ImportError as e:
    sys.exit(f"ERROR: Cannot import partcompat.py from {_HERE}: {e}")

# ---------------------------------------------------------------------------
# TSV reader
# ---------------------------------------------------------------------------

def read_lento_tsv(path, n_tax, tax_order_set):
    """
    Read a partcompat.py / sample_loci.py TSV file.

    Returns:
      support      : {bitmask_int: float}  support_norm values
      conflict     : {bitmask_int: float}  conflict_norm values
      support_raw  : {bitmask_int: float}  support_raw values
      conflict_raw : {bitmask_int: float}  conflict_raw values

    Performs two validation checks:
      1. Every binary string has length == n_tax.
      2. Every taxon name in the grouping column is present in tax_order_set.

    Skips footer rows (rows whose bipartition_binary is not a binary string).
    """
    support      = {}
    conflict     = {}
    support_raw  = {}
    conflict_raw = {}

    len_errors   = []
    taxon_errors = set()  # accumulates all taxon names seen in grouping column
    n_rows       = 0

    try:
        with open(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                bp = row.get('bipartition_binary', '').strip().lstrip("'")
                # Skip non-data rows (footer summary rows, blank lines)
                if not bp or not all(c in '01' for c in bp):
                    continue

                n_rows += 1

                # Check 1: bit string length
                if len(bp) != n_tax:
                    len_errors.append((bp, len(bp)))
                    continue

                bm = int(bp, 2)

                try:
                    sup_n = float(row.get('support_norm',  0) or 0)
                    con_n = float(row.get('conflict_norm', 0) or 0)
                    sup_r = float(row.get('support_raw',   0) or 0)
                    con_r = float(row.get('conflict_raw',  0) or 0)
                except (ValueError, KeyError):
                    continue

                support[bm]      = sup_n
                conflict[bm]     = con_n
                support_raw[bm]  = sup_r
                conflict_raw[bm] = con_r

                # Collect taxon names from grouping column for validation
                grouping = row.get('grouping', '')
                if grouping:
                    names = [n.strip('() ') for n in
                             grouping.replace('|', ',').split(',')]
                    for name in names:
                        if name:
                            taxon_errors.add(name)   # collect ALL names seen

    except FileNotFoundError:
        sys.exit(f"ERROR: TSV file not found: {path}")
    except Exception as e:
        sys.exit(f"ERROR: Cannot read TSV file '{path}': {e}")

    # Report validation issues
    if len_errors:
        print(f"ERROR: {len(len_errors)} bipartition(s) in TSV have wrong "
              f"length (expected {n_tax}):", file=sys.stderr)
        for bp, ln in len_errors[:5]:
            print(f"  '{bp}' has length {ln}", file=sys.stderr)
        if len(len_errors) > 5:
            print(f"  ... and {len(len_errors) - 5} more", file=sys.stderr)
        sys.exit("ERROR: TSV and --tax-order are inconsistent. "
                 "Did you use the correct taxon order file?")

    # Taxon name cross-check using all names collected from grouping column
    tsv_taxa    = taxon_errors   # reuse set; now holds all names from TSV
    in_tsv_not_order = tsv_taxa - tax_order_set
    in_order_not_tsv = tax_order_set - tsv_taxa

    if in_tsv_not_order:
        print(f"ERROR: {len(in_tsv_not_order)} taxon name(s) found in TSV "
              f"grouping column but absent from --tax-order:",
              file=sys.stderr)
        for name in sorted(in_tsv_not_order):
            print(f"  {name}", file=sys.stderr)
        sys.exit("ERROR: --tax-order does not match the file used to produce "
                 "this TSV. Please supply the correct taxon order file.")

    if in_order_not_tsv:
        print(f"NOTE: {len(in_order_not_tsv)} taxon name(s) in --tax-order "
              f"do not appear in any bipartition in this TSV "
              f"(expected if these taxa only appear in trivial bipartitions):",
              file=sys.stderr)
        for name in sorted(in_order_not_tsv):
            print(f"  {name}", file=sys.stderr)

    if not support:
        sys.exit(f"ERROR: No valid bipartition rows found in '{path}'. "
                 f"Check that the file is a partcompat.py / sample_loci.py TSV.")

    return support, conflict, support_raw, conflict_raw

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            'Generate a LentoPlot SVG from a partcompat.py, sample_loci.py, '
            'or analyze_loci.py TSV file without rerunning any analysis.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Accepted TSV sources:
  partcompat.py output (PREFIX.tsv)
  sample_loci.py summary output (summary.tsv)
  analyze_loci.py per-locus output (LOCUS.tsv)

The normalized support and conflict values from the TSV are used directly,
preserving the normalization applied during the original analysis.  The goal
is figure control, not reanalysis.

Example:
  python3 plot_lento.py \\
      -i results.tsv \\
      --tax-order taxa.txt \\
      -o results_fig \\
      -t tree.biparts \\
      -p tree.biparts \\
      --truncate 24
        ''',
    )
    p.add_argument('-i', '--infile',      required=True,
                   help='Input TSV file from partcompat.py, sample_loci.py, '
                        'or analyze_loci.py')
    p.add_argument('--tax-order',         required=True,
                   help='Taxon order file (one name per line); must match '
                        'the file used when producing the TSV')
    p.add_argument('-o', '--outfile',     required=True,
                   help='Output prefix; produces PREFIX.svg '
                        '(and PREFIX.trunc.svg with --truncate)')
    p.add_argument('-t', '--tree-biparts',
                   help='Reference tree bipartitions file (outline bars in SVG)')
    p.add_argument('-p', '--plot-biparts',
                   help='Display filter: only these bipartitions appear in the SVG')
    p.add_argument('-b', '--bipart-order',
                   help='Column ordering file')
    p.add_argument('--omit-biparts',
                   help='File of bipartitions to exclude from the SVG '
                        '(same binary/hex format as --tree-biparts; '
                        'useful for omitting outgroup branches)')
    p.add_argument('--truncate',          type=int, default=None,
                   metavar='N',
                   help='Also write PREFIX.trunc.svg with only the first N columns')
    p.add_argument('--no-hex',            action='store_true',
                   help='Suppress Jakobsen hex labels above the plot')
    return p.parse_args()

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # ── Validate files exist ──────────────────────────────────────────────────
    if not os.path.isfile(args.infile):
        sys.exit(f"ERROR: Input TSV not found: {args.infile}")

    # ── Taxon order ───────────────────────────────────────────────────────────
    order = read_tax_order(args.tax_order)
    n_tax = len(order)
    if n_tax == 0:
        sys.exit(f"ERROR: No taxa found in {args.tax_order}")
    order_set = set(order)
    print(f"Taxon order: {n_tax} taxa from {args.tax_order}", file=sys.stderr)

    # ── Read TSV ──────────────────────────────────────────────────────────────
    support, conflict, support_raw, conflict_raw = read_lento_tsv(
        args.infile, n_tax, order_set)
    n_biparts = len(support)
    print(f"Read {n_biparts} bipartitions from {args.infile}", file=sys.stderr)

    # ── Load bipartition control files ────────────────────────────────────────
    tree_bitmasks  = None
    plot_bitmasks  = None
    order_bitmasks = None
    omit_bitmasks  = set()

    if args.tree_biparts:
        tree_bitmasks = read_tree_biparts(args.tree_biparts, n_tax)
        absent  = tree_bitmasks - set(support.keys())
        present = tree_bitmasks &  set(support.keys())
        if absent:
            print(f"WARNING: {len(absent)} --tree-biparts bipartition(s) absent "
                  f"from TSV (no support; will not appear as outline bars):",
                  file=sys.stderr)
            for bm in sorted(absent, key=lambda b: bitmask_to_binary_str(b, n_tax)):
                print(f"  {bitmask_to_binary_str(bm, n_tax)}, "
                      f"hex {bitmask_to_jakobsen_hex(bm, n_tax)}",
                      file=sys.stderr)
        print(f"Tree bipartitions: {len(tree_bitmasks)} supplied, "
              f"{len(present)} present in TSV", file=sys.stderr)

    if args.plot_biparts:
        plot_bitmasks = read_tree_biparts(args.plot_biparts, n_tax)
        absent  = plot_bitmasks - set(support.keys())
        present = plot_bitmasks &  set(support.keys())
        if absent:
            print(f"Plot filter: {len(absent)} bipartition(s) with zero support "
                  f"in TSV — these will appear as conflict-only bars",
                  file=sys.stderr)
            # Add zero-support entries so they appear in the plot
            for bm in absent:
                support[bm]      = 0.0
                conflict[bm]     = 0.0
                support_raw[bm]  = 0.0
                conflict_raw[bm] = 0.0
        print(f"Plot filter: {len(plot_bitmasks)} supplied, "
              f"{len(present)} with support in TSV", file=sys.stderr)

    if args.bipart_order:
        order_bitmasks = read_bipart_order(args.bipart_order, n_tax)
        absent_order   = [bm for bm in order_bitmasks if bm not in support]
        if absent_order:
            absent_strs = [bitmask_to_binary_str(bm, n_tax) for bm in absent_order]
            print(f"WARNING: {len(absent_order)} --bipart-order bipartition(s) "
                  f"absent from TSV (will be skipped in column ordering): "
                  f"{', '.join(absent_strs)}", file=sys.stderr)
        present_order = [bm for bm in order_bitmasks if bm in support]
        print(f"Bipart order: {len(order_bitmasks)} supplied, "
              f"{len(present_order)} present in TSV "
              f"(these will be the first {len(present_order)} columns)",
              file=sys.stderr)

    if args.omit_biparts:
        omit_bitmasks = read_tree_biparts(args.omit_biparts, n_tax)
        n_omit = sum(1 for bm in omit_bitmasks if bm in support)
        print(f"Omit filter: {len(omit_bitmasks)} supplied, "
              f"{n_omit} present in TSV (will be excluded from SVG)",
              file=sys.stderr)

    # ── Tree bipartition position reporting ───────────────────────────────────
    if tree_bitmasks:
        full_order = build_column_order(
            support, support, conflict,   # support used as both support and norm
            plot_bitmasks=plot_bitmasks,
            order_bitmasks=order_bitmasks,
        )
        pos_map      = {bm: i + 1 for i, bm in enumerate(full_order)}
        tree_in_plot = {bm: pos_map[bm] for bm in tree_bitmasks if bm in pos_map}
        if tree_in_plot:
            max_pos = max(tree_in_plot.values())
            max_bm  = next(bm for bm, p in tree_in_plot.items() if p == max_pos)
            if args.truncate is not None:
                n_in_trunc = sum(1 for p in tree_in_plot.values()
                                 if p <= args.truncate)
                print(f"Tree bipartitions in truncated plot "
                      f"(first {args.truncate} columns): {n_in_trunc}")
            print(f"Column of last tree bipartition = {max_pos} "
                  f"({bitmask_to_binary_str(max_bm, n_tax)}, "
                  f"hex {bitmask_to_jakobsen_hex(max_bm, n_tax)})")

        tbp_path = args.outfile + '.tree_bipart.tsv'
        _write_tree_bipart_tsv(
            tbp_path, tree_bitmasks,
            support_raw, conflict_raw,   # raw counts
            support, conflict,           # normalized (=TSV values)
            n_tax, order, pos_map,
        )
        print(f"Wrote {tbp_path}", file=sys.stderr)

    # ── Build column order then apply omit filter ─────────────────────────────
    full_sorted = build_column_order(
        support, support, conflict,
        plot_bitmasks=plot_bitmasks,
        order_bitmasks=order_bitmasks,
    )
    svg_bm = [bm for bm in full_sorted if bm not in omit_bitmasks]

    # ── Write SVG(s) ──────────────────────────────────────────────────────────
    # write_svg signature: (path, support, conflict_raw, support_norm, conflict_norm, ...)
    # We pass the TSV values directly: support_norm = support_norm column from TSV,
    # conflict_norm = conflict_norm column from TSV, using the already-normalized values.
    svg_path = args.outfile + '.svg'

    # We call write_svg directly but need to pass a pre-filtered bitmask list.
    # write_svg internally calls build_column_order, so we monkey-patch by
    # passing plot_bitmasks that exactly equals our svg_bm set.
    svg_plot_set = set(svg_bm) if omit_bitmasks else plot_bitmasks

    write_svg(svg_path, support, conflict_raw, support, conflict,
              n_tax, order,
              tree_bitmasks=tree_bitmasks,
              plot_bitmasks=svg_plot_set,
              order_bitmasks=order_bitmasks,
              show_hex=(not args.no_hex))
    print(f"Wrote {svg_path}", file=sys.stderr)

    if args.truncate is not None:
        trunc_path = args.outfile + '.trunc.svg'
        write_svg(trunc_path, support, conflict_raw, support, conflict,
                  n_tax, order,
                  tree_bitmasks=tree_bitmasks,
                  plot_bitmasks=svg_plot_set,
                  order_bitmasks=order_bitmasks,
                  show_hex=(not args.no_hex),
                  max_cols=args.truncate)
        print(f"Wrote {trunc_path}", file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
