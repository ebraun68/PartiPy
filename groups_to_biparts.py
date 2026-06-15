#!/usr/bin/env python3
"""
groups_to_biparts.py - Convert manually-specified taxon groupings to
                       bipartition codes (binary string and Jakobsen hex),
                       optionally with support/conflict values from a TSV.

Reads a taxon order file and a file of manually-written bipartitions.  Each
line may specify either a FULL bipartition:

    tax1,tax2,tax3|tax4,tax5,tax6

or a HALF bipartition (no '|'):

    tax1,tax2,tax3

in which case the listed taxa become one side and every other taxon in
--tax-order becomes the other side.  Either way, every named taxon must
appear in --tax-order, with no duplicates; for full bipartitions every
taxon in --tax-order must appear on exactly one side.

For every line, the program echoes the CANONICAL grouping (both sides, in
--tax-order order, derived from the resulting bitmask) alongside the binary
string and Jakobsen hex -- this lets a half-specified line be checked at a
glance, since the complementary side is exactly what gets echoed.

With --tsv FILE, support_raw/conflict_raw/support_norm/conflict_norm values
for each bipartition are looked up from a partcompat.py / sample_loci.py /
analyze_loci.py TSV and added to the output (zeros if the bipartition is
absent from the TSV).

Useful for compiling a focused set of candidate bipartitions (e.g. competing
resolutions of a particular node across several trees) for use with
--tree-biparts, --plot-biparts, --bipart-order, --omit-biparts, or
id_compat.py's --bipart-file -- without having to construct or edit a Newick
tree for each alternative, and for directly inspecting the support/conflict
values associated with each candidate.

Usage:
  python3 groups_to_biparts.py --tax-order FILE GROUPS_FILE [-o PREFIX] [--tsv TSV]
"""

import argparse
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

try:
    from partcompat import (
        read_tax_order,
        _canonicalize,
        _coded_to_bitmask,
        bitmask_to_binary_str,
        bitmask_to_jakobsen_hex,
    )
except ImportError as e:
    sys.exit(f"ERROR: Cannot import partcompat.py from {_HERE}: {e}")

try:
    from plot_lento import read_lento_tsv
except ImportError as e:
    sys.exit(f"ERROR: Cannot import plot_lento.py from {_HERE}: {e}")

# ---------------------------------------------------------------------------
# Groups file parsing
# ---------------------------------------------------------------------------

def read_groups_file(path, order, n_tax):
    """
    Read a groups file: one bipartition per line, either a FULL bipartition
    (tax1,tax2,...|taxN,...) or a HALF bipartition (tax1,tax2,...; the
    complement within `order` becomes the other side).  '#' comments and
    blank lines are ignored.

    Validation (full bipartitions): every taxon in `order` must appear on
    exactly one side (no missing, no duplicated, no unrecognized names).

    Validation (half bipartitions): every listed taxon must be in `order`,
    no duplicates, and the list must be neither empty nor the entire taxon
    set (either of which would leave one side empty).

    Every line is checked, and ALL problems found across the whole file are
    reported together (not just the first); the program then exits with a
    single summary error message listing every problematic line.

    Returns a list of (bitmask_int, line_no, raw_line, is_half) tuples, in
    file order, if and only if no problems were found.
    """
    order_set = set(order)
    results = []
    errors   = []   # list of (line_no, raw_line, [problem messages])

    with open(path) as f:
        for line_no, raw_line in enumerate(f, start=1):
            line = raw_line.strip()
            if not line or line.startswith('#'):
                continue

            problems = []

            if '|' in line:
                # -- Full bipartition --------------------------------------
                parts = line.split('|')
                if len(parts) != 2:
                    problems.append(
                        f"expected exactly one '|' separator, "
                        f"found {len(parts) - 1}"
                    )
                    errors.append((line_no, line, problems))
                    continue

                side0 = [t.strip() for t in parts[0].split(',') if t.strip()]
                side1 = [t.strip() for t in parts[1].split(',') if t.strip()]
                all_names = side0 + side1

                unrecognized = [t for t in all_names if t not in order_set]
                if unrecognized:
                    problems.append(
                        f"taxon name(s) not found in --tax-order: "
                        f"{', '.join(unrecognized)}"
                    )

                seen = set()
                duplicates = set()
                for t in all_names:
                    if t in seen:
                        duplicates.add(t)
                    seen.add(t)
                if duplicates:
                    problems.append(
                        f"taxon name(s) appear more than once: "
                        f"{', '.join(sorted(duplicates))}"
                    )

                missing = order_set - seen
                if missing:
                    problems.append(
                        f"taxon name(s) from --tax-order missing from this "
                        f"line: {', '.join(sorted(missing))}"
                    )

                if problems:
                    errors.append((line_no, line, problems))
                    continue

                side1_set = set(side1)
                raw = [1 if t in side1_set else 0 for t in order]
                is_half = False

            else:
                # -- Half bipartition ---------------------------------------
                side = [t.strip() for t in line.split(',') if t.strip()]

                unrecognized = [t for t in side if t not in order_set]
                if unrecognized:
                    problems.append(
                        f"taxon name(s) not found in --tax-order: "
                        f"{', '.join(unrecognized)}"
                    )

                seen = set()
                duplicates = set()
                for t in side:
                    if t in seen:
                        duplicates.add(t)
                    seen.add(t)
                if duplicates:
                    problems.append(
                        f"taxon name(s) appear more than once: "
                        f"{', '.join(sorted(duplicates))}"
                    )

                if not problems:
                    if len(seen) == n_tax:
                        problems.append(
                            "lists all taxa in --tax-order; the "
                            "complementary side would be empty"
                        )

                if problems:
                    errors.append((line_no, line, problems))
                    continue

                side_set = seen
                raw = [1 if t in side_set else 0 for t in order]
                is_half = True

            # Build canonical bitmask
            canon = _canonicalize(raw)
            bm = _coded_to_bitmask(canon)

            results.append((bm, line_no, line, is_half))

    if errors:
        msg_lines = [
            f"ERROR: {len(errors)} line(s) in {path} have problems:"
        ]
        for line_no, line, problems in errors:
            msg_lines.append(f"  Line {line_no}: '{line}'")
            for problem in problems:
                msg_lines.append(f"    - {problem}")
        sys.exit('\n'.join(msg_lines))

    if not results:
        sys.exit(f"ERROR: No bipartitions found in {path}")

    return results


def canonical_groups(bm, order, n_tax):
    """
    Given a canonical bitmask, return (group0, group1) lists of taxon names
    in --tax-order order, where group0 corresponds to bit '0' and group1 to
    bit '1' in the canonical binary string.
    """
    bin_str = bitmask_to_binary_str(bm, n_tax)
    group0 = [order[i] for i, c in enumerate(bin_str) if c == '0']
    group1 = [order[i] for i, c in enumerate(bin_str) if c == '1']
    return group0, group1

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            'Convert manually-specified taxon groupings to bipartition '
            'codes (binary string and Jakobsen hex), optionally with '
            'support/conflict values from a TSV.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Groups file format: one bipartition per line.  '#' comments and blank lines
are ignored.

  Full bipartition (comma-separated taxa on each side of '|'):
    tax1,tax2,tax3|tax4,tax5,tax6
  Every taxon in --tax-order must appear on exactly one side.

  Half bipartition (no '|'; the complement within --tax-order becomes the
  other side):
    tax1,tax2,tax3
  Every listed taxon must be in --tax-order, with no duplicates, and the
  list must not be the entire taxon set.

For every line, the CANONICAL grouping (both sides, in --tax-order order,
derived from the resulting bitmask) is echoed -- this lets a half-specified
line be checked at a glance.

Example:
  python3 groups_to_biparts.py --tax-order taxa.txt groups.txt -o candidates \\
      --tsv results.tsv
        ''',
    )
    p.add_argument('groups_file',
                   help='File of manually-specified taxon groupings')
    p.add_argument('--tax-order', required=True,
                   help='Taxon order file (one name per line)')
    p.add_argument('-o', '--outfile',
                   help='Optional output prefix; writes PREFIX.biparts '
                        '(binary strings, one per line, with #-comment '
                        'lines giving hex, grouping, and -- if --tsv is '
                        'given -- support/conflict values) for use with '
                        '--tree-biparts, --plot-biparts, --bipart-order, '
                        '--omit-biparts, or id_compat.py --bipart-file')
    p.add_argument('--tsv',
                   help='TSV file from partcompat.py, sample_loci.py '
                        '(summary.tsv), or analyze_loci.py; if given, '
                        'support_raw/conflict_raw/support_norm/conflict_norm '
                        'for each bipartition are looked up and reported '
                        '(zeros if the bipartition is absent from the TSV)')
    return p.parse_args()

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if not os.path.isfile(args.groups_file):
        sys.exit(f"ERROR: Groups file not found: {args.groups_file}")

    order = read_tax_order(args.tax_order)
    n_tax = len(order)
    if n_tax == 0:
        sys.exit(f"ERROR: No taxa found in {args.tax_order}")
    order_set = set(order)
    print(f"Taxon order: {n_tax} taxa", file=sys.stderr)

    results = read_groups_file(args.groups_file, order, n_tax)
    n_half = sum(1 for _, _, _, is_half in results if is_half)
    print(f"Read {len(results)} bipartition(s) from {args.groups_file} "
          f"({n_half} half-specified)", file=sys.stderr)

    # -- Optional TSV lookup --------------------------------------------------
    tsv_data = None
    if args.tsv:
        if not os.path.isfile(args.tsv):
            sys.exit(f"ERROR: TSV file not found: {args.tsv}")
        support, conflict, support_raw, conflict_raw = read_lento_tsv(
            args.tsv, n_tax, order_set)
        tsv_data = (support, conflict, support_raw, conflict_raw)
        print(f"Read {len(support)} bipartition(s) from {args.tsv}",
              file=sys.stderr)

    # -- stdout report ----------------------------------------------------------
    header = ['bipartition_binary', 'jakobsen_hex']
    if tsv_data:
        header += ['support_raw', 'conflict_raw',
                   'support_norm', 'conflict_norm']
    header += ['grouping', 'input_line']
    print('\t'.join(header))

    for bm, line_no, raw_line, is_half in results:
        bin_str = bitmask_to_binary_str(bm, n_tax)
        hex_str = bitmask_to_jakobsen_hex(bm, n_tax)
        group0, group1 = canonical_groups(bm, order, n_tax)
        grouping = f"({','.join(group0)})|({','.join(group1)})"

        row = [bin_str, hex_str]
        if tsv_data:
            support, conflict, support_raw, conflict_raw = tsv_data
            row += [
                f"{support_raw.get(bm, 0.0):.4f}",
                f"{conflict_raw.get(bm, 0.0):.4f}",
                f"{support.get(bm, 0.0):.4f}",
                f"{conflict.get(bm, 0.0):.4f}",
            ]
        row += [grouping, raw_line]
        print('\t'.join(row))

    # -- optional PREFIX.biparts ------------------------------------------------
    if args.outfile:
        out_path = args.outfile + '.biparts'
        with open(out_path, 'w') as f:
            f.write(f"# Bipartitions from {args.groups_file}\n")
            f.write(f"# {n_tax} taxa in bit-string order:\n")
            for i, t in enumerate(order, start=1):
                f.write(f"#   bit {i:2d}: {t}\n")
            f.write("#\n")
            for bm, line_no, raw_line, is_half in results:
                bin_str = bitmask_to_binary_str(bm, n_tax)
                hex_str = bitmask_to_jakobsen_hex(bm, n_tax)
                group0, group1 = canonical_groups(bm, order, n_tax)
                grouping = f"({','.join(group0)})|({','.join(group1)})"
                line = f"#   {bin_str}\thex {hex_str}\t{grouping}"
                if tsv_data:
                    support, conflict, support_raw, conflict_raw = tsv_data
                    line += (
                        f"\tsupport_raw={support_raw.get(bm, 0.0):.4f}"
                        f"\tconflict_raw={conflict_raw.get(bm, 0.0):.4f}"
                        f"\tsupport_norm={support.get(bm, 0.0):.4f}"
                        f"\tconflict_norm={conflict.get(bm, 0.0):.4f}"
                    )
                if is_half:
                    line += "\t[half-specified]"
                f.write(line + "\n")
            f.write("#\n")
            for bm, line_no, raw_line, is_half in results:
                f.write(f"{bitmask_to_binary_str(bm, n_tax)}\n")
        print(f"Wrote {out_path}", file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
