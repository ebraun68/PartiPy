#!/usr/bin/env python3
"""
id_compat.py - Identify alignment sites compatible with specified bipartitions.

For one or more user-specified bipartitions, finds the parsimony-informative
sites (under the chosen coding scheme) that support each bipartition - either
fully (no missing data) or partially (exactly one missing taxon, contributing
the standard +0.5 weight).

Optionally annotates sites using a per-character annotation file, and/or
writes a FASTA file containing the union of all sites supporting any of the
queried bipartitions (as original nucleotides or as recoded binary characters).

Usage:
  python3 id_compat.py -i ALIGNMENT --tax-order FILE -o PREFIX \\
      (--bipart BIPART [--bipart ...] | --bipart-file FILE) [options]
"""

import argparse
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

try:
    from partcompat import (
        detect_and_parse,
        read_tax_order,
        read_tree_biparts,
        is_binary_data,
        process_alignment,
        resolve_mode,
        CODON_POS_MAP,
        _canonicalize,
        _coded_to_bitmask,
        bitmask_to_binary_str,
        bitmask_to_jakobsen_hex,
    )
except ImportError as e:
    sys.exit(f"ERROR: Cannot import partcompat.py from {_HERE}: {e}")

# ---------------------------------------------------------------------------
# Bipartition parsing (single token: binary string or Jakobsen hex)
# ---------------------------------------------------------------------------

def parse_bipart_token(token, n_tax):
    """
    Parse a single bipartition token (binary string or Jakobsen hex) into a
    canonical bitmask integer.  Exits with an error on malformed input.
    """
    token = token.strip()
    full_mask = (1 << n_tax) - 1

    if all(c in '01' for c in token):
        if len(token) != n_tax:
            sys.exit(f"ERROR: bipartition '{token}' has length {len(token)}, "
                     f"expected {n_tax}")
        raw   = [int(c) for c in token]
        canon = _canonicalize(raw)
        return _coded_to_bitmask(canon)
    else:
        hex_tok = token.upper().lstrip('0X') or '0'
        try:
            jak_val = int(hex_tok, 16)
        except ValueError:
            sys.exit(f"ERROR: cannot parse '{token}' as binary or hexadecimal")
        canon_val = jak_val ^ full_mask
        if canon_val.bit_length() > n_tax:
            sys.exit(f"ERROR: hex '{token}' exceeds {n_tax} bits")
        from partcompat import _bitmask_to_coded
        canon_coded = _canonicalize(_bitmask_to_coded(canon_val, n_tax))
        return _coded_to_bitmask(canon_coded)

# ---------------------------------------------------------------------------
# Annotation file parsing
# ---------------------------------------------------------------------------

def read_annotations(path, aln_len):
    """
    Read an annotation file.  Two formats are auto-detected:

      1. One annotation per line, with exactly aln_len lines (1:1 with
         alignment positions, 1-based).
      2. Two tab-separated columns: <1-based position> <tab> <annotation>.
         Positions not listed default to "unannotated".

    Returns a dict {1-based position: annotation_string}.
    """
    with open(path) as f:
        lines = [line.rstrip('\n') for line in f]
    # Drop trailing blank lines for the one-per-line check
    stripped = [l for l in lines if l != '']

    # Detect two-column format: every non-blank line has a tab and the
    # first field is an integer
    is_two_col = True
    for line in stripped:
        parts = line.split('\t', 1)
        if len(parts) != 2:
            is_two_col = False
            break
        try:
            int(parts[0].strip())
        except ValueError:
            is_two_col = False
            break

    annotations = {}
    if is_two_col and stripped:
        for line in stripped:
            pos_str, annot = line.split('\t', 1)
            pos = int(pos_str.strip())
            annotations[pos] = annot.strip()
        # Fill unannotated positions
        for pos in range(1, aln_len + 1):
            if pos not in annotations:
                annotations[pos] = 'unannotated'
    else:
        # One-per-line format: must match alignment length exactly
        if len(lines) != aln_len:
            sys.exit(
                f"ERROR: Annotation file '{path}' has {len(lines)} line(s) "
                f"but the alignment has {aln_len} character(s). "
                f"One-per-line annotation files must have exactly one line "
                f"per alignment position."
            )
        for pos, annot in enumerate(lines, start=1):
            annotations[pos] = annot.strip() if annot.strip() else 'unannotated'

    return annotations

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            'Identify alignment sites compatible with specified bipartitions.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Bipartition specification (at least one required):
  --bipart BIPART      Single bipartition (binary string or Jakobsen hex).
                        May be repeated for multiple bipartitions.
  --bipart-file FILE    File listing bipartitions (binary or hex, one per
                        line; same format as -t/--tree-biparts elsewhere).

Output:
  PREFIX.txt            Taxon list, then one block per bipartition listing
                        site positions (1-based) with full support and
                        partial (1-missing-taxon) support.  Partial-support
                        positions are marked with '*' in annotated lists.
  PREFIX.fasta          (with --fasta-output) union of all sites supporting
                        any queried bipartition, as original nucleotides
                        (default) or recoded binary characters.

Example:
  python3 id_compat.py -i alignment.fasta --tax-order taxa.txt \\
      -o results --bipart-file tree.biparts --mode RY \\
      --annotation positions.txt --fasta-output original
        ''',
    )
    p.add_argument('-i', '--infile',  required=True,
                   help='Input alignment (FASTA, relaxed PHYLIP, or Nexus)')
    p.add_argument('--tax-order',     required=True,
                   help='Taxon order file (one name per line)')
    p.add_argument('-o', '--outfile', required=True,
                   help='Output prefix; produces PREFIX.txt '
                        '(and PREFIX.fasta with --fasta-output)')

    bgrp = p.add_argument_group('bipartition specification')
    bgrp.add_argument('--bipart', action='append', default=[],
                      metavar='BIPART',
                      help='Single bipartition (binary string or Jakobsen '
                           'hex). May be repeated.')
    bgrp.add_argument('--bipart-file',
                      help='File listing bipartitions (binary or hex, one '
                           'per line, # comments allowed)')

    p.add_argument('-m', '--mode', default='RY',
                   help='Site coding mode (default: RY). Same options as '
                        'partcompat.py: RY, KM, WS, perfectRY, perfectKM, '
                        'perfectWS, transitions, transitionsAG, '
                        'transitionsCT, 2state, partimatrix (and reversed '
                        'pairs).')
    p.add_argument('--codon-positions', default='all',
                   choices=['all', '1', '2', '3', '12', '13', '23', '123'],
                   help='Restrict to specific codon positions (default: all)')
    p.add_argument('--annotation',
                   help='Annotation file: either one annotation per line '
                        '(one line per alignment position), or two '
                        'tab-separated columns (1-based position, '
                        'annotation). Unlisted positions in the two-column '
                        'format are labeled "unannotated".')
    p.add_argument('--verbose-outfile', action='store_true',
                   help='Include taxon groupings (tax1,tax2,...|tax3,tax4,...) '
                        'for each bipartition in PREFIX.txt')
    p.add_argument('--fasta-output', choices=['binary', 'original'],
                   default=None,
                   help='Write PREFIX.fasta containing the union of all '
                        'sites supporting any queried bipartition, as '
                        '"original" nucleotides (default representation if '
                        'flag given without further qualification) or '
                        '"binary" recoded characters.')
    return p.parse_args()

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if not args.bipart and not args.bipart_file:
        sys.exit("ERROR: At least one of --bipart or --bipart-file is required.")

    if not os.path.isfile(args.infile):
        sys.exit(f"ERROR: Input file not found: {args.infile}")

    # ── Taxon order ───────────────────────────────────────────────────────────
    order = read_tax_order(args.tax_order)
    n_tax = len(order)
    if n_tax == 0:
        sys.exit(f"ERROR: No taxa found in {args.tax_order}")
    print(f"Taxon order: {n_tax} taxa", file=sys.stderr)

    # ── Read alignment ────────────────────────────────────────────────────────
    with open(args.infile) as f:
        text = f.read()
    file_order, file_seqs = detect_and_parse(text)
    if not file_order:
        sys.exit(f"ERROR: No sequences found in {args.infile}")

    aln_len = len(file_seqs[file_order[0]])
    print(f"Read {len(file_order)} sequences, length {aln_len}", file=sys.stderr)

    # Reorder / pad to tax_order
    missing_taxa = [t for t in order if t not in file_seqs]
    if missing_taxa:
        print(f"WARNING: {len(missing_taxa)} taxa in --tax-order absent from "
              f"alignment (padded with missing data): {missing_taxa}",
              file=sys.stderr)
    seqs = {t: file_seqs.get(t, '?' * aln_len) for t in order}

    binary = is_binary_data(seqs)
    mode   = resolve_mode(args.mode)
    if binary:
        print("Detected binary data; --mode setting ignored.", file=sys.stderr)
    codon_positions = CODON_POS_MAP[args.codon_positions]

    # ── Parse target bipartitions ─────────────────────────────────────────────
    targets = []   # list of (bitmask, source_label) preserving order
    for tok in args.bipart:
        bm = parse_bipart_token(tok, n_tax)
        targets.append(bm)
    if args.bipart_file:
        file_bms = read_tree_biparts(args.bipart_file, n_tax)
        # read_tree_biparts returns a set; preserve file order by re-reading
        # (set is fine for membership but we want deterministic block order)
        for bm in sorted(file_bms,
                         key=lambda b: bitmask_to_binary_str(b, n_tax)):
            if bm not in targets:
                targets.append(bm)

    # Deduplicate while preserving order
    seen = set()
    unique_targets = []
    for bm in targets:
        if bm not in seen:
            seen.add(bm)
            unique_targets.append(bm)
    targets = unique_targets

    print(f"Target bipartitions: {len(targets)}", file=sys.stderr)

    # ── Process alignment ────────────────────────────────────────────────────
    sites = process_alignment(order, seqs, mode, codon_positions, binary,
                              track_positions=True)
    print(f"Parsimony-informative sites: {len(sites)}", file=sys.stderr)

    # ── Annotations ───────────────────────────────────────────────────────────
    annotations = None
    if args.annotation:
        if not os.path.isfile(args.annotation):
            sys.exit(f"ERROR: Annotation file not found: {args.annotation}")
        annotations = read_annotations(args.annotation, aln_len)
        print(f"Loaded annotations for {aln_len} alignment position(s) "
              f"from {args.annotation}", file=sys.stderr)

    # ── For each target bipartition, find full/partial support sites ─────────
    # results[bm] = (list of 1-based full-support positions,
    #                list of 1-based partial-support positions)
    results = {bm: ([], []) for bm in targets}
    target_set = set(targets)

    # union_sites: maps 1-based position -> site record (for FASTA output)
    union_sites = {}

    for site in sites:
        pos1 = site['site_index'] + 1   # 1-based
        for bm, weight in site['bipartitions']:
            if bm in target_set:
                if weight == 1.0:
                    results[bm][0].append(pos1)
                else:
                    results[bm][1].append(pos1)
                if pos1 not in union_sites:
                    union_sites[pos1] = site

    for bm in targets:
        results[bm][0].sort()
        results[bm][1].sort()

    # ── Write text output ────────────────────────────────────────────────────
    txt_path = args.outfile + '.txt'
    with open(txt_path, 'w') as f:
        f.write("Taxa (analysis order):\n")
        for i, t in enumerate(order, start=1):
            f.write(f"  {i}: {t}\n")
        f.write("\n")

        for bm in targets:
            bin_str = bitmask_to_binary_str(bm, n_tax)
            hex_str = bitmask_to_jakobsen_hex(bm, n_tax)
            full_sites, partial_sites = results[bm]

            f.write(f"=== Bipartition: {bin_str}  (hex: {hex_str}) ===\n")

            if args.verbose_outfile:
                group1 = [order[i] for i, c in enumerate(bin_str) if c == '1']
                group0 = [order[i] for i, c in enumerate(bin_str) if c == '0']
                f.write(f"({','.join(group1)})|({','.join(group0)})\n")

            full_str    = ', '.join(str(p) for p in full_sites)
            partial_str = ', '.join(str(p) for p in partial_sites)
            f.write(f"Full support (no missing): {full_str}\n")
            f.write(f"Partial support (1 missing): {partial_str}\n")

            if annotations is not None:
                all_positions = sorted(set(full_sites) | set(partial_sites))
                partial_set   = set(partial_sites)
                for pos in all_positions:
                    marker = '*' if pos in partial_set else ''
                    annot  = annotations.get(pos, 'unannotated')
                    f.write(f"{pos}{marker}\t{annot}\n")

            f.write("\n")

    print(f"Wrote {txt_path}", file=sys.stderr)

    # Summary to stdout
    total_full    = sum(len(results[bm][0]) for bm in targets)
    total_partial = sum(len(results[bm][1]) for bm in targets)
    print(f"\nBipartitions queried: {len(targets)}")
    print(f"Total full-support site instances:    {total_full}")
    print(f"Total partial-support site instances: {total_partial}")
    print(f"Union of sites (for FASTA): {len(union_sites)}")

    # ── FASTA output ──────────────────────────────────────────────────────────
    if args.fasta_output:
        fasta_path = args.outfile + '.fasta'
        positions  = sorted(union_sites.keys())

        seqs_ordered = {t: seqs[t] for t in order}

        with open(fasta_path, 'w') as f:
            for t_idx, t in enumerate(order):
                if args.fasta_output == 'original':
                    chars = [seqs_ordered[t][pos - 1] for pos in positions]
                else:  # binary
                    chars = []
                    for pos in positions:
                        coded_val = union_sites[pos]['coded'][t_idx]
                        chars.append('?' if coded_val is None else str(coded_val))
                f.write(f">{t}\n{''.join(chars)}\n")

        print(f"\nWrote {fasta_path} ({len(positions)} sites, "
              f"{args.fasta_output} representation)", file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
