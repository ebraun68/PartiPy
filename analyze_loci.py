#!/usr/bin/env python3
"""
analyze_loci.py - Run partcompat.py on every alignment in a directory.

Produces one TSV (and optionally a compatibility score file) per locus in a
flat output directory.  No SVG files are produced.  This is the companion to
sample_loci.py: run both programs with the same --tax-order, --mode, and
--plot-biparts arguments, then pass their output directories to
lento_distances.py for comparison.

Usage:
  python3 analyze_loci.py -i INDIR -o OUTDIR --tax-order FILE [options]
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

# ── Import helpers from partcompat.py ─────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
PARTCOMPAT = os.path.join(_HERE, 'partcompat.py')

try:
    from partcompat import (
        detect_and_parse,
        is_binary_data,
        process_alignment,
        read_tax_order,
        resolve_mode,
        CODON_POS_MAP,
    )
except ImportError as e:
    sys.exit(f"ERROR: Cannot import partcompat.py from {_HERE}: {e}")

# ── Alignment file extensions ──────────────────────────────────────────────────
ALIGN_EXTENSIONS = {'.fasta', '.fa', '.fas', '.fna', '.phy', '.phylip',
                    '.nex', '.nexus', '.nxs', '.txt'}


def find_alignment_files(indir):
    files = []
    for name in sorted(os.listdir(indir)):
        path = os.path.join(indir, name)
        if os.path.isfile(path):
            ext = os.path.splitext(name)[1].lower()
            if ext in ALIGN_EXTENSIONS or ext == '':
                files.append(path)
    return files


def count_pi_sites(path, tax_order, mode_str, codon_positions, max_missing_frac):
    """
    Quick check: return (n_pi, n_missing_taxa, skip_reason).
    skip_reason is None if the locus is eligible.
    """
    try:
        with open(path) as f:
            text = f.read()
        file_order, file_seqs = detect_and_parse(text)
    except Exception as exc:
        return 0, 0, f"parse error: {exc}"

    if not file_order:
        return 0, 0, "no sequences found"

    n_missing = sum(1 for t in tax_order if t not in file_seqs)
    missing_frac = n_missing / len(tax_order)
    if missing_frac > max_missing_frac:
        return 0, n_missing, (
            f"too many missing taxa "
            f"({n_missing}/{len(tax_order)} = {missing_frac:.1%} > "
            f"{max_missing_frac:.1%} threshold)"
        )

    aln_len = len(file_seqs[file_order[0]])
    full_seqs = {t: file_seqs.get(t, '?' * aln_len) for t in tax_order}
    binary = is_binary_data(full_seqs)
    sites  = process_alignment(tax_order, full_seqs,
                               resolve_mode(mode_str), codon_positions, binary)
    if not sites:
        return 0, n_missing, "no parsimony-informative sites after coding"

    return len(sites), n_missing, None


def run_partcompat(alignment_path, out_prefix, tax_order_file,
                   mode_str, codon_positions_str,
                   plot_biparts, tree_biparts, bipart_order, no_hex):
    """Run partcompat.py as a subprocess. Returns (success, stderr_text)."""
    cmd = [
        sys.executable, PARTCOMPAT,
        '-i',  alignment_path,
        '-o',  out_prefix,
        '-m',  mode_str,
        '--tax-order', tax_order_file,
        '--no-svg',
        '--score-file',
    ]
    if codon_positions_str and codon_positions_str != 'all':
        cmd += ['--codon-positions', codon_positions_str]
    if plot_biparts:
        cmd += ['-p', plot_biparts]
    if tree_biparts:
        cmd += ['-t', tree_biparts]
    if bipart_order:
        cmd += ['-b', bipart_order]
    if no_hex:
        cmd += ['--no-hex']

    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0, result.stderr


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            'Run partcompat.py on every alignment in a directory, '
            'producing one TSV per locus in a flat output directory. '
            'Run with the same --tax-order, --mode, and --plot-biparts '
            'as sample_loci.py for comparable results.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Output directory contents:
  LOCUS.tsv          bipartition table from partcompat.py
  LOCUS.compat.txt   compatibility score (two lines: score, n_pi_sites)
  run_info.txt       parameters, eligible and skipped loci

Example:
  python3 analyze_loci.py \\
      -i loci_directory/ \\
      -o gene_results/ \\
      --tax-order taxa.txt \\
      -p ref.biparts \\
      --mode RY
        ''',
    )
    p.add_argument('-i', '--indir',      required=True,
                   help='Input directory containing alignment files')
    p.add_argument('-o', '--outdir',     required=True,
                   help='Output directory for TSV files (created if absent)')
    p.add_argument('--tax-order',        required=True,
                   help='Taxon order file (one name per line); required')
    p.add_argument('-m', '--mode',       default='RY',
                   help='Site coding mode (default: RY). '
                        'Options: RY, KM, WS, perfectRY, perfectKM, '
                        'perfectWS, transitions, transitionsAG, '
                        'transitionsCT, 2state, partimatrix (and reversed '
                        'pairs e.g. YR, perfectYR, transitionGA). '
                        'Perfect variants require exactly 2 unambiguous '
                        'nucleotides per column, one per axis side, no gaps '
                        'or IUPAC (Tiley et al. 2020). Transition variants '
                        'exclude sites with any cross-axis character.')
    p.add_argument('--codon-positions',  default='all',
                   choices=['all', '1', '2', '3', '12', '13', '23', '123'],
                   help='Codon positions to include (default: all)')
    p.add_argument('-p', '--plot-biparts',
                   help='Display filter: only these bipartitions appear in TSV. '
                        'Must match the file used with sample_loci.py.')
    p.add_argument('-t', '--tree-biparts',
                   help='Reference tree bipartitions (for TSV grouping column)')
    p.add_argument('-b', '--bipart-order',
                   help='Bipartition column ordering file')
    p.add_argument('--max-missing-taxa', type=float, default=0.5,
                   metavar='FRAC',
                   help='Skip a locus if this fraction of taxa from '
                        '--tax-order are absent (default: 0.5)')
    p.add_argument('--no-hex',           action='store_true',
                   help='Suppress Jakobsen hex labels (passed to partcompat.py)')
    p.add_argument('--redo',             action='store_true',
                   help='Overwrite existing output directory contents')
    return p.parse_args()


def main():
    args = parse_args()

    if not os.path.isfile(PARTCOMPAT):
        sys.exit(f"ERROR: partcompat.py not found at {PARTCOMPAT}")

    # ── Taxon order ───────────────────────────────────────────────────────────
    tax_order = read_tax_order(args.tax_order)
    if not tax_order:
        sys.exit(f"ERROR: No taxa found in {args.tax_order}")
    print(f"Taxon order: {len(tax_order)} taxa", file=sys.stderr)

    # ── Output directory ──────────────────────────────────────────────────────
    os.makedirs(args.outdir, exist_ok=True)
    existing = list(Path(args.outdir).glob('*.tsv'))
    if existing and not args.redo:
        sys.exit(
            f"ERROR: '{args.outdir}' already contains {len(existing)} TSV "
            f"file(s). Use --redo to overwrite."
        )

    # ── Mode and codon positions ──────────────────────────────────────────────
    mode_str         = resolve_mode(args.mode)
    codon_positions  = CODON_POS_MAP[args.codon_positions]

    # ── Find alignment files ──────────────────────────────────────────────────
    if not os.path.isdir(args.indir):
        sys.exit(f"ERROR: Input directory not found: {args.indir}")
    aln_files = find_alignment_files(args.indir)
    if not aln_files:
        sys.exit(f"ERROR: No alignment files found in '{args.indir}'")
    print(f"Found {len(aln_files)} alignment file(s) in {args.indir}",
          file=sys.stderr)

    # ── Process loci ─────────────────────────────────────────────────────────
    n_ok      = 0
    n_failed  = 0
    n_skipped = 0
    skip_log  = []
    fail_log  = []

    for path in aln_files:
        stem = Path(path).stem

        # Pre-screen: check eligibility and PI site count
        n_pi, n_miss, skip_reason = count_pi_sites(
            path, tax_order, mode_str, codon_positions, args.max_missing_taxa)

        if skip_reason:
            n_skipped += 1
            skip_log.append((stem, skip_reason))
            print(f"  SKIP  {stem}: {skip_reason}", file=sys.stderr)
            continue

        out_prefix = os.path.join(args.outdir, stem)
        success, stderr_text = run_partcompat(
            path, out_prefix,
            args.tax_order, mode_str, args.codon_positions,
            args.plot_biparts, args.tree_biparts, args.bipart_order,
            args.no_hex,
        )

        if success:
            n_ok += 1
        else:
            n_failed += 1
            fail_log.append((stem, stderr_text.strip()))
            print(f"  FAIL  {stem}", file=sys.stderr)

        if (n_ok + n_failed + n_skipped) % max(1, len(aln_files) // 10) == 0 \
                or (n_ok + n_failed + n_skipped) == len(aln_files):
            done = n_ok + n_failed + n_skipped
            print(f"  {done}/{len(aln_files)} processed "
                  f"({n_ok} ok, {n_skipped} skipped, {n_failed} failed)",
                  file=sys.stderr)

    # ── run_info.txt ──────────────────────────────────────────────────────────
    info_path = os.path.join(args.outdir, 'run_info.txt')
    with open(info_path, 'w') as f:
        f.write("=== analyze_loci.py run information ===\n\n")
        f.write(f"Command:          {' '.join(sys.argv)}\n")
        f.write(f"Input directory:  {args.indir}\n")
        f.write(f"Output directory: {args.outdir}\n")
        f.write(f"Tax order file:   {args.tax_order}\n")
        f.write(f"Mode:             {mode_str}\n")
        f.write(f"Codon positions:  {args.codon_positions}\n")
        f.write(f"Max missing taxa: {args.max_missing_taxa:.2f}\n")
        if args.plot_biparts:
            f.write(f"Plot biparts:     {args.plot_biparts}\n")
        if args.tree_biparts:
            f.write(f"Tree biparts:     {args.tree_biparts}\n")
        if args.bipart_order:
            f.write(f"Bipart order:     {args.bipart_order}\n")
        f.write(f"\nLoci found:       {len(aln_files)}\n")
        f.write(f"Loci succeeded:   {n_ok}\n")
        f.write(f"Loci skipped:     {n_skipped}\n")
        f.write(f"Loci failed:      {n_failed}\n")
        if skip_log:
            f.write("\nSkipped loci:\n")
            for name, reason in skip_log:
                f.write(f"  {name}: {reason}\n")
        if fail_log:
            f.write("\nFailed loci:\n")
            for name, err in fail_log:
                f.write(f"  {name}: {err[:200]}\n")

    # ── Final report ──────────────────────────────────────────────────────────
    print(f"\nCompleted: {n_ok} succeeded, {n_skipped} skipped, "
          f"{n_failed} failed  (of {len(aln_files)} total)")
    print(f"Output: {args.outdir}/", file=sys.stderr)
    print(f"Wrote {info_path}", file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
