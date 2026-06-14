#!/usr/bin/env python3
"""
paup_pstats.py - Parsimony statistics (CI, RI, RC) for binary character matrices.

Converts binary PHYLIP or Nexus alignments to Nexus format, runs PAUP* to find
the most parsimonious tree(s), and extracts CI, RI, and RC from the pscores output.

Two modes:
  Single-file mode:  -f / --file  (one alignment; optional --compat for compatibility score)
  Directory mode:    -i / --indir (sample_loci.py output directory; compat scores auto-loaded)

Usage:
  python3 paup_pstats.py -f alignment.phy -o results --paup /path/to/paup
  python3 paup_pstats.py -f alignment.phy -o results --compat --paup /path/to/paup
  python3 paup_pstats.py -i sample_loci_outdir/ -o pstats_results/ --paup /path/to/paup
"""

import argparse
import csv
import os
import random
import re
import statistics
import subprocess
import sys
from pathlib import Path

# ── Import helpers from partcompat.py ─────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
PARTCOMPAT = os.path.join(_HERE, 'partcompat.py')

try:
    from partcompat import (
        read_tax_order,
        detect_and_parse,
        is_binary_data,
        process_alignment,
        resolve_mode,
        CODON_POS_MAP,
    )
except ImportError as e:
    sys.exit(f"ERROR: Cannot import partcompat.py from {_HERE}: {e}")

# ---------------------------------------------------------------------------
# Input parsing and PI-site filtering
# ---------------------------------------------------------------------------

def read_raw_alignment(path):
    """
    Auto-detect format (FASTA, PHYLIP, Nexus) and return (taxa_list, seq_dict).
    Uses partcompat's detect_and_parse for consistency with the rest of the suite.
    """
    with open(path) as f:
        text = f.read()
    return detect_and_parse(text)


def process_to_binary_matrix(path, mode_str, codon_positions, tax_order=None):
    """
    Read an alignment file and return a PI-only binary matrix using
    partcompat's coding pipeline (RY, KM, WS, 2state, partimatrix, etc.).

    Returns (taxa_list, binary_seq_dict, n_pi_sites).

    Each character in the output sequences is '0', '1', or '?' (missing).
    Only parsimony-informative sites under the chosen mode are included.
    Single-missing sites contribute their canonical coded pattern (the
    missing position becomes '?').
    """
    taxa, seqs = read_raw_alignment(path)
    if not taxa:
        return [], {}, 0

    # Optional taxon reordering
    if tax_order:
        present = [t for t in tax_order if t in seqs]
        extra   = [t for t in taxa if t not in set(tax_order)]
        taxa    = present + extra

    binary = is_binary_data(seqs)
    mode   = resolve_mode(mode_str)
    sites  = process_alignment(taxa, seqs, mode, codon_positions, binary)

    # Build output matrix from coded arrays
    bin_seqs = {t: [] for t in taxa}
    for site in sites:
        coded = site['coded']
        for i, t in enumerate(taxa):
            v = coded[i]
            bin_seqs[t].append('?' if v is None else str(v))

    bin_seqs = {t: ''.join(chars) for t, chars in bin_seqs.items()}
    return taxa, bin_seqs, len(sites)

# ---------------------------------------------------------------------------
# Nexus output for PAUP*
# ---------------------------------------------------------------------------

def write_nexus(path, taxa, seqs):
    """Write a Nexus file with binary data block."""
    nchar = len(seqs[taxa[0]])
    with open(path, 'w') as f:
        f.write('#NEXUS\n\n')
        f.write('Begin data;\n')
        f.write(f'    Dimensions ntax={len(taxa)} nchar={nchar};\n')
        f.write('    Format datatype=standard symbols="01" missing=? gap=-;\n')
        f.write('    Matrix\n')
        for t in taxa:
            f.write(f'        {t}    {seqs[t]}\n')
        f.write('    ;\nEnd;\n')

# ---------------------------------------------------------------------------
# PAUP* block generation
# ---------------------------------------------------------------------------

def make_paup_block(scores_file, mode, nreps, rseed, n_tax):
    """
    Generate PAUP* command block.

    mode:   'bandb'   - branch and bound (exact, warns if n_tax > 12)
            'thorough' - hsearch increase=auto (default)
            'lazy'    - hsearch maxtrees=1000 increase=no
    """
    scores_file = scores_file.replace('\\', '/')   # PAUP* prefers forward slashes

    if mode == 'bandb':
        search_cmd = 'bandb;'
    elif mode == 'lazy':
        search_cmd = (f'hsearch addseq=random nreps={nreps} '
                      f'maxtrees=1000 increase=no nchuck=10 chuckscore=1 '
                      f'rseed={rseed};')
    else:  # thorough
        search_cmd = (f'hsearch addseq=random nreps={nreps} '
                      f'rseed={rseed};')

    return (
        f'\nBegin PAUP;\n'
        f'    set increase=auto;\n'
        f'    {search_cmd}\n'
        f'    pscores 1 / ci ri rc scorefile={scores_file} replace;\n'
        f'    cleartrees;\n'
        f'    quit;\n'
        f'End;\n'
    )


def append_paup_block(nexus_path, scores_path, mode, nreps, rseed, n_tax):
    """Append PAUP* block to an existing Nexus file."""
    block = make_paup_block(scores_path, mode, nreps, rseed, n_tax)
    with open(nexus_path, 'a') as f:
        f.write(block)

# ---------------------------------------------------------------------------
# PAUP* execution
# ---------------------------------------------------------------------------

def check_paup(paup_path):
    """Return True if PAUP* is reachable at paup_path."""
    try:
        result = subprocess.run([paup_path, '-h'],
                                capture_output=True, text=True, timeout=10)
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False
    except Exception:
        return False


def run_paup(paup_path, nexus_file, cwd=None, quiet=False):
    """Run PAUP* in batch mode. Returns (success, stdout, stderr)."""
    try:
        result = subprocess.run(
            [paup_path, '-n', nexus_file],
            capture_output=True, text=True, cwd=cwd
        )
        success = result.returncode == 0
        if not success and not quiet:
            print(f"  WARNING: PAUP* returned exit code {result.returncode}",
                  file=sys.stderr)
        return success, result.stdout, result.stderr
    except Exception as e:
        if not quiet:
            print(f"  ERROR running PAUP*: {e}", file=sys.stderr)
        return False, '', str(e)

# ---------------------------------------------------------------------------
# pscores file parsing
# ---------------------------------------------------------------------------

def parse_pscores_file(path):
    """
    Parse the pscores output file written by PAUP* scorefile= option.

    Returns dict with keys: n_sites, tree_length, CI, RI, RC
    or None if parsing fails.

    Expected format (whitespace-delimited, header lines start with dashes):
      Tree    Length    CI       RI       RC     ...
      ----  --------  ------  ------  ------
         1      47    0.7021  0.8120  0.5703  ...
    """
    try:
        with open(path) as f:
            lines = f.readlines()
    except FileNotFoundError:
        return None

    for line in lines:
        s = line.strip()
        if not s or s.startswith('-') or s.lower().startswith('tree'):
            continue
        parts = s.split()
        # Expect at least 5 fields: tree_num, length, CI, RI, RC
        if len(parts) >= 5:
            try:
                length = float(parts[1])
                ci     = float(parts[2])
                ri     = float(parts[3])
                rc     = float(parts[4])
                return {'tree_length': length, 'CI': ci, 'RI': ri, 'RC': rc}
            except ValueError:
                continue
    return None

# ---------------------------------------------------------------------------
# Compatibility score reading
# ---------------------------------------------------------------------------

def read_compat_score(prefix):
    """Read compatibility score from OUTFILE.compat.txt line 1."""
    try:
        with open(prefix + '.compat.txt') as f:
            return float(f.readline().strip().split()[0])
    except Exception:
        return None

# ---------------------------------------------------------------------------
# Core per-replicate analysis
# ---------------------------------------------------------------------------

def run_one(phy_path, workdir, prefix, paup_path, mode_str, codon_positions,
            nreps, seed, tax_order=None, quiet=False):
    """
    Process one alignment through the partcompat coding pipeline to extract
    PI sites, write a Nexus file, run PAUP*, and parse CI/RI/RC.

    Returns dict with keys: n_pi_sites, tree_length, CI, RI, RC
    or None on failure.
    """
    taxa, seqs, n_pi = process_to_binary_matrix(
        phy_path, mode_str, codon_positions, tax_order=tax_order)

    if not taxa:
        if not quiet:
            print(f"  WARNING: no sequences in {phy_path}", file=sys.stderr)
        return None
    if n_pi == 0:
        if not quiet:
            print(f"  WARNING: no PI sites in {phy_path} under mode '{mode_str}'",
                  file=sys.stderr)
        return None

    n_tax_actual = len(taxa)

    nex_path    = os.path.join(workdir, prefix + '.nex')
    scores_path = os.path.join(workdir, prefix + '.scores.txt')
    scores_base = prefix + '.scores.txt'

    write_nexus(nex_path, taxa, seqs)

    rseed = seed if seed is not None else random.randint(0, 999999999)
    append_paup_block(nex_path, scores_base, mode_str, nreps, rseed, n_tax_actual)

    success, stdout, stderr = run_paup(paup_path, os.path.basename(nex_path),
                                        cwd=workdir, quiet=quiet)
    if not success:
        return None

    stats = parse_pscores_file(scores_path)
    if stats is None:
        if not quiet:
            print(f"  WARNING: could not parse pscores for {prefix}",
                  file=sys.stderr)
        return None

    stats['n_pi_sites'] = n_pi
    return stats

# ---------------------------------------------------------------------------
# Summary statistics helpers
# ---------------------------------------------------------------------------

def _sd(vals):
    n = len(vals)
    if n <= 1:
        return 0.0
    mu = sum(vals) / n
    return (sum((v - mu)**2 for v in vals) / (n - 1))**0.5


def fmt_summary(label, vals):
    if not vals:
        return f"{label}: no data"
    return (f"{label}: mean={statistics.mean(vals):.4f}  "
            f"SD={_sd(vals):.4f}  "
            f"min={min(vals):.4f}  max={max(vals):.4f}  "
            f"(n={len(vals)})")

# ---------------------------------------------------------------------------
# Write summary TSV
# ---------------------------------------------------------------------------

def write_summary_tsv(path, rows, compat_scores):
    """
    Write per-replicate summary TSV.
    rows: list of dicts with keys replicate, n_sites, tree_length, CI, RI, RC
    compat_scores: dict mapping replicate label -> score (may be empty)
    """
    have_compat = bool(compat_scores)
    fieldnames = ['replicate', 'n_sites', 'tree_length', 'CI', 'RI', 'RC']
    if have_compat:
        fieldnames.append('compat_score')

    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            out = {k: row.get(k, '') for k in fieldnames}
            if have_compat:
                out['compat_score'] = compat_scores.get(row['replicate'], '')
            writer.writerow(out)

        # Summary footer rows
        def stat_row(label, vals):
            if not vals:
                return None
            return {'replicate': label,
                    'n_sites':   '',
                    'tree_length': f"{statistics.mean([r['tree_length'] for r in rows if r.get('tree_length') is not None]):.4f}" if label == 'mean' else '',
                    'CI': f"{statistics.mean(vals):.4f}",
                    'RI': '',
                    'RC': ''}

        # Write blank separator then summary stats
        for metric in ['CI', 'RI', 'RC', 'tree_length']:
            vals = [r[metric] for r in rows if r.get(metric) is not None]
            if not vals:
                continue
            for lbl, fn in [('mean', statistics.mean), ('sd', _sd),
                             ('min', min), ('max', max)]:
                out = {k: '' for k in fieldnames}
                out['replicate'] = f'{lbl}_{metric}'
                out[metric] = f'{fn(vals):.4f}'
                f.write('\t'.join(str(out.get(k, '')) for k in fieldnames) + '\n')
        if have_compat:
            cvals = [v for v in compat_scores.values() if v is not None]
            for lbl, fn in [('mean', statistics.mean), ('sd', _sd),
                             ('min', min), ('max', max)]:
                if not cvals:
                    continue
                out = {k: '' for k in fieldnames}
                out['replicate'] = f'{lbl}_compat_score'
                out['compat_score'] = f'{fn(cvals):.6f}'
                f.write('\t'.join(str(out.get(k, '')) for k in fieldnames) + '\n')

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Parsimony statistics (CI, RI, RC) for binary matrices.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Input modes (mutually exclusive):
  -f / --file      Single alignment file (PHYLIP or Nexus, auto-detected)
  -i / --indir     sample_loci.py output directory (processes all replicate_*.phy)

Search modes (mutually exclusive):
  --bandb          Branch and bound (exact; recommended for ≤12 taxa)
  --lazy           Fast heuristic (maxtrees=1000, nchuck=10)
  (default)        Thorough heuristic (hsearch increase=auto)

Examples:
  python3 paup_pstats.py -f locus.phy -o locus_stats --paup /usr/local/bin/paup
  python3 paup_pstats.py -f locus.phy -o locus_stats --compat --paup paup
  python3 paup_pstats.py -i genome_results/ -o pstats/ --paup paup --lazy
  python3 paup_pstats.py -i genome_results/ -o pstats/ --bandb --paup paup
        ''',
    )

    mode_grp = p.add_mutually_exclusive_group(required=True)
    mode_grp.add_argument('-f', '--file',
                          help='Single input alignment (PHYLIP or Nexus)')
    mode_grp.add_argument('-i', '--indir',
                          help='sample_loci.py output directory')

    p.add_argument('-o', '--outdir', required=True,
                   help='Output directory for working files and summary')
    p.add_argument('--paup', default='paup',
                   help='Path to PAUP* executable (default: paup)')

    search_grp = p.add_mutually_exclusive_group()
    search_grp.add_argument('--bandb', action='store_true',
                             help='Branch-and-bound search (exact; warns if >12 taxa)')
    search_grp.add_argument('--lazy', action='store_true',
                             help='Fast heuristic (maxtrees=1000, increase=no)')

    p.add_argument('--nreps', type=int, default=100,
                   help='Number of hsearch random-addition replicates (default: 100)')
    p.add_argument('--seed', type=int, default=None,
                   help='Base random seed (auto-generated if omitted)')
    p.add_argument('-m', '--mode', default='RY',
                   help='Site coding mode (default: RY). '
                        'Options: RY, KM, WS, perfectRY, perfectKM, '
                        'perfectWS, transitions, transitionsAG, '
                        'transitionsCT, 2state, partimatrix (and reversed '
                        'pairs e.g. YR, perfectYR, transitionGA). '
                        'Perfect variants require exactly 2 unambiguous '
                        'nucleotides per column, one per axis side, no gaps '
                        'or IUPAC (Tiley et al. 2020). Transition variants '
                        'exclude cross-axis characters. Only PI sites under '
                        'the chosen scheme are passed to PAUP*.')
    p.add_argument('--codon-positions', default='all',
                   choices=['all', '1', '2', '3', '12', '13', '23', '123'],
                   help='Restrict to specific codon positions (default: all)')
    p.add_argument('--tax-order',
                   help='Taxon order file (one name per line). Used to reorder '
                        'taxa before processing; ensures consistent ordering '
                        'with partcompat.py outputs.')
    p.add_argument('--compat', action='store_true',
                   help='(Single-file mode only) also run partcompat.py and '
                        'include the compatibility score in output')
    p.add_argument('--redo', action='store_true',
                   help='Overwrite existing output directory contents')
    return p.parse_args()

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # ── Validate PAUP* ───────────────────────────────────────────────────────
    if not check_paup(args.paup):
        sys.exit(f"ERROR: PAUP* not found at '{args.paup}'. "
                 f"Use --paup to specify the path.")

    # ── Search mode ──────────────────────────────────────────────────────────
    if args.bandb:
        search_mode = 'bandb'
    elif args.lazy:
        search_mode = 'lazy'
    else:
        search_mode = 'thorough'

    # ── Base seed ────────────────────────────────────────────────────────────
    base_seed = args.seed if args.seed is not None else random.randrange(2**32)
    rng = random.Random(base_seed)
    print(f"Random seed: {base_seed}", file=sys.stderr)

    # ── Output directory ─────────────────────────────────────────────────────
    os.makedirs(args.outdir, exist_ok=True)
    existing = list(Path(args.outdir).glob('*.scores.txt'))
    if existing and not args.redo:
        sys.exit(
            f"ERROR: '{args.outdir}' already contains {len(existing)} scores "
            f"file(s). Use --redo to overwrite."
        )

    # ── Single-file mode ─────────────────────────────────────────────────────
    # Resolve mode and codon positions (used in both single-file and dir mode)
    mode_str         = args.mode
    codon_positions  = CODON_POS_MAP[args.codon_positions]
    tax_order        = read_tax_order(args.tax_order) if args.tax_order else None

    if args.file:
        if not os.path.isfile(args.file):
            sys.exit(f"ERROR: Input file not found: {args.file}")
        print(f"Single-file mode: {args.file}", file=sys.stderr)
        print(f"  Mode: {mode_str}  Codon positions: {args.codon_positions}",
              file=sys.stderr)

        # Probe taxon count for bandb warning (before full processing)
        _taxa_probe, _, _ = process_to_binary_matrix(
            args.file, mode_str, codon_positions, tax_order=tax_order)
        n_tax = len(_taxa_probe)
        if args.bandb and n_tax > 12:
            print(f"  WARNING: branch-and-bound with {n_tax} taxa may be very slow.",
                  file=sys.stderr)

        prefix = Path(args.file).stem
        rseed  = rng.randint(0, 999999999)
        stats  = run_one(args.file, args.outdir, prefix,
                         args.paup, mode_str, codon_positions,
                         args.nreps, rseed, tax_order=tax_order)

        if stats is None:
            sys.exit("ERROR: PAUP* analysis failed or no PI sites found. "
                     "Check PAUP* output and --mode setting.")

        # Optional compatibility score (run partcompat.py on same file)
        compat = None
        if args.compat:
            if not os.path.isfile(PARTCOMPAT):
                print("WARNING: partcompat.py not found; skipping --compat.",
                      file=sys.stderr)
            else:
                cp_prefix = os.path.join(args.outdir, prefix + '_compat')
                cmd = [sys.executable, PARTCOMPAT,
                       '-i', args.file, '-o', cp_prefix,
                       '--no-svg', '--score-file',
                       '-m', mode_str]
                if args.tax_order:
                    cmd += ['--tax-order', args.tax_order]
                if args.codon_positions != 'all':
                    cmd += ['--codon-positions', args.codon_positions]
                subprocess.run(cmd, capture_output=True)
                compat = read_compat_score(cp_prefix)

        # Print results
        print(f"\nPI sites passed to PAUP*: {stats['n_pi_sites']}")
        print(f"Tree length: {stats['tree_length']:.0f}")
        print(f"CI:          {stats['CI']:.4f}")
        print(f"RI:          {stats['RI']:.4f}")
        print(f"RC:          {stats['RC']:.4f}")
        if compat is not None:
            print(f"Compat:      {compat:.6f}")

        # Write single-file TSV
        tsv_path = os.path.join(args.outdir, prefix + '.pstats.tsv')
        rows = [{'replicate': prefix, 'n_sites': stats['n_pi_sites'],
                 'tree_length': stats['tree_length'],
                 'CI': stats['CI'], 'RI': stats['RI'], 'RC': stats['RC']}]
        compat_scores = {prefix: compat} if compat is not None else {}
        write_summary_tsv(tsv_path, rows, compat_scores)
        print(f"\nWrote {tsv_path}", file=sys.stderr)
        return

    # ── Directory mode ───────────────────────────────────────────────────────
    if not os.path.isdir(args.indir):
        sys.exit(f"ERROR: Input directory not found: {args.indir}")
    phy_files = sorted(Path(args.indir).glob('replicate_*.phy'))
    if not phy_files:
        sys.exit(f"ERROR: No replicate_*.phy files found in '{args.indir}'.")

    n_reps = len(phy_files)
    print(f"Directory mode: {n_reps} replicate(s) in {args.indir}",
          file=sys.stderr)
    print(f"  Mode: {mode_str}  Codon positions: {args.codon_positions}",
          file=sys.stderr)

    # Check n_tax from first file for bandb warning
    try:
        taxa0, _, _ = process_to_binary_matrix(
            str(phy_files[0]), mode_str, codon_positions, tax_order=tax_order)
        if args.bandb and len(taxa0) > 12:
            print(f"  WARNING: branch-and-bound with {len(taxa0)} taxa may be "
                  f"very slow.", file=sys.stderr)
    except Exception:
        pass

    rows          = []
    compat_scores = {}
    n_failed      = 0

    for idx, phy_path in enumerate(phy_files, 1):
        prefix = phy_path.stem          # e.g. replicate_001
        rseed  = rng.randint(0, 999999999)

        stats = run_one(str(phy_path), args.outdir, prefix,
                        args.paup, mode_str, codon_positions,
                        args.nreps, rseed, tax_order=tax_order, quiet=False)

        if stats is None:
            n_failed += 1
            print(f"  WARNING: {prefix} failed", file=sys.stderr)
        else:
            rows.append({'replicate': prefix,
                         'n_sites':     stats['n_pi_sites'],
                         'tree_length': stats['tree_length'],
                         'CI': stats['CI'],
                         'RI': stats['RI'],
                         'RC': stats['RC']})

        # Load compat score from sample_loci output
        compat_prefix = os.path.join(args.indir, prefix)
        sc = read_compat_score(compat_prefix)
        if sc is not None:
            compat_scores[prefix] = sc

        if idx % max(1, n_reps // 10) == 0 or idx == n_reps:
            print(f"  {idx}/{n_reps} replicates processed", file=sys.stderr)

    n_ok = len(rows)
    print(f"\nCompleted: {n_ok} succeeded, {n_failed} failed", file=sys.stderr)

    if n_ok == 0:
        sys.exit("ERROR: All replicates failed.")

    # Summary to stdout
    for metric in ['CI', 'RI', 'RC']:
        vals = [r[metric] for r in rows]
        print(fmt_summary(metric, vals))
    tl_vals = [r['tree_length'] for r in rows]
    print(fmt_summary('Tree length', tl_vals))
    if compat_scores:
        cvals = [v for v in compat_scores.values() if v is not None]
        print(fmt_summary('Compat', cvals))

    # Write summary TSV
    tsv_path = os.path.join(args.outdir, 'pstats_summary.tsv')
    write_summary_tsv(tsv_path, rows, compat_scores)
    print(f"\nWrote {tsv_path}", file=sys.stderr)

    # Write run info
    info_path = os.path.join(args.outdir, 'pstats_run_info.txt')
    with open(info_path, 'w') as f:
        f.write("=== paup_pstats.py run information ===\n\n")
        f.write(f"Command:       {' '.join(sys.argv)}\n")
        f.write(f"Seed:          {base_seed}\n")
        f.write(f"Search mode:   {search_mode}\n")
        if search_mode != 'bandb':
            f.write(f"hsearch nreps: {args.nreps}\n")
        f.write(f"PAUP* path:    {args.paup}\n")
        f.write(f"Replicates:    {n_ok} succeeded, {n_failed} failed\n\n")
        for metric in ['CI', 'RI', 'RC', 'tree_length']:
            vals = [r[metric] for r in rows]
            if vals:
                f.write(f"{metric}: mean={statistics.mean(vals):.4f}  "
                        f"sd={_sd(vals):.4f}  "
                        f"min={min(vals):.4f}  max={max(vals):.4f}\n")
        if compat_scores:
            cvals = [v for v in compat_scores.values() if v is not None]
            if cvals:
                f.write(f"Compat: mean={statistics.mean(cvals):.6f}  "
                        f"sd={_sd(cvals):.6f}  "
                        f"min={min(cvals):.6f}  max={max(cvals):.6f}\n")
    print(f"Wrote {info_path}", file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
