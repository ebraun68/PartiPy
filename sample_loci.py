#!/usr/bin/env python3
"""
sample_loci.py - Phylogenomic compatibility analysis via per-locus site sampling.

For each replicate:
  1. Randomly sample one parsimony-informative site from each eligible locus.
  2. Write a binary relaxed PHYLIP file (one column per locus).
  3. Run partcompat.py on the replicate to produce a TSV and compatibility score.

Then summarise across replicates:
  - Mean ± SD compatibility score
  - Mean support/conflict per bipartition (Lento-normalised)
  - Summary TSV and SVG

Usage:
  python3 sample_loci.py -i INDIR -o OUTDIR --tax-order ORDER.txt [options]
"""

import argparse
import csv
import math
import os
import random
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# ── Import shared functions from partcompat.py ─────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
PARTCOMPAT = os.path.join(_HERE, 'partcompat.py')

try:
    from partcompat import (
        detect_and_parse,
        is_binary_data,
        process_alignment,
        read_tax_order,
        read_tree_biparts,
        read_bipart_order,
        normalize_lento,
        build_column_order,
        write_tsv,
        write_svg,
        bitmask_to_binary_str,
        bitmask_to_jakobsen_hex,
        _write_tree_bipart_tsv,
        CODON_POS_MAP,
        resolve_mode,
    )
except ImportError as e:
    sys.exit(f"ERROR: Cannot import partcompat.py from {_HERE}: {e}")

# ── Locus scanning ─────────────────────────────────────────────────────────────

ALIGN_EXTENSIONS = {'.fasta', '.fa', '.fas', '.fna', '.phy', '.phylip',
                    '.nex', '.nexus', '.nxs', '.txt'}


def find_alignment_files(indir):
    """Return sorted list of alignment file paths in indir."""
    if not os.path.isdir(indir):
        sys.exit(f"ERROR: Input directory not found: {indir}")
    files = []
    for name in sorted(os.listdir(indir)):
        path = os.path.join(indir, name)
        if os.path.isfile(path):
            ext = os.path.splitext(name)[1].lower()
            if ext in ALIGN_EXTENSIONS or ext == '':
                files.append(path)
    return files


def load_locus(filepath, tax_order, mode, codon_positions, max_missing_frac):
    """
    Parse one alignment file and return its valid PI site records.

    Returns (sites, skip_reason):
      sites       - list of site dicts (from process_alignment) if eligible
      skip_reason - None if ok, else a string explaining why the locus was skipped
    """
    try:
        with open(filepath) as f:
            text = f.read()
        file_order, file_seqs = detect_and_parse(text)
    except Exception as exc:
        return None, f"parse error: {exc}"

    if not file_order:
        return None, "no sequences found"

    # Sequence length from the file
    aln_len = len(file_seqs[file_order[0]])

    # Count taxa absent from this locus (present in tax_order but not in file)
    n_present = sum(1 for t in tax_order if t in file_seqs)
    n_missing = len(tax_order) - n_present
    missing_frac = n_missing / len(tax_order)

    if missing_frac > max_missing_frac:
        return None, (
            f"too many missing taxa "
            f"({n_missing}/{len(tax_order)} = {missing_frac:.1%} > "
            f"{max_missing_frac:.1%} threshold)"
        )

    # Pad absent taxa with all-? sequences
    full_seqs = {}
    for t in tax_order:
        full_seqs[t] = file_seqs.get(t, '?' * aln_len)

    binary = is_binary_data(full_seqs)
    sites  = process_alignment(tax_order, full_seqs, mode, codon_positions, binary)

    if not sites:
        return None, "no parsimony-informative sites after coding"

    return sites, None


def scan_loci(indir, tax_order, mode, codon_positions, max_missing_frac):
    """
    Scan all alignment files in indir.

    Returns:
      eligible  - list of (filepath, sites) for usable loci
      skipped   - list of (filepath, reason) for skipped loci
    """
    files   = find_alignment_files(indir)
    eligible = []
    skipped  = []

    for path in files:
        sites, reason = load_locus(path, tax_order, mode,
                                   codon_positions, max_missing_frac)
        if sites is not None:
            eligible.append((path, sites))
        else:
            skipped.append((path, reason))

    return eligible, skipped

# ── Replicate sampling ─────────────────────────────────────────────────────────

def parse_dist(dist_str):
    """
    Parse a distribution specification string.
    Accepted formats:
      uniform:MIN:MAX  - draw uniformly from integers in [MIN, MAX]
      fixed:N          - always return N
    Returns a callable(rng) -> int.
    """
    parts = dist_str.split(':')
    kind  = parts[0].lower()
    if kind == 'fixed':
        if len(parts) != 2:
            sys.exit(f"ERROR: --variable-sites fixed:N requires one value, "
                     f"got '{dist_str}'")
        try:
            n = int(parts[1])
            if n < 1:
                raise ValueError
        except ValueError:
            sys.exit(f"ERROR: --variable-sites fixed:N requires a positive "
                     f"integer, got '{parts[1]}'")
        return lambda rng, _n=n: _n
    elif kind == 'uniform':
        if len(parts) != 3:
            sys.exit(f"ERROR: --variable-sites uniform:MIN:MAX requires two "
                     f"values, got '{dist_str}'")
        try:
            lo, hi = int(parts[1]), int(parts[2])
            if lo < 1 or hi < lo:
                raise ValueError
        except ValueError:
            sys.exit(f"ERROR: --variable-sites uniform:MIN:MAX requires MIN>=1 "
                     f"and MAX>=MIN, got '{dist_str}'")
        return lambda rng, _lo=lo, _hi=hi: rng.randint(_lo, _hi)
    else:
        sys.exit(f"ERROR: Unknown distribution '{parts[0]}'. "
                 f"Use uniform:MIN:MAX or fixed:N")


def sample_replicate(eligible_loci, boot_locus, rng,
                     dist_fn=None, downsample=False):
    """
    Sample PI sites from eligible loci, returning a list of coded vectors.

    Default (dist_fn=None, downsample=False):
      Each eligible locus contributes exactly one site (or multiple copies
      with boot_locus=True via multinomial resampling).

    With dist_fn supplied (--variable-sites):
      N is drawn from dist_fn(rng).

      downsample=False:
        N sites drawn i.i.d. from the weighted locus pool:
        select a locus (uniformly, or by bootstrap weight if boot_locus=True),
        draw one PI site from it, repeat N times.
        N may exceed n_loci; loci may contribute multiple sites.

      downsample=True (--downsample):
        N loci selected WITHOUT replacement (N must be <= n_eligible).
        One site drawn from each selected locus.
        Incompatible with boot_locus (enforced in main()).

    Returns list of coded vectors (list of 0/1/None, length = n_taxa).
    """
    n_loci = len(eligible_loci)

    if dist_fn is None:
        # Original behaviour: one site per locus (or per bootstrap instance)
        if boot_locus:
            weights = defaultdict(int)
            for idx in rng.choices(range(n_loci), k=n_loci):
                weights[idx] += 1
        else:
            weights = {i: 1 for i in range(n_loci)}
        coded_vectors = []
        for locus_idx, weight in sorted(weights.items()):
            _, sites = eligible_loci[locus_idx]
            site = rng.choice(sites)
            for _ in range(weight):
                coded_vectors.append(site['coded'])
        return coded_vectors

    # Variable-sites mode
    n_sites = dist_fn(rng)

    if downsample:
        # Without-replacement locus selection; one site per selected locus
        selected = rng.sample(range(n_loci), k=n_sites)
        return [rng.choice(eligible_loci[idx][1])['coded'] for idx in selected]
    else:
        # i.i.d. draw from weighted locus pool
        if boot_locus:
            boot_counts = defaultdict(int)
            for idx in rng.choices(range(n_loci), k=n_loci):
                boot_counts[idx] += 1
            locus_indices = list(boot_counts.keys())
            locus_weights = [boot_counts[i] for i in locus_indices]
        else:
            locus_indices = list(range(n_loci))
            locus_weights = [1] * n_loci

        coded_vectors = []
        for _ in range(n_sites):
            idx = rng.choices(locus_indices, weights=locus_weights, k=1)[0]
            _, sites = eligible_loci[idx]
            coded_vectors.append(rng.choice(sites)['coded'])
        return coded_vectors

# ── PHYLIP output ──────────────────────────────────────────────────────────────

def write_binary_phylip(coded_vectors, tax_order, outpath):
    """
    Write a relaxed PHYLIP file with binary data.
    Rows = taxa (in tax_order), columns = loci (one per coded vector).
    """
    n_tax   = len(tax_order)
    n_sites = len(coded_vectors)
    with open(outpath, 'w') as f:
        f.write(f"{n_tax} {n_sites}\n")
        for t_idx, taxon in enumerate(tax_order):
            chars = []
            for vec in coded_vectors:
                v = vec[t_idx]
                chars.append('?' if v is None else str(v))
            f.write(f"{taxon} {''.join(chars)}\n")

# ── partcompat.py subprocess ───────────────────────────────────────────────────

def run_partcompat(phy_path, prefix, tax_order_file, mode,
                   tree_biparts, plot_biparts, bipart_order,
                   codon_positions, no_hex, truncate=None):
    """
    Run partcompat.py as a subprocess (no SVG, score file always written).

    Returns (success: bool, stderr: str).
    """
    cmd = [
        sys.executable, PARTCOMPAT,
        '-i',  phy_path,
        '-o',  prefix,
        '-m',  mode,
        '--tax-order', tax_order_file,
        '--no-svg',
        '--score-file',
    ]
    if tree_biparts:
        cmd += ['-t', tree_biparts]
    if plot_biparts:
        cmd += ['-p', plot_biparts]
    if bipart_order:
        cmd += ['-b', bipart_order]
    if codon_positions and codon_positions != 'all':
        cmd += ['--codon-positions', codon_positions]
    if no_hex:
        cmd += ['--no-hex']
    if truncate is not None:
        cmd += ['--truncate', str(truncate)]

    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0, result.stderr

# ── Result readers ─────────────────────────────────────────────────────────────

def read_compat_score(prefix):
    """Read OUTFILE.compat.txt; return float or None on failure."""
    path = prefix + '.compat.txt'
    try:
        with open(path) as f:
            return float(f.read().strip().split()[0])
    except Exception:
        return None


def read_compat_n_biparts(prefix):
    """Read n_supported from line 2 of OUTFILE.compat.txt; return int or None."""
    path = prefix + '.compat.txt'
    try:
        with open(path) as f:
            lines = f.read().splitlines()
        if len(lines) >= 2:
            return int(lines[1].strip())
    except Exception:
        pass
    return None


def read_tsv_biparts(prefix):
    """
    Read OUTFILE.tsv; return dict mapping bitmask_int → (support_raw, conflict_raw).
    """
    path = prefix + '.tsv'
    result = {}
    try:
        with open(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                bin_str = row['bipartition_binary'].lstrip("'")
                bm      = int(bin_str, 2)
                result[bm] = (float(row['support_raw']),
                              float(row['conflict_raw']))
    except Exception:
        pass
    return result

# ── Descriptive statistics helpers ────────────────────────────────────────────

def _median(vals):
    """Median of a list of floats."""
    s = sorted(vals)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else (s[mid - 1] + s[mid]) / 2.0


def _sd(vals):
    """Population-corrected (n-1) standard deviation; 0.0 for n<=1."""
    n = len(vals)
    if n <= 1:
        return 0.0
    mu  = sum(vals) / n
    return (sum((v - mu) ** 2 for v in vals) / (n - 1)) ** 0.5


def _scaled_mad(vals):
    """
    Scaled median absolute deviation (MAD × 1.4826).
    Consistent estimator of the standard deviation under normality.
    Returns 0.0 for n<=1.
    """
    if len(vals) <= 1:
        return 0.0
    med  = _median(vals)
    mads = sorted(abs(v - med) for v in vals)
    raw_mad = _median(mads)
    return raw_mad * 1.4826


# ── Summarisation ──────────────────────────────────────────────────────────────

def summarise(outdir, rep_prefixes, n_tax, order,
              tree_biparts_file, plot_biparts_file, bipart_order_file,
              no_hex, robust=False, truncate=None):
    """
    Read all replicate outputs and write summary files to outdir.

    robust=True: report median and scaled MAD instead of mean and SD.
    Variability columns are always raw (un-normalised) values.
    """
    n_reps = len(rep_prefixes)
    summary_prefix = os.path.join(outdir, 'summary')

    # ── Compatibility scores ──────────────────────────────────────────────────
    scores = [s for s in (read_compat_score(p) for p in rep_prefixes)
              if s is not None]
    if scores:
        mean_score = statistics.mean(scores)
        sd_score   = statistics.stdev(scores) if len(scores) > 1 else 0.0
        min_score  = min(scores)
        max_score  = max(scores)
    else:
        mean_score = sd_score = min_score = max_score = float('nan')

    # summary.compat.txt written below after all_bm is known (for union count)

    # ── Bipartition count statistics ─────────────────────────────────────────
    n_bipart = [b for b in (read_compat_n_biparts(p) for p in rep_prefixes)
                if b is not None]

    if n_bipart:
        mean_nb = sum(n_bipart) / len(n_bipart)
        sd_nb   = _sd(n_bipart)
        min_nb  = float(min(n_bipart))
        max_nb  = float(max(n_bipart))
        def _pct(vals, p):
            s = sorted(vals)
            idx = (len(s) - 1) * p / 100.0
            lo, hi = int(idx), min(int(idx) + 1, len(s) - 1)
            return s[lo] + (s[hi] - s[lo]) * (idx - lo)
        pct05 = _pct(n_bipart,  5)
        pct25 = _pct(n_bipart, 25)
        pct50 = _pct(n_bipart, 50)
        pct75 = _pct(n_bipart, 75)
        pct95 = _pct(n_bipart, 95)
    else:
        mean_nb = sd_nb = min_nb = max_nb = float('nan')
        pct05 = pct25 = pct50 = pct75 = pct95 = float('nan')

    n_possible = 2 ** (n_tax - 1) - n_tax - 1

    # ── Bipartition averages and variability ─────────────────────────────────
    rep_data = [read_tsv_biparts(p) for p in rep_prefixes]
    all_bm   = set().union(*(d.keys() for d in rep_data))

    union_n_biparts = len(all_bm)

    with open(summary_prefix + '.compat.txt', 'w') as f:
        f.write(f"mean_score\t{mean_score:.6f}\n")
        f.write(f"sd_score\t{sd_score:.6f}\n")
        f.write(f"min_score\t{min_score:.6f}\n")
        f.write(f"max_score\t{max_score:.6f}\n")
        f.write(f"n_replicates\t{len(scores)}\n")
        f.write(f"mean_n_biparts\t{mean_nb:.2f}\n")
        f.write(f"sd_n_biparts\t{sd_nb:.2f}\n")
        f.write(f"min_n_biparts\t{int(min_nb) if not __import__('math').isnan(min_nb) else 'nan'}\n")
        f.write(f"max_n_biparts\t{int(max_nb) if not __import__('math').isnan(max_nb) else 'nan'}\n")
        f.write(f"pct05_n_biparts\t{pct05:.1f}\n")
        f.write(f"pct25_n_biparts\t{pct25:.1f}\n")
        f.write(f"pct50_n_biparts\t{pct50:.1f}\n")
        f.write(f"pct75_n_biparts\t{pct75:.1f}\n")
        f.write(f"pct95_n_biparts\t{pct95:.1f}\n")
        f.write(f"union_n_biparts\t{union_n_biparts}\n")
        f.write(f"total_possible\t{n_possible}\n")

    center_support  = {}
    center_conflict = {}
    spread_support  = {}
    spread_conflict = {}

    for bm in all_bm:
        sups = [rep_data[i].get(bm, (0.0, 0.0))[0] for i in range(n_reps)]
        cons = [rep_data[i].get(bm, (0.0, 0.0))[1] for i in range(n_reps)]
        if robust:
            center_support[bm]  = _median(sups)
            center_conflict[bm] = _median(cons)
            spread_support[bm]  = _scaled_mad(sups)
            spread_conflict[bm] = _scaled_mad(cons)
        else:
            center_support[bm]  = sum(sups) / n_reps
            center_conflict[bm] = sum(cons) / n_reps
            spread_support[bm]  = _sd(sups)
            spread_conflict[bm] = _sd(cons)

    support_dd  = defaultdict(float, center_support)
    conflict_dd = defaultdict(float, center_conflict)

    # Apply plot-biparts filter: ensure all requested bipartitions appear,
    # including those with zero mean support across all replicates
    plot_bitmasks = None
    if plot_biparts_file:
        plot_bitmasks = read_tree_biparts(plot_biparts_file, n_tax)
        for bm in plot_bitmasks:
            if bm not in support_dd:
                support_dd[bm]  = 0.0
                conflict_dd[bm] = 0.0

    support_norm, conflict_norm = normalize_lento(support_dd, conflict_dd)

    order_bitmasks = None
    if bipart_order_file:
        order_bitmasks = read_bipart_order(bipart_order_file, n_tax)

    sd_label = 'mad' if robust else 'sd'
    write_tsv(summary_prefix + '.tsv',
              support_dd, conflict_dd, support_norm, conflict_norm,
              n_tax, order,
              plot_bitmasks=plot_bitmasks,
              order_bitmasks=order_bitmasks,
              sd_support=spread_support,
              sd_conflict=spread_conflict,
              sd_label=sd_label)

    tree_bitmasks = None
    if tree_biparts_file:
        tree_bitmasks = read_tree_biparts(tree_biparts_file, n_tax)

    write_svg(summary_prefix + '.svg',
              support_dd, conflict_dd, support_norm, conflict_norm,
              n_tax, order,
              tree_bitmasks=tree_bitmasks,
              plot_bitmasks=plot_bitmasks,
              order_bitmasks=order_bitmasks,
              show_hex=(not no_hex))

    if truncate is not None:
        write_svg(summary_prefix + '.trunc.svg',
                  support_dd, conflict_dd, support_norm, conflict_norm,
                  n_tax, order,
                  tree_bitmasks=tree_bitmasks,
                  plot_bitmasks=plot_bitmasks,
                  order_bitmasks=order_bitmasks,
                  show_hex=(not no_hex),
                  max_cols=truncate)

    # Tree bipartition position reporting and tree_bipart.tsv
    if tree_bitmasks:
        full_order = build_column_order(
            support_dd, support_norm, conflict_norm,
            plot_bitmasks=plot_bitmasks,
            order_bitmasks=order_bitmasks,
        )
        pos_map       = {bm: i + 1 for i, bm in enumerate(full_order)}
        tree_in_plot  = {bm: pos_map[bm] for bm in tree_bitmasks if bm in pos_map}
        if tree_in_plot:
            max_pos = max(tree_in_plot.values())
            max_bm  = next(bm for bm, p in tree_in_plot.items() if p == max_pos)
            if truncate is not None:
                n_in_trunc = sum(1 for p in tree_in_plot.values()
                                 if p <= truncate)
                print(
                    f"Tree bipartitions in truncated summary plot "
                    f"(first {truncate} columns): {n_in_trunc}",
                    file=sys.stderr,
                )
            print(
                f"Column of last tree bipartition (summary.svg) = {max_pos} "
                f"({bitmask_to_binary_str(max_bm, n_tax)}, "
                f"hex {bitmask_to_jakobsen_hex(max_bm, n_tax)})",
                file=sys.stderr,
            )
        tbp_path = summary_prefix + '.tree_bipart.tsv'
        _write_tree_bipart_tsv(
            tbp_path, tree_bitmasks, support_dd, conflict_dd,
            support_norm, conflict_norm, n_tax, order, pos_map,
        )
        print(f"Wrote {tbp_path}", file=sys.stderr)

    return mean_score, sd_score, len(scores), union_n_biparts, n_possible

# ── Argument parsing ───────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description='Phylogenomic compatibility analysis via per-locus site sampling.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output directory contents:
  replicate_NNNN.phy        binary relaxed PHYLIP for each replicate
  replicate_NNNN.tsv        bipartition table from partcompat.py
  replicate_NNNN.compat.txt compatibility score from partcompat.py
  summary.tsv               mean support/conflict per bipartition
  summary.svg               LentoPlot of bipartition averages
  summary.compat.txt        mean ± SD of compatibility scores
  run_info.txt              seed, parameters, loci used/skipped

Coding modes (--mode) are the same as partcompat.py:
  RY / YR  KM / MK  WS / SW  2state  partimatrix  (default: RY)
        """,
    )
    p.add_argument('-i', '--indir',      required=True,
                   help='Input directory containing alignment files')
    p.add_argument('-o', '--outdir',     required=True,
                   help='Output directory (created if absent)')
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
    p.add_argument('-n', '--replicates', type=int, default=100,
                   help='Number of replicates (default: 100)')
    p.add_argument('--seed',             type=int, default=None,
                   help='Random seed (auto-generated and reported if omitted)')
    p.add_argument('--boot-locus',       action='store_true',
                   help='Locus bootstrapping: multinomial resampling of loci '
                        '(default: one copy of each eligible locus per replicate)')
    p.add_argument('--variable-sites',   default=None,
                   metavar='DIST',
                   help='Draw a variable number of sites per replicate from '
                        'DIST instead of one site per locus. '
                        'DIST formats: uniform:MIN:MAX or fixed:N. '
                        'Without --downsample, sites are drawn i.i.d. from '
                        'the locus pool (N may exceed n_loci). '
                        'With --downsample, N loci are selected without '
                        'replacement and one site is drawn from each '
                        '(requires N <= n_eligible loci).')
    p.add_argument('--downsample',       action='store_true',
                   help='Used with --variable-sites: select N loci without '
                        'replacement and draw one site per locus. '
                        'N must be <= number of eligible loci. '
                        'Incompatible with --boot-locus.')
    p.add_argument('--max-missing-taxa', type=float, default=0.5,
                   metavar='FRAC',
                   help='Skip a locus if the fraction of taxa from --tax-order '
                        'that are absent from that alignment exceeds this value '
                        '(0.0–1.0, default: 0.5)')
    p.add_argument('--codon-positions',  default='all',
                   choices=['all', '1', '2', '3', '12', '13', '23', '123'],
                   help='Codon positions to consider (default: all)')
    p.add_argument('-t', '--tree-biparts',
                   help='Tree bipartitions file (passed to partcompat.py; '
                        'matching bipartitions shown as outline bars)')
    p.add_argument('-p', '--plot-biparts',
                   help='Plot filter file (passed to partcompat.py; '
                        'only these bipartitions appear in SVG/TSV outputs)')
    p.add_argument('-b', '--bipart-order',
                   help='Bipartition column ordering file (passed to partcompat.py '
                        'and applied to summary outputs; same format as '
                        '--tree-biparts)')
    p.add_argument('--robust',           action='store_true',
                   help='Report median and scaled MAD instead of mean and SD '
                        'in summary TSV variability columns')
    p.add_argument('--truncate',          type=int, default=None,
                   metavar='N',
                   help='Also write a truncated Lento plot (first N columns) '
                        'for each replicate and the summary SVG. '
                        'Passed through to partcompat.py.')
    p.add_argument('--no-hex',           action='store_true',
                   help='Suppress Jakobsen hex labels in SVG outputs')
    p.add_argument('--redo',             action='store_true',
                   help='Overwrite existing output directory contents')
    return p.parse_args()

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # ── Validate partcompat.py is available ──────────────────────────────────
    if not os.path.isfile(PARTCOMPAT):
        sys.exit(f"ERROR: partcompat.py not found at {PARTCOMPAT}")

    # ── Seed ─────────────────────────────────────────────────────────────────
    if args.seed is None:
        seed = random.randrange(2 ** 32)
    else:
        seed = args.seed
    rng = random.Random(seed)

    # ── Mode ─────────────────────────────────────────────────────────────────
    mode = resolve_mode(args.mode)

    # ── Codon positions ───────────────────────────────────────────────────────
    codon_positions = CODON_POS_MAP[args.codon_positions]

    # ── Taxon order ───────────────────────────────────────────────────────────
    tax_order = read_tax_order(args.tax_order)
    n_tax     = len(tax_order)
    if n_tax == 0:
        sys.exit(f"ERROR: No taxa found in {args.tax_order}")
    print(f"Taxon order: {n_tax} taxa from {args.tax_order}", file=sys.stderr)

    # ── Output directory ──────────────────────────────────────────────────────
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # Check for existing replicate files
    existing = list(Path(outdir).glob('replicate_*.tsv'))
    if existing and not args.redo:
        sys.exit(
            f"ERROR: Output directory '{outdir}' already contains "
            f"{len(existing)} replicate TSV file(s). "
            f"Use --redo to overwrite."
        )

    # ── Scan loci ─────────────────────────────────────────────────────────────
    print(f"Scanning loci in {args.indir} ...", file=sys.stderr)
    eligible, skipped = scan_loci(
        args.indir, tax_order, mode, codon_positions,
        args.max_missing_taxa,
    )
    n_loci     = len(eligible) + len(skipped)
    n_eligible = len(eligible)
    n_skipped  = len(skipped)

    print(f"  Found {n_loci} file(s): "
          f"{n_eligible} eligible, {n_skipped} skipped",
          file=sys.stderr)

    if n_eligible == 0:
        sys.exit("ERROR: No eligible loci found. "
                 "Check --max-missing-taxa, --mode, and input files.")

    # ── Validate --variable-sites / --downsample ──────────────────────────────
    dist_fn    = None
    downsample = False

    if args.downsample and not args.variable_sites:
        sys.exit("ERROR: --downsample requires --variable-sites.")
    if args.downsample and args.boot_locus:
        sys.exit("ERROR: --downsample and --boot-locus are mutually exclusive.")

    if args.variable_sites:
        dist_fn    = parse_dist(args.variable_sites)
        downsample = args.downsample

        if downsample:
            # Validate that max of distribution <= n_eligible
            # We check by sampling the distribution 1000 times
            test_vals = [dist_fn(rng) for _ in range(1000)]
            max_test  = max(test_vals)
            if max_test > n_eligible:
                sys.exit(
                    f"ERROR: --downsample requires N <= n_eligible ({n_eligible}), "
                    f"but --variable-sites '{args.variable_sites}' can produce "
                    f"N={max_test} which exceeds this. "
                    f"Reduce the MAX value in your distribution."
                )
            # Reset rng to avoid consuming seed values in validation
            rng = random.Random(seed)

    # ── run_info.txt (header) ─────────────────────────────────────────────────
    info_path = os.path.join(outdir, 'run_info.txt')
    with open(info_path, 'w') as info:
        info.write("=== sample_loci.py run information ===\n\n")
        info.write(f"Command:          {' '.join(sys.argv)}\n")
        info.write(f"Seed:             {seed}\n")
        info.write(f"Mode:             {mode}\n")
        info.write(f"Replicates:       {args.replicates}\n")
        info.write(f"Boot-locus:       {args.boot_locus}\n")
        if args.variable_sites:
            info.write(f"Variable-sites:   {args.variable_sites}\n")
            info.write(f"Downsample:       {args.downsample}\n")
        info.write(f"Max missing taxa: {args.max_missing_taxa:.2f}\n")
        info.write(f"Codon positions:  {args.codon_positions}\n")
        info.write(f"Tax order file:   {args.tax_order}\n")
        if args.tree_biparts:
            info.write(f"Tree biparts:     {args.tree_biparts}\n")
        if args.plot_biparts:
            info.write(f"Plot biparts:     {args.plot_biparts}\n")
        info.write(f"\nLoci found:    {n_loci}\n")
        info.write(f"Loci eligible: {n_eligible}\n")
        info.write(f"Loci skipped:  {n_skipped}\n")
        if skipped:
            info.write("\nSkipped loci:\n")
            for path, reason in skipped:
                info.write(f"  {os.path.basename(path)}: {reason}\n")
        info.write(f"\nEligible loci:\n")
        for path, sites in eligible:
            info.write(f"  {os.path.basename(path)}: {len(sites)} PI sites\n")
        info.write("\n")

    # ── Replicates ────────────────────────────────────────────────────────────
    n_digits      = len(str(args.replicates))
    rep_prefixes  = []
    n_failed      = 0

    if dist_fn is not None:
        if downsample:
            _mode_desc = f"downsample {args.variable_sites}"
        elif args.boot_locus:
            _mode_desc = f"locus bootstrap + variable sites {args.variable_sites}"
        else:
            _mode_desc = f"variable sites {args.variable_sites}"
    elif args.boot_locus:
        _mode_desc = "locus bootstrap"
    else:
        _mode_desc = "uniform sampling"

    print(f"Running {args.replicates} replicate(s) ({_mode_desc}) ...",
          file=sys.stderr)

    for rep in range(1, args.replicates + 1):
        rep_label  = str(rep).zfill(n_digits)
        phy_path   = os.path.join(outdir, f"replicate_{rep_label}.phy")
        rep_prefix = os.path.join(outdir, f"replicate_{rep_label}")

        # Sample sites
        coded_vectors = sample_replicate(
            eligible, args.boot_locus, rng,
            dist_fn=dist_fn, downsample=downsample,
        )

        # Write binary PHYLIP
        write_binary_phylip(coded_vectors, tax_order, phy_path)

        # Run partcompat.py
        success, stderr_text = run_partcompat(
            phy_path, rep_prefix, args.tax_order, mode,
            args.tree_biparts, args.plot_biparts, args.bipart_order,
            args.codon_positions, args.no_hex, truncate=args.truncate,
        )

        if success:
            rep_prefixes.append(rep_prefix)
        else:
            n_failed += 1
            print(f"  WARNING: replicate {rep_label} failed:\n{stderr_text.strip()}",
                  file=sys.stderr)

        # Progress
        if rep % max(1, args.replicates // 10) == 0 or rep == args.replicates:
            print(f"  {rep}/{args.replicates} replicates complete",
                  file=sys.stderr)

    n_ok = len(rep_prefixes)
    print(f"Replicates: {n_ok} succeeded, {n_failed} failed", file=sys.stderr)

    if n_ok == 0:
        sys.exit("ERROR: All replicates failed. "
                 "Check that the alignment has parsimony-informative sites.")

    # ── Summarise ─────────────────────────────────────────────────────────────
    print("Summarising replicates ...", file=sys.stderr)
    mean_score, sd_score, n_scored, union_nb, n_possible = summarise(
        outdir, rep_prefixes, n_tax, tax_order,
        args.tree_biparts, args.plot_biparts, args.bipart_order,
        args.no_hex, robust=args.robust, truncate=args.truncate,
    )

    # ── run_info.txt (footer) ─────────────────────────────────────────────────
    with open(info_path, 'a') as info:
        info.write(f"Replicates completed: {n_ok}/{args.replicates}\n")
        info.write(f"Replicates failed:    {n_failed}\n")
        if not math.isnan(mean_score):
            info.write(f"Mean compatibility:   {mean_score:.6f}\n")
            info.write(f"SD compatibility:     {sd_score:.6f}\n")
            info.write(f"Union bipartitions:   {union_nb}\n")
            info.write(f"Total possible:       {n_possible}\n")

    # ── Final report ──────────────────────────────────────────────────────────
    summary_prefix = os.path.join(outdir, 'summary')
    print(f"\nCompatibility: mean={mean_score:.6f}  SD={sd_score:.6f}  "
          f"(n={n_scored})")
    print(f"Bipartitions:  union={union_nb} across all replicates  "
          f"of {n_possible} possible")
    print(f"Wrote {summary_prefix}.tsv", file=sys.stderr)
    print(f"Wrote {summary_prefix}.svg", file=sys.stderr)
    print(f"Wrote {summary_prefix}.compat.txt", file=sys.stderr)
    print(f"Wrote {info_path}", file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
