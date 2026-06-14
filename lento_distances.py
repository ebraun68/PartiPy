#!/usr/bin/env python3
"""
lento_distances.py - Compare Lento spectra of individual loci to the
                     genome-wide average from sample_loci.py.

Reads TSV files produced by sample_loci.py (per-replicate) and
analyze_loci.py (per-gene), computes distances from each dataset's
Lento vector to the genome-wide centroid, and writes a summary TSV
and strip-plot SVG.

Distance measures:
  euclidean  - L2 distance on sum-to-1 normalised vectors
  manhattan  - L1 distance on sum-to-1 normalised vectors
  cosine     - 1 - cosine similarity (scale-invariant)

Vector options:
  support    - support_norm values only
  both       - concatenated [support_norm, conflict_norm]

Usage:
  python3 lento_distances.py -s SAMPLE_DIR -g GENE_DIR -o PREFIX [options]
"""

import argparse
import csv
import math
import os
import statistics
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# TSV reading
# ---------------------------------------------------------------------------

def read_lento_tsv(path):
    """
    Read a partcompat.py TSV file.

    Returns a dict:
      bipartition_binary -> {'support_norm': float, 'conflict_norm': float,
                             'support_raw': float, 'conflict_raw': float}
    Skips footer rows (rows where bipartition_binary starts with a letter,
    i.e. summary stat rows written by sample_loci.py).
    Returns None on failure.
    """
    result = {}
    try:
        with open(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                bp = row.get('bipartition_binary', '').strip().lstrip("'")
                # Skip summary footer rows
                if not bp or not all(c in '01' for c in bp):
                    continue
                try:
                    result[bp] = {
                        'support_norm':  float(row.get('support_norm',  0)),
                        'conflict_norm': float(row.get('conflict_norm', 0)),
                        'support_raw':   float(row.get('support_raw',   0)),
                        'conflict_raw':  float(row.get('conflict_raw',  0)),
                    }
                except (ValueError, KeyError):
                    continue
    except Exception:
        return None
    return result if result else None


def read_compat_score(prefix):
    """Read compatibility score from OUTFILE.compat.txt line 1."""
    try:
        with open(prefix + '.compat.txt') as f:
            return float(f.readline().strip().split()[0])
    except Exception:
        return None

# ---------------------------------------------------------------------------
# Vector construction and normalisation
# ---------------------------------------------------------------------------

def build_vector(data, bipart_keys, vector_type):
    """
    Build a numeric vector from a bipartition data dict.

    data        : {bipartition_binary -> {'support_norm', 'conflict_norm', ...}}
    bipart_keys : ordered list of bipartition binary strings (common set)
    vector_type : 'support' or 'both'

    Returns list of floats (absent bipartitions contribute 0).
    """
    if vector_type == 'support':
        return [data.get(bp, {}).get('support_norm', 0.0) for bp in bipart_keys]
    else:  # 'both'
        sup = [data.get(bp, {}).get('support_norm',  0.0) for bp in bipart_keys]
        con = [data.get(bp, {}).get('conflict_norm', 0.0) for bp in bipart_keys]
        return sup + con


def sum_to_one(vec):
    """
    Scale a vector so the first half (support values) sums to 1.
    For 'support' vectors the whole vector sums to 1.
    For 'both' vectors the support half sums to 1; the conflict half is scaled
    by the same factor (preserving the Lento relationship).
    Returns the scaled vector and the original support sum.
    """
    n = len(vec)
    # Support is always the first n//2 elements for 'both', or all for 'support'
    sup_sum = sum(vec[:n // 2]) if n > 1 else sum(vec)
    if sup_sum == 0:
        return vec[:], 0.0
    scale = 1.0 / sup_sum
    return [v * scale for v in vec], sup_sum


def raw_support_sum(data):
    """Sum of support_raw across all bipartitions in a dataset."""
    return sum(v['support_raw'] for v in data.values())

# ---------------------------------------------------------------------------
# Distance functions
# ---------------------------------------------------------------------------

def euclidean(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def manhattan(a, b):
    return sum(abs(x - y) for x, y in zip(a, b))


def cosine(a, b):
    dot  = sum(x * y for x, y in zip(a, b))
    na   = math.sqrt(sum(x * x for x in a))
    nb   = math.sqrt(sum(y * y for y in b))
    if na == 0 or nb == 0:
        return float('nan')
    return 1.0 - dot / (na * nb)


DIST_FUNCS = {
    'euclidean': euclidean,
    'manhattan': manhattan,
    'cosine':    cosine,
}

# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

def _sd(vals):
    n = len(vals)
    if n <= 1:
        return 0.0
    mu = sum(vals) / n
    return math.sqrt(sum((v - mu) ** 2 for v in vals) / (n - 1))


def _pct(vals, p):
    if not vals:
        return float('nan')
    s   = sorted(v for v in vals if not math.isnan(v))
    idx = (len(s) - 1) * p / 100.0
    lo  = int(idx)
    hi  = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (idx - lo)

# ---------------------------------------------------------------------------
# SVG strip plot
# ---------------------------------------------------------------------------

def write_svg(path, rows, dist_cols, show_hex=True):
    """
    Strip plot SVG: one horizontal panel per distance column.
    Each panel has two rows of dots: replicates (blue) and genes (orange).
    Dots are jittered vertically within each row.
    Mean ± SD annotations are added for each population.
    """
    import random
    rng = random.Random(42)

    # Separate replicates and genes
    rep_rows  = [r for r in rows if r['type'] == 'replicate']
    gene_rows = [r for r in rows if r['type'] == 'gene']

    n_panels     = len(dist_cols)
    PANEL_H      = 120
    PANEL_GAP    = 30
    MARGIN_LEFT  = 160
    MARGIN_RIGHT = 30
    MARGIN_TOP   = 30
    MARGIN_BOT   = 20
    STRIP_W      = max(500, 6 * (len(rep_rows) + len(gene_rows)))
    DOT_R        = 4
    JITTER       = 12    # vertical jitter range (px)
    LABEL_FONT   = 10
    ANNOT_FONT   = 9

    total_h = MARGIN_TOP + n_panels * PANEL_H + (n_panels - 1) * PANEL_GAP + MARGIN_BOT
    total_w = MARGIN_LEFT + STRIP_W + MARGIN_RIGHT

    svg = []
    svg.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{total_w}" height="{total_h}" '
        f'font-family="Helvetica,Arial,sans-serif">'
    )
    svg.append(
        f'<rect width="{total_w}" height="{total_h}" fill="white" stroke="none"/>'
    )

    # Legend
    lx = MARGIN_LEFT
    ly = MARGIN_TOP - 10
    svg.append(
        f'<circle cx="{lx}" cy="{ly}" r="{DOT_R}" fill="#1565c0"/>'
        f'<text x="{lx + 8}" y="{ly + 4}" font-size="{LABEL_FONT}" fill="#333">'
        f'Replicates (n={len(rep_rows)})</text>'
    )
    svg.append(
        f'<circle cx="{lx + 160}" cy="{ly}" r="{DOT_R}" fill="#e65100"/>'
        f'<text x="{lx + 168}" y="{ly + 4}" font-size="{LABEL_FONT}" fill="#333">'
        f'Genes (n={len(gene_rows)})</text>'
    )

    for pi, col in enumerate(dist_cols):
        panel_y = MARGIN_TOP + pi * (PANEL_H + PANEL_GAP)
        axis_y  = panel_y + PANEL_H * 0.75   # axis line at 75% of panel height

        # Collect values for this column (skip nan)
        rep_vals  = [r[col] for r in rep_rows  if not math.isnan(r.get(col, float('nan')))]
        gene_vals = [r[col] for r in gene_rows if not math.isnan(r.get(col, float('nan')))]
        all_vals  = rep_vals + gene_vals
        if not all_vals:
            continue

        x_min = 0.0
        x_max = max(all_vals) * 1.05 if all_vals else 1.0
        x_range = x_max - x_min if x_max > x_min else 1.0

        def to_x(v):
            return MARGIN_LEFT + (v - x_min) / x_range * STRIP_W

        # Axis line
        svg.append(
            f'<line x1="{MARGIN_LEFT}" y1="{axis_y:.1f}" '
            f'x2="{MARGIN_LEFT + STRIP_W}" y2="{axis_y:.1f}" '
            f'stroke="#888" stroke-width="1"/>'
        )

        # Axis label
        svg.append(
            f'<text x="{MARGIN_LEFT - 8}" y="{axis_y + 4:.1f}" '
            f'text-anchor="end" font-size="{LABEL_FONT}" fill="#333">'
            f'{_xml_escape(col)}</text>'
        )

        # X-axis ticks (at 0, 25%, 50%, 75%, 100% of range)
        for frac in [0, 0.25, 0.5, 0.75, 1.0]:
            tv  = x_min + frac * x_range
            tx  = to_x(tv)
            svg.append(
                f'<line x1="{tx:.1f}" y1="{axis_y:.1f}" '
                f'x2="{tx:.1f}" y2="{axis_y + 4:.1f}" '
                f'stroke="#888" stroke-width="0.8"/>'
            )
            svg.append(
                f'<text x="{tx:.1f}" y="{axis_y + 14:.1f}" '
                f'text-anchor="middle" font-size="7" fill="#666">'
                f'{tv:.3f}</text>'
            )

        # Dots: replicates above axis, genes below
        for row_set, color, y_base, sign in [
            (rep_rows,  '#1565c0', axis_y - 18, -1),
            (gene_rows, '#e65100', axis_y + 18,  1),
        ]:
            for row in row_set:
                v = row.get(col, float('nan'))
                if math.isnan(v):
                    continue
                cx = to_x(v)
                cy = y_base + rng.uniform(-JITTER / 2, JITTER / 2)
                svg.append(
                    f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{DOT_R}" '
                    f'fill="{color}" fill-opacity="0.55" '
                    f'stroke="{color}" stroke-width="0.5"/>'
                )

        # Mean ± SD annotations
        for vals, color, y_base, label in [
            (rep_vals,  '#1565c0', axis_y - 38, 'R'),
            (gene_vals, '#e65100', axis_y + 34, 'G'),
        ]:
            if not vals:
                continue
            mu  = sum(vals) / len(vals)
            sd  = _sd(vals)
            mx  = to_x(mu)
            # Mean tick
            svg.append(
                f'<line x1="{mx:.1f}" y1="{y_base - 6:.1f}" '
                f'x2="{mx:.1f}" y2="{y_base + 6:.1f}" '
                f'stroke="{color}" stroke-width="2"/>'
            )
            # SD bar
            x1 = to_x(max(x_min, mu - sd))
            x2 = to_x(min(x_max, mu + sd))
            svg.append(
                f'<line x1="{x1:.1f}" y1="{y_base:.1f}" '
                f'x2="{x2:.1f}" y2="{y_base:.1f}" '
                f'stroke="{color}" stroke-width="1.5"/>'
            )
            svg.append(
                f'<text x="{MARGIN_LEFT - 8}" y="{y_base + 3:.1f}" '
                f'text-anchor="end" font-size="{ANNOT_FONT}" fill="{color}">'
                f'{label}: {mu:.3f}±{sd:.3f}</text>'
            )

    svg.append('</svg>')
    with open(path, 'w') as f:
        f.write('\n'.join(svg))


def _xml_escape(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Compute distances between individual-locus and '
                    'genome-wide Lento spectra.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Both -s/--sample-dir and -g/--gene-dir must have been produced using the
same --tax-order, --mode, and --plot-biparts settings.

Distance columns in output TSV are named as DISTANCE_VECTOR, e.g.:
  euclidean_support, euclidean_both, manhattan_support, cosine_both, ...

Example:
  python3 lento_distances.py \\
      -s genome_wide_results/ \\
      -g gene_results/ \\
      -o distances/comparison
        ''',
    )
    p.add_argument('-s', '--sample-dir', required=True,
                   help='sample_loci.py output directory '
                        '(contains summary.tsv and replicate_*.tsv)')
    p.add_argument('-g', '--gene-dir',   required=True,
                   help='analyze_loci.py output directory '
                        '(contains one TSV per gene)')
    p.add_argument('-o', '--outfile',    required=True,
                   help='Output prefix; produces PREFIX.tsv and PREFIX.svg')
    p.add_argument('--distances', default='all',
                   help="Comma-separated distances to compute: euclidean, "
                        "manhattan, cosine, or 'all' (default: all)")
    p.add_argument('--vectors', default='both',
                   choices=['support', 'both'],
                   help="Vector composition: 'support' (support_norm only) or "
                        "'both' (support_norm + conflict_norm) (default: both)")
    p.add_argument('--no-svg', action='store_true',
                   help='Suppress SVG output')
    return p.parse_args()

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    # ── Parse distance selection ──────────────────────────────────────────────
    if args.distances.lower() == 'all':
        dist_names = ['euclidean', 'manhattan', 'cosine']
    else:
        dist_names = [d.strip().lower() for d in args.distances.split(',')]
        unknown = [d for d in dist_names if d not in DIST_FUNCS]
        if unknown:
            sys.exit(f"ERROR: Unknown distance(s): {unknown}. "
                     f"Choose from: euclidean, manhattan, cosine")

    vec_types = ['support', 'both'] if args.vectors == 'both' else ['support']

    dist_cols = [f"{d}_{v}" for d in dist_names for v in vec_types]

    # ── Load summary TSV (centroid) ───────────────────────────────────────────
    summary_path = os.path.join(args.sample_dir, 'summary.tsv')
    if not os.path.isfile(summary_path):
        sys.exit(f"ERROR: summary.tsv not found in '{args.sample_dir}'. "
                 f"Run sample_loci.py first.")

    centroid_data = read_lento_tsv(summary_path)
    if not centroid_data:
        sys.exit(f"ERROR: Could not read bipartition data from {summary_path}")

    bipart_keys = sorted(centroid_data.keys())
    n_biparts   = len(bipart_keys)
    print(f"Centroid: {n_biparts} bipartitions from {summary_path}",
          file=sys.stderr)

    # Build centroid vectors (one per vec_type; sum-to-1 normalised)
    centroid_vecs = {}
    for vt in vec_types:
        raw_vec, _ = sum_to_one(build_vector(centroid_data, bipart_keys, vt))
        centroid_vecs[vt] = raw_vec

    # ── Load replicate TSVs ───────────────────────────────────────────────────
    rep_files = sorted(Path(args.sample_dir).glob('replicate_*.tsv'))
    if not rep_files:
        print(f"WARNING: No replicate_*.tsv files found in '{args.sample_dir}'",
              file=sys.stderr)

    print(f"Loading {len(rep_files)} replicate TSV(s) ...", file=sys.stderr)
    rep_rows = []
    for rf in rep_files:
        data = read_lento_tsv(str(rf))
        if data is None:
            print(f"  WARNING: could not read {rf.name}", file=sys.stderr)
            continue

        # Warn if bipartition set differs significantly
        rep_bps = set(data.keys())
        diff    = len(rep_bps.symmetric_difference(set(bipart_keys)))
        if diff > n_biparts * 0.5:
            print(f"  WARNING: {rf.name} has very different bipartitions "
                  f"({diff} differing); check --plot-biparts consistency.",
                  file=sys.stderr)

        row = {'dataset': rf.stem, 'type': 'replicate',
               'n_biparts':     len(data),
               'sum_support_raw': raw_support_sum(data)}

        # Compat score
        compat_prefix = str(rf).replace('.tsv', '')
        row['compat_score'] = read_compat_score(compat_prefix)

        for vt in vec_types:
            vec, _ = sum_to_one(build_vector(data, bipart_keys, vt))
            for dn in dist_names:
                fn  = DIST_FUNCS[dn]
                key = f"{dn}_{vt}"
                row[key] = fn(vec, centroid_vecs[vt])
        rep_rows.append(row)

    # ── Load gene TSVs ────────────────────────────────────────────────────────
    if not os.path.isdir(args.gene_dir):
        sys.exit(f"ERROR: Gene directory not found: {args.gene_dir}")
    gene_files = sorted(
        f for f in Path(args.gene_dir).glob('*.tsv')
        if not f.name.startswith('run_info')
    )
    if not gene_files:
        sys.exit(f"ERROR: No TSV files found in '{args.gene_dir}'.")

    print(f"Loading {len(gene_files)} gene TSV(s) ...", file=sys.stderr)
    gene_rows = []
    for gf in gene_files:
        data = read_lento_tsv(str(gf))
        if data is None:
            print(f"  WARNING: could not read {gf.name}", file=sys.stderr)
            continue

        gene_bps = set(data.keys())
        diff     = len(gene_bps.symmetric_difference(set(bipart_keys)))
        if diff > n_biparts * 0.5:
            print(f"  WARNING: {gf.name} has very different bipartitions; "
                  f"check --plot-biparts consistency.", file=sys.stderr)

        row = {'dataset': gf.stem, 'type': 'gene',
               'n_biparts':       len(data),
               'sum_support_raw': raw_support_sum(data)}

        compat_prefix = str(gf).replace('.tsv', '')
        row['compat_score'] = read_compat_score(compat_prefix)

        for vt in vec_types:
            vec, _ = sum_to_one(build_vector(data, bipart_keys, vt))
            for dn in dist_names:
                fn  = DIST_FUNCS[dn]
                key = f"{dn}_{vt}"
                row[key] = fn(vec, centroid_vecs[vt])
        gene_rows.append(row)

    all_rows = rep_rows + gene_rows
    print(f"Distances computed for {len(rep_rows)} replicates and "
          f"{len(gene_rows)} genes.", file=sys.stderr)

    # ── Summary to stdout ─────────────────────────────────────────────────────
    print()
    for col in dist_cols:
        rv = [r[col] for r in rep_rows  if not math.isnan(r.get(col, float('nan')))]
        gv = [r[col] for r in gene_rows if not math.isnan(r.get(col, float('nan')))]
        print(f"{col}:")
        if rv:
            print(f"  replicates: mean={sum(rv)/len(rv):.4f}  "
                  f"SD={_sd(rv):.4f}  "
                  f"median={_pct(rv,50):.4f}  (n={len(rv)})")
        if gv:
            print(f"  genes:      mean={sum(gv)/len(gv):.4f}  "
                  f"SD={_sd(gv):.4f}  "
                  f"median={_pct(gv,50):.4f}  (n={len(gv)})")

    # ── Write TSV ─────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(args.outfile + '.tsv')),
                exist_ok=True)
    tsv_path = args.outfile + '.tsv'
    fieldnames = (['dataset', 'type', 'n_biparts', 'sum_support_raw',
                   'compat_score'] + dist_cols)

    with open(tsv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        for row in all_rows:
            out = {k: (f"{row[k]:.6f}" if isinstance(row[k], float)
                       else row[k])
                   for k in fieldnames}
            writer.writerow(out)

        # Summary footer
        f.write('\n')
        for label, subset in [('mean_replicates', rep_rows),
                               ('sd_replicates',  rep_rows),
                               ('mean_genes',      gene_rows),
                               ('sd_genes',        gene_rows)]:
            fn = (statistics.mean if label.startswith('mean')
                  else _sd)
            out = {k: '' for k in fieldnames}
            out['dataset'] = label
            out['type']    = 'summary'
            for col in dist_cols:
                vals = [r[col] for r in subset
                        if not math.isnan(r.get(col, float('nan')))]
                if vals:
                    out[col] = f"{fn(vals):.6f}"
            f.write('\t'.join(str(out.get(k, '')) for k in fieldnames) + '\n')

    print(f"\nWrote {tsv_path}", file=sys.stderr)

    # ── Write SVG ─────────────────────────────────────────────────────────────
    if not args.no_svg:
        svg_path = args.outfile + '.svg'
        write_svg(svg_path, all_rows, dist_cols)
        print(f"Wrote {svg_path}", file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        sys.exit(0)
