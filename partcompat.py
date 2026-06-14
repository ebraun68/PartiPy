#!/usr/bin/env python3
"""
partcompat.py - Compatibility and partition analysis of aligned sequences.

Implements ideas from:
  Jakobsen & Easteal (1996) CABIOS 12:291-295  [reticulate / compatibility matrix]
  Jakobsen, Wilson & Easteal (1997) MBE 14:474-484  [PartiMatrix / partition matrix]

Outputs:
  OUTFILE.tsv  - bipartition support/conflict table
  OUTFILE.svg  - LentoPlot (support up, conflict down per bipartition)
  stdout       - compatibility score
"""

import argparse
import os
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# IUPAC nucleotide ambiguity codes -> sets of unambiguous bases
IUPAC = {
    'A': frozenset('A'),
    'C': frozenset('C'),
    'G': frozenset('G'),
    'T': frozenset('T'),
    'U': frozenset('T'),
    'R': frozenset('AG'),
    'Y': frozenset('CT'),
    'S': frozenset('CG'),
    'W': frozenset('AT'),
    'K': frozenset('GT'),
    'M': frozenset('AC'),
    'B': frozenset('CGT'),
    'D': frozenset('AGT'),
    'H': frozenset('ACT'),
    'V': frozenset('ACG'),
    'N': frozenset('ACGT'),
    '-': frozenset(),
    '?': frozenset(),
    '.': frozenset(),
}

# Two-state IUPAC codes (exactly two bases)
TWO_STATE_IUPAC = frozenset('RYSWKM')

# Mode definitions: (group0_bases, group1_bases,
#                    iupac_codes_to_0, iupac_codes_to_1)
# partimatrix and 2state are handled separately (column-dependent coding).
MODE_AXES = {
    'RY': (frozenset('AG'),  frozenset('CT'),  {'R': 0}, {'Y': 1}),
    'KM': (frozenset('AC'),  frozenset('GT'),  {'M': 0}, {'K': 1}),
    'WS': (frozenset('AT'),  frozenset('CG'),  {'W': 0}, {'S': 1}),
}

# Codon position sets (1-based)
CODON_POS_MAP = {
    '1':   {1},
    '2':   {2},
    '3':   {3},
    '12':  {1, 2},
    '13':  {1, 3},
    '23':  {2, 3},
    '123': {1, 2, 3},
    'all': None,
}

# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def parse_fasta(text):
    seqs = {}
    order = []
    name = None
    buf = []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('>'):
            if name is not None:
                seqs[name] = ''.join(buf).upper()
                order.append(name)
            name = line[1:].split()[0]
            buf = []
        elif line:
            buf.append(line)
    if name is not None:
        seqs[name] = ''.join(buf).upper()
        order.append(name)
    return order, seqs


def parse_phylip(text):
    """Relaxed PHYLIP: name then sequence, separated by whitespace."""
    lines = [l for l in text.splitlines() if l.strip()]
    header = lines[0].split()
    ntax = int(header[0])
    seqs = {}
    order = []
    i = 1
    while len(order) < ntax and i < len(lines):
        parts = lines[i].split(None, 1)
        if len(parts) == 2:
            name, seq = parts
            seqs[name] = seq.replace(' ', '').upper()
            order.append(name)
        i += 1
    # interleaved blocks
    while i < len(lines):
        line = lines[i].strip()
        if line:
            parts = line.split()
            if parts[0] in seqs:
                seqs[parts[0]] += ''.join(parts[1:]).upper()
            else:
                idx = (i - ntax - 1) % ntax
                if idx < len(order):
                    seqs[order[idx]] += ''.join(parts).upper()
        i += 1
    return order, seqs


def parse_nexus(text):
    seqs = {}
    order = []
    in_matrix = False
    for line in text.splitlines():
        stripped = line.strip()
        upper = stripped.upper()
        if upper == 'MATRIX':
            in_matrix = True
            continue
        if in_matrix:
            if upper == ';' or stripped == ';':
                in_matrix = False
                continue
            if stripped and not stripped.startswith('['):
                parts = stripped.split(None, 1)
                if len(parts) == 2:
                    name, seq = parts
                    seq = seq.replace(' ', '').replace(';', '').upper()
                    if name not in seqs:
                        order.append(name)
                        seqs[name] = seq
                    else:
                        seqs[name] += seq
    return order, seqs


def detect_and_parse(text):
    stripped = text.strip().upper()
    if stripped.startswith('#NEXUS'):
        return parse_nexus(text)
    elif stripped.startswith('>'):
        return parse_fasta(text)
    else:
        return parse_phylip(text)


def read_alignment(path):
    if not os.path.isfile(path):
        sys.exit(f"ERROR: Input file not found: {path}")
    with open(path) as f:
        text = f.read()
    order, seqs = detect_and_parse(text)
    if not order:
        sys.exit(f"ERROR: No sequences found in {path}")
    lengths = set(len(s) for s in seqs.values())
    if len(lengths) > 1:
        sys.exit(f"ERROR: Sequences have different lengths: {lengths}")
    return order, seqs


def read_tax_order(path):
    if not os.path.isfile(path):
        sys.exit(f"ERROR: Taxon order file not found: {path}")
    with open(path) as f:
        return [l.strip() for l in f if l.strip()]


def is_binary_data(seqs):
    """Return True if all characters are 0/1 (plus gaps/missing)."""
    valid = set('01-?. \t\n')
    for seq in seqs.values():
        if not all(c in valid for c in seq):
            return False
    return True


def read_tree_biparts(path, n_tax):
    """
    Read a plain-text file of bipartitions, one per line.
    Each line is either:
      - a binary string of length n_tax (e.g. '1100010')
      - a hex string in Jakobsen convention (e.g. '6', '0x6', '6A')
        (Jakobsen hex = canonical bitmask XOR full_mask)
    Returns a set of canonical bitmask integers.
    Lines beginning with '#' and blank lines are ignored.
    """
    if not os.path.isfile(path):
        sys.exit(f"ERROR: Bipartition file not found: {path}")
    bitmasks = set()
    full_mask = (1 << n_tax) - 1
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            token = line.strip()
            if not token or token.startswith('#'):
                continue
            if all(c in '01' for c in token):
                # Binary string
                if len(token) != n_tax:
                    sys.exit(
                        f"ERROR: {path} line {lineno}: binary string '{token}' "
                        f"has length {len(token)}, expected {n_tax}"
                    )
                raw   = [int(c) for c in token]
                canon = _canonicalize(raw)
                bitmasks.add(_coded_to_bitmask(canon))
            else:
                # Hex string (Jakobsen convention)
                hex_tok = token.upper().lstrip('0X') or '0'
                try:
                    jak_val   = int(hex_tok, 16)
                    canon_val = jak_val ^ full_mask
                    if canon_val.bit_length() > n_tax:
                        sys.exit(
                            f"ERROR: {path} line {lineno}: hex '{token}' "
                            f"exceeds {n_tax} bits"
                        )
                    canon_coded = _bitmask_to_coded(canon_val, n_tax)
                    canon_coded = _canonicalize(canon_coded)
                    bitmasks.add(_coded_to_bitmask(canon_coded))
                except ValueError:
                    sys.exit(
                        f"ERROR: {path} line {lineno}: cannot parse '{token}' "
                        f"as binary or hexadecimal"
                    )
    return bitmasks

def read_bipart_order(path, n_tax):
    """
    Read a plain-text ordering file for bipartitions.
    Same format as --tree-biparts (binary strings or Jakobsen hex, # comments,
    blank lines ignored), but each bipartition must appear exactly once.
    Returns an ordered list of canonical bitmask integers.
    """
    if not os.path.isfile(path):
        sys.exit(f"ERROR: Bipartition order file not found: {path}")
    bitmasks = []
    seen     = {}         # bitmask -> line number (for duplicate detection)
    full_mask = (1 << n_tax) - 1
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            token = line.strip()
            if not token or token.startswith('#'):
                continue
            if all(c in '01' for c in token):
                if len(token) != n_tax:
                    sys.exit(
                        f"ERROR: {path} line {lineno}: binary string '{token}' "
                        f"has length {len(token)}, expected {n_tax}"
                    )
                raw   = [int(c) for c in token]
                canon = _canonicalize(raw)
                bm    = _coded_to_bitmask(canon)
            else:
                hex_tok = token.upper().lstrip('0X') or '0'
                try:
                    jak_val     = int(hex_tok, 16)
                    canon_val   = jak_val ^ full_mask
                    if canon_val.bit_length() > n_tax:
                        sys.exit(
                            f"ERROR: {path} line {lineno}: hex '{token}' "
                            f"exceeds {n_tax} bits"
                        )
                    canon_coded = _bitmask_to_coded(canon_val, n_tax)
                    canon_coded = _canonicalize(canon_coded)
                    bm          = _coded_to_bitmask(canon_coded)
                except ValueError:
                    sys.exit(
                        f"ERROR: {path} line {lineno}: cannot parse '{token}' "
                        f"as binary or hexadecimal"
                    )
            if bm in seen:
                sys.exit(
                    f"ERROR: {path} line {lineno}: duplicate bipartition "
                    f"'{bitmask_to_binary_str(bm, n_tax)}' "
                    f"(first seen at line {seen[bm]})"
                )
            seen[bm] = lineno
            bitmasks.append(bm)
    return bitmasks


# ---------------------------------------------------------------------------
# Site coding
# ---------------------------------------------------------------------------

def code_site_partimatrix(col_chars):
    """
    PartiMatrix/historical mode:
      - IUPAC two-state ambiguities (RYSWKM) and all other non-ACGT -> missing
      - Exactly 2 distinct unambiguous bases: treat as binary
      - 3 or 4 distinct bases: RY-collapse then treat as binary
    """
    resolved = []
    for c in col_chars:
        if c in ('A', 'C', 'G', 'T'):
            resolved.append(c)
        else:
            resolved.append(None)

    bases = set(b for b in resolved if b is not None)

    if len(bases) <= 1:
        return [0 if b is not None else None for b in resolved]

    if len(bases) == 2:
        b0, b1 = sorted(bases)
        raw = [None if b is None else (0 if b == b0 else 1) for b in resolved]
    else:
        # 3 or 4 bases: RY collapse
        def ry(b):
            if b in 'AG': return 0
            if b in 'CT': return 1
        raw = [None if b is None else ry(b) for b in resolved]

    return _canonicalize(raw)


def code_site_axis(col_chars, group0, group1, iupac_to_0, iupac_to_1):
    """Axis-based modes (RY, KM, WS)."""
    raw = []
    for c in col_chars:
        if c in group0:
            raw.append(0)
        elif c in group1:
            raw.append(1)
        elif c in iupac_to_0:
            raw.append(0)
        elif c in iupac_to_1:
            raw.append(1)
        else:
            raw.append(None)
    return _canonicalize(raw)


def code_site_2state(col_chars):
    """
    2state mode: only columns with exactly 2 distinct unambiguous bases.
    Returns None if the column should be excluded.
    """
    resolved = [c if c in ('A', 'C', 'G', 'T') else None for c in col_chars]
    bases = set(b for b in resolved if b is not None)
    if len(bases) != 2:
        return None
    b0, b1 = sorted(bases)
    raw = [None if b is None else (0 if b == b0 else 1) for b in resolved]
    return _canonicalize(raw)


def code_site_binary(col_chars):
    """Binary data (0/1 input). All other characters -> missing."""
    raw = []
    for c in col_chars:
        if c == '0':
            raw.append(0)
        elif c == '1':
            raw.append(1)
        else:
            raw.append(None)
    return _canonicalize(raw)



def is_perfect_site(col_chars, group0, group1):
    """
    Return True if the raw nucleotide column qualifies as a perfect site
    for the given axis (group0/group1).

    A perfect site has:
      - No IUPAC ambiguity codes, gaps, or missing characters of any kind
      - Exactly two distinct nucleotide states (from A, C, G, T only)
      - The two states fall on opposite sides of the axis
        (one in group0, one in group1)

    Matches the definition in Tiley et al. (2020): parsimony-informative
    sites with exactly two states, no gaps, differing by a transversion
    (for perfectRY) or the equivalent for SW/KM axes.
    """
    unambig = frozenset("ACGT")
    seen = set()
    for c in col_chars:
        if c not in unambig:
            return False   # any ambiguity, gap, or missing disqualifies
        seen.add(c)
        if len(seen) > 2:
            return False   # more than two nucleotide states
    if len(seen) < 2:
        return False       # invariant
    return bool(seen & group0) and bool(seen & group1)


_TRANSAG_EXCLUDE = frozenset({'C','T','Y','S','W','K','M','B','D','H','V'})
_TRANSCT_EXCLUDE = frozenset({'A','G','R','S','W','K','M','B','D','H','V'})
_UNAMBIG_BASES   = frozenset('ACGT')


def code_site_transitions(col, require_ag=False, require_ct=False):
    """
    Code a site for transition-only analysis.

    require_ag : True for transitionsAG (only AG columns)
    require_ct : True for transitionsCT (only CT columns)
    Default    : general transitions - axis inferred from unambiguous bases

    Any character definitively indicating a transversion excludes the site.
    R in an AG column and Y in a CT column become missing (None), triggering
    the standard +0.5 logic when they are the single missing taxon.
    N, ?, - are always missing.  Returns list of 0/1/None or None to exclude.
    """
    if require_ag:
        if any(c in _TRANSAG_EXCLUDE for c in col):
            return None
        return [0 if c == 'A' else (1 if c == 'G' else None) for c in col]

    if require_ct:
        if any(c in _TRANSCT_EXCLUDE for c in col):
            return None
        return [0 if c == 'C' else (1 if c == 'T' else None) for c in col]

    unambig = {c for c in col if c in _UNAMBIG_BASES}
    has_ag  = bool(unambig & {'A', 'G'})
    has_ct  = bool(unambig & {'C', 'T'})

    if has_ag and has_ct:
        return None
    if not has_ag and not has_ct:
        return None
    if has_ag:
        if any(c in _TRANSAG_EXCLUDE for c in col):
            return None
        return [0 if c == 'A' else (1 if c == 'G' else None) for c in col]
    else:
        if any(c in _TRANSCT_EXCLUDE for c in col):
            return None
        return [0 if c == 'C' else (1 if c == 'T' else None) for c in col]


def _canonicalize(raw):
    """
    Canonical form: first non-None value = 1.
    If first non-None = 0, flip all non-None bits.
    """
    first = next((v for v in raw if v is not None), None)
    if first is None:
        return list(raw)
    if first == 0:
        return [None if v is None else (1 - v) for v in raw]
    return list(raw)


def _coded_to_bitmask(coded):
    """Convert list of 0/1 (no None allowed) to integer bitmask."""
    bm = 0
    for v in coded:
        bm = (bm << 1) | v
    return bm


def _bitmask_to_coded(bm, n_tax):
    """Convert integer bitmask to list of 0/1, left=first taxon."""
    return [(bm >> (n_tax - 1 - i)) & 1 for i in range(n_tax)]

# ---------------------------------------------------------------------------
# Site processing
# ---------------------------------------------------------------------------

def process_alignment(order, seqs, mode, codon_positions, binary,
                      track_positions=False):
    """
    Process all alignment columns and return list of PI site records.
    Each record is a dict:
      'coded'        : list of 0/1/None per taxon (canonical)
      'n_missing'    : 0 or 1 (sites with >=2 missing are excluded)
      'bipartitions' : list of (bitmask_int, weight) tuples

    If track_positions=True, each record also includes:
      'site_index'   : 0-based position in the alignment
      'raw_col'      : list of raw (uppercase) characters per taxon,
                       in `order`, before coding
    """
    seqs_ordered = [seqs[name] for name in order]
    aln_len = len(seqs_ordered[0])
    sites = []

    for i in range(aln_len):
        if codon_positions is not None:
            if ((i % 3) + 1) not in codon_positions:
                continue

        col = [seq[i] for seq in seqs_ordered]

        # Detect perfect modes (e.g. 'perfectRY' -> base 'RY', perfect=True)
        perfect_mode = mode.startswith('perfect')
        base_mode    = mode[len('perfect'):] if perfect_mode else mode

        if binary:
            coded = code_site_binary(col)
        elif base_mode in ('transitions', 'transitionsAG', 'transitionsCT'):
            coded = code_site_transitions(
                col,
                require_ag=(base_mode == 'transitionsAG'),
                require_ct=(base_mode == 'transitionsCT'),
            )
            if coded is None:
                continue
        elif base_mode == 'partimatrix':
            coded = code_site_partimatrix(col)
        elif base_mode == '2state':
            coded = code_site_2state(col)
            if coded is None:
                continue
        else:
            g0, g1, iu0, iu1 = MODE_AXES[base_mode]
            # Perfect check: any ambiguity/gap disqualifies
            if perfect_mode and not is_perfect_site(col, g0, g1):
                continue
            coded = code_site_axis(col, g0, g1, iu0, iu1)

        n_missing = sum(1 for v in coded if v is None)
        if n_missing >= 2:
            continue

        n0 = sum(1 for v in coded if v == 0)
        n1 = sum(1 for v in coded if v == 1)
        if not (n0 >= 2 and n1 >= 2):
            continue

        if n_missing == 0:
            bipartitions = [(_coded_to_bitmask(coded), 1.0)]
        else:
            miss_idx = next(k for k, v in enumerate(coded) if v is None)
            bipartitions = []
            for fill in (0, 1):
                c2 = coded[:]
                c2[miss_idx] = fill
                c2 = _canonicalize(c2)
                if sum(1 for v in c2 if v == 0) >= 2 and sum(1 for v in c2 if v == 1) >= 2:
                    bipartitions.append((_coded_to_bitmask(c2), 0.5))
            if not bipartitions:
                continue

        record = {
            'coded':        coded,
            'n_missing':    n_missing,
            'bipartitions': bipartitions,
        }
        if track_positions:
            record['site_index'] = i
            record['raw_col']    = col
        sites.append(record)

    return sites

# ---------------------------------------------------------------------------
# Compatibility calculation
# ---------------------------------------------------------------------------

def sites_compatible(coded_a, coded_b):
    """
    Four-gamete test. Returns True if fewer than all 4 gametes {00,01,10,11}
    are observed across non-missing position pairs.
    """
    gametes = set()
    for a, b in zip(coded_a, coded_b):
        if a is None or b is None:
            continue
        gametes.add((a, b))
        if len(gametes) == 4:
            return False
    return True


def compute_compatibility_score(sites):
    """Fraction of compatible pairs among all C(n_PI, 2) site pairs."""
    n = len(sites)
    if n < 2:
        return None
    total  = n * (n - 1) // 2
    compat = 0
    for i in range(n):
        for j in range(i + 1, n):
            if sites_compatible(sites[i]['coded'], sites[j]['coded']):
                compat += 1
    return compat / total

# ---------------------------------------------------------------------------
# Partition support / conflict
# ---------------------------------------------------------------------------

def compute_partition_support(sites, n_tax):
    """
    For each observed bipartition, compute:
      support     : sum of weights of sites whose bitmask == this bipartition
      conflict_raw: sum of weights of sites four-gamete-incompatible with this
                    bipartition's induced binary pattern
    """
    support = defaultdict(float)
    for site in sites:
        for bm, weight in site['bipartitions']:
            support[bm] += weight

    part_coded   = {bm: _bitmask_to_coded(bm, n_tax) for bm in support}
    conflict_raw = defaultdict(float)

    for site in sites:
        for bm_site, weight in site['bipartitions']:
            sc = _bitmask_to_coded(bm_site, n_tax)
            for bm_part, pc in part_coded.items():
                if bm_site == bm_part:
                    continue
                if not sites_compatible(sc, pc):
                    conflict_raw[bm_part] += weight

    return support, conflict_raw


def normalize_lento(support, conflict_raw):
    """
    Lento normalization: scale conflict so sum(conflict_norm) == sum(support).
    support_norm is identical to support (returned as a plain dict for clarity).
    """
    total_support  = sum(support.values())
    total_conflict = sum(conflict_raw.values())
    if total_conflict == 0:
        conflict_norm = {bm: 0.0 for bm in conflict_raw}
    else:
        scale = total_support / total_conflict
        conflict_norm = {bm: v * scale for bm, v in conflict_raw.items()}
    return dict(support), conflict_norm

# ---------------------------------------------------------------------------
# Bipartition formatting helpers
# ---------------------------------------------------------------------------

def bitmask_to_binary_str(bm, n_tax):
    """Left-to-right binary string; first taxon leftmost. Always starts with 1."""
    return format(bm, f'0{n_tax}b')


def bitmask_to_jakobsen_hex(bm, n_tax):
    """
    Jakobsen hex = canonical bitmask XOR full_mask, formatted as uppercase hex.
    Inverting our canonical form recovers the original program's convention.
    """
    full_mask = (1 << n_tax) - 1
    return format(bm ^ full_mask, 'X')


def sort_key(bm, support_norm, conflict_norm):
    """Descending sort key: support_norm - conflict_norm."""
    return support_norm.get(bm, 0.0) - conflict_norm.get(bm, 0.0)

def build_column_order(support, support_norm, conflict_norm,
                       plot_bitmasks=None, order_bitmasks=None):
    """
    Determine the final left-to-right column order for TSV and SVG.

    1. Start with the set of bipartitions to display:
       - if plot_bitmasks supplied: restrict to that set
       - otherwise: all keys in support

    2. Order them:
       - bipartitions present in order_bitmasks come first, in file order,
         skipping any not in the display set
       - remaining display bipartitions follow, sorted by
         support_norm - conflict_norm descending
    """
    if plot_bitmasks is not None:
        display = set(support.keys()) & plot_bitmasks
        # zero-support plot-biparts entries are already in support (added earlier)
        display |= (plot_bitmasks & set(support.keys()))
    else:
        display = set(support.keys())

    if order_bitmasks:
        ordered_head = [bm for bm in order_bitmasks if bm in display]
        ordered_tail = sorted(
            display - set(ordered_head),
            key=lambda bm: sort_key(bm, support_norm, conflict_norm),
            reverse=True,
        )
        return ordered_head + ordered_tail
    else:
        return sorted(
            display,
            key=lambda bm: sort_key(bm, support_norm, conflict_norm),
            reverse=True,
        )


def _write_tree_bipart_tsv(path, tree_bitmasks, support, conflict_raw,
                           support_norm, conflict_norm, n_tax, order, pos_map):
    """
    Write a TSV restricted to tree bipartitions, sorted by support_norm -
    conflict_norm descending (same ordering as the main TSV).  Includes a
    plot_position column showing the 1-based column index in the full Lento
    plot (empty if the bipartition is absent from the plot).
    """
    # Include all tree bipartitions: those in data (with support) and
    # zero-support ones added via --plot-biparts (support == 0.0)
    tree_bms = sorted(
        tree_bitmasks,
        key=lambda bm: sort_key(bm, support_norm, conflict_norm),
        reverse=True,
    )
    with open(path, 'w') as f:
        f.write('\t'.join([
            'bipartition_binary', 'jakobsen_hex',
            'support_raw', 'conflict_raw',
            'support_norm', 'conflict_norm',
            'plot_position',
            'grouping',
        ]) + '\n')
        for bm in tree_bms:
            bin_str  = bitmask_to_binary_str(bm, n_tax)
            hex_str  = bitmask_to_jakobsen_hex(bm, n_tax)
            sup_r    = support.get(bm, 0.0)
            con_r    = conflict_raw.get(bm, 0.0)
            sup_n    = support_norm.get(bm, 0.0)
            con_n    = conflict_norm.get(bm, 0.0)
            pos      = pos_map.get(bm, '')   # empty if not in plot
            group1   = [order[i] for i, c in enumerate(bin_str) if c == '1']
            group0   = [order[i] for i, c in enumerate(bin_str) if c == '0']
            grouping = f"({','.join(group1)})|({','.join(group0)})"
            f.write('\t'.join([
                "'" + bin_str, hex_str,
                f'{sup_r:.4f}', f'{con_r:.4f}',
                f'{sup_n:.4f}', f'{con_n:.4f}',
                str(pos),
                grouping,
            ]) + '\n')


# ---------------------------------------------------------------------------
# TSV output
# ---------------------------------------------------------------------------

def write_tsv(path, support, conflict_raw, support_norm, conflict_norm,
              n_tax, order, plot_bitmasks=None, order_bitmasks=None,
              sd_support=None, sd_conflict=None, sd_label='sd'):
    """
    Write bipartition TSV.

    Optional sd_support / sd_conflict dicts add per-bipartition variability
    columns (support_{sd_label} and conflict_{sd_label}).  These are raw
    (un-normalised) values, on the same scale as support_raw / conflict_raw.
    sd_label is 'sd' for standard deviation or 'mad' for scaled MAD.
    """
    bipartitions = build_column_order(
        support, support_norm, conflict_norm,
        plot_bitmasks=plot_bitmasks,
        order_bitmasks=order_bitmasks,
    )
    have_sd = sd_support is not None and sd_conflict is not None
    header = [
        'bipartition_binary', 'jakobsen_hex',
        'support_raw', 'conflict_raw',
        'support_norm', 'conflict_norm',
    ]
    if have_sd:
        header += [f'support_{sd_label}', f'conflict_{sd_label}']
    header += ['grouping']

    with open(path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for bm in bipartitions:
            bin_str  = bitmask_to_binary_str(bm, n_tax)
            hex_str  = bitmask_to_jakobsen_hex(bm, n_tax)
            sup_r    = support.get(bm, 0.0)
            con_r    = conflict_raw.get(bm, 0.0)
            sup_n    = support_norm.get(bm, 0.0)
            con_n    = conflict_norm.get(bm, 0.0)
            group1   = [order[i] for i, c in enumerate(bin_str) if c == '1']
            group0   = [order[i] for i, c in enumerate(bin_str) if c == '0']
            grouping = f"({','.join(group1)})|({','.join(group0)})"
            row = [
                "'" + bin_str, hex_str,
                f'{sup_r:.4f}', f'{con_r:.4f}',
                f'{sup_n:.4f}', f'{con_n:.4f}',
            ]
            if have_sd:
                row += [
                    f'{sd_support.get(bm, 0.0):.4f}',
                    f'{sd_conflict.get(bm, 0.0):.4f}',
                ]
            row += [grouping]
            f.write('\t'.join(row) + '\n')

# ---------------------------------------------------------------------------
# SVG LentoPlot
# ---------------------------------------------------------------------------

def write_svg(path, support, conflict_raw, support_norm, conflict_norm,
              n_tax, order, tree_bitmasks=None, plot_bitmasks=None,
              order_bitmasks=None, show_hex=True, max_cols=None):
    """
    LentoPlot SVG.

    Column order is determined by build_column_order():
      - order_bitmasks (from --bipart-order) pins leading columns in file order
      - remaining columns follow, sorted by support_norm - conflict_norm desc
    Bipartitions present in tree_bitmasks are drawn with outline bars.
    If plot_bitmasks is supplied, only those bipartitions are plotted.
    Jakobsen hex labels are drawn above the plot unless show_hex=False.
    If max_cols is supplied, only the first max_cols columns are plotted.
    """
    if tree_bitmasks is None:
        tree_bitmasks = set()

    bipartitions = build_column_order(
        support, support_norm, conflict_norm,
        plot_bitmasks=plot_bitmasks,
        order_bitmasks=order_bitmasks,
    )

    # Truncate to max_cols if requested
    if max_cols is not None and max_cols < len(bipartitions):
        bipartitions = bipartitions[:max_cols]

    n_parts = len(bipartitions)
    if n_parts == 0:
        print("WARNING: No bipartitions to plot.", file=sys.stderr)
        return

    # Layout constants
    BAR_W           = 12
    BAR_GAP         = 2
    COL_W           = BAR_W + BAR_GAP
    MAX_BAR_H       = 120
    SCHEMATIC_ROW_H = 10
    SCHEMATIC_GAP   = 8
    HEX_FONT        = 7
    HEX_GAP         = 6     # gap between bottom of hex labels and top of bars
    MARGIN_BASE     = 10    # blank margin above hex labels (or above bars if no hex)
    MARGIN_LEFT     = 120
    MARGIN_RIGHT    = 20
    AXIS_H          = 2
    LABEL_FONT      = 10
    TITLE_FONT      = 11

    # Hex labels sit ABOVE the Lento plot, reading upward (rotate -90).
    # text-anchor="end" anchors the bottom of the upward-reading text at
    # hex_anchor_y = top of support bars; labels extend into the top margin.
    if show_hex:
        max_hex_len = max(len(bitmask_to_jakobsen_hex(bm, n_tax))
                         for bm in bipartitions)
        hex_label_h = max_hex_len * HEX_FONT
    else:
        hex_label_h = 0

    # Top margin = blank gap above + hex label height + gap below labels
    MARGIN_TOP  = MARGIN_BASE + hex_label_h + (HEX_GAP if show_hex else 0)

    schematic_h = n_tax * SCHEMATIC_ROW_H
    total_w     = MARGIN_LEFT + n_parts * COL_W + MARGIN_RIGHT
    total_h     = (MARGIN_TOP
                   + MAX_BAR_H + AXIS_H + MAX_BAR_H
                   + SCHEMATIC_GAP + schematic_h
                   + MARGIN_BASE)

    max_val = max(
        max((support_norm.get(bm, 0) for bm in bipartitions), default=0.0),
        max((conflict_norm.get(bm, 0) for bm in bipartitions), default=0.0),
        0.001,
    )

    def scale(v):
        return (v / max_val) * MAX_BAR_H

    axis_y       = MARGIN_TOP + MAX_BAR_H
    hex_anchor_y = MARGIN_TOP   # bottom of hex label text aligns with top of bars
    schematic_y0 = axis_y + AXIS_H + MAX_BAR_H + SCHEMATIC_GAP
    dot_r        = SCHEMATIC_ROW_H * 0.35

    svg = []
    svg.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{total_w}" height="{total_h}" '
        f'font-family="Helvetica,Arial,sans-serif">'
    )
    svg.append(
        f'<rect width="{total_w}" height="{total_h}" '
        f'fill="white" stroke="none"/>'
    )

    # Axis line
    svg.append(
        f'<line x1="{MARGIN_LEFT}" y1="{axis_y}" '
        f'x2="{MARGIN_LEFT + n_parts * COL_W}" y2="{axis_y}" '
        f'stroke="#333" stroke-width="1.5"/>'
    )

    # Bars
    for i, bm in enumerate(bipartitions):
        x       = MARGIN_LEFT + i * COL_W
        sup_h   = scale(support_norm.get(bm, 0.0))
        con_h   = scale(conflict_norm.get(bm, 0.0))
        in_tree = bm in tree_bitmasks

        if in_tree:
            # Outline only for tree bipartitions
            if sup_h > 0:
                svg.append(
                    f'<rect x="{x}" y="{axis_y - sup_h:.1f}" '
                    f'width="{BAR_W}" height="{sup_h:.1f}" '
                    f'fill="none" stroke="#2e7d32" stroke-width="1.2"/>'
                )
            if con_h > 0:
                svg.append(
                    f'<rect x="{x}" y="{axis_y + AXIS_H}" '
                    f'width="{BAR_W}" height="{con_h:.1f}" '
                    f'fill="none" stroke="#c62828" stroke-width="1.2"/>'
                )
        else:
            # Solid fill for all other bipartitions
            if sup_h > 0:
                svg.append(
                    f'<rect x="{x}" y="{axis_y - sup_h:.1f}" '
                    f'width="{BAR_W}" height="{sup_h:.1f}" '
                    f'fill="#2e7d32" stroke="none"/>'
                )
            if con_h > 0:
                svg.append(
                    f'<rect x="{x}" y="{axis_y + AXIS_H}" '
                    f'width="{BAR_W}" height="{con_h:.1f}" '
                    f'fill="#c62828" stroke="none"/>'
                )

    # Taxon labels (left margin)
    for t_idx, name in enumerate(order):
        cy = schematic_y0 + t_idx * SCHEMATIC_ROW_H + SCHEMATIC_ROW_H // 2
        svg.append(
            f'<text x="{MARGIN_LEFT - 4}" y="{cy + 3}" '
            f'text-anchor="end" font-size="{LABEL_FONT}" fill="#222">'
            f'{_xml_escape(name)}</text>'
        )

    # Dot schematic
    for i, bm in enumerate(bipartitions):
        bin_str = bitmask_to_binary_str(bm, n_tax)
        cx = MARGIN_LEFT + i * COL_W + BAR_W // 2
        for t_idx, bit in enumerate(bin_str):
            cy   = schematic_y0 + t_idx * SCHEMATIC_ROW_H + SCHEMATIC_ROW_H // 2
            fill = '#222' if bit == '1' else 'white'
            svg.append(
                f'<circle cx="{cx}" cy="{cy}" r="{dot_r:.1f}" '
                f'fill="{fill}" stroke="#555" stroke-width="0.8"/>'
            )

    # Hex labels: rotated -90deg above the Lento plot, reading upward.
    # text-anchor="end" means the bottom of the upward-reading label sits at
    # hex_anchor_y (= top of support bars); text extends upward into top margin.
    if show_hex:
        for i, bm in enumerate(bipartitions):
            hex_str = bitmask_to_jakobsen_hex(bm, n_tax)
            cx = MARGIN_LEFT + i * COL_W + BAR_W // 2
            svg.append(
                f'<text transform="rotate(-90,{cx},{hex_anchor_y})" '
                f'x="{cx}" y="{hex_anchor_y}" '
                f'text-anchor="end" font-size="{HEX_FONT}" fill="#444">'
                f'{_xml_escape(hex_str)}</text>'
            )

    # Y-axis labels (rotated, centred on the support and conflict bar regions)
    Y_LABEL_X = MARGIN_LEFT - 60
    for label, my in [('Support',  axis_y - MAX_BAR_H // 2),
                      ('Conflict', axis_y + AXIS_H + MAX_BAR_H // 2)]:
        svg.append(
            f'<text transform="rotate(-90,{Y_LABEL_X},{my})" '
            f'x="{Y_LABEL_X}" y="{my + 4}" '
            f'text-anchor="middle" font-size="{TITLE_FONT}" fill="#333">'
            f'{label}</text>'
        )

    # Y-axis scale ticks (at 50% and 100% of max)
    for frac in [0.5, 1.0]:
        val   = max_val * frac
        ty_up = axis_y - scale(val)
        ty_dn = axis_y + AXIS_H + scale(val)
        for ty in [ty_up, ty_dn]:
            svg.append(
                f'<line x1="{MARGIN_LEFT - 4}" y1="{ty:.1f}" '
                f'x2="{MARGIN_LEFT}" y2="{ty:.1f}" '
                f'stroke="#666" stroke-width="1"/>'
            )
            svg.append(
                f'<text x="{MARGIN_LEFT - 6}" y="{ty + 3:.1f}" '
                f'text-anchor="end" font-size="8" fill="#555">'
                f'{val:.2f}</text>'
            )

    svg.append('</svg>')

    with open(path, 'w') as f:
        f.write('\n'.join(svg))


def _xml_escape(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Compatibility and partition analysis of aligned sequences.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Modes (--mode, case-insensitive; reversed pairs accepted e.g. YR == RY):
  RY / YR        purines (A,G)=1 vs pyrimidines (C,T)=0  [default]
  KM / MK        amino (A,C)=1 vs keto (G,T)=0
  WS / SW        weak (A,T)=1 vs strong (C,G)=0
  perfectRY / perfectYR
                 like RY but requires exactly 2 unambiguous nucleotides
                 (one purine, one pyrimidine) with no gaps or IUPAC codes
                 (perfect transversions; Tiley et al. 2020)
  perfectKM / perfectMK
                 like KM but requires exactly 2 unambiguous nucleotides
                 (one amino, one keto) with no gaps or IUPAC codes
  perfectWS / perfectSW
                 like WS but requires exactly 2 unambiguous nucleotides
                 (one weak, one strong) with no gaps or IUPAC codes
  2state         only sites with exactly 2 distinct unambiguous nucleotides
  partimatrix    historical behaviour: 2-state columns kept as binary;
                 3/4-state columns RY-collapsed

Tree bipartitions file (-t / --tree-biparts):
  Plain text, one bipartition per line; blank lines and # comments ignored.
  Each line is a binary string of length n_tax (e.g. 1100010)
  or a Jakobsen hex string (e.g. 6, 0x6, 6A).
  Matching bipartitions are drawn as outline bars in the LentoPlot SVG.
        """,
    )
    p.add_argument('-i', '--infile',  required=True,
                   help='Input alignment (FASTA, relaxed PHYLIP, or Nexus)')
    p.add_argument('-o', '--outfile', required=True,
                   help='Output file prefix (produces PREFIX.tsv and PREFIX.svg)')
    p.add_argument('-m', '--mode',    default='RY',
                   help='Site coding mode (default: RY). '
                        'Perfect variants (perfectRY, perfectKM, perfectWS) '
                        'require exactly 2 unambiguous nucleotides, one per '
                        'side of axis, no gaps or IUPAC (Tiley et al. 2020). '
                        'Transition variants (transitions, transitionsAG, '
                        'transitionsCT) restrict to transition-only sites; '
                        'cross-axis IUPAC codes exclude the site; within-axis '
                        'IUPAC (R for AG, Y for CT) is treated as missing.')
    p.add_argument('-t', '--tree-biparts',
                   help='Optional file of tree bipartitions (binary or Jakobsen hex)')
    p.add_argument('--tax-order',
                   help='File with taxon names one per line, '
                        'setting the bit-string ordering')
    p.add_argument('--codon-positions', default='all',
                   choices=['all', '1', '2', '3', '12', '13', '23', '123'],
                   help='Codon positions to include (default: all)')
    p.add_argument('-p', '--plot-biparts',
                   help='Optional file of bipartitions to plot (binary or Jakobsen hex); '
                        'same format as --tree-biparts. Only these bipartitions appear '
                        'in the SVG; all computation and the TSV are unaffected. '
                        'Compatible with --tree-biparts: a bipartition in both files '
                        'will be plotted with an outline bar.')
    p.add_argument('--no-svg', action='store_true',
                   help='Suppress SVG output (useful for batch/replicate runs)')
    p.add_argument('--score-file', action='store_true',
                   help='Write compatibility score to OUTFILE.compat.txt')
    p.add_argument('-b', '--bipart-order',
                   help='Optional file specifying bipartition column order '
                        '(binary or Jakobsen hex, one per line, no duplicates). '
                        'Listed bipartitions are plotted first (left) in file '
                        'order; remaining bipartitions follow sorted by '
                        'support_norm - conflict_norm. Bipartitions in this '
                        'file absent from the data trigger a warning.')
    p.add_argument('--truncate', type=int, default=None,
                   metavar='N',
                   help='Also write PREFIX.trunc.svg showing only the first N '
                        'columns of the Lento plot (same ordering as the full '
                        'plot; full PREFIX.svg is unchanged)')
    p.add_argument('--no-hex', action='store_true',
                   help='Suppress Jakobsen hex labels above the Lento plot')
    return p.parse_args()


def resolve_mode(mode_str):
    m = mode_str.upper()
    if m in ('RY', 'YR'):               return 'RY'
    if m in ('KM', 'MK'):               return 'KM'
    if m in ('WS', 'SW'):               return 'WS'
    if m in ('PERFECTRY', 'PERFECTYR'): return 'perfectRY'
    if m in ('PERFECTKM', 'PERFECTMK'): return 'perfectKM'
    if m in ('PERFECTWS', 'PERFECTSW'): return 'perfectWS'
    if m in ('TRANSITIONS', 'TRANSITION'): return 'transitions'
    if m in ('TRANSITIONSAG','TRANSITIONAG','TRANSITIONSGA','TRANSITIONGA'): return 'transitionsAG'
    if m in ('TRANSITIONSCT','TRANSITIONCT','TRANSITIONSTC','TRANSITIONTC'): return 'transitionsCT'
    if m == '2STATE':                   return '2state'
    if m == 'PARTIMATRIX':              return 'partimatrix'
    sys.exit(
        f"ERROR: Unknown mode '{mode_str}'. "
        f"Choose from: RY, KM, WS, perfectRY, perfectKM, perfectWS, "
        f"transitions, transitionsAG, transitionsCT, 2state, partimatrix"
    )


def main():
    args  = parse_args()
    mode  = resolve_mode(args.mode)

    # Read alignment
    order, seqs = read_alignment(args.infile)
    n_tax = len(order)
    print(f"Read {n_tax} sequences, length {len(seqs[order[0]])}", file=sys.stderr)

    # Taxon reordering
    if args.tax_order:
        requested = read_tax_order(args.tax_order)
        missing   = [n for n in requested if n not in seqs]
        if missing:
            sys.exit(f"ERROR: Taxa in --tax-order not in alignment: {missing}")
        extra = [n for n in order if n not in set(requested)]
        if extra:
            print(f"WARNING: taxa not in --tax-order, appending at end: {extra}",
                  file=sys.stderr)
        order = [n for n in requested if n in seqs] + extra
        n_tax = len(order)

    # Binary detection
    binary = is_binary_data(seqs)
    if binary:
        print("Detected binary data; --mode setting ignored.", file=sys.stderr)

    # Codon positions
    codon_positions = CODON_POS_MAP[args.codon_positions]

    # Process alignment
    sites = process_alignment(order, seqs, mode, codon_positions, binary)
    n_pi  = len(sites)
    print(f"Parsimony-informative sites: {n_pi}", file=sys.stderr)
    if n_pi == 0:
        sys.exit("ERROR: No parsimony-informative sites found.")

    # Compatibility score
    print("Computing compatibility score...", file=sys.stderr)
    compat_score = compute_compatibility_score(sites)
    total_pairs  = n_pi * (n_pi - 1) // 2
    compat_pairs = round(compat_score * total_pairs)
    print(
        f"\nCompatibility score: {compat_score:.6f}  "
        f"({compat_pairs} compatible pairs / {total_pairs} total pairs; "
        f"{n_pi} PI sites)"
    )

    # Partition support / conflict (observed bipartitions only)
    print("Computing partition support and conflict...", file=sys.stderr)
    support, conflict_raw = compute_partition_support(sites, n_tax)
    print(f"Distinct bipartitions observed in data: {len(support)}", file=sys.stderr)

    # Tree bipartitions (controls outline-bar styling; does not affect support/conflict)
    tree_bitmasks = None
    if args.tree_biparts:
        tree_bitmasks = read_tree_biparts(args.tree_biparts, n_tax)
        absent  = tree_bitmasks - set(support.keys())
        present = tree_bitmasks &  set(support.keys())
        if absent:
            print(
                f"WARNING: {len(absent)} --tree-biparts bipartition(s) absent from "
                f"data (no support; will not appear as outline bars):",
                file=sys.stderr,
            )
            for bm in sorted(absent, key=lambda b: bitmask_to_binary_str(b, n_tax)):
                print(
                    f"  {bitmask_to_binary_str(bm, n_tax)}, "
                    f"hex {bitmask_to_jakobsen_hex(bm, n_tax)}",
                    file=sys.stderr,
                )
        print(
            f"Tree bipartitions: {len(tree_bitmasks)} supplied, "
            f"{len(present)} present in data",
            file=sys.stderr,
        )

    # Plot bipartition filter.
    # Zero-support bipartitions in --plot-biparts are fully included: their conflict
    # is computed from the data, they appear in both the TSV and SVG (sorted to the
    # right because sort_key = 0 - conflict_norm <= 0).
    plot_bitmasks = None
    if args.plot_biparts:
        plot_bitmasks = read_tree_biparts(args.plot_biparts, n_tax)
        absent  = plot_bitmasks - set(support.keys())
        present = plot_bitmasks &  set(support.keys())
        if absent:
            print(
                f"Plot filter: {len(absent)} bipartition(s) with zero support in "
                f"data; computing conflict-only entries...",
                file=sys.stderr,
            )
            for bm in absent:
                pc           = _bitmask_to_coded(bm, n_tax)
                raw_conflict = 0.0
                for site in sites:
                    for bm_site, weight in site['bipartitions']:
                        if not sites_compatible(_bitmask_to_coded(bm_site, n_tax), pc):
                            raw_conflict += weight
                support[bm]      = 0.0
                conflict_raw[bm] = raw_conflict
        n_conflict_only = len(absent)
        print(
            f"Plot filter: {len(plot_bitmasks)} supplied, "
            f"{len(present)} with support, "
            f"{n_conflict_only} conflict-only "
            f"(SVG will show {len(plot_bitmasks)} columns)",
            file=sys.stderr,
        )

    # Lento normalization (runs after any zero-support bipartitions are added)
    support_norm, conflict_norm = normalize_lento(support, conflict_raw)

    # Bipartition count reporting (after support is fully populated)
    n_zero_plot = 0
    if args.plot_biparts and plot_bitmasks:
        n_zero_plot = sum(1 for bm in plot_bitmasks
                         if support.get(bm, -1) == 0.0)
    n_supported = len(support) - n_zero_plot
    n_possible  = 2 ** (n_tax - 1) - n_tax - 1
    pct         = 100.0 * n_supported / n_possible if n_possible > 0 else 0.0
    print(f"Bipartitions with support: {n_supported} of {n_possible} possible "
          f"({pct:.6f}%)")
    if n_zero_plot > 0:
        print(f"Bipartitions added via --plot-biparts with zero support: "
              f"{n_zero_plot}")

    # Bipartition ordering
    order_bitmasks = None
    if args.bipart_order:
        order_bitmasks = read_bipart_order(args.bipart_order, n_tax)
        # Warn about entries absent from data (not in support)
        absent_order = [bm for bm in order_bitmasks if bm not in support]
        if absent_order:
            absent_strs = [bitmask_to_binary_str(bm, n_tax)
                           for bm in absent_order]
            print(
                f"WARNING: {len(absent_order)} --bipart-order bipartition(s) "
                f"absent from data (will be skipped in column ordering): "
                f"{', '.join(absent_strs)}",
                file=sys.stderr,
            )
        present_order = [bm for bm in order_bitmasks if bm in support]
        print(
            f"Bipart order: {len(order_bitmasks)} supplied, "
            f"{len(present_order)} present in data "
            f"(these will be the first {len(present_order)} columns)",
            file=sys.stderr,
        )

    # Write outputs
    tsv_path   = args.outfile + '.tsv'
    svg_path   = args.outfile + '.svg'
    score_path = args.outfile + '.compat.txt'

    write_tsv(tsv_path, support, conflict_raw, support_norm, conflict_norm,
              n_tax, order,
              plot_bitmasks=plot_bitmasks,
              order_bitmasks=order_bitmasks)
    print(f"Wrote {tsv_path}", file=sys.stderr)

    # Tree-bipartition TSV and position reporting
    if tree_bitmasks:
        # Build the full ordered column list (same as SVG; plot_bitmasks applied)
        full_order = build_column_order(
            support, support_norm, conflict_norm,
            plot_bitmasks=plot_bitmasks,
            order_bitmasks=order_bitmasks,
        )
        # Map each tree bipartition to its 1-based position in the full plot
        pos_map = {bm: i + 1 for i, bm in enumerate(full_order)}
        tree_in_plot = {bm: pos_map[bm] for bm in tree_bitmasks if bm in pos_map}
        if tree_in_plot:
            max_pos = max(tree_in_plot.values())
            max_bm  = next(bm for bm, p in tree_in_plot.items() if p == max_pos)
            # Only report truncated count when --truncate is active;
            # without truncation this would just repeat "N present in data"
            if args.truncate is not None:
                n_in_trunc = sum(1 for p in tree_in_plot.values()
                                 if p <= args.truncate)
                print(
                    f"Tree bipartitions in truncated plot "
                    f"(first {args.truncate} columns): {n_in_trunc}"
                )
            print(
                f"Column of last tree bipartition = {max_pos} "
                f"({bitmask_to_binary_str(max_bm, n_tax)}, "
                f"hex {bitmask_to_jakobsen_hex(max_bm, n_tax)})"
            )
        # Write PREFIX.tree_bipart.tsv
        tbp_path = args.outfile + '.tree_bipart.tsv'
        _write_tree_bipart_tsv(
            tbp_path, tree_bitmasks, support, conflict_raw,
            support_norm, conflict_norm, n_tax, order, pos_map,
        )
        print(f"Wrote {tbp_path}", file=sys.stderr)

    if args.score_file:
        with open(score_path, 'w') as _sf:
            _sf.write(f"{compat_score:.6f}\n")
            _sf.write(f"{n_supported}\n")
        print(f"Wrote {score_path}", file=sys.stderr)

    if not args.no_svg:
        write_svg(svg_path, support, conflict_raw, support_norm, conflict_norm,
                  n_tax, order,
                  tree_bitmasks=tree_bitmasks,
                  plot_bitmasks=plot_bitmasks,
                  order_bitmasks=order_bitmasks,
                  show_hex=(not args.no_hex))
        print(f"Wrote {svg_path}", file=sys.stderr)

        if args.truncate is not None:
            trunc_path = args.outfile + '.trunc.svg'
            write_svg(trunc_path, support, conflict_raw, support_norm, conflict_norm,
                      n_tax, order,
                      tree_bitmasks=tree_bitmasks,
                      plot_bitmasks=plot_bitmasks,
                      order_bitmasks=order_bitmasks,
                      show_hex=(not args.no_hex),
                      max_cols=args.truncate)
            print(f"Wrote {trunc_path}", file=sys.stderr)


if __name__ == '__main__':
    main()
