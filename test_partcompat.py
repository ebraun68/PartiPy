#!/usr/bin/env python3
"""
Tests for partcompat.py

Run with: python3 test_partcompat.py
"""

import sys
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from partcompat import (
    _canonicalize,
    _coded_to_bitmask,
    bitmask_to_binary_str,
    bitmask_to_jakobsen_hex,
    code_site_partimatrix,
    code_site_axis,
    code_site_2state,
    code_site_binary,
    sites_compatible,
    compute_compatibility_score,
    MODE_AXES,
    parse_fasta,
    parse_phylip,
    process_alignment,
    compute_partition_support,
    normalize_lento,
)

PASS = 0
FAIL = 0

def check(name, got, expected):
    global PASS, FAIL
    if got == expected:
        print(f"  PASS  {name}")
        PASS += 1
    else:
        print(f"  FAIL  {name}")
        print(f"        got:      {got}")
        print(f"        expected: {expected}")
        FAIL += 1

def check_close(name, got, expected, tol=1e-9):
    global PASS, FAIL
    if got is None and expected is None:
        print(f"  PASS  {name}")
        PASS += 1
        return
    if got is not None and expected is not None and abs(got - expected) <= tol:
        print(f"  PASS  {name}")
        PASS += 1
    else:
        print(f"  FAIL  {name}")
        print(f"        got:      {got}")
        print(f"        expected: {expected}")
        FAIL += 1

# ---------------------------------------------------------------------------
print("\n=== Canonicalization ===")

# First non-None is 0 → flip
check("canon flip 0->1", _canonicalize([0, 0, 1, 1]), [1, 1, 0, 0])
# First non-None is 1 → no flip
check("canon no-flip", _canonicalize([1, 1, 0, 0]), [1, 1, 0, 0])
# With None values
check("canon with None", _canonicalize([0, None, 1, 0]), [1, None, 0, 1])
# First non-None is None: all None
check("canon all None", _canonicalize([None, None]), [None, None])
# First is None then 0 → flip
check("canon None then 0", _canonicalize([None, 0, 1]), [None, 1, 0])

# ---------------------------------------------------------------------------
print("\n=== Bitmask conversion ===")

# 4 taxa: coded [1,1,0,0] → binary 1100 = 12
check("coded_to_bitmask 1100", _coded_to_bitmask([1,1,0,0]), 0b1100)
# 4 taxa: [1,0,0,1] → 1001 = 9
check("coded_to_bitmask 1001", _coded_to_bitmask([1,0,0,1]), 0b1001)

# bitmask_to_binary_str
check("bitmask 12 n=4", bitmask_to_binary_str(0b1100, 4), '1100')
check("bitmask 9 n=4",  bitmask_to_binary_str(0b1001, 4), '1001')

# Jakobsen hex: flip all bits of canonical
# canonical 1100 (n=4), flip → 0011 = 3 → hex '3'
check("jakobsen hex 1100 n=4", bitmask_to_jakobsen_hex(0b1100, 4), '3')
# canonical 1001, flip → 0110 = 6 → hex '6'
check("jakobsen hex 1001 n=4", bitmask_to_jakobsen_hex(0b1001, 4), '6')

# ---------------------------------------------------------------------------
print("\n=== Bipartition canonicalization examples from spec ===")

# ABCDE|FGH with 8 taxa → groups: ABCDE=1, FGH=0
# canonical (first bit=1): 11111000
bm = _coded_to_bitmask([1,1,1,1,1,0,0,0])
check("ABCDE|FGH bitmask", bm, 0b11111000)
check("ABCDE|FGH binary str", bitmask_to_binary_str(bm, 8), '11111000')

# 00000111 has first bit=0 → must be flipped to 11111000
# i.e., site pattern FGHIJ|ABCDE encodes as 00000111
# after canonicalize → 11111000
raw = [0,0,0,0,0,1,1,1]
canon = _canonicalize(raw)
check("00000111 canon → 11111000", canon, [1,1,1,1,1,0,0,0])

# ---------------------------------------------------------------------------
print("\n=== code_site_partimatrix ===")

# 2-state: AAAAGG → coded binary
# A→group0, G→group1 (two-state rule: nucleotides in alphabetical order)
# bases = {A,G}, b0=A, b1=G
# raw = [0,0,0,0,1,1] → first is 0 → flip → [1,1,1,1,0,0]
col = list('AAAAGG')
coded = code_site_partimatrix(col)
check("partimatrix AAAAGG", coded, [1,1,1,1,0,0])

# AAAACC → same bipartition structure
col = list('AAAACC')
coded = code_site_partimatrix(col)
check("partimatrix AAAACC", coded, [1,1,1,1,0,0])

# Both map to same bipartition 1111|00
bm1 = _coded_to_bitmask(code_site_partimatrix(list('AAAAGG')))
bm2 = _coded_to_bitmask(code_site_partimatrix(list('AAAACC')))
check("partimatrix AAAAGG == AAAACC bitmask", bm1, bm2)

# AAGGCC → 3-state → RY collapse: A,A→0 G,G→0 C,C→1 → [0,0,0,0,1,1] → flip → [1,1,1,1,0,0]
col = list('AAGGCC')
coded = code_site_partimatrix(col)
check("partimatrix AAGGCC (3-state RY)", coded, [1,1,1,1,0,0])

# IUPAC R in partimatrix → treated as ?
col = list('AAAAGR')  # R is ?
coded = code_site_partimatrix(col)
# 5 non-missing: AAAAG → 2-state, A→0, G→1 → [0,0,0,0,1,?] → flip → [1,1,1,1,0,?]
check("partimatrix R→missing", coded, [1,1,1,1,0,None])

# ---------------------------------------------------------------------------
print("\n=== code_site_axis (RY mode) ===")

g0, g1, iu0, iu1 = MODE_AXES['RY']

# AAAAGG: A→0(purine), G→0(purine) → not PI (only one state)
col = list('AAAAGG')
coded = code_site_axis(col, g0, g1, iu0, iu1)
# All purines → all 0 after coding → not PI (handled upstream)
# But coded should be [0,0,0,0,0,0] → flip → [1,1,1,1,1,1]  (all same)
check("RY AAAAGG all purines", coded, [1,1,1,1,1,1])

# AAAACC: A→purine=0, C→pyrimidine=1 → [0,0,0,0,1,1] → flip → [1,1,1,1,0,0]
col = list('AAAACC')
coded = code_site_axis(col, g0, g1, iu0, iu1)
check("RY AAAACC", coded, [1,1,1,1,0,0])

# R→0 (purine), Y→1 (pyrimidine)
col = list('AAAAYY')
coded = code_site_axis(col, g0, g1, iu0, iu1)
check("RY AAAAYY (Y→1)", coded, [1,1,1,1,0,0])

col = list('AAAARY')  # R→purine=0, Y→pyrimidine=1
# raw=[0,0,0,0,0,1] → first=0 → flip → [1,1,1,1,1,0]
coded = code_site_axis(col, g0, g1, iu0, iu1)
coded_expected = _canonicalize([0,0,0,0,0,1])
check("RY AAAARY", coded, coded_expected)  # [1,1,1,1,1,0]

# S (CG) → ? under RY
col = list('AAAASS')
coded = code_site_axis(col, g0, g1, iu0, iu1)
check("RY AAAASS (S→?)", coded, [1,1,1,1,None,None])  # [0,0,0,0,?,?] → flip [1,1,1,1,?,?]

# ---------------------------------------------------------------------------
print("\n=== code_site_2state ===")

# AAAACC → exactly 2 states → binary
col = list('AAAACC')
coded = code_site_2state(col)
check("2state AAAACC", coded, [1,1,1,1,0,0])

# AAGGCC → 3 states → excluded (returns None)
col = list('AAGGCC')
coded = code_site_2state(col)
check("2state AAGGCC excluded", coded, None)

# AAAAGG → exactly 2 states
col = list('AAAAGG')
coded = code_site_2state(col)
# A,G → sorted [A,G], A→0, G→1 → [0,0,0,0,1,1] → flip → [1,1,1,1,0,0]
check("2state AAAAGG", coded, [1,1,1,1,0,0])

# ---------------------------------------------------------------------------
print("\n=== code_site_binary ===")

# Binary data
col = list('000011')
coded = code_site_binary(col)
# raw [0,0,0,0,1,1] → first=0 → flip → [1,1,1,1,0,0]
check("binary 000011", coded, [1,1,1,1,0,0])

col = list('110000')
coded = code_site_binary(col)
check("binary 110000", coded, [1,1,0,0,0,0])

col = list('0000?1')
coded = code_site_binary(col)
# [0,0,0,0,None,1] → flip → [1,1,1,1,None,0]
check("binary 0000?1", coded, [1,1,1,1,None,0])

# ---------------------------------------------------------------------------
print("\n=== Four-gamete compatibility test ===")

# Compatible: only 3 gametes
a = [1,1,0,0,1,0]
b = [1,1,0,0,1,0]
check("compat identical", sites_compatible(a, b), True)

# Compatible: 00, 01, 10 (no 11)
a = [1,1,0,0]
b = [1,0,1,0]
# pairs: (1,1),(1,0),(0,1),(0,0) → 4 gametes → INCOMPATIBLE
check("incompat all 4 gametes", sites_compatible(a, b), False)

# Missing data handled: skip None pairs
a = [1,1,0,None]
b = [1,0,1,0]
# pairs: (1,1),(1,0),(0,1) → 3 gametes → compatible
check("compat with None", sites_compatible(a, b), True)

# Classic incompatible 4-taxon example
a = [1,1,0,0]
b = [1,0,1,0]
check("classic 4-gamete incompat", sites_compatible(a, b), False)

# ---------------------------------------------------------------------------
print("\n=== Compatibility score ===")

# All 4 sites on same tree → fully compatible
sites_data = [
    {'coded': [1,1,0,0], 'bipartitions': [(_coded_to_bitmask([1,1,0,0]), 1.0)]},
    {'coded': [1,0,1,0], 'bipartitions': [(_coded_to_bitmask([1,0,1,0]), 1.0)]},
]
# (1100) vs (1010): gametes (1,1),(1,0),(0,1),(0,0) → incompatible
score = compute_compatibility_score(sites_data)
check_close("2 sites incompatible → score=0", score, 0.0)

sites_data2 = [
    {'coded': [1,1,0,0], 'bipartitions': [(_coded_to_bitmask([1,1,0,0]), 1.0)]},
    {'coded': [1,1,0,0], 'bipartitions': [(_coded_to_bitmask([1,1,0,0]), 1.0)]},
]
score2 = compute_compatibility_score(sites_data2)
check_close("2 identical sites → score=1", score2, 1.0)

# ---------------------------------------------------------------------------
print("\n=== Single-missing bipartition resolution ===")

# Site pattern: AAACC? (6 taxa)
# coded (post-coding, say RY or 2state): [1,1,1,0,0,None] 
# (first=1 so already canonical for non-missing part)
# Resolution A: missing→0: [1,1,1,0,0,0] → first=1, canonical → bitmask 111000
# Resolution B: missing→1: [1,1,1,0,0,1] → first=1, canonical → bitmask 111001

# We simulate this via process_alignment with known data
fasta_text = """>t1
AAACC-
>t2
AAACC-
>t3
AAACC-
>t4
CCCAA-
>t5
CCCAA-
>t6
ACGTA-
"""
# Actually let's use a cleaner example with explicit missing
fasta_text2 = """>t1
AAAC
>t2
AAAC
>t3
AAAC
>t4
CCCA
>t5
CCCA
>t6
CCN-
"""
# Column 0: AAACCC → after RY: all pyrimidines/purines... let's keep it simple
# Use binary data for the missing test
fasta_bin = """>t1
0001
>t2
0001
>t3
0001
>t4
1110
>t5
1110
>t6
111?
"""
from partcompat import parse_fasta, process_alignment
order, seqs = parse_fasta(fasta_bin)
sites = process_alignment(order, seqs, 'RY', None, binary=True)
# col 3: 1,1,1,0,0,? → n_missing=1
# coded after binary: raw=[0,0,0,1,1,?] ... wait let me think through
# raw col3: t1=1, t2=1, t3=1, t4=0, t5=0, t6=?
# first non-None = 1 → no flip needed
# Resolution A (t6=0): [1,1,1,0,0,0] → bitmask 111000
# Resolution B (t6=1): [1,1,1,0,0,1] → bitmask 111001
# Both PI (3 ones, 3 zeros; and 3 ones, 2 zeros 1 one...)
# Actually 111001 has three 1s and three 0s... wait: t1=1,t2=1,t3=1,t4=0,t5=0,t6=1
# n1=4, n0=2 → still PI

# check we get 2 bipartitions for the last site
last_site = sites[-1]
bitmasks = [bm for bm, w in last_site['bipartitions']]
weights   = [w  for bm, w in last_site['bipartitions']]

bm_111000 = _coded_to_bitmask([1,1,1,0,0,0])
bm_111001 = _coded_to_bitmask([1,1,1,0,0,1])

check("single-missing gives 2 bipartitions",
      len(last_site['bipartitions']), 2)
check("single-missing weights sum to 1.0",
      abs(sum(weights) - 1.0) < 1e-9, True)
bitmasks_set = set(bitmasks)
check("single-missing correct bipartitions",
      bitmasks_set, {bm_111000, bm_111001})

# ---------------------------------------------------------------------------
print("\n=== FASTA parsing ===")

fasta = """>seq1
ACGT
>seq2
TGCA
"""
order, seqs = parse_fasta(fasta)
check("fasta order", order, ['seq1', 'seq2'])
check("fasta seq1", seqs['seq1'], 'ACGT')
check("fasta seq2", seqs['seq2'], 'TGCA')

# ---------------------------------------------------------------------------
print("\n=== Full pipeline: small known example ===")

# 5 taxa, 4 sites, designed so we know which bipartitions should appear
# Site 1: 11100 → bipartition 11100
# Site 2: 11100 → same bipartition (full support)
# Site 3: 11000 → bipartition 11000
# Site 4: 10100 → bipartition 10100

fasta_known = """>t1
1111
>t2
1110
>t3
1101
>t4
0010
>t5
0000
"""
order5, seqs5 = parse_fasta(fasta_known)
sites5 = process_alignment(order5, seqs5, 'RY', None, binary=True)
# Check all sites are PI
check("known example: 4 PI sites", len(sites5), 4)

support5, conflict5 = compute_partition_support(sites5, 5)
# Site col0: t1=1,t2=1,t3=1,t4=0,t5=0 → coded=[1,1,1,0,0] → bitmask 11100=28
# Site col1: t1=1,t2=1,t3=1,t4=1,t5=0 ... wait let me re-read
# fasta_known col by col:
# col0: 1,1,1,0,0 → coded [1,1,1,0,0] canonical (first=1) → bitmask = 0b11100=28
# col1: 1,1,0,0,0 → coded [1,1,0,0,0] → bitmask 0b11000=24
# col2: 1,1,0,1,0 ... wait:
#   t1 col2=1, t2 col2=1, t3 col2=0, t4 col2=1, t5 col2=0
#   raw=[1,1,0,1,0] → first=1 → bitmask 0b11010=26
# col3: 1,0,1,0,0 → first=1 → bitmask 0b10100=20

bm_11100 = 0b11100  # 28
bm_11000 = 0b11000  # 24
bm_11010 = 0b11010  # 26
bm_10100 = 0b10100  # 20

check("known: bm 11100 has support", support5.get(bm_11100, 0) > 0, True)

# ---------------------------------------------------------------------------# ---------------------------------------------------------------------------
# label_tree.py tests
# ---------------------------------------------------------------------------
from label_tree import (
    parse_newick, get_leaves, find_leaf_path, reroot, label_nodes, newick_str
)

print("\n--- label_tree.py: Newick parsing ---")

nwk_simple = "((t1:0.1,t2:0.2):0.3,(t3:0.4,(t4:0.5,(t5:0.6,t6:0.7):0.8):0.9):1.0);"
tree_simple = parse_newick(nwk_simple)
check("parse: 6 leaves", len(get_leaves(tree_simple)), 6)
check("parse: root has 2 children", len(tree_simple['children']), 2)

nwk_nolen = "((t1,t2),(t3,(t4,(t5,t6))));"
tree_nolen = parse_newick(nwk_nolen)
check("parse no branch lengths: 6 leaves", len(get_leaves(tree_nolen)), 6)

nwk_trifurc = "((t1,t2),(t3,t4),(t5,t6));"
tree_tri = parse_newick(nwk_trifurc)
check("parse trifurcation: root has 3 children", len(tree_tri['children']), 3)

print("\n--- label_tree.py: re-rooting ---")

order6 = ['t1','t2','t3','t4','t5','t6']

# Root on t1: path = [root, (t1,t2), t1]; old root has 2 children → absorb
rerooted = reroot(parse_newick(nwk_simple), 't1')
leaves_r = sorted(get_leaves(rerooted))
check("reroot: all 6 leaves present", leaves_r, ['t1','t2','t3','t4','t5','t6'])
check("reroot: new root has 2 children", len(rerooted['children']), 2)
check("reroot: first child is outgroup t1", rerooted['children'][0]['name'], 't1')
check("reroot: t1 branch length preserved", rerooted['children'][0]['branch_length'], 0.1)

# Branch length preservation: total distance t1→t3 should be unchanged
# Original: t1(0.1) + (t1t2 node 0.3) + (t3t4t5t6 node 1.0) + t3(0.4) = 1.8
# After reroot+absorb: t1(0.1) + (1.3) + t3(0.4) = 1.8
subtree = rerooted['children'][1]   # the ingroup
# Find the branch to t3's containing clade
def find_bl(node, target, acc):
    if not node['children'] and node['name'] == target:
        return acc
    for child in node['children']:
        bl = child['branch_length'] or 0.0
        result = find_bl(child, target, acc + bl)
        if result is not None:
            return result
    return None

dist_t1_to_t3 = (rerooted['children'][0]['branch_length'] or 0.0) + \
                find_bl(rerooted['children'][1], 't3',
                        rerooted['children'][1]['branch_length'] or 0.0)
check("reroot: t1→t3 distance preserved", round(dist_t1_to_t3, 6), round(1.8, 6))

# Reroot when outgroup is direct child of root (n=2 path)
nwk_direct = "(t1:0.5,(t2:0.1,(t3:0.2,(t4:0.3,(t5:0.4,t6:0.5):0.6):0.7):0.8):0.9);"
rerooted_direct = reroot(parse_newick(nwk_direct), 't1')
check("reroot direct: 6 leaves", sorted(get_leaves(rerooted_direct)),
      ['t1','t2','t3','t4','t5','t6'])
check("reroot direct: first child is t1",
      rerooted_direct['children'][0]['name'], 't1')

# Trifurcating root: old root has 3 children, 1 removed → 2 remain → NOT absorbed
rerooted_tri = reroot(parse_newick(nwk_trifurc), 't1')
leaves_tri = sorted(get_leaves(rerooted_tri))
check("reroot trifurc: all 6 leaves", leaves_tri, ['t1','t2','t3','t4','t5','t6'])
check("reroot trifurc: 2 top-level children", len(rerooted_tri['children']), 2)

print("\n--- label_tree.py: hex labeling ---")

# Known 6-taxon tree with tax-order [t1..t6]
# After reroot on t1, internal nodes should get predictable hex codes.
# Clade {t5,t6}: last 2 taxa, jak_val has bits 0,1 set = 3
# Clade {t4,t5,t6}: last 3 taxa, bits 0,1,2 set = 7
# Clade {t3,t4,t5,t6}: last 4, bits 0,1,2,3 = F
rerooted6 = reroot(parse_newick(nwk_nolen), 't1')
label_nodes(rerooted6, order6, 6)

def find_label(node, clade_leaves):
    """Find the label of the internal node whose clade is exactly clade_leaves."""
    if not node['children']:
        return None
    leaves = frozenset(get_leaves(node))
    if leaves == frozenset(clade_leaves):
        return node.get('name', '')
    for child in node['children']:
        result = find_label(child, clade_leaves)
        if result is not None:
            return result
    return None

check("label: {t5,t6} = '3'",       find_label(rerooted6, ['t5','t6']), '3')
check("label: {t4,t5,t6} = '7'",    find_label(rerooted6, ['t4','t5','t6']), '7')
check("label: {t3,t4,t5,t6} = 'F'", find_label(rerooted6, ['t3','t4','t5','t6']), 'F')
check("label: {t2,t3,t4,t5,t6} trivial (not labeled)",
      find_label(rerooted6, ['t2','t3','t4','t5','t6']), '')

# Test --label-prefix: should prepend to every non-empty label
rerooted6p = reroot(parse_newick(nwk_nolen), 't1')
label_nodes(rerooted6p, order6, 6, prefix='b_')
check("label prefix: {t5,t6} = 'b_3'", find_label(rerooted6p, ['t5','t6']), 'b_3')
check("label prefix: {t3,t4,t5,t6} = 'b_F'",
      find_label(rerooted6p, ['t3','t4','t5','t6']), 'b_F')

print("\n--- label_tree.py: Newick serialization ---")

rerooted6_nolen = reroot(parse_newick(nwk_nolen), 't1')
label_nodes(rerooted6_nolen, order6, 6)
serialized = newick_str(rerooted6_nolen, is_root=True)
check("serialize: ends with ';'", serialized.endswith(';'), True)
check("serialize: 6 leaves parseable",
      sorted(get_leaves(parse_newick(serialized))),
      ['t1','t2','t3','t4','t5','t6'])
# Labels should appear in the serialized string
check("serialize: label '3' present", '3' in serialized, True)
check("serialize: label 'F' present", 'F' in serialized, True)

print("\n--- label_tree.py: validation (mismatch detection) ---")
import subprocess, sys as _sys
result = subprocess.run(
    [_sys.executable, 'label_tree.py',
     '/tmp/test6.nwk', '--tax-order', '/tmp/taxorder6.txt'],
    capture_output=True, text=True, cwd=_HERE
)
# Write the test tree
import os as _os
with open('/tmp/test6.nwk', 'w') as _f:
    _f.write("((t1,t2),(t3,(t4,(t5,t6))));\n")
with open('/tmp/taxorder6.txt', 'w') as _f:
    _f.write('\n'.join(['t1','t2','t3','t4','t5','t6']) + '\n')

result_ok = subprocess.run(
    [_sys.executable, 'label_tree.py',
     '/tmp/test6.nwk', '--tax-order', '/tmp/taxorder6.txt'],
    capture_output=True, text=True, cwd=_HERE
)
check("validation: valid tree exits 0", result_ok.returncode, 0)
check("validation: output contains semicolon", ';' in result_ok.stdout, True)

# Bad tax-order (extra taxon not in tree)
with open('/tmp/bad_order.txt', 'w') as _f:
    _f.write('\n'.join(['t1','t2','t3','t4','t5','t7']) + '\n')
result_bad = subprocess.run(
    [_sys.executable, 'label_tree.py',
     '/tmp/test6.nwk', '--tax-order', '/tmp/bad_order.txt'],
    capture_output=True, text=True, cwd=_HERE
)
check("validation: mismatch exits non-zero", result_bad.returncode != 0, True)
check("validation: mismatch reports both directions",
      'In tree but not in --tax-order' in result_bad.stderr and
      'In --tax-order but not in tree' in result_bad.stderr, True)

# ---------------------------------------------------------------------------
print(f"\n=== Results: {PASS} passed, {FAIL} failed ===")
if FAIL > 0:
    sys.exit(1)
