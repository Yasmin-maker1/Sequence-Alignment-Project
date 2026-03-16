# experiments/biological_experiments.py
# Three biological use cases:
#   1. Protein comparison       → Needleman-Wunsch (global)
#   2. Mutation detection       → Smith-Waterman (local)
#   3. Short motif matching     → Smith-Waterman (local)
#
# Run from project root:
#   python experiments/biological_experiments.py

import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.needleman_wunsch import needleman_wunsch
from src.smith_waterman   import smith_waterman

def print_alignment(a1, a2, label1="Seq1", label2="Seq2", max_width=60):
    middle = ""
    for c1, c2 in zip(a1, a2):
        if c1 == c2:          middle += "|"
        elif "-" in (c1, c2): middle += " "
        else:                 middle += "."
    for start in range(0, len(a1), max_width):
        end = start + max_width
        print(f"  {label1:<10} {a1[start:end]}")
        print(f"  {'':10} {middle[start:end]}")
        print(f"  {label2:<10} {a2[start:end]}")
        print()

def identity_percent(a1, a2):
    matches = sum(1 for c1, c2 in zip(a1, a2) if c1 == c2 and c1 != "-")
    return 100.0 * matches / len(a1)

# ══════════════════════════════════════════════════════════════
# EXPERIMENT 1: Protein Comparison — Needleman-Wunsch
# Human hemoglobin alpha (UniProt P69905) vs
# Horse hemoglobin alpha (UniProt P01958)
# ══════════════════════════════════════════════════════════════
print("=" * 65)
print("EXPERIMENT 1: Human vs Horse Hemoglobin Alpha")
print("Algorithm: Needleman-Wunsch (Global Alignment)")
print("Source: UniProt P69905 (human), P01958 (horse)")
print("=" * 65)

human_hba = (
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"
    "KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP"
    "AVHASLDKFLASVSTVLTSKYR"
)

horse_hba = (
    "MVLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFPTTKTYFPHFDLSHGSAQVKAHG"
    "KKVGDALTLAVGHLDDLPGALSNLSDLHAHKLRVDPVNFKLLSHCLLSTLAVHLPNDFTPAVHASLDKFLSSVSTVLTSKYR"
)

a1, a2, sc = needleman_wunsch(human_hba, horse_hba)
pct_id = identity_percent(a1, a2)

print(f"\n  Alignment score    : {sc}")
print(f"  Alignment length   : {len(a1)} positions")
print(f"  Percent identity   : {pct_id:.1f}%")
print(f"  Gaps in human seq  : {a1.count('-')}")
print(f"  Gaps in horse seq  : {a2.count('-')}")
print()
print_alignment(a1, a2, "Human", "Horse")
print("→ High identity confirms evolutionary conservation.")
print("→ NW is appropriate: both sequences are full-length homologs.\n")


# ══════════════════════════════════════════════════════════════
# EXPERIMENT 2: Mutation Detection — Smith-Waterman
# BRCA1 wild-type vs 185delAG mutant (pathogenic frameshift)
# ══════════════════════════════════════════════════════════════
print("=" * 65)
print("EXPERIMENT 2: BRCA1 Wild-Type vs 185delAG Mutant")
print("Algorithm: Smith-Waterman (Local Alignment)")
print("Source: NCBI NM_007294, known pathogenic variant 185delAG")
print("=" * 65)

brca1_wildtype = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAAC"

# Simulate 185delAG: delete 2 bases at position 39
brca1_mutant = brca1_wildtype[:39] + brca1_wildtype[41:]

a1, a2, sc, pos = smith_waterman(brca1_wildtype, brca1_mutant)

print(f"\n  Best local score   : {sc}")
print(f"  Alignment end pos  : row={pos[0]}, col={pos[1]}")
print(f"  Aligned length     : {len(a1)} positions")
print(f"  Gaps detected      : {a1.count('-') + a2.count('-')} (marks deletion site)")
print()
print_alignment(a1[:80], a2[:80], "WT", "Mutant")
print("→ SW identifies the local region of maximum similarity.")
print("→ Gaps mark the exact 2-base deletion site.\n")


# ══════════════════════════════════════════════════════════════
# EXPERIMENT 3: Motif Matching — Smith-Waterman
# TATA box consensus motif vs synthetic promoter sequence
# ══════════════════════════════════════════════════════════════
print("=" * 65)
print("EXPERIMENT 3: TATA Box Motif in Promoter Sequence")
print("Algorithm: Smith-Waterman (Local Alignment)")
print("Source: Consensus TATA box vs synthetic promoter region")
print("=" * 65)

tata_motif = "TATAAAA"
promoter        = "GCGCGCGCGCGCATGCATGCATGCATGCATGCTATAAAAGCGCGCGCGCGCATGCATGCATGCATGCATGC"
promoter_no_tata = "GCGCGCGCGCGCATGCATGCATGCATGCATGCGGGGGGCGCGCGCGCGCATGCATGCATGCATGCATGC"

a1, a2, sc, pos      = smith_waterman(tata_motif, promoter)
_, _,  sc_no, _      = smith_waterman(tata_motif, promoter_no_tata)

print(f"\n  Motif                  : {tata_motif}")
print(f"  Score WITH TATA box    : {sc}")
print(f"  Score WITHOUT TATA box : {sc_no}")
print(f"  Match end position     : col={pos[1]} in promoter")
print()
print_alignment(a1, a2, "Motif", "Promoter")
print(f"→ Score drops from {sc} to {sc_no} when TATA box is absent.")
print("→ SW reliably locates the motif and quantifies match quality.\n")

print("=" * 65)
print("All biological experiments complete.")