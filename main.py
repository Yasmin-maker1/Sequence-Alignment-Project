# main.py
from src.needleman_wunsch import needleman_wunsch
from src.smith_waterman   import smith_waterman

def display_alignment(seq1, seq2):
    middle = ""
    for a, b in zip(seq1, seq2):
        if a == b:       middle += "|"
        elif "-" in (a, b): middle += " "
        else:            middle += "."
    print(f"  Seq1: {seq1}")
    print(f"        {middle}")
    print(f"  Seq2: {seq2}")

# ── TEST 1: Perfect match ──────────────────────────
print("=" * 50)
print("TEST 1: Perfect match")
print("=" * 50)
a1, a2, sc = needleman_wunsch("ATCG", "ATCG")
print(f"\nNeedleman-Wunsch | Score: {sc}")
display_alignment(a1, a2)

a1, a2, sc, pos = smith_waterman("ATCG", "ATCG")
print(f"\nSmith-Waterman   | Score: {sc}")
display_alignment(a1, a2)

# ── TEST 2: Needle in a haystack ───────────────────
print("\n" + "=" * 50)
print("TEST 2: Needle in a haystack")
print("=" * 50)
a1, a2, sc = needleman_wunsch("ATCG", "GGGGGATCGTTTT")
print(f"\nNeedleman-Wunsch | Score: {sc}")
display_alignment(a1, a2)

a1, a2, sc, pos = smith_waterman("ATCG", "GGGGGATCGTTTT")
print(f"\nSmith-Waterman   | Score: {sc}")
display_alignment(a1, a2)
print("\n→ SW finds the exact match. NW score is low (forced full alignment).")