# tests/test_algorithms.py
# ─────────────────────────────────────────────────────────────
# Formal validation suite for Needleman-Wunsch and Smith-Waterman
# Scoring: match = +1, mismatch = -1, gap = -2
#
# Run from project root with:
#   python tests/test_algorithms.py
# ─────────────────────────────────────────────────────────────

import sys
import os

# This line lets Python find the src/ folder from the tests/ folder
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.needleman_wunsch import needleman_wunsch
from src.smith_waterman   import smith_waterman

# ── Simple pass/fail tracker ──────────────────────────────────
PASS = 0
FAIL = 0

def check(name, condition, detail=""):
    """Prints PASS or FAIL for each test and tracks the count."""
    global PASS, FAIL
    if condition:
        print(f"  PASS ✓  {name}")
        PASS += 1
    else:
        print(f"  FAIL ✗  {name}")
        print(f"          Expected: {detail}")
        FAIL += 1


# ══════════════════════════════════════════════════════════════
# NEEDLEMAN-WUNSCH TESTS
#
# NW always aligns the FULL length of both sequences.
# It uses gap penalties in the first row and column.
# The final score is always in the bottom-right cell.
# ══════════════════════════════════════════════════════════════
print("\n" + "=" * 55)
print("  NEEDLEMAN-WUNSCH TESTS")
print("=" * 55)

# ── NW Test 1: Perfect match ──────────────────────────────────
# ATCG vs ATCG → every base matches → 4 matches x (+1) = 4
# Expected alignment:
#   ATCG
#   ||||
#   ATCG
a1, a2, sc = needleman_wunsch("ATCG", "ATCG")
check("NW-1 | Perfect match → score = 4",
      sc == 4, f"score 4, got {sc}")
check("NW-1 | Perfect match → no gaps inserted",
      "-" not in a1 and "-" not in a2,
      f"no dashes, got a1={a1}, a2={a2}")

# ── NW Test 2: One mismatch ───────────────────────────────────
# ATCG vs ATCC → last base G != C → 3 matches + 1 mismatch
# Score: 3(+1) + 1(-1) = 2
a1, a2, sc = needleman_wunsch("ATCG", "ATCC")
check("NW-2 | One mismatch → score = 2",
      sc == 2, f"score 2, got {sc}")

# ── NW Test 3: One gap required ───────────────────────────────
# ATCG vs ACG → NW must insert one gap to align fully
# Best alignment:  A T C G
#                  A - C G
# Score: 3 matches + 1 gap = 3(+1) + 1(-2) = 1
a1, a2, sc = needleman_wunsch("ATCG", "ACG")
check("NW-3 | One gap required → score = 1",
      sc == 1, f"score 1, got {sc}")
check("NW-3 | One gap required → gap appears in alignment",
      "-" in a1 or "-" in a2,
      f"a dash expected, got a1={a1}, a2={a2}")

# ── NW Test 4: All mismatches ─────────────────────────────────
# AT vs GC → both bases mismatch → 2(-1) = -2
a1, a2, sc = needleman_wunsch("AT", "GC")
check("NW-4 | All mismatches → score = -2",
      sc == -2, f"score -2, got {sc}")

# ── NW Test 5: Empty sequence ─────────────────────────────────
# "" vs "ATG" → NW must gap-fill the entire first sequence
# Score: 3 gaps x (-2) = -6
a1, a2, sc = needleman_wunsch("", "ATG")
check("NW-5 | Empty vs ATG → score = -6",
      sc == -6, f"score -6, got {sc}")

# ── NW Test 6: Classic textbook example ───────────────────────
# GATTACA vs GCATGCU is a standard bioinformatics validation case
# No matter how many gaps are added, original bases must all appear
a1, a2, sc = needleman_wunsch("GATTACA", "GCATGCU")
check("NW-6 | GATTACA/GCATGCU → aligned lengths are equal",
      len(a1) == len(a2),
      f"equal lengths, got {len(a1)} vs {len(a2)}")
check("NW-6 | GATTACA/GCATGCU → all original bases preserved",
      a1.replace("-", "") == "GATTACA" and a2.replace("-", "") == "GCATGCU",
      f"bases preserved, got a1={a1}, a2={a2}")

# ── NW Test 7: Single character inputs ───────────────────────
check("NW-7 | Single match 'A' vs 'A' → score = 1",
      needleman_wunsch("A", "A")[2] == 1, "score 1")
check("NW-7 | Single mismatch 'A' vs 'T' → score = -1",
      needleman_wunsch("A", "T")[2] == -1, "score -1")


# ══════════════════════════════════════════════════════════════
# SMITH-WATERMAN TESTS
#
# SW finds the BEST LOCAL alignment anywhere in both sequences.
# Key differences from NW:
#   - All cells have a zero floor (never go below 0)
#   - Traceback starts at the MAXIMUM cell, not bottom-right
#   - Traceback STOPS when it hits a zero
# ══════════════════════════════════════════════════════════════
print("\n" + "=" * 55)
print("  SMITH-WATERMAN TESTS")
print("=" * 55)

# ── SW Test 1: Perfect match ──────────────────────────────────
# Full match is also the best local match → score = 4
a1, a2, sc, pos = smith_waterman("ATCG", "ATCG")
check("SW-1 | Perfect match → score = 4",
      sc == 4, f"score 4, got {sc}")

# ── SW Test 2: Needle in a haystack ──────────────────────────
# ATCG embedded inside GGGGGATCGTTTT
# SW ignores flanking Gs and Ts and finds exact match → score = 4
a1, a2, sc, pos = smith_waterman("ATCG", "GGGGGATCGTTTT")
check("SW-2 | Needle in haystack → score = 4",
      sc == 4, f"score 4, got {sc}")
check("SW-2 | Needle in haystack → finds correct subsequence",
      a1 == "ATCG" and a2 == "ATCG",
      f"both aligned = ATCG, got a1={a1}, a2={a2}")

# ── SW Test 3: Score is NEVER negative ───────────────────────
# Zero floor means SW can always start fresh → score >= 0
_, _, sc, _ = smith_waterman("AAAA", "TTTT")
check("SW-3 | Score never negative (zero floor)",
      sc >= 0, f"score >= 0, got {sc}")

# ── SW Test 4: No similarity → score = 0 ─────────────────────
# AAAA vs TTTT → no matching bases → all cells reset to 0
_, _, sc, _ = smith_waterman("AAAA", "TTTT")
check("SW-4 | No similarity → score = 0",
      sc == 0, f"score 0, got {sc}")

# ── SW Test 5: Finds best local region, ignores flanking ──────
# TTTACGTTT has ACGT buried in the middle → SW finds it → score = 4
a1, a2, sc, pos = smith_waterman("TTTACGTTT", "ACGT")
check("SW-5 | Finds best local region → score = 4",
      sc == 4, f"score 4, got {sc}")
check("SW-5 | Aligned region is ACGT",
      "ACGT" in a1 or "ACGT" in a2,
      f"ACGT expected, got a1={a1}, a2={a2}")

# ── SW Test 6: Aligned strings always same length ────────────
# Gaps are added to make both aligned strings equal length
a1, a2, sc, _ = smith_waterman("ACACACTA", "AGCACACA")
check("SW-6 | Aligned strings have equal length",
      len(a1) == len(a2),
      f"equal lengths, got {len(a1)} vs {len(a2)}")

# ── SW Test 7: SW score always >= NW score ────────────────────
# SW is more flexible → for same input SW score >= NW score
_, _, nw_sc    = needleman_wunsch("ATCG", "GGGGGATCGTTTT")
_, _, sw_sc, _ = smith_waterman("ATCG",  "GGGGGATCGTTTT")
check("SW-7 | SW score >= NW score for needle-in-haystack",
      sw_sc >= nw_sc, f"SW({sw_sc}) >= NW({nw_sc})")

# ── SW Test 8: Single character match ────────────────────────
_, _, sc, _ = smith_waterman("A", "A")
check("SW-8 | Single match 'A' vs 'A' → score = 1",
      sc == 1, f"score 1, got {sc}")

# ── SW Test 9: Single mismatch hits zero floor ────────────────
# A vs T → would be -1 but zero floor applies → score = 0
_, _, sc, _ = smith_waterman("A", "T")
check("SW-9 | Single mismatch 'A' vs 'T' → score = 0 (zero floor)",
      sc == 0, f"score 0, got {sc}")


# ══════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════
print("\n" + "=" * 55)
total = PASS + FAIL
print(f"  Results: {PASS}/{total} tests passed", "!" if FAIL == 0 else "check failures above")
print("=" * 55 + "\n")