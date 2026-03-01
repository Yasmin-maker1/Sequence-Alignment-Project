# src/scoring.py

MATCH_SCORE    =  1
MISMATCH_SCORE = -1
GAP_PENALTY    = -2

def score(base1, base2):
    if base1 == base2:
        return MATCH_SCORE
    else:
        return MISMATCH_SCORE