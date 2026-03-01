# src/smith_waterman.py
from src.scoring import score, GAP_PENALTY

def smith_waterman(seq1, seq2):
    m = len(seq1)
    n = len(seq2)

    # Create matrix (all zeros)
    matrix = [[0] * (n + 1) for _ in range(m + 1)]

    max_score = 0 # Tracks the highest score anywhere in the matrix (starts at 0)
    max_pos   = (0, 0) # Tracks the position (row, column) of that highest score (starts at origin)

    # Fill the matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diagonal = matrix[i-1][j-1] + score(seq1[i-1], seq2[j-1])
            up       = matrix[i-1][j]   + GAP_PENALTY
            left     = matrix[i][j-1]   + GAP_PENALTY

            matrix[i][j] = max(0, diagonal, up, left)

            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos   = (i, j)

    # Traceback
    aligned1 = ""
    aligned2 = ""
    i, j = max_pos

    while i > 0 and j > 0 and matrix[i][j] != 0:
        current = matrix[i][j]

        if current == matrix[i-1][j-1] + score(seq1[i-1], seq2[j-1]):
            aligned1 = seq1[i-1] + aligned1
            aligned2 = seq2[j-1] + aligned2
            i -= 1
            j -= 1
        elif current == matrix[i-1][j] + GAP_PENALTY:
            aligned1 = seq1[i-1] + aligned1
            aligned2 = "-"       + aligned2
            i -= 1
        else:
            aligned1 = "-"       + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1

    return aligned1, aligned2, max_score, max_pos