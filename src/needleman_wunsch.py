# src/needleman_wunsch.py
from src.scoring import score, GAP_PENALTY # Bring "score" function and "GAP_PENALTY" variable into this file

# Definition of the main function
def needleman_wunsch(seq1, seq2):   
    m = len(seq1) # "len()" counts how many characters are in a string
    n = len(seq2) # O(mn)

    # Create scoring matrix
    matrix = [[0] * (n + 1) for _ in range(m + 1)] # The result is a grid of zeros with (m+1) rows and (n+1) columns


    # Initialize first column (fills column 0 to m)
    for i in range(m + 1):
        matrix[i][0] = i * GAP_PENALTY # Example: Row 0: 0 * -2 = 0, Row 1: 1 * -2 = -2, etc.

    # Initialize first row (fills row 0 to n))
    for j in range(n + 1):
        matrix[0][j] = j * GAP_PENALTY # Example: Row 0: 0 * -2 = 0, Row 1: 1 * -2 = -2, etc.


    # Fill the matrix 
    """ (Outer loop --> for i in range(1, m + 1): goes thorugh every row, 
    the innner loop --> for j in range(1, n + 1):goes through every column, 
    together they visit every single cell in the matrix exactly once),
    That is the O(mn) part --> m * n total cell visits."""
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diagonal = matrix[i-1][j-1] + score(seq1[i-1], seq2[j-1]) # Match or mismatch
            up       = matrix[i-1][j]   + GAP_PENALTY # Up/Top gap
            left     = matrix[i][j-1]   + GAP_PENALTY # Left/Horizontal gap
            matrix[i][j] = max(diagonal, up, left) # Selection of the largest of the three values above

    # Traceback
    aligned1 = ""   # Start with two empty strings that will become the final aligned sequence
    aligned2 = ""
    i, j = m, n  # Sets the starting position to the "Bottom-right corner"

    while i > 0 or j > 0:
        current = matrix[i][j]

        # Checks: "Did this cell's value come from the diagonal?""
        if i > 0 and j > 0 and current == matrix[i-1][j-1] + score(seq1[i-1], seq2[j-1]):
            aligned1 = seq1[i-1] + aligned1
            aligned2 = seq2[j-1] + aligned2
            i -= 1
            j -= 1
        # Check if the value come from above (up)
        elif i > 0 and current == matrix[i-1][j] + GAP_PENALTY:
            aligned1 = seq1[i-1] + aligned1
            aligned2 = "-"       + aligned2
            i -= 1
        # If neither diagonal nor up, so the value must come from left
        else:
            aligned1 = "-"       + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1

    return aligned1, aligned2, matrix[m][n]
    # aligned1: The first sequence with gaps inserted.
    # aligned2: The Second sequence with gaps inserted.
    # matrix[m][n]: The final alignment score (Always the the bottom-right cell in NW)