def visualize_matrix(F, T, seq1, seq2):
    print("      " + "  ".join(" " + c for c in "-" + seq2))
    for i, row in enumerate(F):
        line = f"{'-' if i == 0 else seq1[i-1]} "
        for j, val in enumerate(row):
            if i == 0 and j == 0:
                line += "  0"
            elif T[i][j] is None:
                line += f" {val:2}"
            else:
                arrow = {'diag': '\\', 'up': '^', 'left': '<'}[T[i][j]]
                line += f" {arrow}{val:2}"
        print(line)


def align(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Needleman-Wunsch DP algo
    """
    m, n = len(seq1), len(seq2)
    
    # init score matrix F and traceback matrix T
    F = [[0] * (n + 1) for _ in range(m + 1)]
    T = [[None] * (n + 1) for _ in range(m + 1)]  # 'diag', 'up', 'left'

    # init first row and column
    for i in range(1, m + 1):
        F[i][0] = i * gap
        T[i][0] = 'up'
    for j in range(1, n + 1):
        F[0][j] = j * gap
        T[0][j] = 'left'

    # fill in the matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diag = F[i - 1][j - 1] + s
            up = F[i - 1][j] + gap
            left = F[i][j - 1] + gap
            F[i][j] = max(diag, up, left)

            # traceback
            if F[i][j] == diag:
                T[i][j] = 'diag'
            elif F[i][j] == up:
                T[i][j] = 'up'
            else:
                T[i][j] = 'left'

    # get the alignment
    aligned1, aligned2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if T[i][j] == 'diag':
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif T[i][j] == 'up':
            aligned1.append(seq1[i - 1])
            aligned2.append('-')
            i -= 1
        else:  # 'left'
            aligned1.append('-')
            aligned2.append(seq2[j - 1])
            j -= 1

    # visualize_matrix(F,T,seq1,seq2)

    return F[m][n], ''.join(reversed(aligned1)), ''.join(reversed(aligned2))

seq1 = "ACGTGCTAGCTAGCTAGGCTAGCTACGTAGCTAGCTAGCA"
seq2 = "ACGTCGATGCTAGATAGGATAGCTCGGCTAAGGTAGCTAGTTAGCA"

score, a1, a2 = align(seq1, seq2)
print("Score:", score)
print(a1)
print(a2)


