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


def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-2):
    m, n = len(seq1), len(seq2)
    
    # init
    F = [[0] * (n + 1) for _ in range(m + 1)]
    T = [[None] * (n + 1) for _ in range(m + 1)]  # 'diag', 'up', 'left', or None

    max_score = 0
    max_pos = (0, 0)

    # fill DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diag = F[i - 1][j - 1] + s
            up = F[i - 1][j] + gap
            left = F[i][j - 1] + gap
            F[i][j] = max(0, diag, up, left)

            # track pointers for traceback
            if F[i][j] == 0:
                T[i][j] = None
            elif F[i][j] == diag:
                T[i][j] = 'diag'
            elif F[i][j] == up:
                T[i][j] = 'up'
            else:
                T[i][j] = 'left'

            # update score and pos
            if F[i][j] > max_score:
                max_score = F[i][j]
                max_pos = (i, j)

    # traceback
    aligned1, aligned2 = [], []
    i, j = max_pos
    while T[i][j] and F[i][j] > 0:
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
    return max_score, ''.join(reversed(aligned1)), ''.join(reversed(aligned2))


seq1 = "TTTACGTGCTAGCTAGCTTTTGGG"
seq2 = "GGGACGTCGATGCTAGATAGCCC"
print("Input:", "\n", seq1,"\n", seq2)

score, a1, a2 = smith_waterman(seq1, seq2)
print("Score:", score)
print(a1)
print(a2)

