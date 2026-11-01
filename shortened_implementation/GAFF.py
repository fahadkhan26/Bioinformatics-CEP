#!/usr/bin/env python3
"""Rosalind GAFF - Global alignment with BLOSUM62 and affine gap penalties (open 11, extend 1)."""

AA = "ARNDCQEGHILKMFPSTWYV"
_BLOSUM = [
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1, 1,-1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
]
B = {a: {b: _BLOSUM[i][j] for j, b in enumerate(AA)} for i, a in enumerate(AA)}
NEG = -10**12
MOVES = [(1, 1), (1, 0), (0, 1)]          # state 0 = M (pair), 1 = X (gap in t), 2 = Y (gap in s)


def clean_sequence(raw):
    """Strip FASTA header lines and spaces; validate the 20 amino-acid letters."""
    lines = [l.strip() for l in raw.strip().splitlines() if l.strip() and not l.strip().startswith(">")]
    seq = "".join("".join(lines).split()).upper()
    if not seq:
        raise ValueError("Sequence cannot be empty.")
    bad = sorted(set(seq) - set(AA))
    if bad:
        raise ValueError(f"Invalid amino acid symbols: {bad}")
    return seq


def affine_global_alignment(s, t, go=11, ge=1):
    """Three-state DP (M, X, Y) with predecessor pointers, then traceback from (n, m)."""
    n, m = len(s), len(t)
    grid = lambda v: [[v] * (m + 1) for _ in range(n + 1)]
    M, X, Y = grid(NEG), grid(NEG), grid(NEG)
    P = [grid(None), grid(None), grid(None)]          # predecessor state for M, X, Y
    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0], P[1][i][0] = -go - (i - 1) * ge, (0 if i == 1 else 1)
    for j in range(1, m + 1):
        Y[0][j], P[2][0][j] = -go - (j - 1) * ge, (0 if j == 1 else 2)

    best = lambda open_, ext, k: (open_, 0) if open_ >= ext else (ext, k)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            prev = (M[i-1][j-1], X[i-1][j-1], Y[i-1][j-1])
            k = max(range(3), key=prev.__getitem__)
            M[i][j], P[0][i][j] = prev[k] + B[s[i-1]][t[j-1]], k
            X[i][j], P[1][i][j] = best(M[i-1][j] - go, X[i-1][j] - ge, 1)
            Y[i][j], P[2][i][j] = best(M[i][j-1] - go, Y[i][j-1] - ge, 2)

    end = (M[n][m], X[n][m], Y[n][m])
    k = max(range(3), key=end.__getitem__)
    score, i, j, a, b, path = end[k], n, m, "", "", set()
    while i or j:
        path.add((i, j, k))
        di, dj = MOVES[k]
        a, b = (s[i-1] if di else "-") + a, (t[j-1] if dj else "-") + b
        k, i, j = P[k][i][j], i - di, j - dj
    path.add((0, 0, k))
    return {"score": score, "a": a, "b": b, "M": M, "X": X, "Y": Y, "P": P, "path": path}


def rescore(a, b, go, ge):
    """Independent check: add BLOSUM62 pairs and charge go + ge*(L-1) per gap run."""
    total, in_gap = 0, False
    for x, y in zip(a, b):
        if "-" in (x, y):
            total -= ge if in_gap else go
            in_gap = True
        else:
            total += B[x][y]
            in_gap = False
    return total


def show_matrices(s, t, r):
    if len(s) > 22 or len(t) > 22:
        print("  DP matrices omitted (a sequence is longer than 22 aa).")
        return
    w = 7
    cols = ["-(0)"] + [f"{t[j-1]}({j})" for j in range(1, len(t) + 1)]
    head = "i\\j".center(w) + "".join(f"{c:^{w}}" for c in cols)
    for k, (name, what, arrow) in enumerate((("M", "paired residues", "\\"),
                                             ("X", "gap in t", "^"), ("Y", "gap in s", "<"))):
        print(f"\n  State {name}: {what}\n{head}")
        for i in range(len(s) + 1):
            row = f"{('-' if i == 0 else s[i-1]) + f'({i})':^{w}}"
            for j in range(len(t) + 1):
                v = r[name][i][j]
                cell = "-INF" if v < NEG // 2 else str(v) + ("" if r["P"][k][i][j] is None else arrow)
                if v > NEG // 2 and (i, j, k) in r["path"]:
                    cell = f"[{cell}]"
                row += f"{cell:^{w}}"
            print(row)
    print("\n  [cell] = traceback path;  -INF = unreachable;  \\ pair  ^ gap in t  < gap in s")


def run_case(s, t, go=11, ge=1):
    r = affine_global_alignment(s, t, go, ge)
    a, b = r["a"], r["b"]
    pairs = [(x, y) for x, y in zip(a, b) if "-" not in (x, y)]
    ident = sum(x == y for x, y in pairs)
    mid = "".join(" " if "-" in (x, y) else "|" if x == y else ":" if B[x][y] > 0
                  else "." if B[x][y] == 0 else " " for x, y in zip(a, b))

    print(f"\n[STEP 1] INPUT\n  s = {s}  ({len(s)} aa)\n  t = {t}  ({len(t)} aa)")
    print(f"  Matrix = BLOSUM62, gap opening = {go}, gap extension = {ge}")
    print(f"  A gap of length L costs {go} + {ge} x (L - 1)")
    print(f"\n[STEP 2] AFFINE-GAP DP MATRICES (optimal score {r['score']})")
    show_matrices(s, t, r)
    print(f"\n[STEP 3] OPTIMAL GLOBAL ALIGNMENT ({len(a)} columns)")
    print(f"  s : {a}\n      {mid}\n  t : {b}")
    print(f"  Identity {ident}/{len(pairs)} paired residues, gap symbols {len(a) - len(pairs)}")
    print("\n[STEP 4] FINAL ROSALIND 'GAFF' OUTPUT")
    print(r["score"], a, b, sep="\n")
    check = rescore(a, b, go, ge)
    print("\n[STEP 5] VALIDATION")
    print(f"  Recalculated score = {check} -> {'MATCH' if check == r['score'] else 'MISMATCH'}")
    print(f"  Gaps removed give back s and t: {a.replace('-', '') == s and b.replace('-', '') == t}")


def read_input():
    seqs = []
    for k in (1, 2):
        while True:
            try:
                seqs.append(clean_sequence(input(f"  Enter Sequence {k}: ")))
                break
            except ValueError as e:
                print("  Invalid:", e)
    return seqs


def main():
    sample = ("PRTEINS", "PRTWPSEIN")
    print("\n=== GLOBAL ALIGNMENT, AFFINE GAP (ROSALIND: GAFF) ===")
    print("[1] Sample Case  [2] Custom Case  [3] Sample + Custom  [4] Exit")
    choice = input("Select an option [1-4]: ").strip()
    if choice in ("1", "3"):
        print("\nRunning the sample case...")
        run_case(*sample)
    if choice in ("2", "3"):
        run_case(*read_input())
    if choice not in ("1", "2", "3", "4"):
        print("Invalid choice.")


if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, EOFError):
        print("\nProgram stopped by user.")
