#!/usr/bin/env python3
"""Rosalind LOCA - Local Alignment with Scoring Matrix (Smith-Waterman, PAM250)."""

AA = "ACDEFGHIKLMNPQRSTVWY"
_PAM250 = [
    [ 2,-2, 0, 0,-3, 1,-1,-1,-1,-2,-1, 0, 1, 0,-2, 1, 1, 0,-6,-3],
    [-2,12,-5,-5,-4,-3,-3,-2,-5,-6,-5,-4,-3,-5,-4, 0,-2,-2,-8, 0],
    [ 0,-5, 4, 3,-6, 1, 1,-2, 0,-4,-3, 2,-1, 2,-1, 0, 0,-2,-7,-4],
    [ 0,-5, 3, 4,-5, 0, 1,-2, 0,-3,-2, 1,-1, 2,-1, 0, 0,-2,-7,-4],
    [-3,-4,-6,-5, 9,-5,-2, 1,-5, 2, 0,-3,-5,-5,-4,-3,-3,-1, 0, 7],
    [ 1,-3, 1, 0,-5, 5,-2,-3,-2,-4,-3, 0, 0,-1,-3, 1, 0,-1,-7,-5],
    [-1,-3, 1, 1,-2,-2, 6,-2, 0,-2,-2, 2, 0, 3, 2,-1,-1,-2,-3, 0],
    [-1,-2,-2,-2, 1,-3,-2, 5,-2, 2, 2,-2,-2,-2,-2,-1, 0, 4,-5,-1],
    [-1,-5, 0, 0,-5,-2, 0,-2, 5,-3, 0, 1,-1, 1, 3, 0, 0,-2,-3,-4],
    [-2,-6,-4,-3, 2,-4,-2, 2,-3, 6, 4,-3,-3,-2,-3,-3,-2, 2,-2,-1],
    [-1,-5,-3,-2, 0,-3,-2, 2, 0, 4, 6,-2,-2,-1, 0,-2,-1, 2,-4,-2],
    [ 0,-4, 2, 1,-3, 0, 2,-2, 1,-3,-2, 2, 0, 1, 0, 1, 0,-2,-4,-2],
    [ 1,-3,-1,-1,-5, 0, 0,-2,-1,-3,-2, 0, 6, 0, 0, 1, 0,-1,-6,-5],
    [ 0,-5, 2, 2,-5,-1, 3,-2, 1,-2,-1, 1, 0, 4, 1,-1,-1,-2,-5,-4],
    [-2,-4,-1,-1,-4,-3, 2,-2, 3,-3, 0, 0, 0, 1, 6, 0,-1,-2, 2,-4],
    [ 1, 0, 0, 0,-3, 1,-1,-1, 0,-3,-2, 1, 1,-1, 0, 2, 1,-1,-2,-3],
    [ 1,-2, 0, 0,-3, 0,-1, 0, 0,-2,-1, 0, 0,-1,-1, 1, 3, 0,-5,-3],
    [ 0,-2,-2,-2,-1,-1,-2, 4,-2, 2, 2,-2,-1,-2,-2,-1, 0, 4,-6,-2],
    [-6,-8,-7,-7, 0,-7,-3,-5,-3,-2,-4,-4,-6,-5, 2,-2,-5,-6,17, 0],
    [-3, 0,-4,-4, 7,-5, 0,-1,-4,-1,-2,-2,-5,-4,-4,-3,-3,-2, 0,10],
]
PAM = {a: {b: _PAM250[i][j] for j, b in enumerate(AA)} for i, a in enumerate(AA)}
ARROW = {(1, 1): "\\", (1, 0): "^", (0, 1): "<"}      # diagonal, up, left


def clean_sequence(raw):
    """Strip FASTA header lines and spaces; validate the 20 amino-acid letters."""
    lines = [l.strip() for l in raw.strip().splitlines() if not l.strip().startswith(">")]
    seq = "".join("".join(lines).split()).upper()
    if not seq:
        raise ValueError("Sequence cannot be empty.")
    bad = sorted(set(seq) - set(AA))
    if bad:
        raise ValueError(f"Invalid amino acid symbol(s): {bad}")
    return seq


def smith_waterman(s, t, gap):
    """Fill H, keep every tied direction, then trace back from every peak cell."""
    n, m = len(s), len(t)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    T = [[[] for _ in range(m + 1)] for _ in range(n + 1)]
    best, peaks = 0, []
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            moves = {(1, 1): H[i-1][j-1] + PAM[s[i-1]][t[j-1]],
                     (1, 0): H[i-1][j] - gap, (0, 1): H[i][j-1] - gap}
            H[i][j] = v = max(0, *moves.values())
            if v > 0:
                T[i][j] = [mv for mv, sc in moves.items() if sc == v]
            if v > best:
                best, peaks = v, [(i, j)]
            elif v == best and v > 0:
                peaks.append((i, j))

    found, seen = [], set()

    def trace(i, j, a, b, path, end):
        if H[i][j] == 0:                       # reached a zero cell: alignment complete
            key = (i + 1, j + 1, end, a, b)
            if key not in seen:
                seen.add(key)
                found.append({"start": (i + 1, j + 1), "end": end, "a": a, "b": b, "path": set(path)})
            return
        for di, dj in T[i][j]:
            trace(i - di, j - dj, (s[i-1] if di else "-") + a,
                  (t[j-1] if dj else "-") + b, path + [(i - di, j - dj)], end)

    for pi, pj in peaks:
        trace(pi, pj, "", "", [(pi, pj)], (pi, pj))
    return {"H": H, "T": T, "score": best, "all": found,
            "path": set().union(*(f["path"] for f in found)) if found else set()}


def show_matrix(s, t, r):
    """Print H[i][j]; cells on a traceback path are shown as [value]."""
    n, m, H, T = len(s), len(t), r["H"], r["T"]
    ei, ej = r["all"][0]["end"] if r["all"] else (n, m)
    i0, i1 = (0, n) if n <= 22 else (max(0, ei - 14), ei)
    j0, j1 = (0, m) if m <= 22 else (max(0, ej - 14), ej)
    print("      " + "".join(f"{('-' if j == 0 else t[j-1]) + f'({j})':^7}" for j in range(j0, j1 + 1)))
    for i in range(i0, i1 + 1):
        row = f"{('-' if i == 0 else s[i-1]) + f'({i})':<6}"
        for j in range(j0, j1 + 1):
            cell = f"{H[i][j]}{''.join(ARROW[d] for d in T[i][j])}"
            row += f"{'[' + cell + ']' if (i, j) in r['path'] else cell:^7}"
        print(row)
    print("  [value] = traceback path;  \\ diagonal  ^ up (gap in t)  < left (gap in s)")


def run_case(s, t, gap):
    r = smith_waterman(s, t, gap)
    first = r["all"][0] if r["all"] else {"a": "", "b": "", "start": (1, 1), "end": (0, 0)}
    a, b = first["a"], first["b"]
    mid = "".join(" " if "-" in (x, y) else "|" if x == y else ":" if PAM[x][y] > 0
                  else "." if PAM[x][y] == 0 else " " for x, y in zip(a, b))
    ident, pos = mid.count("|"), mid.count("|") + mid.count(":")
    L = len(a) or 1

    print(f"\n[STEP 1] INPUT\n  s = {s}  ({len(s)} aa)\n  t = {t}  ({len(t)} aa)")
    print(f"  Matrix = PAM250, linear gap penalty = {gap}")
    print(f"\n[STEP 2] SMITH-WATERMAN MATRIX H (peak score {r['score']})")
    show_matrix(s, t, r)
    print(f"\n[STEP 3] OPTIMAL LOCAL ALIGNMENT  (s[{first['start'][0]}..{first['end'][0]}], "
          f"t[{first['start'][1]}..{first['end'][1]}])")
    print(f"  s : {a}\n      {mid}\n  t : {b}")
    print(f"  Identity {ident}/{L} ({ident / L:.1%}), similarity {pos}/{L} ({pos / L:.1%}), "
          f"gaps {a.count('-') + b.count('-')}/{L}")
    print("\n[STEP 4] FINAL ROSALIND 'LOCA' OUTPUT")
    print(r["score"], a.replace("-", ""), b.replace("-", ""), sep="\n")
    print(f"\n[STEP 5] ALL DISTINCT OPTIMAL ALIGNMENTS: {len(r['all'])}")
    for k, f in enumerate(r["all"], 1):
        print(f"  #{k}  s[{f['start'][0]}..{f['end'][0]}] t[{f['start'][1]}..{f['end'][1]}]\n"
              f"      {f['a']}\n      {f['b']}")


def read_input():
    """Read two validated sequences and a gap penalty (default 5)."""
    seqs = []
    for k in (1, 2):
        while True:
            try:
                seqs.append(clean_sequence(input(f"  Enter Sequence {k}: ")))
                break
            except ValueError as e:
                print("  Invalid:", e)
    while True:
        raw = input("  Enter gap penalty [5]: ").strip() or "5"
        if raw.isdigit():
            return seqs[0], seqs[1], int(raw)
        print("  Please enter a whole number.")


def main():
    sample = ("MEANLYPRTEINSTRING", "PLEASANTLYEINSTEIN", 5)
    print("\n=== LOCAL ALIGNMENT (ROSALIND: LOCA) ===")
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