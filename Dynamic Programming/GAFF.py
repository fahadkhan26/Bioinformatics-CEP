#!/usr/bin/env python3
"""CS-434 Bioinformatics PBL — Rosalind GAFF
Global alignment using BLOSUM62 and affine gap penalties.
Terminal-only report; no output files are created.
"""

import sys
from pathlib import Path

# BLOSUM62 values for the 20 standard amino acids, in this fixed order.
AA = "ARNDCQEGHILKMFPSTWYV"
_BLOSUM_ROWS = [
    [4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    [-1,5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2,0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    [-2,-2,1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1,1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    [-1,0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    [0,-2,0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    [-2,0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
    [-1,2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1, 1,-1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
    [1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
    [0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
    [0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
]
assert len(_BLOSUM_ROWS) == len(AA) == 20, "BLOSUM62 must contain 20 rows and columns."
for row_index, row in enumerate(_BLOSUM_ROWS):
    assert len(row) == len(AA), (
        f"BLOSUM62 row {row_index} ({AA[row_index]}) has "
        f"{len(row)} values instead of {len(AA)}."
    )

BLOSUM62 = {a: {b: _BLOSUM_ROWS[i][j] for j, b in enumerate(AA)}
            for i, a in enumerate(AA)}
NEG_INF = -10**12


def clean_sequence(raw):
    lines = [line.strip() for line in raw.strip().splitlines()
             if line.strip() and not line.strip().startswith(">")]
    seq = "".join("".join(lines).split()).upper()
    invalid = sorted(set(seq) - set(AA))
    if not seq:
        raise ValueError("Sequence cannot be empty.")
    if invalid:
        raise ValueError(f"Invalid amino acid symbols: {invalid}. Use the 20 standard amino acids.")
    return seq


def read_fasta_file(path):
    text = Path(path).read_text(encoding="utf-8")
    records, current = [], []
    for line in text.splitlines():
        if line.startswith(">"):
            if current:
                records.append("".join(current))
                current = []
        elif line.strip():
            current.append(line.strip())
    if current:
        records.append("".join(current))
    if len(records) < 2:
        raise ValueError("FASTA file must contain at least two sequence records.")
    return clean_sequence(records[0]), clean_sequence(records[1])


def get_sequences():
    sample_s, sample_t = "PRTEINS", "PRTWPSEIN"
    if len(sys.argv) > 1 and Path(sys.argv[1]).is_file():
        return (*read_fasta_file(sys.argv[1]), "file")

    print("=" * 96)
    print("GLOBAL ALIGNMENT WITH AFFINE GAP PENALTY (ROSALIND: GAFF)")
    print("=" * 96)
    print("  [1] Rosalind Sample Case")
    print("  [2] Custom Case")
    print("  [3] Sample + Custom Case")
    print("  [4] Exit")
    print("-" * 96)
    choice = input("  Select an option [1-4]: ").strip()
    if choice == "4":
        return None
    if choice == "1":
        return sample_s, sample_t, "sample"
    if choice == "2":
        s = clean_sequence(input("  Enter Sequence 1 (s): "))
        t = clean_sequence(input("  Enter Sequence 2 (t): "))
        return s, t, "custom"
    if choice == "3":
        s = clean_sequence(input("  Enter Sequence 1 (s): "))
        t = clean_sequence(input("  Enter Sequence 2 (t): "))
        return sample_s, sample_t, "sample+custom", s, t
    raise ValueError("Please select 1, 2, 3, or 4.")


def affine_global_alignment(s, t, gap_open=11, gap_extend=1):
    """Three-state global DP: M=pair, X=gap in t, Y=gap in s."""
    n, m = len(s), len(t)
    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    X = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    # Traceback predecessor state: 0=M, 1=X, 2=Y.
    pM = [[None] * (m + 1) for _ in range(n + 1)]
    pX = [[None] * (m + 1) for _ in range(n + 1)]
    pY = [[None] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = -gap_open - (i - 1) * gap_extend
        pX[i][0] = 0 if i == 1 else 1
    for j in range(1, m + 1):
        Y[0][j] = -gap_open - (j - 1) * gap_extend
        pY[0][j] = 0 if j == 1 else 2

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            prev = (M[i-1][j-1], X[i-1][j-1], Y[i-1][j-1])
            state = max(range(3), key=lambda k: prev[k])
            M[i][j] = prev[state] + BLOSUM62[s[i-1]][t[j-1]]
            pM[i][j] = state

            open_x = M[i-1][j] - gap_open
            extend_x = X[i-1][j] - gap_extend
            if open_x >= extend_x:
                X[i][j], pX[i][j] = open_x, 0
            else:
                X[i][j], pX[i][j] = extend_x, 1

            open_y = M[i][j-1] - gap_open
            extend_y = Y[i][j-1] - gap_extend
            if open_y >= extend_y:
                Y[i][j], pY[i][j] = open_y, 0
            else:
                Y[i][j], pY[i][j] = extend_y, 2

    finals = (M[n][m], X[n][m], Y[n][m])
    state = max(range(3), key=lambda k: finals[k])
    score = finals[state]
    i, j = n, m
    a_s, a_t, path = [], [], []
    while i > 0 or j > 0:
        path.append((i, j, state))
        if state == 0:
            if i == 0 or j == 0:
                raise RuntimeError("Invalid diagonal traceback at matrix boundary.")
            prev_state = pM[i][j]
            a_s.append(s[i-1]); a_t.append(t[j-1])
            i -= 1; j -= 1; state = prev_state
        elif state == 1:
            if i == 0:
                raise RuntimeError("Invalid vertical traceback at matrix boundary.")
            prev_state = pX[i][j]
            a_s.append(s[i-1]); a_t.append("-")
            i -= 1; state = prev_state
        else:
            if j == 0:
                raise RuntimeError("Invalid horizontal traceback at matrix boundary.")
            prev_state = pY[i][j]
            a_s.append("-"); a_t.append(t[j-1])
            j -= 1; state = prev_state
    path.append((0, 0, state))
    aligned_s, aligned_t = "".join(reversed(a_s)), "".join(reversed(a_t))
    return {"score": score, "aligned_s": aligned_s, "aligned_t": aligned_t,
            "M": M, "X": X, "Y": Y, "pM": pM, "pX": pX, "pY": pY,
            "path": set(path), "gap_open": gap_open, "gap_extend": gap_extend}


def alignment_score(a, b, gap_open, gap_extend):
    total, i, gap_runs, gap_symbols = 0, 0, 0, 0
    while i < len(a):
        if a[i] != "-" and b[i] != "-":
            total += BLOSUM62[a[i]][b[i]]
            i += 1
        else:
            gap_runs += 1
            start = i
            while i < len(a) and (a[i] == "-" or b[i] == "-"):
                if a[i] == "-" and b[i] == "-":
                    raise ValueError("Alignment column cannot contain two gaps.")
                gap_symbols += 1
                i += 1
            length = i - start
            total -= gap_open + gap_extend * (length - 1)
    return total, gap_runs, gap_symbols


def active_matrix_table(s, t):
    """Render the active BLOSUM62 sub-matrix in the LOCA table style."""
    rows, cols = sorted(set(s)), sorted(set(t))
    cell_w = 4
    header = ["s \\ t"] + cols
    border = "+" + "+".join("-" * (cell_w + 2) for _ in header) + "+"
    out = [border,
           "|" + "|".join(f" {x:^{cell_w}} " for x in header) + "|",
           border]
    for a in rows:
        out.append("|" + "|".join(
            [f" {a:^{cell_w}} "] +
            [f" {BLOSUM62[a][b]:^{cell_w}d} " for b in cols]) + "|")
    out.append(border)
    return "\n".join(out)


def matrix_table(s, t, res):
    """Render M, X and Y DP matrices as aligned bordered tables."""
    n, m = len(s), len(t)
    if n > 22 or m > 22:
        return "  DP matrices omitted because at least one sequence is longer than 22 aa."

    out = []
    states = (("M", res["M"], res["pM"], 0, "paired residues"),
              ("X", res["X"], res["pX"], 1, "gap in sequence 2"),
              ("Y", res["Y"], res["pY"], 2, "gap in sequence 1"))
    path = res["path"]
    w = 7
    headers = ["i\\j"] + ["-(0)" if j == 0 else f"{t[j-1]}({j})" for j in range(m+1)]
    border = "+" + "+".join("-" * w for _ in headers) + "+"

    for name, mat, ptr, state_id, description in states:
        out.append(f"  State {name}: {description}")
        out.append(border)
        out.append("|" + "|".join(f"{h:^{w}}" for h in headers) + "|")
        out.append(border)
        for i in range(n+1):
            label = "-(0)" if i == 0 else f"{s[i-1]}({i})"
            cells = []
            for j in range(m+1):
                v = mat[i][j]
                if v <= NEG_INF // 2:
                    marker = "-INF"
                else:
                    # Movement direction for this DP state. The predecessor
                    # state (M/X/Y) remains stored separately in ptr for traceback.
                    arrow = "\\" if state_id == 0 else ("^" if state_id == 1 else "<")
                    marker = str(v) if ptr[i][j] is None else f"{v}{arrow}"
                    if (i, j, state_id) in path:
                        marker = "[" + marker + "]"
                cells.append(f"{marker:^{w}}")
            out.append("|" + "|".join([f"{label:^{w}}"] + cells) + "|")
        out.append(border)
    out.append("  Legend: \\ = diagonal (paired residues), ^ = up (gap in sequence 2), < = left (gap in sequence 1).")
    out.append("          [cell] marks a state/cell visited by the selected optimal traceback; -INF = unreachable state.")
    return "\n".join(out)


def render_report(s, t, res, case_label):
    a, b = res["aligned_s"], res["aligned_t"]
    checked, gap_runs, gap_symbols = alignment_score(a, b, res["gap_open"], res["gap_extend"])
    align_len = len(a)
    identities = sum(x == y for x, y in zip(a, b) if x != "-" and y != "-")
    paired = sum(x != "-" and y != "-" for x, y in zip(a, b))
    positives = sum(BLOSUM62[x][y] > 0 for x, y in zip(a, b) if x != "-" and y != "-")
    midline = "".join(
        " " if x == "-" or y == "-" else
        "|" if x == y else
        ":" if BLOSUM62[x][y] > 0 else
        "." if BLOSUM62[x][y] == 0 else " "
        for x, y in zip(a, b))
    pct_id = identities / align_len * 100 if align_len else 0
    pct_pos = positives / align_len * 100 if align_len else 0
    pct_gap = gap_symbols / align_len * 100 if align_len else 0

    out = []
    banner = "=" * 96
    out += [banner,
            "  GLOBAL ALIGNMENT WITH AFFINE GAP PENALTY REPORT (ROSALIND: GAFF)",
            banner,
            "\n[STEP 1] INPUT SEQUENCES & SUMMARY TELEMETRY"]
    s_disp = s if len(s) <= 65 else s[:62] + "..."
    t_disp = t if len(t) <= 65 else t[:62] + "..."
    summary_rows = [
        ("Sequence 1 (s)", s_disp),
        ("Sequence 2 (t)", t_disp),
        ("Sequence Lengths", f"|s| = {len(s)} aa,  |t| = {len(t)} aa,  DP Matrix Size = {len(s)+1} x {len(t)+1}"),
        ("Scoring Parameters", f"Matrix = BLOSUM62,  Gap Opening = {res['gap_open']},  Gap Extension = {res['gap_extend']}"),
        ("Optimal Global Score", str(res["score"])),
        ("Alignment Length", f"{align_len} columns; Paired residues = {paired}"),
        ("Identity & Positive Scores", f"Identity: {identities}/{paired} ({identities/paired*100 if paired else 0:.1f}%),  Positive-scoring pairs: {positives}/{paired} ({positives/paired*100 if paired else 0:.1f}%)"),
        ("Gap Statistics", f"Gap runs = {gap_runs},  Gap symbols = {gap_symbols},  Gap columns = {pct_gap:.1f}% of alignment"),
    ]
    tbl_border = "+" + "-" * 23 + "+" + "-" * 71 + "+"
    out.append(tbl_border)
    for label, value in summary_rows:
        out.append(f"| {label:<21s} | {value[:69]:<69s} |")
    out.append(tbl_border)

    out += ["\n[STEP 2] GENERATED BLOSUM62 SCORING SUB-MATRIX (For Amino Acids Present in s and t)",
            active_matrix_table(s, t),
            "\n[STEP 3] GENERATED AFFINE-GAP DYNAMIC PROGRAMMING MATRICES",
            "  M = residue pair | X = gap in sequence 2 | Y = gap in sequence 1",
            "  Gap of length L costs: gap opening + gap extension x (L - 1)",
            matrix_table(s, t, res),
            "\n[STEP 4] OPTIMAL GLOBAL ALIGNMENT & SCORE BREAKDOWN",
            "  Legend: '*' Identical residue   ':' Positive substitution   '.' Neutral (0)   ' ' Negative score or gap"]
    box = "+" + "-" * 94 + "+"
    out.append(box)
    chunk = 68
    pos_s = pos_t = 1
    for idx in range(0, align_len, chunk):
        aa, mm, bb = a[idx:idx+chunk], midline[idx:idx+chunk], b[idx:idx+chunk]
        count_s = sum(ch != "-" for ch in aa)
        count_t = sum(ch != "-" for ch in bb)
        end_s = pos_s + count_s - 1 if count_s else pos_s
        end_t = pos_t + count_t - 1 if count_t else pos_t
        line_s = f"Seq s ({pos_s:>3d}..{end_s:<3d}) :  {aa}"
        line_m = f"Conservation     :  {mm}"
        line_t = f"Seq t ({pos_t:>3d}..{end_t:<3d}) :  {bb}"
        out.extend([f"|  {line_s:<90s}  |", f"|  {line_m:<90s}  |", f"|  {line_t:<90s}  |"])
        pos_s += count_s
        pos_t += count_t
    out.append(box)

    if align_len <= 25:
        cw = 4
        sep = "+" + "-"*10 + "+" + "+".join(["-"*cw for _ in range(align_len)]) + "+-------+"
        out.append("\n  Column-by-Column Score Verification Table:")
        out.append(sep)
        out.append("| Column   |" + "|".join(f"{k:^{cw}d}" for k in range(1, align_len+1)) + "| Total |")
        out.append(sep)
        out.append("| s_i      |" + "|".join(f"{x:^{cw}s}" for x in a) + "|       |")
        score_marks = [
            "*" if x == y and x != "-" else
            ":" if x != "-" and y != "-" and BLOSUM62[x][y] > 0 else
            "." if x != "-" and y != "-" and BLOSUM62[x][y] == 0 else
            " "
            for x, y in zip(a, b)
        ]
        out.append("| Match    |" + "|".join(f"{x:^{cw}s}" for x in score_marks) + "|       |")
        out.append("| t_j      |" + "|".join(f"{x:^{cw}s}" for x in b) + "|       |")
        out.append(sep)
        contributions = []
        i = 0
        while i < align_len:
            x, y = a[i], b[i]
            if x != "-" and y != "-":
                contributions.append(BLOSUM62[x][y]); i += 1
            else:
                start = i
                while i < align_len and (a[i] == "-" or b[i] == "-"):
                    i += 1
                length = i - start
                contributions.append(-res["gap_open"])
                contributions.extend([-res["gap_extend"]] * (length - 1))
        total_score = sum(contributions)
        formatted_scores = "|".join(f"{v:>{cw-1}d} " if v == 0 else f"{v:+{cw-1}d} " for v in contributions)
        formatted_total = f"{total_score:4d}" if total_score == 0 else f"{total_score:+4d}"
        out.append(f"| Score    |{formatted_scores}| {formatted_total}  |")
        out.append(sep)

    out += ["\n[STEP 5] FINAL ROSALIND 'GAFF' OUTPUT", box,
            f"| {str(res['score']):<92s} |"]
    for idx in range(0, len(a), 90):
        out.append(f"| {a[idx:idx+90]:<92s} |")
    for idx in range(0, len(b), 90):
        out.append(f"| {b[idx:idx+90]:<92s} |")
    out += [box,
            "\n[STEP 6] ALIGNMENT VALIDATION",
            f"  DP score equals independently recalculated alignment score: {'YES' if checked == res['score'] else 'NO'}",
            f"  Original sequence 1 recovered after removing gaps: {'YES' if a.replace('-', '') == s else 'NO'}",
            f"  Original sequence 2 recovered after removing gaps: {'YES' if b.replace('-', '') == t else 'NO'}",
            f"  Recalculated affine-gap score: {checked}",
            "=" * 96,
            "  RUN COMPLETED SUCCESSFULLY",
            "=" * 96]
    return "\n".join(out)


def show_menu():
    print("\n" + "=" * 78)
    print("         GLOBAL ALIGNMENT WITH AFFINE GAP PENALTY (ROSALIND: GAFF)")
    print("=" * 78)
    print("  [1] Sample Case")
    print("  [2] Custom Case")
    print("  [3] Sample + Custom Case")
    print("  [4] Exit")
    print("-" * 78)


def read_custom_sequences():
    while True:
        try:
            s = clean_sequence(input("\n  Enter Sequence 1 (s): "))
            break
        except ValueError as exc:
            print(f"  Invalid Sequence 1: {exc}")
    while True:
        try:
            t = clean_sequence(input("  Enter Sequence 2 (t): "))
            break
        except ValueError as exc:
            print(f"  Invalid Sequence 2: {exc}")
    return s, t

def run_case(s, t, label):
    result = affine_global_alignment(s, t, 11, 1)
    print("\n" + render_report(s, t, result, label))


def main():
    sample_s, sample_t = "PRTEINS", "PRTWPSEIN"
    while True:
        show_menu()
        choice = input("  Select an option [1-4]: ").strip()
        if choice == "1":
            print("\n  Running the sample case...")
            print(f"  Sequence 1: {sample_s}")
            print(f"  Sequence 2: {sample_t}")
            print("  Gap opening penalty: 11")
            print("  Gap extension penalty: 1")
            run_case(sample_s, sample_t, "Rosalind sample")
            break
        elif choice == "2":
            s, t = read_custom_sequences()
            run_case(s, t, "Custom")
            break
        elif choice == "3":
            print("\n  Running the sample case first...")
            print(f"  Sequence 1: {sample_s}")
            print(f"  Sequence 2: {sample_t}")
            print("  Gap opening penalty: 11")
            print("  Gap extension penalty: 1")
            run_case(sample_s, sample_t, "Rosalind sample")
            print("\n\n  Now enter your custom case.")
            s, t = read_custom_sequences()
            run_case(s, t, "Custom")
            break
        elif choice == "4":
            print("\n  Exiting. Goodbye!")
            break
        else:
            print("\n  Invalid choice. Please enter 1, 2, 3, or 4.")

if __name__ == "__main__":
    try:
        main()
    except (KeyboardInterrupt, EOFError):
        print("\n\n  Program stopped by user.")
