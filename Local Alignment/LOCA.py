#!/usr/bin/env python3
"""
CS-434 Bioinformatics — Problem Based Learning (PBL)
Topic 1: Alignment | Local Alignment with Scoring Matrix (Rosalind ID: LOCA)
Algorithm: Smith-Waterman Local Alignment (Interactive Input, Matrix Generation & Traceback)
"""

import sys
from pathlib import Path

# Standard PAM250 Evolutionary Substitution Matrix (20 Standard Amino Acids)
PAM250_DATA = {
    "A": {"A":  2, "C": -2, "D":  0, "E":  0, "F": -3, "G":  1, "H": -1, "I": -1, "K": -1, "L": -2, "M": -1, "N":  0, "P":  1, "Q":  0, "R": -2, "S":  1, "T":  1, "V":  0, "W": -6, "Y": -3},
    "C": {"A": -2, "C": 12, "D": -5, "E": -5, "F": -4, "G": -3, "H": -3, "I": -2, "K": -5, "L": -6, "M": -5, "N": -4, "P": -3, "Q": -5, "R": -4, "S":  0, "T": -2, "V": -2, "W": -8, "Y":  0},
    "D": {"A":  0, "C": -5, "D":  4, "E":  3, "F": -6, "G":  1, "H":  1, "I": -2, "K":  0, "L": -4, "M": -3, "N":  2, "P": -1, "Q":  2, "R": -1, "S":  0, "T":  0, "V": -2, "W": -7, "Y": -4},
    "E": {"A":  0, "C": -5, "D":  3, "E":  4, "F": -5, "G":  0, "H":  1, "I": -2, "K":  0, "L": -3, "M": -2, "N":  1, "P": -1, "Q":  2, "R": -1, "S":  0, "T":  0, "V": -2, "W": -7, "Y": -4},
    "F": {"A": -3, "C": -4, "D": -6, "E": -5, "F":  9, "G": -5, "H": -2, "I":  1, "K": -5, "L":  2, "M":  0, "N": -3, "P": -5, "Q": -5, "R": -4, "S": -3, "T": -3, "V": -1, "W":  0, "Y":  7},
    "G": {"A":  1, "C": -3, "D":  1, "E":  0, "F": -5, "G":  5, "H": -2, "I": -3, "K": -2, "L": -4, "M": -3, "N":  0, "P":  0, "Q": -1, "R": -3, "S":  1, "T":  0, "V": -1, "W": -7, "Y": -5},
    "H": {"A": -1, "C": -3, "D":  1, "E":  1, "F": -2, "G": -2, "H":  6, "I": -2, "K":  0, "L": -2, "M": -2, "N":  2, "P":  0, "Q":  3, "R":  2, "S": -1, "T": -1, "V": -2, "W": -3, "Y":  0},
    "I": {"A": -1, "C": -2, "D": -2, "E": -2, "F":  1, "G": -3, "H": -2, "I":  5, "K": -2, "L":  2, "M":  2, "N": -2, "P": -2, "Q": -2, "R": -2, "S": -1, "T":  0, "V":  4, "W": -5, "Y": -1},
    "K": {"A": -1, "C": -5, "D":  0, "E":  0, "F": -5, "G": -2, "H":  0, "I": -2, "K":  5, "L": -3, "M":  0, "N":  1, "P": -1, "Q":  1, "R":  3, "S":  0, "T":  0, "V": -2, "W": -3, "Y": -4},
    "L": {"A": -2, "C": -6, "D": -4, "E": -3, "F":  2, "G": -4, "H": -2, "I":  2, "K": -3, "L":  6, "M":  4, "N": -3, "P": -3, "Q": -2, "R": -3, "S": -3, "T": -2, "V":  2, "W": -2, "Y": -1},
    "M": {"A": -1, "C": -5, "D": -3, "E": -2, "F":  0, "G": -3, "H": -2, "I":  2, "K":  0, "L":  4, "M":  6, "N": -2, "P": -2, "Q": -1, "R":  0, "S": -2, "T": -1, "V":  2, "W": -4, "Y": -2},
    "N": {"A":  0, "C": -4, "D":  2, "E":  1, "F": -3, "G":  0, "H":  2, "I": -2, "K":  1, "L": -3, "M": -2, "N":  2, "P":  0, "Q":  1, "R":  0, "S":  1, "T":  0, "V": -2, "W": -4, "Y": -2},
    "P": {"A":  1, "C": -3, "D": -1, "E": -1, "F": -5, "G":  0, "H":  0, "I": -2, "K": -1, "L": -3, "M": -2, "N":  0, "P":  6, "Q":  0, "R":  0, "S":  1, "T":  0, "V": -1, "W": -6, "Y": -5},
    "Q": {"A":  0, "C": -5, "D":  2, "E":  2, "F": -5, "G": -1, "H":  3, "I": -2, "K":  1, "L": -2, "M": -1, "N":  1, "P":  0, "Q":  4, "R":  1, "S": -1, "T": -1, "V": -2, "W": -5, "Y": -4},
    "R": {"A": -2, "C": -4, "D": -1, "E": -1, "F": -4, "G": -3, "H":  2, "I": -2, "K":  3, "L": -3, "M":  0, "N":  0, "P":  0, "Q":  1, "R":  6, "S":  0, "T": -1, "V": -2, "W":  2, "Y": -4},
    "S": {"A":  1, "C":  0, "D":  0, "E":  0, "F": -3, "G":  1, "H": -1, "I": -1, "K":  0, "L": -3, "M": -2, "N":  1, "P":  1, "Q": -1, "R":  0, "S":  2, "T":  1, "V": -1, "W": -2, "Y": -3},
    "T": {"A":  1, "C": -2, "D":  0, "E":  0, "F": -3, "G":  0, "H": -1, "I":  0, "K":  0, "L": -2, "M": -1, "N":  0, "P":  0, "Q": -1, "R": -1, "S":  1, "T":  3, "V":  0, "W": -5, "Y": -3},
    "V": {"A":  0, "C": -2, "D": -2, "E": -2, "F": -1, "G": -1, "H": -2, "I":  4, "K": -2, "L":  2, "M":  2, "N": -2, "P": -1, "Q": -2, "R": -2, "S": -1, "T":  0, "V":  4, "W": -6, "Y": -2},
    "W": {"A": -6, "C": -8, "D": -7, "E": -7, "F":  0, "G": -7, "H": -3, "I": -5, "K": -3, "L": -2, "M": -4, "N": -4, "P": -6, "Q": -5, "R":  2, "S": -2, "T": -5, "V": -6, "W": 17, "Y":  0},
    "Y": {"A": -3, "C":  0, "D": -4, "E": -4, "F":  7, "G": -5, "H":  0, "I": -1, "K": -4, "L": -1, "M": -2, "N": -2, "P": -5, "Q": -4, "R": -4, "S": -3, "T": -3, "V": -2, "W":  0, "Y": 10},
}


def clean_sequence(raw_seq: str) -> str:
    """Removes FASTA headers, whitespace, and validates amino acid characters."""
    lines = [
        line.strip()
        for line in raw_seq.strip().splitlines()
        if not line.strip().startswith(">")
    ]
    cleaned = "".join("".join(lines).split()).upper()
    invalid = [ch for ch in cleaned if ch not in PAM250_DATA]
    if invalid:
        raise ValueError(
            f"Invalid amino acid symbol(s): {sorted(set(invalid))}. "
            f"Allowed symbols: {sorted(PAM250_DATA.keys())}"
        )
    return cleaned


def get_user_sequences() -> tuple[str, str, int]:
    """
    Prompts the user interactively for Sequence 1, Sequence 2, and Gap Penalty.
    Also supports passing a FASTA file path via command-line argument.
    """
    if len(sys.argv) > 1 and Path(sys.argv[1]).is_file():
        content = Path(sys.argv[1]).read_text(encoding="utf-8")
        blocks = [b.strip() for b in content.split(">") if b.strip()]
        if len(blocks) >= 2:
            seq1 = clean_sequence(">" + blocks[0])
            seq2 = clean_sequence(">" + blocks[1])
            return seq1, seq2, 5

    print("=" * 96)
    print("  SMITH-WATERMAN LOCAL ALIGNMENT GENERATOR (ROSALIND: LOCA)")
    print("  Enter your two protein sequences below (or press [Enter] to run the Rosalind sample).")
    print("=" * 96)

    try:
        raw_s = input("  Enter Sequence 1 (s) [default: MEANLYPRTEINSTRING] : ").strip()
    except EOFError:
        raw_s = ""
    seq1 = clean_sequence(raw_s) if raw_s else "MEANLYPRTEINSTRING"

    try:
        raw_t = input("  Enter Sequence 2 (t) [default: PLEASANTLYEINSTEIN] : ").strip()
    except EOFError:
        raw_t = ""
    seq2 = clean_sequence(raw_t) if raw_t else "PLEASANTLYEINSTEIN"

    try:
        raw_gap = input("  Enter Linear Gap Penalty sigma [default: 5]        : ").strip()
    except EOFError:
        raw_gap = ""
    gap_penalty = int(raw_gap) if raw_gap.isdigit() else 5

    print()
    return seq1, seq2, gap_penalty


def generate_local_alignment(s: str, t: str, gap_penalty: int = 5) -> dict:
    """Build the Smith-Waterman matrices and find every distinct best alignment."""
    n, m = len(s), len(t)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    # Each cell stores every direction tied for its best positive score.
    # 1 = diagonal, 2 = up, 3 = left.
    T = [[[] for _ in range(m + 1)] for _ in range(n + 1)]
    max_score = 0
    peak_cells = []

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = H[i - 1][j - 1] + PAM250_DATA[s[i - 1]][t[j - 1]]
            up = H[i - 1][j] - gap_penalty
            left = H[i][j - 1] - gap_penalty
            best = max(0, diag, up, left)
            H[i][j] = best

            if best > 0:
                if diag == best:
                    T[i][j].append(1)
                if up == best:
                    T[i][j].append(2)
                if left == best:
                    T[i][j].append(3)

            if best > max_score:
                max_score = best
                peak_cells = [(i, j)]
            elif best == max_score and best > 0:
                peak_cells.append((i, j))

    all_alignments = []
    seen = set()

    def traceback(i, j, aligned_s_rev, aligned_t_rev, path, end_i, end_j):
        if H[i][j] == 0:
            aligned_s = "".join(reversed(aligned_s_rev))
            aligned_t = "".join(reversed(aligned_t_rev))
            start_i = i + 1
            start_j = j + 1
            key = (start_i, end_i, start_j, end_j, aligned_s, aligned_t)
            if key not in seen:
                seen.add(key)
                all_alignments.append({
                    "start_i": start_i,
                    "end_i": end_i,
                    "start_j": start_j,
                    "end_j": end_j,
                    "aligned_s": aligned_s,
                    "aligned_t": aligned_t,
                    "motif_s": aligned_s.replace("-", ""),
                    "motif_t": aligned_t.replace("-", ""),
                    "traceback_cells": set(path),
                    "score": max_score,
                })
            return

        for direction in T[i][j]:
            if direction == 1:
                traceback(i - 1, j - 1,
                          aligned_s_rev + [s[i - 1]],
                          aligned_t_rev + [t[j - 1]],
                          path + [(i - 1, j - 1)], end_i, end_j)
            elif direction == 2:
                traceback(i - 1, j,
                          aligned_s_rev + [s[i - 1]],
                          aligned_t_rev + ["-"],
                          path + [(i - 1, j)], end_i, end_j)
            elif direction == 3:
                traceback(i, j - 1,
                          aligned_s_rev + ["-"],
                          aligned_t_rev + [t[j - 1]],
                          path + [(i, j - 1)], end_i, end_j)

    for peak_i, peak_j in peak_cells:
        traceback(peak_i, peak_j, [], [], [(peak_i, peak_j)], peak_i, peak_j)

    if not all_alignments:
        all_alignments.append({
            "start_i": 1, "end_i": 0, "start_j": 1, "end_j": 0,
            "aligned_s": "", "aligned_t": "", "motif_s": "", "motif_t": "",
            "traceback_cells": set(), "score": 0,
        })

    # Keep the first alignment in the original result fields for compatibility.
    result = dict(all_alignments[0])
    result.update({
        "H": H,
        "T": T,
        "max_score": max_score,
        "gap_penalty": gap_penalty,
        "all_alignments": all_alignments,
        "traceback_cells": set().union(
            *(alignment["traceback_cells"] for alignment in all_alignments)
        ),
    })
    return result


def render_active_pam250_table(s: str, t: str) -> str:
    """Generates the PAM250 scoring sub-matrix for the unique amino acids present in s and t."""
    unique_s = sorted(set(s))
    unique_t = sorted(set(t))

    col_w = 4
    header_cells = [" s \\ t "] + [f"{aa:^{col_w}s}" for aa in unique_t]
    header_line = "|" + "|".join(header_cells) + "|"
    border_line = "+" + "+".join(["-" * len(c) for c in header_cells]) + "+"

    rows = [border_line, header_line, border_line]
    for aa_s in unique_s:
        row_cells = [f"   {aa_s}   "] + [
            f"{PAM250_DATA[aa_s][aa_t]:>{col_w - 1}d} " for aa_t in unique_t
        ]
        rows.append("|" + "|".join(row_cells) + "|")
    rows.append(border_line)
    return "\n".join(rows)


def render_dp_matrix_table(s: str, t: str, res: dict) -> str:
    """
    Renders the generated (n+1) x (m+1) Smith-Waterman Scoring Matrix H.
    Cells along the optimal local alignment traceback path are bracketed [val].
    """
    H = res["H"]
    T = res["T"]
    path = res["traceback_cells"]
    n, m = len(s), len(t)

    # Display full matrix for sequences <= 22 aa; otherwise show 15x15 window ending at the peak
    r_end = n if n <= 22 else res["end_i"]
    r_start = 0 if n <= 22 else max(0, r_end - 14)
    c_end = m if m <= 22 else res["end_j"]
    c_start = 0 if m <= 22 else max(0, c_end - 14)

    dir_symbols = {0: " ", 1: "\\", 2: "^", 3: "<"}
    cell_w = 6

    res_header = ["  i\\j  "]
    for j in range(c_start, c_end + 1):
        aa_t = "-" if j == 0 else t[j - 1]
        res_header.append(f"{aa_t + '(' + str(j) + ')':^{cell_w}s}")

    border = "+" + "+".join(["-" * len(c) for c in res_header]) + "+"
    lines = [border, "|" + "|".join(res_header) + "|", border]

    for i in range(r_start, r_end + 1):
        aa_s = "-" if i == 0 else s[i - 1]
        row_label = f" {aa_s}({i:>2d}) "
        row_cells = [row_label]
        for j in range(c_start, c_end + 1):
            val = H[i][j]
            d_sym = "".join(dir_symbols[d] for d in T[i][j]) if val > 0 else " "
            if (i, j) in path:
                cell_str = f"[{val:>2d}{d_sym}]"
                row_cells.append(f"{cell_str:^{cell_w}s}")
            else:
                cell_str = f" {val:>2d}{d_sym} "
                row_cells.append(f"{cell_str:^{cell_w}s}")
        lines.append("|" + "|".join(row_cells) + "|")

    lines.append(border)
    return "\n".join(lines)


def display_complete_results(s: str, t: str, res: dict) -> str:
    """Formats the complete 5-step output neatly for terminal and lab report."""
    g_s = res["aligned_s"]
    g_t = res["aligned_t"]
    gap_p = res["gap_penalty"]
    align_len = len(g_s)

    midline = []
    col_scores = []
    identities = 0
    positives = 0
    gaps = 0

    for a, b in zip(g_s, g_t):
        if a == "-" or b == "-":
            midline.append(" ")
            col_scores.append(-gap_p)
            gaps += 1
        else:
            sc = PAM250_DATA[a][b]
            col_scores.append(sc)
            if a == b:
                midline.append("|")
                identities += 1
                positives += 1
            elif sc > 0:
                midline.append(":")
                positives += 1
            elif sc == 0:
                midline.append(".")
            else:
                midline.append(" ")

    mid_str = "".join(midline)
    pct_id = (identities / align_len * 100.0) if align_len else 0.0
    pct_pos = (positives / align_len * 100.0) if align_len else 0.0
    pct_gap = (gaps / align_len * 100.0) if align_len else 0.0

    out = []
    banner = "=" * 96

    out.append(banner)
    out.append("  CS-434 BIOINFORMATICS | SMITH-WATERMAN LOCAL ALIGNMENT REPORT (ROSALIND: LOCA)")
    out.append(banner)

    # SECTION 1: Input Sequences & Summary Table
    out.append("\n[STEP 1] INPUT SEQUENCES & SUMMARY TELEMETRY")
    s_disp = s if len(s) <= 65 else s[:62] + "..."
    t_disp = t if len(t) <= 65 else t[:62] + "..."
    r_disp = f"s[{res['start_i']}..{res['end_i']}] = {res['motif_s']} ({len(res['motif_s'])} aa)"
    u_disp = f"t[{res['start_j']}..{res['end_j']}] = {res['motif_t']} ({len(res['motif_t'])} aa)"
    if len(r_disp) > 69:
        r_disp = r_disp[:66] + "..."
    if len(u_disp) > 69:
        u_disp = u_disp[:66] + "..."

    summary_rows = [
        ("Sequence 1 (s)", s_disp),
        ("Sequence 2 (t)", t_disp),
        ("Sequence Lengths", f"|s| = {len(s)} aa,  |t| = {len(t)} aa,  DP Matrix Size = {len(s)+1} x {len(t)+1}"),
        ("Scoring Parameters", f"Matrix = PAM250,  Linear Gap Penalty (sigma) = {gap_p}"),
        ("Optimal Local Score", f"{res['max_score']}"),
        ("Optimal Peak -> Zero", f"Peak H[{res['end_i']}][{res['end_j']}] = {res['max_score']}  --->  Origin H[{res['start_i']-1}][{res['start_j']-1}] = 0"),
        ("Aligned Substring r", r_disp),
        ("Aligned Substring u", u_disp),
        ("Identity & Similarity", f"Identity: {identities}/{align_len} ({pct_id:.1f}%),  Similarity: {positives}/{align_len} ({pct_pos:.1f}%),  Gaps: {gaps}/{align_len} ({pct_gap:.1f}%)"),
    ]

    tbl_border = "+" + "-" * 23 + "+" + "-" * 71 + "+"
    out.append(tbl_border)
    for label, val in summary_rows:
        out.append(f"| {label:<21s} | {val:<69s} |")
    out.append(tbl_border)

    # SECTION 2: Active PAM250 Lookup Sub-Matrix
    out.append("\n[STEP 2] GENERATED PAM250 SCORING SUB-MATRIX (For Amino Acids Present in s and t)")
    out.append(render_active_pam250_table(s, t))

    # SECTION 3: Dynamically Generated Smith-Waterman DP Matrix
    out.append("\n[STEP 3] GENERATED SMITH-WATERMAN DP MATRIX H[i][j]")
    out.append("  Legend: [val] = Optimal Traceback Path | '\\' = Diag (Match/Sub) | '^' = Up (Del) | '<' = Left (Ins)")
    out.append(render_dp_matrix_table(s, t, res))

    # SECTION 4: Visual Pairwise Local Alignment & Per-Column Score Breakdown
    out.append("\n[STEP 4] OPTIMAL LOCAL ALIGNMENT & SCORE BREAKDOWN")
    out.append("  Legend: '|' Identical residue   ':' Conservative (PAM250 > 0)   '.' Neutral (0)   ' ' Gap/Mismatch")
    box_border = "+" + "-" * 94 + "+"
    out.append(box_border)

    chunk_size = 68
    cur_i, cur_j = res["start_i"], res["start_j"]
    for idx in range(0, align_len, chunk_size):
        c_s = g_s[idx : idx + chunk_size]
        c_m = mid_str[idx : idx + chunk_size]
        c_t = g_t[idx : idx + chunk_size]
        ng_s = sum(1 for ch in c_s if ch != "-")
        ng_t = sum(1 for ch in c_t if ch != "-")
        end_i = cur_i + ng_s - 1 if ng_s else cur_i
        end_j = cur_j + ng_t - 1 if ng_t else cur_j

        line_s = f"Seq s ({cur_i:>3d}..{end_i:<3d}) :  {c_s}"
        line_m = f"Conservation     :  {c_m}"
        line_t = f"Seq t ({cur_j:>3d}..{end_j:<3d}) :  {c_t}"
        out.append(f"|  {line_s:<90s}  |")
        out.append(f"|  {line_m:<90s}  |")
        out.append(f"|  {line_t:<90s}  |")
        cur_i += ng_s
        cur_j += ng_t
    out.append(box_border)

    if align_len <= 25:
        col_w = 4
        hdr_idx = "| Column   |" + "|".join(f"{k:^{col_w}d}" for k in range(1, align_len + 1)) + "| Total |"
        sep_row = "+" + "-" * 10 + "+" + "+".join(["-" * col_w for _ in range(align_len)]) + "+-------+"
        row_s   = "| s_i      |" + "|".join(f"{c:^{col_w}s}" for c in g_s) + "|       |"
        row_m   = "| Match    |" + "|".join(f"{m:^{col_w}s}" for m in mid_str) + "|       |"
        row_t   = "| t_j      |" + "|".join(f"{c:^{col_w}s}" for c in g_t) + "|       |"
        row_sc  = "| Score    |" + "|".join(f"{sc:+{col_w - 1}d} " for sc in col_scores) + f"| {sum(col_scores):+4d}  |"
        out.append("\n  Column-by-Column Score Verification Table:")
        out.extend([sep_row, hdr_idx, sep_row, row_s, row_m, row_t, sep_row, row_sc, sep_row])

    # SECTION 5: Clean Rosalind Output
    out.append("\n[STEP 5] FINAL ROSALIND 'LOCA' OUTPUT (Also saved to 'loca_submission.txt')")
    out.append(box_border)
    out.append(f"| {str(res['max_score']):<92s} |")
    for idx in range(0, len(res["motif_s"]), 90):
        out.append(f"| {res['motif_s'][idx:idx+90]:<92s} |")
    for idx in range(0, len(res["motif_t"]), 90):
        out.append(f"| {res['motif_t'][idx:idx+90]:<92s} |")
    out.append(box_border)

    # SECTION 6: Every distinct optimal local alignment
    out.append("\n[STEP 6] ALL DISTINCT OPTIMAL LOCAL ALIGNMENTS")
    out.append(f"Maximum score: {res['max_score']}")
    out.append(f"Number of distinct optimal alignments: {len(res['all_alignments'])}")
    for number, item in enumerate(res["all_alignments"], 1):
        out.append(f"\nAlignment {number}")
        out.append(f"Sequence 1 positions: {item['start_i']}..{item['end_i']}")
        out.append(f"Sequence 2 positions: {item['start_j']}..{item['end_j']}")
        out.append(f"Score: {item['score']}")
        out.append(item["aligned_s"])
        out.append(item["aligned_t"])

    return "\n".join(out)


def show_menu():
    """Displays the interactive menu in the same style as the motif program."""
    print("\n" + "=" * 78)
    print("  CS-434 BIOINFORMATICS | SMITH-WATERMAN LOCAL ALIGNMENT (ROSALIND: LOCA)")
    print("=" * 78)
    print("  [1] Sample Case")
    print("  [2] Custom Case")
    print("  [3] Sample + Custom Case")
    print("  [4] Exit")
    print("-" * 78)


def read_custom_sequences():
    """Reads and validates two protein sequences and a gap penalty."""
    while True:
        raw_s = input("\n  Enter Sequence 1 (s): ").strip()
        try:
            seq_s = clean_sequence(raw_s)
            if not seq_s:
                print("  Sequence 1 cannot be empty. Please try again.")
                continue
            break
        except ValueError as error:
            print(f"  Invalid Sequence 1: {error}")

    while True:
        raw_t = input("  Enter Sequence 2 (t): ").strip()
        try:
            seq_t = clean_sequence(raw_t)
            if not seq_t:
                print("  Sequence 2 cannot be empty. Please try again.")
                continue
            break
        except ValueError as error:
            print(f"  Invalid Sequence 2: {error}")

    while True:
        raw_gap = input("  Enter gap penalty [default: 5]: ").strip()
        if not raw_gap:
            gap_pen = 5
            break
        try:
            gap_pen = int(raw_gap)
            if gap_pen < 0:
                print("  Gap penalty must be zero or a positive integer.")
                continue
            break
        except ValueError:
            print("  Please enter a whole number, such as 5.")

    return seq_s, seq_t, gap_pen


def save_and_show_results(seq_s, seq_t, gap_pen):
    """Runs the alignment, displays the report, and saves both text files."""
    alignment_data = generate_local_alignment(seq_s, seq_t, gap_penalty=gap_pen)
    full_report = display_complete_results(seq_s, seq_t, alignment_data)

    submission_lines = [
        "ALL DISTINCT OPTIMAL LOCAL ALIGNMENTS",
        f"Maximum score: {alignment_data['max_score']}",
        f"Number of alignments: {len(alignment_data['all_alignments'])}",
    ]
    for number, item in enumerate(alignment_data["all_alignments"], 1):
        submission_lines.extend([
            "",
            f"Alignment {number}",
            f"Sequence 1 positions: {item['start_i']}..{item['end_i']}",
            f"Sequence 2 positions: {item['start_j']}..{item['end_j']}",
            f"Score: {item['score']}",
            item["aligned_s"],
            item["aligned_t"],
        ])

    # Always save the files beside this Python script, not the terminal's folder.
    script_folder = Path(__file__).resolve().parent
    submission_path = script_folder / "loca_submission.txt"
    report_path = script_folder / "loca_full_report.txt"

    # Write mode replaces old results, so the files contain this run's output.
    submission_path.write_text("\n".join(submission_lines) + "\n", encoding="utf-8")
    report_path.write_text(full_report + "\n", encoding="utf-8")

    print("\n" + full_report)
    print("\n" + "=" * 78)
    print("  RUN COMPLETED SUCCESSFULLY")
    print("=" * 78)
    print(f"  Full report saved to: {report_path}")
    print(f"  All optimal alignments saved to: {submission_path}")
    print("  Previous contents of these files were overwritten.")
    print("=" * 78)


def main():
    sample_s = "MEANLYPRTEINSTRING"
    sample_t = "PLEASANTLYEINSTEIN"
    sample_gap = 5

    while True:
        show_menu()
        choice = input("  Select an option [1-4]: ").strip()

        if choice == "1":
            print("\n  Running the sample case...")
            print(f"  Sequence 1: {sample_s}")
            print(f"  Sequence 2: {sample_t}")
            print(f"  Gap penalty: {sample_gap}")
            save_and_show_results(sample_s, sample_t, sample_gap)
            break
        elif choice == "2":
            seq_s, seq_t, gap_pen = read_custom_sequences()
            save_and_show_results(seq_s, seq_t, gap_pen)
            break
        elif choice == "3":
            print("\n  Running the sample case first...")
            print(f"  Sequence 1: {sample_s}")
            print(f"  Sequence 2: {sample_t}")
            print(f"  Gap penalty: {sample_gap}")
            save_and_show_results(sample_s, sample_t, sample_gap)

            print("\n\n  Now enter your custom case.")
            seq_s, seq_t, gap_pen = read_custom_sequences()
            save_and_show_results(seq_s, seq_t, gap_pen)
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
