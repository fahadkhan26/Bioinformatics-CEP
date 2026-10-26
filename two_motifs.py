
#!/usr/bin/env python3
"""
CS-434 Bioinformatics
Problem: Shortest Common Supersequence (SCSP)
Algorithm: Dynamic Programming + Backtracking
"""

import time
from pathlib import Path


# Report is saved beside this Python file.
REPORT_FILE = Path(__file__).resolve().parent / "motif_report.txt"


def make_line(char="-", width=64):
    return char * width


def show_title():
    print("\n" + make_line("="))
    print("       BIOINFORMATICS | SCSP SOLVER")
    print("       Shortest Common Supersequence")
    print(make_line("="))


def is_subsequence(small, big):
    """Check whether small appears in big in the same order."""
    pos = 0

    for ch in big:
        if pos < len(small) and small[pos] == ch:
            pos += 1

    return pos == len(small)


def find_shortest(s, t):
    """Find one shortest common supersequence using DP."""

    m = len(s)
    n = len(t)

    # dp[i][j] stores the shortest length for two prefixes.
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize the first row and column.
    for i in range(m + 1):
        dp[i][0] = i

    for j in range(n + 1):
        dp[0][j] = j

    # Fill the dynamic programming table.
    for i in range(1, m + 1):
        for j in range(1, n + 1):

            if s[i - 1] == t[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1

            else:
                dp[i][j] = min(
                    dp[i - 1][j] + 1,
                    dp[i][j - 1] + 1
                )

    # Backtrack to construct the sequence.
    i = m
    j = n
    result = []

    while i > 0 and j > 0:

        if s[i - 1] == t[j - 1]:
            result.append(s[i - 1])
            i -= 1
            j -= 1

        elif dp[i - 1][j] < dp[i][j - 1]:
            result.append(s[i - 1])
            i -= 1

        else:
            result.append(t[j - 1])
            j -= 1

    # Include any remaining characters.
    while i > 0:
        result.append(s[i - 1])
        i -= 1

    while j > 0:
        result.append(t[j - 1])
        j -= 1

    result.reverse()
    answer = "".join(result)

    return answer, dp[m][n]


def lcs_length(s, t):
    """Return the longest common subsequence length."""

    m = len(s)
    n = len(t)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(1, m + 1):
        for j in range(1, n + 1):

            if s[i - 1] == t[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1

            else:
                dp[i][j] = max(
                    dp[i - 1][j],
                    dp[i][j - 1]
                )

    return dp[m][n]


def create_report(s, t, answer, expected, elapsed):
    """Create a formatted report as a string."""

    first_ok = is_subsequence(s, answer)
    second_ok = is_subsequence(t, answer)

    # A shortest common supersequence has length m + n - LCS.
    minimum_ok = (
        len(answer) == len(s) + len(t) - lcs_length(s, t)
    )

    lines = [
        "",
        "=" * 64,
        "             BIOINFORMATICS - SCSP RESULT REPORT",
        "=" * 64,
        "",
        f"First sequence  : {s}",
        f"Second sequence : {t}",
        "",
        "-" * 64,
        "SHORTEST COMMON SUPERSEQUENCE",
        "-" * 64,
        f"Result          : {answer}",
        f"Result length   : {len(answer)}",
        "",
        "VERIFICATION",
        "-" * 64,
        f"First sequence preserved  : {'PASS' if first_ok else 'FAIL'}",
        f"Second sequence preserved : {'PASS' if second_ok else 'FAIL'}",
        f"Minimum length verified   : {'PASS' if minimum_ok else 'FAIL'}",
    ]

    if expected is not None:
        lines.append(
            f"Matches sample output      : "
            f"{'YES' if answer == expected else 'NO'}"
        )
        lines.append(
            "Note: another shortest valid sequence may also exist."
        )

    lines.extend([
        "",
        "ALGORITHM INFORMATION",
        "-" * 64,
        "Algorithm        : Dynamic Programming + Backtracking",
        "Time complexity  : O(m * n)",
        "Space complexity : O(m * n)",
        f"Execution time   : {elapsed:.6f} seconds",
        "=" * 64,
        ""
    ])

    return "\n".join(lines)


def save_report(s, t, answer, expected, elapsed):
    """Append the latest result to motif_report.txt."""

    report = create_report(s, t, answer, expected, elapsed)

    try:
        with open(REPORT_FILE, "w", encoding="utf-8") as file:
            file.write(report)

        print(f"\nReport saved successfully to:\n{REPORT_FILE}")

    except OSError as error:
        print(f"\nCould not save report: {error}")


def show_result(s, t, answer, expected=None, elapsed=0):
    """Display a formatted result in the terminal."""

    print("\n" + make_line("="))
    print("                 RESULT REPORT")
    print(make_line("="))

    print(f"\n  First sequence  : {s}")
    print(f"  Second sequence : {t}")

    print("\n" + make_line())
    print("  SHORTEST COMMON SUPERSEQUENCE")
    print(make_line())

    print(f"\n  {answer}\n")

    print(f"  Length of first sequence  : {len(s)}")
    print(f"  Length of second sequence : {len(t)}")
    print(f"  Length of result          : {len(answer)}")

    first_ok = is_subsequence(s, answer)
    second_ok = is_subsequence(t, answer)

    minimum_ok = (
        len(answer) == len(s) + len(t) - lcs_length(s, t)
    )

    print("\n" + make_line())
    print("  VERIFICATION")
    print(make_line())

    print(f"  First sequence preserved  : {'PASS' if first_ok else 'FAIL'}")
    print(f"  Second sequence preserved : {'PASS' if second_ok else 'FAIL'}")
    print(f"  Minimum length verified   : {'PASS' if minimum_ok else 'FAIL'}")

    if expected is not None:
        print(
            f"  Matches sample output     : "
            f"{'YES' if answer == expected else 'NO'}"
        )
        print("  Note: another shortest answer may also be valid.")

    print("\n" + make_line())
    print("  ALGORITHM INFORMATION")
    print(make_line())

    print("  Algorithm        : Dynamic Programming + Backtracking")
    print("  Time complexity  : O(m * n)")
    print("  Space complexity : O(m * n)")
    print(f"  Execution time   : {elapsed:.6f} seconds")

    print("\n" + make_line("="))


def run_solver(s, t, expected=None):
    """Solve, display, and save a test result."""

    start = time.perf_counter()
    answer, _ = find_shortest(s, t)
    elapsed = time.perf_counter() - start

    show_result(s, t, answer, expected, elapsed)

    # This call ensures every test result is saved.
    save_report(s, t, answer, expected, elapsed)


def get_custom_sequences():
    """Read and validate custom DNA sequences."""

    print("\n" + make_line())
    print("  CUSTOM SEQUENCE INPUT")
    print(make_line())

    s = input("  Enter first sequence : ").strip().upper()
    t = input("  Enter second sequence: ").strip().upper()

    if not s or not t:
        print("\n  Error: Both sequences must be non-empty.")
        return None, None

    if any(ch.isspace() for ch in s + t):
        print("\n  Error: Do not include spaces inside sequences.")
        return None, None

    return s, t


def main():
    # Sample case from Rosalind SCSP.
    sample_s = "ATCTGAT"
    sample_t = "TGCATA"
    sample_answer = "ATGCATGAT"

    while True:
        show_title()

        print("\n  Choose an option:\n")
        print("  1. Run Rosalind sample test case")
        print("  2. Enter custom sequences")
        print("  3. Run sample and custom test cases")
        print("  4. Exit")

        choice = input("\n  Your choice: ").strip()

        if choice == "1":
            run_solver(sample_s, sample_t, sample_answer)

        elif choice == "2":
            s, t = get_custom_sequences()

            if s is not None:
                run_solver(s, t)

        elif choice == "3":
            print("\n  Running sample test case...")
            run_solver(sample_s, sample_t, sample_answer)

            print("\n  Enter a custom test case.")
            s, t = get_custom_sequences()

            if s is not None:
                run_solver(s, t)

        elif choice == "4":
            print("\n  Thank you for using the SCSP Solver!")
            print(f"  Report location: {REPORT_FILE}")
            print(make_line("="))
            break

        else:
            print("\n  Invalid choice. Select 1, 2, 3, or 4.")

        input("\n  Press Enter to return to the menu...")


if __name__ == "__main__":
    main()
