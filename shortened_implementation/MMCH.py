#!/usr/bin/env python3
"""Rosalind MMCH - Maximum Matchings and RNA Secondary Structures (A-U and C-G pairs)."""
from math import factorial

SAMPLE_RNA = "AUGCUUC"


def pair_count(x, y):
    """Ways to form the maximum number of pairs between two groups: max! / (max - min)!"""
    big, small = max(x, y), min(x, y)
    return factorial(big) // factorial(big - small)


def run_case(rna, name):
    rna = "".join(rna.split()).upper()
    print(f"\n===== {name} =====")
    invalid = sorted(set(rna) - set("AUGC"))
    if not rna or invalid:
        print("Invalid input:", "empty sequence" if not rna else "bad symbol(s) " + ", ".join(invalid))
        return
    c = {b: rna.count(b) for b in "AUGC"}
    print(f"RNA sequence : {rna}  (length {len(rna)})")
    print("Base counts  : " + "  ".join(f"{b}={c[b]}" for b in "AUGC"))

    total = 1
    for x, y in ("AU", "CG"):
        big, small = max(c[x], c[y]), min(c[x], c[y])
        ways = pair_count(c[x], c[y])
        print(f"{x}-{y}: {small} pairs, {big}!/{big - small}! = {ways}")
        total *= ways
    print(f"\nMaximum matchings (A-U x C-G) = {total}")
    print(total)                                   # Rosalind output


def main():
    while True:
        print("\n[1] Sample Case  [2] Custom Case  [3] Sample + Custom  [4] Exit")
        choice = input("Select an option (1-4): ").strip()
        if choice in ("1", "3"):
            run_case(SAMPLE_RNA, "SAMPLE CASE")
        if choice in ("2", "3"):
            run_case(input("Enter RNA sequence (A, U, G, C): "), "CUSTOM CASE")
        if choice == "4":
            break
        if choice not in ("1", "2", "3"):
            print("Invalid option.")


if __name__ == "__main__":
    try:
        main()
    except (EOFError, KeyboardInterrupt):
        print("\nProgram terminated.")
