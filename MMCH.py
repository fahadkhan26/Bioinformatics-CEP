#!/usr/bin/env python3
"""
Rosalind MMCH - Maximum Matchings and RNA Secondary Structures

Terminal-only interactive report.
Counts maximum noncrossing-unrestricted Watson-Crick base-pair matchings
for an RNA string, using A-U and C-G pairs.
"""

from math import factorial
import sys


SAMPLE_RNA = "AUGCUUC"


def count_bases(rna):
    return {base: rna.count(base) for base in "AUGC"}


def pair_matching_count(first_count, second_count):
    """
    Number of ways to form the maximum number of pairs between two
    complementary nucleotide groups:
        max(x, y)! / (max(x, y) - min(x, y))!
    """
    larger = max(first_count, second_count)
    smaller = min(first_count, second_count)
    return factorial(larger) // factorial(larger - smaller)


def validate_rna(rna):
    if not rna:
        return False, "RNA sequence cannot be empty."
    invalid = sorted(set(rna) - set("AUGC"))
    if invalid:
        return False, "Invalid RNA symbol(s): " + ", ".join(invalid)
    return True, "RNA sequence contains only valid bases (A, U, G, C)."


def print_header(title):
    print("\n" + "=" * 78)
    print(title.center(78))
    print("=" * 78)


def print_input_report(rna):
    counts = count_bases(rna)
    print_header("STEP 1 - INPUT RNA SEQUENCE")
    print(f"RNA Sequence : {rna}")
    print(f"Length       : {len(rna)} nucleotides")
    print("\nNucleotide composition:")
    for base in "AUGC":
        percentage = counts[base] / len(rna) * 100
        print(f"  {base} : {counts[base]:>3}  ({percentage:6.2f}%)")
    print(f"\nTotal bases  : {sum(counts.values())}")


def print_pair_report(label, base1, base2, count1, count2):
    larger = max(count1, count2)
    smaller = min(count1, count2)
    difference = larger - smaller
    result = pair_matching_count(count1, count2)

    print(f"\n{label} PAIRING ({base1}-{base2})")
    print("-" * 58)
    print(f"{base1} count                    : {count1}")
    print(f"{base2} count                    : {count2}")
    print(f"Maximum possible pairs       : {smaller}")
    print(f"Unpaired excess bases        : {difference}")
    print(f"Number of maximum matchings  : {larger}! / {difference}!")
    print(f"                              = {factorial(larger)} / {factorial(difference)}")
    print(f"                              = {result}")
    return result, smaller, difference


def run_case(rna, case_name):
    rna = "".join(rna.split()).upper()
    print_header(f"{case_name} - MMCH ANALYSIS")

    valid, message = validate_rna(rna)
    if not valid:
        print("\nVALIDATION ERROR")
        print(message)
        return

    print_input_report(rna)
    counts = count_bases(rna)

    print_header("STEP 2 - COMPLEMENTARY BASE COUNTS")
    print("Only Watson-Crick RNA pairs are considered:")
    print("  A pairs with U")
    print("  C pairs with G")
    print("\nPair group | First count | Second count | Maximum pairs")
    print("-" * 58)
    print(f"A-U       | {counts['A']:>11} | {counts['U']:>12} | {min(counts['A'], counts['U']):>13}")
    print(f"C-G       | {counts['C']:>11} | {counts['G']:>12} | {min(counts['C'], counts['G']):>13}")

    print_header("STEP 3 - MAXIMUM MATCHING CALCULATIONS")
    au_result, au_pairs, au_unpaired = print_pair_report(
        "ADENINE-URACIL", "A", "U", counts["A"], counts["U"]
    )
    cg_result, cg_pairs, cg_unpaired = print_pair_report(
        "CYTOSINE-GUANINE", "C", "G", counts["C"], counts["G"]
    )

    total_pairs = au_pairs + cg_pairs
    final_answer = au_result * cg_result

    print_header("STEP 4 - COMBINED RESULT")
    print(f"Maximum A-U pairs : {au_pairs}")
    print(f"Maximum C-G pairs : {cg_pairs}")
    print(f"Total pairs       : {total_pairs}")
    print("\nIndependent pair-group choices are multiplied:")
    print(f"  A-U matchings   : {au_result}")
    print(f"  C-G matchings   : {cg_result}")
    print(f"  Total matchings : {au_result} x {cg_result} = {final_answer}")

    print_header("STEP 5 - FINAL ROSALIND MMCH OUTPUT")
    print(final_answer)

    print_header("STEP 6 - VALIDATION")
    print(f"Input validation       : {message}")
    print(f"RNA length             : {len(rna)}")
    print(f"Maximum paired bases   : {total_pairs * 2}")
    print(f"Unpaired A/U bases     : {au_unpaired}")
    print(f"Unpaired C/G bases     : {cg_unpaired}")
    print("Pairing rule           : A-U and C-G only")
    print("Counting method        : factorial permutation formula")
    print("Result status          : Successfully calculated")


def get_custom_rna():
    print("\nEnter an RNA sequence using only A, U, G, C.")
    return input("RNA sequence: ").strip()


def main():
    while True:
        print_header("ROSALIND MMCH - MAXIMUM MATCHINGS")
        print("[1] Sample Case")
        print("[2] Custom Case")
        print("[3] Sample + Custom Case")
        print("[4] Exit")

        choice = input("\nSelect an option (1-4): ").strip()

        if choice == "1":
            run_case(SAMPLE_RNA, "SAMPLE CASE")
        elif choice == "2":
            custom = get_custom_rna()
            run_case(custom, "CUSTOM CASE")
        elif choice == "3":
            run_case(SAMPLE_RNA, "SAMPLE CASE")
            custom = get_custom_rna()
            run_case(custom, "CUSTOM CASE")
        elif choice == "4":
            print("\nExiting Rosalind MMCH program. Goodbye!")
            break
        else:
            print("\nInvalid option. Please enter 1, 2, 3, or 4.")

        input("\nPress Enter to return to the main menu...")


if __name__ == "__main__":
    try:
        main()
    except (EOFError, KeyboardInterrupt):
        print("\n\nProgram terminated.")
        sys.exit(0)
