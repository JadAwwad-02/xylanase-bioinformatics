"""mytools.py — my own sequence functions.

Written during Lessons 1–3. Import them in any notebook in this folder with:

    from mytools import complement, reverse_complement, my_translate
"""

from Bio.Data import CodonTable

# The standard genetic code (table 1).
# table.forward_table = dict of 61 codons -> amino acid
# table.stop_codons   = list of the 3 stop codons
table = CodonTable.unambiguous_dna_by_id[1]


def complement(seq):
    """Complement of a DNA sequence, NOT reversed."""
    seq = seq.upper()
    pairs = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    comp = ""
    for base in seq:
        comp = comp + pairs[base]
    return comp


def reverse_complement(seq):
    """Reverse complement: complement first, then reverse."""
    return complement(seq)[::-1]


def base_counts(seq):
    """How many times each base appears."""
    seq = seq.upper()
    counts = {"A": seq.count("A"),
              "T": seq.count("T"),
              "G": seq.count("G"),
              "C": seq.count("C")}
    return counts


def is_palindrome(seq):
    """True if the sequence equals its own reverse complement."""
    seq = seq.upper()
    return seq == reverse_complement(seq)


def find_residue(seq, target):
    """Positions (counting from 1) where target appears in seq."""
    positions = []
    for position, aa in enumerate(seq, start=1):
        if aa == target:
            positions.append(position)
    return positions


def my_translate(seq):
    """Translate DNA to protein. Stops at a stop codon or a partial codon."""
    seq = seq.upper()
    protein = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            break
        if codon in table.stop_codons:
            break
        protein = protein + table.forward_table[codon]
    return protein

def six_frames(seq):
    seq = seq.upper()
    rc = reverse_complement(seq)
    frames = []
    for i in range(3):
        frames.append(my_translate(seq[i:]))
    for i in range(3):
        frames.append(my_translate(rc[i:]))
    return frames