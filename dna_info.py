from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, nt_search, molecular_weight

#counts the number and percentage of A, T, G, C nucleotides
def count_nucleotides(sequence):
    seq = Seq(sequence)
    length = len(seq)

    nucleotide_counts = {base: seq.count(base) for base in "ATGN"}
    nucleotide_percentage = {base: round((seq.count(base) / length) * 100, 2) for base in "ATGCN"}

    return nucleotide_counts, nucleotide_percentage


#gc content is calculated taking into the account ambiguous bases like "N"
def gc_content(sequence):
    gc_percentage = gc_fraction(sequence.upper()) *100
    return round(gc_percentage, 2)


#generates the RNA transcript of the DNA sequence
def transcribe(sequence):
    return str(Seq(sequence.upper()).transcribe())


#generates the reverse complement of the sequence
def reverse_complement(sequence):
    return str(Seq(sequence.upper()).reverse_complement())


#finds the position of the motif in a DNA sequence
def motif_search(sequence,  motif):
    output = nt_search(str(sequence.upper(), motif.upper()))
    motif_positions = [position + 1 for position in output[1:]]
    return motif_positions


#Calculates the molecular weight of a DNA sequence
def calc_molecular_weight(sequence, type = "DNA"):
    seq = Seq(sequence.upper())
    weight = molecular_weight(seq, type = type)
    return round(weight, 2)
