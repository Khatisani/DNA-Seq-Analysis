import math 
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
def calc_molecular_weight(sequence, seq_type = "DNA"):
    seq = Seq(sequence.upper())
    weight = molecular_weight(seq, seq_type = seq_type)
    return round(weight, 2)

#Calculates the entropy of a sequence to measure randmness
def cal_entropy (sequence):
    s = sequence.upper()
    tot_length = len(s)

    if tot_length == 0:
        return 0.0 
    
    entropy = 0.0

    for base in "ATGC":
        count = s.count(base)
        if count > 0:
            probability = count / tot_length
            entropy -= probability * math.log2(probability)

    return round(entropy, 3)

#Calculates the gc skew of a sequence. 
def cal_gc_skew (sequence):

    s = sequence.upper()
    g_count = s.count("G")
    c_count = s.count("C")

    total_gc = g_count + c_count

    if total_gc == 0:
        return 0.0
    
    skew = (g_count - c_count) / total_gc
    return round(skew, 4)