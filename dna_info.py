from Bio.Seq import Seq


#counts the number and percentage of A, T, G, C nucleotides
def count_nucleotides(sequence):
    seq = Seq(sequence)
    length = len(seq)

    nucleotide_counts = {base: seq.count(base) for base in "ATGN"}
    nucleotide_percentage = {base: round((seq.count(base) / length) * 100, 2) for base in "ATGCN"}

    return nucleotide_counts, nucleotide_percentage

#gc content is calculated taking into the account ambiguous bases like "N"
def gc_content(a, t, c, g):
    total_nucleotides = a + t + c + g
    gc = round(((c + g )/ total_nucleotides) * 100, 2)

    return gc

#generates the RNA transcript of the DNA sequence
def transcribe(sequence):
    mRNA = sequence.replace("T", "U")

    return mRNA

#generates the reverse complement of the sequence
def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())

#finds the position of the motif in a DNA sequence
def motif_search(sequence,  motif):

    sequence = sequence.upper()
    motif = motif.upper()
    
    positions = []
        
    i = 0
    while True:
        i = sequence.find(motif, i)
        if i == -1:
            break
        positions.append(i + 1)
        i += 1 
        
    return positions
    

