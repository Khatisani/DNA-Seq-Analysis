
#counts the number of A, T, G, C nucleotides
def count_nucleotides(sequence):
    a = sequence.count("A")
    t = sequence.count("T")
    c = sequence.count("C")
    g = sequence.count("G")
    
    return a, t, c, g

#gc content is calculated taking into the account ambiguous bases like "N"
def gc_content(a, t, c, g):
    total_nucleotides = a + t + c + g
    gc= round(((c + g )/ total_nucleotides) * 100, 2)

    return gc

#generates the RNA transcript of the DNA sequence
def transcribe(sequence):
    mRNA = sequence.replace("T", "U")

    return mRNA

#generates the reverse complement of the sequence
def reverse_complement(sequence):
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N" : "N"}

    reversed_seq = reversed(sequence)
    base_pairs = []

    for nucleotide in reversed_seq:
        complementary_base= complement[nucleotide]
        base_pairs.append(complementary_base)

    rev_comp = "".join(base_pairs)

    return rev_comp

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
    

