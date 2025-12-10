
sequence = "ATGCGTAAACGTTNAAGA"
sequence= sequence.strip().upper()

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
    gc_content = round(((c + g )/ total_nucleotides) * 100, 2)

    return gc_content

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
        "G": "C"}

    reversed_seq = reversed(sequence)
    base_pairs = []

    for nucleotide in reversed_seq:
        base_pairs = complement[nucleotide]
        base_pairs.append(base_pairs)

    reverse_complement = "".join(base_pairs)
    
    return reverse_complement

    

