
sequence = "ATGCGTAAACGTTAAGA"
#function counts the number of A, T, G, C nucleotides
def count_nucleotides(sequence):
    a = sequence.count("A")
    t = sequence.count("T")
    c = sequence.count("C")
    g = sequence.count("G")
    
    return a, t, c, g

#gc content is calculated taking into the account ambiguous bases like "N"
def gc_content(a, t, c, g):
    total_nucleotides = a + t + c + g
    gc_content = ((c + g )/ total_nucleotides) * 100

    return gc_content

def transcribe(sequence):
    ...

def reverse_compliment(sequence):
    ...

