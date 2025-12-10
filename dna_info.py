
sequence = "ATGCGTAAACGTTAAGA"
def count_nucleotides(sequence):
    A = sequence.count("A")
    T = sequence.count("T")
    C = sequence.count("C")
    G = sequence.count("G")
    
    return A, T, C, G



print(count_nucleotides(sequence))

def gc_content(sequence):
    ...

def transcribe(sequence):
    ...

def reverse_compliment(sequence):
    ...

