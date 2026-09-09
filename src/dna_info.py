import math 
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, nt_search, molecular_weight

def validate_sequence(sequence):
    """
    Validates a DNA sequence string to ensure it is not empty and 
    contains only valid nucleotide characters (A, T, G, C, N).

    Args:
        sequence (str): The raw DNA sequence string.

    Returns:
        tuple: (bool, str) indicating whether the sequence is valid and a descriptive message.
    """

    s = sequence.upper().strip()

    if len(s) == 0:
        return False, "Sequence is empty."
    
    valid_bases = set("ATGCN")
    sequence_set = set(s)

    invalid = sequence_set - valid_bases

    if invalid:
        return False, f"Contains invalid characters: {', '.join(invalid)}"
    return True, "Valid sequence."


def count_nucleotides(sequence):
    """
    Counts the absolute number and percentage of each nucleotide (A, T, G, C, N) 
    in a given DNA sequence.

    Args:
        sequence (str): The DNA sequence string.

    Returns:
        tuple: Two dictionaries containing the absolute counts and percentage breakdown respectively.
    """

    seq = Seq(sequence.upper().strip())
    length = len(seq)

    if length == 0:
        empty_dict = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
        return empty_dict, empty_dict
    
    nucleotide_counts = {base: seq.count(base) for base in "ATGCN"}
    nucleotide_percentage = {base: round((seq.count(base) / length) * 100, 2) for base in "ATGCN"}

    return nucleotide_counts, nucleotide_percentage


def gc_content(sequence):
    """
    Calculates the GC content percentage of a DNA sequence, 
    accounting for ambiguous bases like 'N'.

    Args:
        sequence (str): The DNA sequence string.

    Returns:
        float: The GC percentage rounded to 2 decimal places.
    """
    gc_percentage = gc_fraction(sequence.upper()) *100
    return round(gc_percentage, 2)


def transcribe(sequence):
    """
    Generates the messenger RNA (mRNA) transcript of a DNA sequence 
    by replacing thymine (T) with uracil (U).

    Args:
        sequence (str): The DNA sequence string.

    Returns:
        str: The transcribed RNA sequence string.
    """
    return str(Seq(sequence.upper()).transcribe())


def reverse_complement(sequence):
    """
    Generates the reverse complement of a DNA sequence.

    Args:
        sequence (str): The DNA sequence string.

    Returns:
        str: The reverse complement DNA sequence string.
    """
    return str(Seq(sequence.upper()).reverse_complement())


def motif_search(sequence,  motif):
    """
    Searches for the starting positions of a given motif within a DNA sequence.

    Args:
        sequence (str): The target DNA sequence string.
        motif (str): The nucleotide pattern/motif to search for.

    Returns:
        tuple: A list of integer start positions and the total match count.
    """

    s_clean = str(sequence).upper().strip()
    m_clean = str(motif).upper().strip()

    if not s_clean or not m_clean or len(m_clean) > len(s_clean):
        return [], 0

    try:
        output = nt_search(s_clean, m_clean)
        
        if len(output) > 1:
            positions = output[1:]
            count = len(positions)
            return positions, count
        else:
            return [], 0
            
    except Exception:
        return [], 0

    
def validate_motif (motif):
    """
    Validates a search motif string to ensure it is not empty and 
    contains only valid nucleotide characters.

    Args:
        motif (str): The motif string to validate.

    Returns:
        tuple: (bool, str) indicating whether the motif is valid and a descriptive message.
    """
    
    m = motif.upper().strip()

    if len(m) == 0:
        return False, "Motif is empty"
    
    valid_bases = set("ATGCN")
    motif_set = set(m)

    invalid = motif_set - valid_bases

    if invalid:
        return False, f"Contains invalid characters: {', '.join(invalid)}"
    return True, "Valid motif."


def calc_molecular_weight(sequence, seq_type = "DNA"):
    """
    Calculates the molecular weight of a nucleic acid sequence in Daltons (Da).

    Args:
        sequence (str): The nucleotide sequence string.
        seq_type (str): The type of sequence (default is "DNA").

    Returns:
        float: The molecular weight rounded to 2 decimal places.
    """

    seq = Seq(sequence).upper().strip()

    if not seq:
        return 0.0 
    
    clean_seq = seq.replace("N", "A")
    weight = molecular_weight(clean_seq, seq_type = seq_type)
    return round(weight, 2)

def cal_entropy (sequence):
    """
    Calculates the Shannon sequence entropy to measure the randomness 
    and complexity of nucleotide distribution.

    Args:
        sequence (str): The DNA sequence string.

    Returns:
        float: The entropy value rounded to 3 decimal places.
    """

    seq = sequence.upper()
    tot_length = len(seq)

    if tot_length == 0:
        return 0.0 
    
    entropy = 0.0

    for base in "ATGC":
        count = seq.count(base)
        if count > 0:
            probability = count / tot_length
            entropy -= probability * math.log2(probability)

    return round(entropy, 3)


def cal_gc_skew (sequence):
    """
    Calculates the GC skew ((G - C) / (G + C)) of a sequence, 
    useful for locating leading and lagging replication strands.

    Args:
        sequence (str): The DNA sequence string.

    Returns:
        float: The GC skew value rounded to 4 decimal places.
    """

    seq = sequence.upper()
    g_count = seq.count("G")
    c_count = seq.count("C")

    total_gc = g_count + c_count

    if total_gc == 0:
        return 0.0
    
    skew = (g_count - c_count) / total_gc
    return round(skew, 4)