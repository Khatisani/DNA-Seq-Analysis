import pytest
import sys
import os
import math 

from dna_info import count_nucleotides, gc_content, transcribe, reverse_complement, motif_search, calc_molecular_weight, cal_entropy, cal_gc_skew, validate_sequence, validate_motif
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def test_valid_dna_sequence():
    is_valid, msg = validate_sequence("ATCGN")
    assert is_valid == True
    assert msg == "Valid sequence."

def test_empty_dna_sequence():
    is_valid, msg = validate_sequence("")
    assert is_valid == False
    assert msg == "Sequence is empty."

def test_invalid_dna_sequence():
    is_valid, msg = validate_sequence("XYZ")
    assert is_valid == False
    assert "Contains invalid characters" in msg

def test_valid_motif():
    is_valid, msg = validate_motif("ATG")
    assert is_valid == True

def test_invalid_motif():
    is_valid, msg = validate_motif("X!")
    assert is_valid == False

def test_count_nucleotides():
    """Test that nucleotide counting works accurately."""
    seq = "ATCGATCG"
    counts, percentages = count_nucleotides(seq)
    
    assert counts["A"] == 2
    assert counts["T"] == 2
    assert counts["C"] == 2
    assert counts["G"] == 2
    assert counts["N"] == 0
    
    assert percentages["A"] == 25.0
    assert percentages["T"] == 25.0
    assert percentages["C"] == 25.0
    assert percentages["G"] == 25.0

def test_gc_content():
    seq = "ATCG"
    gc = gc_content(seq)
    assert gc == 50.0

def test_gc_content_case_insensitivity():
    seq = "atcg"
    gc = gc_content(seq)
    assert gc == 50.0

def test_transcribe_standard():
    dna = "ATCG"
    rna = transcribe(dna)
    assert rna == "AUCG"

def test_transcribe_lowercase():
    dna = "atcg"
    rna = transcribe(dna)
    assert rna == "AUCG"

def test_transcribe_empty():
    dna = ""
    rna = transcribe(dna)
    assert rna == ""

def test_reverse_complement_standard():
    dna = "ATCG"
    rev_comp = reverse_complement(dna)
    assert rev_comp == "CGAT"

def test_reverse_complement_with_n():
    dna = "ATN"
    rev_comp = reverse_complement(dna)
    assert rev_comp == "NAT"

def test_reverse_complement_empty():
    dna = ""
    rev_comp = reverse_complement(dna)
    assert rev_comp == ""

def test_molecular_weight_standard():
    dna = "AAAA"
    mw = calc_molecular_weight(dna)
    assert mw > 0

def test_molecular_weight_empty():
    mw = calc_molecular_weight("")
    assert mw == 0.0

def test_entropy_uniform():
    dna = "ATCG"
    entropy = cal_entropy(dna)
    assert math.isclose(entropy, 2.0, rel_tol=1e-9)

def test_entropy_single_base():
    dna = "AAAA"
    entropy = cal_entropy(dna)
    assert math.isclose(entropy, 0.0, rel_tol=1e-9)

def test_entropy_empty():
    entropy = cal_entropy("")
    assert entropy == 0.0

def test_motif_search_standard():
    dna = "ATCGATCG"
    motif = "ATC"
    positions, count = motif_search(dna, motif)
    
    assert count == 2
    assert positions == [0, 4]

def test_motif_search_case_insensitive():
    dna = "atcgatcg"
    motif = "ATC"
    positions, count = motif_search(dna, motif)
    
    assert count == 2
    assert positions == [0, 4]

def test_motif_search_missing():
    dna = "ATCG"
    motif = "GGG"
    positions, count = motif_search(dna, motif)
    
    assert count == 0
    assert positions == []

def test_motif_search_empty_sequence():
    dna = ""
    motif = "ATC"
    positions, count = motif_search(dna, motif)
    
    assert count == 0
    assert positions == []

def test_motif_search_empty_motif():
    dna = "ATCG"
    motif = ""
    positions, count = motif_search(dna, motif)
    
    assert count == 0
    assert positions == []