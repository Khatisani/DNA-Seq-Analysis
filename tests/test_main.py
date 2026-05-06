import pytest
import sys
import os

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