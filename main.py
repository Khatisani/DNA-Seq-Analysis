import csv
from Bio import SeqIO
from dna_info import count_nucleotides, gc_content, transcribe, reverse_complement, motif_search, calc_molecular_weight, cal_entropy, cal_gc_skew, validate_sequence, validate_motif
from visualize import create_visualizations
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from blast_integration import run_blast

def main():
    input_file = "example.fasta"
    output_file = "results.csv"
    output = []

    try:
        dna_seqs = list(SeqIO.parse(input_file, "fasta"))
    except FileNotFoundError:
        print("Error: {input_file} not found")
        return
    
    if not dna_seqs:
        print("Error: No sequences found")
        return
    
    motif = input("Enter motif to search (press Enter to skip): ").strip()

    motif_valid = True
    if motif:
        is_valid, error_message = validate_motif(motif)
        if not is_valid:
            print(f"Skipping motif search '{motif}': {error_message}")
            motif_valid = False
            motif = " "

    for record in dna_seqs:
        dna_seq_str = str(record.seq)

        is_valid, error_message = validate_sequence(dna_seq_str)
        if not is_valid:
            print(f"Skipping sequence '{record.id}': {error_message}")
            continue

        nucleotide_counts, nucleotide_percentages = count_nucleotides(dna_seq_str) 
        gc = gc_content(dna_seq_str)
        mRNA = transcribe(dna_seq_str)
        reverse_comp = reverse_complement(dna_seq_str)
        mol_weight = calc_molecular_weight (dna_seq_str, seq_type = "DNA")
        entropy = cal_entropy(dna_seq_str)
        gc_skew = cal_gc_skew(dna_seq_str)

        if motif and motif_valid:
            motif_positions = motif_search(dna_seq_str, motif)
            motif_counts = len(motif_positions)
            motif_positions_str = " ".join(map(str, motif_positions)) if motif_positions else "None"

        else:
            motif_counts = 0
            motif_positions_str = "N/A"

        run_search = input(f"Do you want to run a BLAST search for sequence {record.id}? (y/n): ").strip().lower()
        if run_search == 'y':
            matches = run_blast(dna_seq_str)
            if matches:
                print(f"\nTop matches for {record.id}:")
                for i, match in enumerate(matches, 1):
                    print(f"{i}. {match['title']} (E-value: {match['e_value']})")
            else:
                print("No matches found or connection failed.")


        row = {
            "id": record.id,
            "length": len(dna_seq_str),
            "GC_Percent": gc,
            "GC_Skew": gc_skew,
            "Entropy": entropy,
            "Mol_Weight_(Da)": mol_weight,
            "A%": nucleotide_percentages["A"],
            "T%": nucleotide_percentages["T"],
            "C%": nucleotide_percentages["C"],
            "G%": nucleotide_percentages["G"],
            "N%": nucleotide_percentages["N"],
            "RNA": mRNA,
            "Rev_Comp": reverse_comp,
            "Motif_Count": motif_counts,
            "Motif_Positions": motif_positions_str
        }

        output.append(row)

    if output:
        df = pd.DataFrame(output)

        df.to_csv(output_file, index = False)
        print(f"Data saved to {output_file}. Preview: ")
        print(df[["id", "length", "GC_Percent", "Mol_Weight_(Da)"]].head())

        create_visualizations(df)
    else:
        print("No valid sequences were processed.")

if __name__ == "__main__":
    main()