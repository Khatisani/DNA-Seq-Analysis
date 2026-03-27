import csv
from Bio import SeqIO
from dna_info import count_nucleotides, gc_content, transcribe, reverse_complement, motif_search

def write_to_csv(output, output_file):
    output_file = "results.csv"

    if not output:
        print("No output to write.")
        return

    keys = output[0].keys()

    with open(output_file, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=keys)
        writer.writeheader()
        writer.writerows(output)
    print(f"Success! Analyzed {len(output)} sequences. See {output_file}.")


def main():
    input_file = "example.fasta"
    output = []

    try:
        dna_seqs = list(SeqIO.parse(input_file, "fasta"))
    except FileNotFoundError:
        print("Error: {input_file} not found")
        return
    
    motif = input("Enter motif to search (press Enter to skip): ").strip()
    
    for dna_seq in dna_seqs:
        dna_seq_str = str(dna_seq.sequence)

        nucleotide_counts, nucleotide_percentages = count_nucleotides(dna_seq_str) 
        gc = gc_content(dna_seq_str)
        mRNA = transcribe(dna_seq_str)
        reverse_comp = reverse_complement(dna_seq_str)

        if motif:
            motif_positions = motif_search(dna_seq_str, motif)
            motif_counts = len(motif_positions)
            motif_positions_str = " ".join(map(str, motif_positions)) if motif_positions else "None"

        else:
            motif_counts = 0
            motif_positions_str = "None"

        row = {
            "id": dna_seq.id,
            "length": len(dna_seq_str),
            "GC%": gc,
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

if __name__ == "__main__":
    main()