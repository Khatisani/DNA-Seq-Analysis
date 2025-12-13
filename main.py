import csv
from dna_info import count_nucleotides, gc_content, transcribe, reverse_complement, motif_search

def parse_fasta(input_file):

    sequences = {}
    id = None
    lines = []

    with open(input_file, "r") as file:
        for line in file:
            line = line.strip()  
            if not line:
                    continue 
            if line.startswith(">"):
                if id:
                    sequences[id] = "".join(lines).upper()
                       
                id = line[1:].split()[0]
                lines = []
            else:
                lines.append(line)
  
        if id:
            sequences[id] = "".join(lines).upper()

    return sequences


def write_to_csv():
    if not output:
        print("No output to write.")
        return

    headers = output[0].keys()

    with open(output_file, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=headers)
        writer.writeheader()
        writer.writerows(output)

def main():
    motif = None
    input_file = "example.fasta"
    sequences = parse_fasta("example.fasta")
    output = {}
    for id, sequence in sequences.items():
        a, t, c, g = count_nucleotides(sequence)
        gc = gc_content(a, t, c, g)
        mrna = transcribe(sequence)
        rev_comp = reverse_complement(sequence)

        if motif:
            motif_positions = motif_search(sequence, motif)
            motif_positions = ";".join(map(str, motif_positions)) if motif_positions else "None"
        else:
            motif_positions = "N/A"

        result = {
            "id": id,
            "length": len(sequences),
            "A": a,
            "T": t,
            "C": c,
            "G": g,
            "GC_content": gc_content,
            "RNA_transcript": mRNA,
            "Reverse_complement": reverse_complement,
            "Motif_positions": positions
        }

        output.append(result)

    write_to_csv(output, "results.csv")


if __name__ == "__main__":
    main()