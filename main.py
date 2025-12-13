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


def write_to_csv(output, output_file):
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
    output = []
    motif = input("Enter motif to search (press Enter to skip): ").strip()

    if motif == "":
        motif = None

    for id, sequence in sequences.items():
        
        percentage_a, a, percentage_t, t, percentage_c, c, percentage_g, g = count_nucleotides(sequence)
        gc = gc_content(a, t, c, g)
        mrna = transcribe(sequence)
        rev_comp = reverse_complement(sequence)

        if motif:
            positions = motif_search(sequence, motif)
            if positions:
                positions = " ".join(map(str,positions))
                pos = positions.split(";")
                number = len(pos)
            else: 
                positions = "None"
                number = 0
        else:
            positions = "N/A"
            number = 0

        result = {
            "id": id,
            "length": len(sequence),
            "A%": percentage_a,
            "T%": percentage_t,
            "C%": percentage_c,
            "G%": percentage_g,
            "GC_content(%)": gc,
            "RNA_transcript": mrna,
            "Reverse_complement": rev_comp,
            "Motif_Count": number,
            "Motif_positions": positions
        }
        

        output.append(result)

    write_to_csv(output, "results.csv")


if __name__ == "__main__":
    main()