import csv

input_file = "example.fasta"
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
    ...