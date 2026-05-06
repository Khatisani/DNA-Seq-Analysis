# DNA Sequence Analysis

Parse FASTA files, calculate DNA statistics, transcribe to RNA, reverse complement sequences, and optionally search for DNA motifs. 

--- 

## **Features**

- **File Processing & Validation**: Parses FASTA files and validates both sequence integrity and search motifs.

- **Nucleotide Analysis**: Calculates sequence length, nucleotide percentages (A, T, C, G, N), and GC content.

- **Genomic Signatures**: Calculates GC Skew for locating leading/lagging replication strands.

- **Biochemical Properties**: Computes sequence molecular weight (in Daltons) and Shannon sequence entropy.

- **Sequence Manipulation**: Transcribes DNA → RNA and generates reverse complements.

- **Motif Searching**: Locates specific nucleotide patterns within the sequence.

- **Data Integration & Dashboard**: Exports analysis metrics to a Pandas DataFrame/CSV and saves a visualization dashboard.
---

## ** Requirements & Installation **

- Install Dependencies 
`pip install biopython pandas matplotlib seaborn` 

## ** Next Steps **
- ✔ Manual data handling
- ✔ Biopython Integration 

- ✔ Molecular Weight 
- ✔ Sequence Entropy
- ✔ GC skew

- ✔ Error Handling and Validation
- ✔ Pandas and Matplotlib
- BLAST Integration
- Translation (start stop codon, ORF)

- Argument Parsing 
- File upload (Tkinter)
- Web App (Streamlit)




