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

- **Data Integration & Dashboard**: Exports analysis metrics to CSV, provides interactive visualizations, and runs live via Streamlit.

---

## ** Requirements & Installation **

Install Dependencies 
`pip install -r requirements.txt` 

## ** Running the Web App ** 

Launch the interactive Streamlit dashboard locally:
`streamlit run app.py`

## ** Next Steps **

- Translation (start stop codon, ORF)





