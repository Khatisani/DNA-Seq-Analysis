import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import io
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Import custom modules from the src package
from dna_info import (
    validate_sequence, count_nucleotides, gc_content, 
    transcribe, reverse_complement, motif_search, 
    calc_molecular_weight, cal_entropy, cal_gc_skew, validate_motif
)
from blast_integration import run_blast

# --- Page Configuration ---
st.set_page_config(
    page_title="DNA Sequence Analysis Dashboard",
    page_icon="🧬",
    layout="wide"
)

st.title("DNA Sequence Analysis Dashboard")
st.markdown("Upload multiple FASTA files, paste multi-sequence FASTA text, or analyze sequences interactively.")

# --- Sidebar Controls ---
st.sidebar.header("Configuration & Input")

upload_option = st.sidebar.radio("Choose Input Method", ["Upload FASTA File", "Manual Sequence Input"])

records_to_process = []

if upload_option == "Upload FASTA File":
    uploaded_files = st.sidebar.file_uploader(
        "Upload .fasta or .fa file(s)", 
        type=["fasta", "fa", "txt"], 
        accept_multiple_files=True
    )
    if uploaded_files:
        for uploaded_file in uploaded_files:
            stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
            records_to_process.extend(list(SeqIO.parse(stringio, "fasta")))
else:
    st.sidebar.markdown("Paste multiple sequences in **FASTA format** below:")
    default_fasta = ">Seq_1\nATGCGATCGATCGATCGATCGATAGCTAGCTA\n>Seq_2\nCGATCGATCGATCGATCGATAGCTAGCTAGC"
    manual_fasta_input = st.sidebar.text_area("FASTA Text Area", value=default_fasta, height=180)
    
    if manual_fasta_input:
        try:
            stringio = io.StringIO(manual_fasta_input)
            records_to_process = list(SeqIO.parse(stringio, "fasta"))
        except Exception as e:
            st.sidebar.error(f"Error parsing FASTA text: {e}")

# Motif search input

motif_input = st.sidebar.text_input("Search Motif (Optional)", value="").strip()
motif_valid = True
if motif_input:
    is_valid, err_msg = validate_motif(motif_input)
    if not is_valid:
        st.sidebar.error(f"Invalid motif: {err_msg}")
        motif_valid = False
    else:
        st.sidebar.success("Motif is valid.")

# Run Analysis Button
run_button = st.sidebar.button("Run Analysis", type="primary")

# --- Main Logic ---
if run_button:
    if not records_to_process:
        st.warning("Please upload valid FASTA file(s) or provide valid FASTA text.")
    else:
        output_data = []
        
        with st.spinner("Analyzing sequences..."):
            for record in records_to_process:
                dna_seq_str = str(record.seq)
                

                is_valid, err_msg = validate_sequence(dna_seq_str)
                if not is_valid:
                    st.error(f"Skipping '{record.id}': {err_msg}")
                    continue
                

                _, nucleotide_percentages = count_nucleotides(dna_seq_str)
                gc = gc_content(dna_seq_str)
                mRNA = transcribe(dna_seq_str)
                reverse_comp = reverse_complement(dna_seq_str)
                mol_weight = calc_molecular_weight(dna_seq_str, seq_type="DNA")
                entropy = cal_entropy(dna_seq_str)
                gc_skew = cal_gc_skew(dna_seq_str)
                
                if motif_input and motif_valid:
                    positions, m_count = motif_search(dna_seq_str, motif_input)
                    motif_positions_str = ", ".join(map(str, positions)) if positions else "None"
                else:
                    m_count = 0
                    motif_positions_str = "N/A"
                
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
                    "Motif_Count": m_count,
                    "Motif_Positions": motif_positions_str
                }
                output_data.append(row)
        
        if output_data:
            df = pd.DataFrame(output_data)
            os.makedirs("outputs", exist_ok=True)
            output_csv_path = "outputs/results.csv"
            df.to_csv(output_csv_path, index=False)

            st.session_state["df"] = df
            st.session_state["records_to_process"] = records_to_process
            st.success(f"Successfully processed {len(output_data)} sequence(s)!")

# --- Display Results if Available in Session State ---

if "df" in st.session_state:
    df = st.session_state["df"]
    records_to_process = st.session_state.get("records_to_process", [])
    
    st.subheader("Summary Metrics Table")
    st.dataframe(df.drop(columns=["RNA", "Rev_Comp"]), use_container_width=True)

    # Download CSV button

    csv = df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="Download Results as CSV",
        data=csv,
        file_name="sequence_analysis_results.csv",
        mime="text/csv",
    )

    # --- Visualizations Section ---

    st.subheader("Visualizations")
    col1, col2 = st.columns(2)
    
    sns.set_theme(style="whitegrid")
    
    with col1:
        st.markdown("### GC Percentage per Sequence")
        fig1, ax1 = plt.subplots(figsize=(6, 4))
        sns.barplot(data=df, x="id", y="GC_Percent", ax=ax1, palette="viridis", legend=False)
        ax1.set_ylabel("GC Content (%)")
        ax1.set_xticklabels(df["id"], rotation=45, ha="right")
        plt.tight_layout()
        st.pyplot(fig1)
        plt.close(fig1)
        
    with col2:
        st.markdown("### Length vs. Molecular Weight")
        fig2, ax2 = plt.subplots(figsize=(6, 4))
        sns.scatterplot(data=df, x="length", y="Mol_Weight_(Da)", ax=ax2, hue="GC_Percent", palette="magma", s=100)
        ax2.set_xlabel("Sequence Length (bp)")
        ax2.set_ylabel("Molecular Weight (Da)")
        plt.tight_layout()
        st.pyplot(fig2)
        plt.close(fig2)

# --- Individual Sequence Explorer ---
    st.markdown("---")
    st.subheader("Deep Dive Explorer")
    selected_id = st.selectbox("Select sequence for detailed sequence inspection:", df["id"].tolist())
    
    if selected_id:
        seq_row = df[df["id"] == selected_id].iloc[0]
        
        m_col1, m_col2, m_col3, m_col4 = st.columns(4)
        m_col1.metric("Length", f"{seq_row['length']} bp")
        m_col2.metric("GC Content", f"{seq_row['GC_Percent']}%")
        m_col3.metric("Molecular Weight", f"{seq_row['Mol_Weight_(Da)']} Da")
        m_col4.metric("Entropy", seq_row['Entropy'])
        
        with st.expander("View RNA Transcript & Reverse Complement"):
            st.text_area("mRNA Transcript:", seq_row["RNA"], height=100)
            st.text_area("Reverse Complement:", seq_row["Rev_Comp"], height=100)
            
        if motif_input:
            st.info(f"**Motif Search Results for '{motif_input}'**: Found {seq_row['Motif_Count']} time(s) at position(s): {seq_row['Motif_Positions']}")

# Optional BLAST Integration

        if st.button(f"Run Online NCBI BLAST for {selected_id}"):
            with st.spinner("Connecting to NCBI BLAST server (this may take a minute)..."):
                raw_seq_str = next((str(r.seq) for r in records_to_process if r.id == selected_id), None)
                if raw_seq_str:
                    matches = run_blast(raw_seq_str)
                    if matches:
                        st.success("Top BLAST Matches:")
                        for idx, match in enumerate(matches, 1):
                            st.write(f"**{idx}.** {match['title']} *(E-value: {match['e_value']}, Score: {match['score']})*")
                    else:
                        st.warning("No matches found or connection failed.")
                else:
                    st.error("Could not retrieve raw sequence data for BLAST.")