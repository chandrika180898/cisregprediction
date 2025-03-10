import streamlit as st
import pandas as pd
from Bio import SeqIO
from io import StringIO
import re
import plotly.express as px
from concurrent.futures import ProcessPoolExecutor
from reportlab.pdfgen import canvas
from Bio.Seq import Seq

# Sidebar Navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "About", "Contact"])

if page == "About":
    st.title("About DNA Motif Analysis")
    st.write("""
        This tool analyzes **Non-B DNA motifs** in DNA sequences, including:
        - **Slipped DNA**
        - **Z-DNA**
        - **Short Tandem Repeats**
        - **I-Motif**
        - **R-Loop**
        - **Cruciform**
        - **G-Quadruplex**
        - **Hairpin**
        - **Triplex**
        - **H-DNA**
        - **Triplex-forming oligonucleotide (TFO)**
    
        Upload a FASTA file and get detailed motif predictions!
    """)

elif page == "Contact":
    st.title("Contact")
    st.write("""
        **Dr. Y V Rajesh**  
        📧 Email: yvrajesh_bt@kluniversity.in 
        
        **G. Aruna Sesha Chandrika**  
        📧 Email: chandrikagummadi1@gmail.com  
    """)

elif page == "Home":
    # Home Page Content
    st.title('Advanced DNA Promoter Prediction and Non-B DNA Motif Analysis')
    st.write('Upload multiple FASTA files to analyze DNA motifs, predict promoter regions, and visualize results.')

    st.image('https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/utr%20image.jpg', 
             caption='DNA Structure')

    uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)

    # Define motifs dictionary
    motifs = {
        "Slipped DNA": re.compile(r'([ATGC]{2,6})\1{1,}'),
        "Z-DNA": re.compile(r'(CG){6,}'),
        "Short Tandem Repeat": re.compile(r'([ATGC]{2,6})\1{2,}'),
        "I-Motif": re.compile(r'((C[A,T]C){3,})'),
        "R-Loop": re.compile(r'(A{4,}[CG]{2,}A{4,})'),
        "Cruciform": re.compile(r'([ATGC]{4,})\1{2,}'),
        "G-Quadruplex": re.compile(r'(G{3,}[ATGC]{1,5}G{3,}[ATGC]{1,5}G{3,}[ATGC]{1,5}G{3,})'),
        "Hairpin": re.compile(r'([ATGC]{4,})\1{1,}'),
        "Triplex": re.compile(r'(A{3,}[ATGC]{1,}A{3,})'),
        "H-DNA": re.compile(r'([AG]{4,}[CT]{4,}[AG]{4,})'),
        "Triplex-forming oligonucleotide (TFO)": re.compile(r'([GATC]{6,}[AG]{4,}[CT]{4,})')
    }

    # Function to find motifs in a DNA sequence
    def find_motifs(sequence):
        results = []
        for motif_name, motif_pattern in motifs.items():
            for match in motif_pattern.finditer(str(sequence)):
                results.append({
                    "Motif": motif_name,
                    "Start": match.start() + 1,
                    "End": match.end(),
                    "Matched Sequence": sequence[match.start():match.end()]
                })
        return results

    # Process uploaded files
    def process_uploaded_files(uploaded_files):
        all_results = pd.DataFrame()
        for uploaded_file in uploaded_files:
            fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.getvalue().decode('utf-8')), 'fasta'))
            results_data = []
            for record in fasta_sequences:
                motif_results = find_motifs(record.seq)
                for motif in motif_results:
                    results_data.append({
                        "Sequence ID": record.id,
                        **motif,
                        "Length": len(record.seq)
                    })
            results_df = pd.DataFrame(results_data)
            all_results = pd.concat([all_results, results_df], ignore_index=True)
        return all_results

    if uploaded_files:
        try:
            results_df = process_uploaded_files(uploaded_files)

            if not results_df.empty:
                st.write("### Motif Analysis Results")
                st.dataframe(results_df)

                # Download as CSV
                csv = results_df.to_csv(index=False)
                st.download_button(
                    label="Download CSV",
                    data=csv,
                    file_name="motif_analysis_results.csv",
                    mime="text/csv"
                )
            else:
                st.warning("No motifs found in the uploaded sequences!")

        except Exception as e:
            st.error(f"An error occurred: {e}")
