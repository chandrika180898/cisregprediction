import streamlit as st
import pandas as pd
import plotly.express as px
from Bio import SeqIO
from io import StringIO
import re
from concurrent.futures import ProcessPoolExecutor
from reportlab.pdfgen import canvas
from Bio.Seq import Seq

# Sidebar Navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "Visualization", "Download Report", "About", "Contact"])

# ----------- HOME PAGE -----------
if page == "Home":
    st.title("Welcome to NON-B DNA Motif Analysis Tool")
    
    # Fixed: Corrected GitHub Image Path (using the raw URL)
    st.image("https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/New%20Microsoft%20PowerPoint%20Presentation.jpg")

    st.write("Upload or paste DNA sequences to analyze Non-B DNA motifs.")

# ----------- UPLOAD & ANALYZE PAGE -----------
elif page == "Upload & Analyze":
    st.title("Upload and Analyze NON-B DNA Sequences")

    uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)
    pasted_sequence = st.text_area("Or paste your DNA sequence here:")

    # Define Motif Patterns
    motifs = {
        "Slipped DNA": re.compile(r'([ATGC]{2,6})\1{1,}'),
        "Z-DNA": re.compile(r'(CG){6,}'),
        "Short Tandem Repeat": re.compile(r'([ATGC]{2,6})\1{2,}'),
        "I-Motif": re.compile(r'((C[A,T]C){3,})'),
        "R-Loop": re.compile(r'(A{4,}[CG]{2,}A{4,})'),
        "Cruciform": re.compile(r'([ATGC]{4,})\1{2,}'),
        "G-Quadruplex": re.compile(r'(G{3,7})([ATCG]{1,7})(G{3,7})\1{2,}'),
        "Bipartite G-Quadruplex": re.compile(r'(G{3}N{1,3}G{3}N{1,3}G{3})N{1,7}(G{3}N{1,3}G{3}N{1,3}G{3})'),
        "G-Triplex DNA (G3-DNA)": re.compile(r'(G{3}N{1,7}){2}G{3}'),
        "G-Hairpin": re.compile(r'(G{3,})N{1,7}(G{3,})'),
        "G-Guanine Slip-Strand DNA": re.compile(r'(GGG){3,}'),
        "Hairpin": re.compile(r'([ATGC]{4,})\1{1,}'),
        "Triplex": re.compile(r'(A{3,}[ATGC]{1,}A{3,})'),
        "H-DNA": re.compile(r'([AG]{4,}[CT]{4,}[AG]{4,})'),
       
    }

    def find_motifs(sequence, seq_id="Pasted Sequence"):
        results = []
        for motif_name, motif_pattern in motifs.items():
            for match in motif_pattern.finditer(str(sequence)):
                results.append({
                    "Sequence ID": seq_id,
                    "Motif": motif_name,
                    "Start": match.start() + 1,
                    "End": match.end(),
                    "Matched Sequence": match.group()
                })
        return results

    def process_uploaded_files(uploaded_files):
        all_results = []
        for uploaded_file in uploaded_files:
            fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.getvalue().decode('utf-8')), 'fasta'))
            for record in fasta_sequences:
                all_results.extend(find_motifs(record.seq, record.id))
        return pd.DataFrame(all_results)

    results_df = pd.DataFrame()

    if uploaded_files:
        results_df = process_uploaded_files(uploaded_files)
    elif pasted_sequence:
        results_df = pd.DataFrame(find_motifs(pasted_sequence))

    if not results_df.empty:
        st.session_state["results_df"] = results_df
        st.success("Analysis completed! Go to 'Results' or 'Visualization'.")
