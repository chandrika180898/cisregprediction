import streamlit as st 
import pandas as pd
from Bio import SeqIO
from io import StringIO
import re
import plotly.express as px
from concurrent.futures import ProcessPoolExecutor
from reportlab.pdfgen import canvas
from Bio.Seq import Seq

# Page Configuration
st.set_page_config(
    page_title="DNA Motif & Promoter Analysis",
    layout="wide",
)

# Sidebar for File Upload
st.sidebar.image('https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/utr%20image.jpg', 
                 use_column_width=True)
st.sidebar.title("Upload FASTA Files")
uploaded_files = st.sidebar.file_uploader("Select FASTA Files", type=['fasta'], accept_multiple_files=True)

# Page Title
st.title("🔬 DNA Motif & Promoter Prediction")
st.markdown("This tool identifies **non-B DNA motifs, promoter regions**, and visualizes results.")

# Motif Definitions
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

# Function to Find Motifs
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

# Function for Parallel Analysis
def analyze_sequences_parallel(sequences):
    data = []
    with ProcessPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(find_motifs, [record.seq for record in sequences]))
        for record, motif_results in zip(sequences, results):
            for motif in motif_results:
                data.append({
                    "Sequence ID": record.id,
                    **motif,
                    "Length": len(record.seq)
                })
    return pd.DataFrame(data)

# Visualization
def visualize_motifs(df):
    fig = px.scatter(df, x='Start', y='Sequence ID', color='Motif',
                     hover_data=['Matched Sequence'], title="Motif Distribution")
    st.plotly_chart(fig)

# Process Uploaded Files
if uploaded_files:
    st.subheader("📊 Processing Uploaded Files...")
    all_results = pd.DataFrame()
    
    for uploaded_file in uploaded_files:
        fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.getvalue().decode('utf-8')), 'fasta'))
        results_df = analyze_sequences_parallel(fasta_sequences)
        all_results = pd.concat([all_results, results_df], ignore_index=True)

    if not all_results.empty:
        # Show Summary
        motif_counts = all_results["Motif"].value_counts().reset_index()
        motif_counts.columns = ["Motif", "Count"]
        st.sidebar.write("### 📌 Summary Statistics")
        st.sidebar.dataframe(motif_counts)

        # Display Results
        st.write("### 🔍 Motif Analysis Results")
        st.dataframe(all_results)

        # Visualization
        visualize_motifs(all_results)

        # Download CSV
        csv = all_results.to_csv(index=False).encode('utf-8')
        st.download_button(label="📥 Download Results (CSV)", data=csv, file_name="motif_analysis_results.csv", mime="text/csv")
