import streamlit as st
import pandas as pd
from Bio import SeqIO
from io import StringIO
import re
import plotly.express as px
from concurrent.futures import ProcessPoolExecutor
from reportlab.pdfgen import canvas
from Bio.Seq import Seq

# Set Page Configuration
st.set_page_config(page_title="DNA Motif Analysis", layout="wide")

# Add Custom Navbar with ZSeeker Style
st.markdown("""
    <style>
        .navbar {
            background-color: #2C2F36;
            padding: 10px;
            display: flex;
            align-items: center;
            justify-content: space-between;
            position: fixed;
            width: 100%;
            top: 0;
            left: 0;
            z-index: 1000;
        }
        .navbar img {
            height: 50px;
            margin-right: 15px;
        }
        .navbar a {
            color: white;
            text-decoration: none;
            font-size: 18px;
            margin-left: 20px;
        }
        .navbar a:hover {
            color: #00bfff;
        }
        .content {
            margin-top: 80px;
        }
    </style>
""", unsafe_allow_html=True)

st.markdown("""
    <div class="navbar">
        <div style="display: flex; align-items: center;">
            <img src="https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/utr%20image.jpg">
            <span style="color:white; font-size:22px; font-weight:bold;">DNA Motif Analysis</span>
        </div>
        <div>
            <a href="#">Job Submission</a>
            <a href="#">About/Contact</a>
            <a href="#">Help</a>
        </div>
    </div>
    <div class="content">
""", unsafe_allow_html=True)

# Page Title
st.title('🔬 Advanced DNA Promoter Prediction and Non-B DNA Motif Analysis')
st.write('Upload multiple FASTA files to analyze DNA motifs, predict promoter regions, and visualize results.')

# File Upload Section
uploaded_files = st.file_uploader("📂 Upload FASTA Files", type=['fasta'], accept_multiple_files=True)

# DNA Motifs Dictionary (Including H-DNA and R-Loop)
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

# Function to Find Inverted Repeats
def find_inverted_repeats(sequence):
    inverted_repeat_results = []
    pattern = r'([ATGC]{3,})[ATGC]{0,10}([ATGC]{3,})'
    for match in re.finditer(pattern, str(sequence)):
        part1 = match.group(1)
        part2 = match.group(2)[::-1]  # Reverse the second part
        if part1 == part2:
            inverted_repeat_results.append({
                "Motif": "Inverted Repeat",
                "Start": match.start() + 1,
                "End": match.end(),
                "Matched Sequence": sequence[match.start():match.end()]
            })
    return inverted_repeat_results

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
    results.extend(find_inverted_repeats(sequence))
    return results

# Parallel Processing for Faster Analysis
def analyze_sequences_parallel(sequences):
    data = []
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(find_motifs, [record.seq for record in sequences]))
        for record, motif_results in zip(sequences, results):
            for motif in motif_results:
                data.append({
                    "Sequence ID": record.id,
                    **motif,
                    "Length": len(record.seq)
                })
    return pd.DataFrame(data)

# Visualization Function
def visualize_motifs(df):
    fig = px.scatter(df, x='Start', y='Sequence ID', color='Motif',
                     hover_data=['Matched Sequence'],
                     title="📊 Motif Distribution Across Sequences")
    st.plotly_chart(fig)

# PDF Generation Function
def generate_pdf(df):
    c = canvas.Canvas("motif_report.pdf")
    c.drawString(100, 800, "DNA Motif Analysis Report")
    y = 780
    for i, row in df.iterrows():
        c.drawString(100, y, f"{row['Sequence ID']} | {row['Motif']} | Start: {row['Start']} | End: {row['End']}")
        y -= 20
    c.save()

# Processing Uploaded FASTA Files
def process_uploaded_files(uploaded_files):
    all_results = pd.DataFrame()
    for uploaded_file in uploaded_files:
        fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.getvalue().decode('utf-8')), 'fasta'))
        results_df = analyze_sequences_parallel(fasta_sequences)
        all_results = pd.concat([all_results, results_df], ignore_index=True)
    return all_results

# Handling File Uploads and Displaying Results
if uploaded_files:
    try:
        results_df = process_uploaded_files(uploaded_files)

        if 'Matched Sequence' in results_df.columns:
            results_df['Matched Sequence'] = results_df['Matched Sequence'].astype(str)
        else:
            st.error("⚠️ No motifs found or 'Matched Sequence' column is missing!")

        st.write("### 📑 Motif Analysis Results")
        st.dataframe(results_df)
        visualize_motifs(results_df)

        if st.button("📄 Generate PDF Report"):
            generate_pdf(results_df)
            with open("motif_report.pdf", "rb") as pdf:
                st.download_button("📥 Download PDF Report", pdf, file_name="motif_analysis_report.pdf")

        csv = results_df.to_csv(index=False)
        st.download_button("📥 Download CSV", data=csv, file_name="motif_analysis_results.csv", mime="text/csv")

    except Exception as e:
        st.error(f"❌ An error occurred: {e}")
