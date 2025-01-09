import streamlit as st
import pandas as pd
from Bio import SeqIO
from io import StringIO
import re
import plotly.express as px
from concurrent.futures import ProcessPoolExecutor
from reportlab.pdfgen import canvas
from Bio.Seq import Seq

# Title and description
st.title('Advanced DNA Promoter Prediction and Non-B DNA Motif Analysis')
st.write('Upload multiple FASTA files to analyze DNA motifs, predict promoter regions, and visualize results.')

# Displaying images
st.image('images/k.png', caption='NON-B-DNA STRUCTURES')

# Uploading files
uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)

# Dictionary of motifs
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

# Function to find inverted repeats
def find_inverted_repeats(sequence):
    inverted_repeat_results = []
    pattern = r'([ATGC]{3,})[ATGC]{0,10}([ATGC]{3,})'
    for match in re.finditer(pattern, str(sequence)):
        part1 = match.group(1)
        part2 = match.group(2)[::-1]  # Reverse complement
        if part1 == part2:
            inverted_repeat_results.append({
                "Motif": "Inverted Repeat",
                "Start": match.start() + 1,
                "End": match.end(),
                "Matched Sequence": sequence[match.start():match.end()]
            })
    return inverted_repeat_results

# Main logic
if uploaded_files:
    results = []
    for uploaded_file in uploaded_files:
        fasta_content = uploaded_file.read().decode('utf-8')
        sequences = SeqIO.parse(StringIO(fasta_content), 'fasta')
        for record in sequences:
            sequence = str(record.seq)
            sequence_results = []
            # Analyze motifs
            for motif_name, motif_pattern in motifs.items():
                for match in re.finditer(motif_pattern, sequence):
                    sequence_results.append({
                        "Motif": motif_name,
                        "Start": match.start() + 1,
                        "End": match.end(),
                        "Matched Sequence": match.group()
                    })
            # Analyze inverted repeats
            sequence_results.extend(find_inverted_repeats(sequence))
            # Add results
            for result in sequence_results:
                result.update({"Sequence ID": record.id})
                results.append(result)
    
    # Create a DataFrame from results
    results_df = pd.DataFrame(results)
    st.write("Results Overview:")
    st.dataframe(results_df)

    # Download results
    csv = results_df.to_csv(index=False)
    st.download_button(
        label="Download Results as CSV",
        data=csv,
        file_name="dna_motif_analysis_results.csv",
        mime="text/csv"
    )
else:
    st.info("Please upload FASTA files to proceed.")
