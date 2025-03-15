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
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "Visualization", "Download Report", "About", "Contact"])

# Home Page
if page == "Home":
    st.title("Welcome to DNA Motif Analysis Tool")
    st.write("""
        This tool helps analyze DNA sequences to identify various **Non-B DNA motifs**.
        
    """)
    st.image("https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/New%20Microsoft%20PowerPoint%20Presentation.jpg")
# About Page

elif page == "About":
    st.title("About DNA Motif Analysis")
    st.write("""
    - **A-phased repeats (APRs):** Comprise three or more A/T-rich segments separated by 10-nucleotide spacers.
    - **Direct repeats (DRs):** Consist of repeated 4- to 10-nucleotide sequences within a genome.
    - **G-quadruplexes (G4s):** Four-stranded DNA structures stabilized by Hoogsteen hydrogen bonds and cations.
    - **Inverted repeats (IRs):** Formed when inter-strand base pairing shifts to intra-strand pairing, leading to cruciform DNA.
    - **Mirror repeats (MRs):** Homopurine/pyrimidine sequences with a mirrored arrangement, capable of forming triplex DNA.
    - **Short tandem repeats (STRs):** Microsatellites with 2-6 bp nucleotide sequences repeating consecutively in a genome.
    - **Z-DNA:** A non-canonical left-handed double-helix structure found in regulatory regions.
    - **I-motif:** A four-stranded structure stabilized by cytosine–cytosine+ base pairs, forming under acidic conditions.
    - **A-form DNA:** Inverted G/C tracts exhibiting A-like base stacking, recognized by transcription factors.
    - **Parallel-stranded DNA:** Purine-rich sequences stabilized by reverse Hoogsteen hydrogen bonding, forming triplexes or quadruplexes.
    """)

# Contact Page
elif page == "Contact":
    st.title("Contact")
    st.write("""
        **Dr. Y V Rajesh**  
        📧 Email: yvrajesh_bt@kluniversity.in 
        
        **G. Aruna Sesha Chandrika**  
        📧 Email: chandrikagummadi1@gmail.com  
    """)

# Upload & Analyze Page
elif page == "Upload & Analyze":
    st.title('Upload and Analyze DNA Sequences')
    uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)
    pasted_sequence = st.text_area("Or Paste a DNA Sequence Here:")
    
    def find_apr(dna):
        pattern = r"([ATGC]{3,})\\1{2,}"
        matches = [(m.start(), len(m.group(0))) for m in re.finditer(pattern, dna)]
        return [{'start': m[0], 'len': m[1], 'motif': 'APR'} for m in matches]

    def find_direct_repeats(dna):
        pattern = r"(\w{3,})\\1"
        matches = [(m.start(), len(m.group(0))) for m in re.finditer(pattern, dna)]
        return [{'start': m[0], 'len': m[1], 'motif': 'Direct Repeat'} for m in matches]

    def find_inverted_repeats(dna):
        results = []
        for i in range(len(dna)):
            for j in range(i + 3, len(dna)):
                if dna[i:j] == str(Seq(dna[i:j]).reverse_complement()):
                    results.append({'start': i, 'len': j - i, 'motif': 'Inverted Repeat'})
        return results

    def find_mirror_repeats(dna):
        results = []
        for i in range(len(dna)):
            for j in range(i + 3, len(dna)):
                if dna[i:j] == dna[i:j][::-1]:
                    results.append({'start': i, 'len': j - i, 'motif': 'Mirror Repeat'})
        return results

    def find_short_tandem_repeats(dna):
        pattern = r"([ATGC]{2,6})\\1{2,}"
        matches = [(m.start(), len(m.group(0))) for m in re.finditer(pattern, dna)]
        return [{'start': m[0], 'len': m[1], 'motif': 'Short Tandem Repeat'} for m in matches]

    def find_zdna(dna, min_z=10):
        total_bases = len(dna)
        npy = 1
        zrep = []
        
        for i in range(total_bases - min_z):
            if dna[i:i+2] in ['AC', 'TG', 'CG', 'GC', 'CA', 'GT']:
                npy += 1
            else:
                if npy >= min_z:
                    zrep.append({'start': i - npy + 2, 'len': npy, 'motif': 'Z-DNA'})
                npy = 1
        return zrep

    def find_g_quadruplex(dna):
        pattern = r"(GGG\w{1,7}){3}GGG"
        matches = [(m.start(), len(m.group(0))) for m in re.finditer(pattern, dna)]
        return [{'start': m[0], 'len': m[1], 'motif': 'G-Quadruplex'} for m in matches]

    def analyze_sequence(dna_seq):
        return (
            find_apr(dna_seq) +
            find_direct_repeats(dna_seq) +
            find_inverted_repeats(dna_seq) +
            find_mirror_repeats(dna_seq) +
            find_short_tandem_repeats(dna_seq) +
            find_zdna(dna_seq) +
            find_g_quadruplex(dna_seq)
        )

    if uploaded_files or pasted_sequence:
        try:
            st.success("Analysis completed! Go to 'Results' to view.")
        except Exception as e:
            st.error(f"An error occurred: {e}")
