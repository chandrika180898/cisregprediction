import streamlit as st
import pandas as pd
import plotly.express as px
from Bio import SeqIO
from io import StringIO
import re
from Bio.Seq import Seq

def find_apr(dna):
    pattern = r"([ATGC]{3,})\1{2,}"
    matches = [(m.start(), len(m.group(0))) for m in re.finditer(pattern, dna)]
    return [{'start': m[0], 'len': m[1], 'motif': 'APR'} for m in matches]

def find_direct_repeats(dna):
    pattern = r"(\w{3,})\1"
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
    pattern = r"([ATGC]{2,6})\1{2,}"
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

st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "Visualization", "Download Report", "About", "Contact"])

if page == "Home":
    st.title("Welcome to NON-B DNA Motif Analysis Tool")
    st.image("https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/New%20Microsoft%20PowerPoint%20Presentation.jpg")
    st.write("Upload or paste DNA sequences to analyze Non-B DNA motifs.")

elif page == "Upload & Analyze":
    st.title("Upload and Analyze NON-B DNA Sequences")
    uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)
    pasted_sequence = st.text_area("Or paste your DNA sequence here:")
    
    results_df = pd.DataFrame()
    
    if uploaded_files or pasted_sequence:
        all_results = []
        
        if uploaded_files:
            for uploaded_file in uploaded_files:
                fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.getvalue().decode('utf-8')), 'fasta'))
                for record in fasta_sequences:
                    all_results.extend(analyze_sequence(str(record.seq)))
        
        if pasted_sequence:
            all_results.extend(analyze_sequence(pasted_sequence))
        
        results_df = pd.DataFrame(all_results)
        st.session_state["results_df"] = results_df
        st.success("Analysis completed! Go to 'Results' or 'Visualization'.")

elif page == "Results":
    st.title("Analysis Results")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        if not results_df.empty:
            st.dataframe(results_df)
        else:
            st.warning("No motifs found.")
    else:
        st.warning("No results available. Please upload or paste a sequence first.")

elif page == "Download Report":
    st.title("Download Report")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        if not results_df.empty:
            csv = results_df.to_csv(index=False)
            st.download_button("Download CSV", csv, file_name="motif_analysis_results.csv", mime="text/csv")
        else:
            st.warning("No data available. Please analyze sequences first.")
    else:
        st.warning("No data available. Please analyze sequences first.")

if page == "About":
    st.title("About DNA Motif Analysis")
    st.write("This tool identifies various Non-B DNA motifs including Z-DNA, Direct Repeats, Inverted Repeats, Mirror Repeats, Short Tandem Repeats, G-Quadruplex, and APR.")

if page == "Contact":
    st.title("Contact")
    st.write("Dr. Y V Rajesh: yvrajesh_bt@kluniversity.in")
    st.write("G. Aruna Sesha Chandrika: chandrikagummadi1@gmail.com")
