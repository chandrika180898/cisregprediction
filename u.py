import streamlit as st
import pandas as pd
import plotly.express as px
from Bio import SeqIO
from io import StringIO
import re
from reportlab.pdfgen import canvas
from Bio.Seq import Seq

# Sidebar Navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "Visualization", "Download Report", "About", "Contact"])

# ----------- HOME PAGE -----------
if page == "Home":
    st.title("Welcome to NON-B DNA Motif Analysis Tool")
    
    st.image("https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/New%20Microsoft%20PowerPoint%20Presentation.jpg")

    st.write("Upload or paste DNA sequences to analyze Non-B DNA motifs.")

# ----------- FUNCTION DEFINITIONS -----------
def find_direct_repeats(dna):
    pattern = r"(\w{3,})\1"
    matches = [(m.start(), len(m.group(0))) for m in re.finditer(pattern, dna)]
    return [{'start': m[0], 'len': m[1], 'motif': 'Direct Repeat'} for m in matches]


def find_inverted_repeats(dna):
    results = []
    dna_seq = Seq(dna)
    rev_comp = str(dna_seq.reverse_complement())
    for i in range(len(dna) - 5):  # Minimum repeat length 6
        for j in range(i + 6, len(dna)):
            if dna[i:j] in rev_comp:
                results.append({'start': i, 'len': j - i, 'motif': 'Inverted Repeat'})
    return results


def find_zdna(dna, min_z):
    total_bases = len(dna)
    npy = 1
    kvsum = 0
    zrep = []

    i = 0
    while i < (total_bases - min_z):
        if dna[i:i+2] in ['AC', 'TG', 'CG', 'GC', 'CA', 'GT']:
            npy += 1
            kvsum += 3
        else:
            if npy >= min_z:
                zrep.append({'start': i - npy + 2, 'len': npy, 'motif': 'Z-DNA'})
            npy = 1
            kvsum = 0
        i += 1
    return zrep

# ----------- UPLOAD & ANALYZE PAGE -----------
if page == "Upload & Analyze":
    st.title("Upload and Analyze NON-B DNA Sequences")
    uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)
    pasted_sequence = st.text_area("Or paste your DNA sequence here:")

    def process_uploaded_files(uploaded_files):
        all_results = []
        for uploaded_file in uploaded_files:
            fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.getvalue().decode('utf-8')), 'fasta'))
            for record in fasta_sequences:
                dna_seq = str(record.seq)
                all_results.extend(find_zdna(dna_seq, 10))
                all_results.extend(find_direct_repeats(dna_seq))
                all_results.extend(find_inverted_repeats(dna_seq))
        return pd.DataFrame(all_results)

    results_df = pd.DataFrame()

    if uploaded_files:
        results_df = process_uploaded_files(uploaded_files)
    elif pasted_sequence:
        results_df = pd.DataFrame(find_zdna(pasted_sequence, 10) + find_direct_repeats(pasted_sequence) + find_inverted_repeats(pasted_sequence))

    if not results_df.empty:
        st.session_state["results_df"] = results_df
        st.success("Analysis completed! Go to 'Results' or 'Visualization'.")

# ----------- RESULTS PAGE -----------
elif page == "Results":
    st.title("Analysis Results")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        st.dataframe(results_df)
    else:
        st.warning("No results available. Please upload or paste a sequence first.")

# ----------- VISUALIZATION PAGE -----------
elif page == "Visualization":
    st.title("Visualization of Motif Analysis")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        motif_counts = results_df["motif"].value_counts().reset_index()
        motif_counts.columns = ["Motif Type", "Count"]
        
        st.subheader("Motif Frequency Bar Chart")
        fig_bar = px.bar(motif_counts, x="Motif Type", y="Count", title="Frequency of Motifs", color="Motif Type")
        st.plotly_chart(fig_bar)
    else:
        st.warning("No data available for visualization.")

# ----------- DOWNLOAD REPORT PAGE -----------
elif page == "Download Report":
    st.title("Download Report")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        csv = results_df.to_csv(index=False)
        st.download_button("Download CSV", csv, file_name="motif_analysis_results.csv", mime="text/csv")
    else:
        st.warning("No data available for download.")

# ----------- ABOUT PAGE -----------
elif page == "About":
    st.title("About DNA Motif Analysis")
    st.write("This tool identifies Non-B DNA motifs such as Z-DNA, Direct Repeats, and Inverted Repeats in uploaded or pasted sequences.")

# ----------- CONTACT PAGE -----------
elif page == "Contact":
    st.title("Contact")
    st.write("Dr. Y V Rajesh: yvrajesh_bt@kluniversity.in")
    st.write("G. Aruna Sesha Chandrika: chandrikagummadi1@gmail.com")
