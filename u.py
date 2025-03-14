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
    
    st.image("https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/New%20Microsoft%20PowerPoint%20Presentation.jpg")

    st.write("Upload or paste DNA sequences to analyze Non-B DNA motifs.")

# ----------- FUNCTION DEFINITIONS -----------
def pupy(dna, pos):
    is_ppy = 0
    if dna[pos] == 'a':
        if dna[pos + 1] == 'c':
            is_ppy = 3
    elif dna[pos] == 't':
        if dna[pos + 1] == 'g':
            is_ppy = 3
    elif dna[pos] == 'c':
        if dna[pos + 1] == 'g':
            is_ppy = 25
        elif dna[pos + 1] == 'a':
            is_ppy = 3
    elif dna[pos] == 'g':
        if dna[pos + 1] == 'c':
            is_ppy = 25
        elif dna[pos + 1] == 't':
            is_ppy = 3
    return is_ppy

def find_zdna(dna, min_z):
    total_bases = len(dna)
    npy = 1
    kvsum = 0
    zrep = []

    i = 0
    while i < (total_bases - min_z):
        tmp_ppy = pupy(dna, i)
        if tmp_ppy > 0:
            npy += 1
            kvsum += tmp_ppy
        else:
            if npy >= min_z:
                zrep.append({
                    'start': i - npy + 2,
                    'len': npy,
                    'loop': kvsum // 2,
                    'num': 0,
                    'end': i + 1,
                    'sub': 0,
                    'strand': 0
                })
            npy = 1
            kvsum = 0
        i += 1
    return zrep

# ----------- UPLOAD & ANALYZE PAGE -----------
elif page == "Upload & Analyze":
    st.title("Upload and Analyze NON-B DNA Sequences")
    uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)
    pasted_sequence = st.text_area("Or paste your DNA sequence here:")

    def process_uploaded_files(uploaded_files):
        all_results = []
        for uploaded_file in uploaded_files:
            fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.getvalue().decode('utf-8')), 'fasta'))
            for record in fasta_sequences:
                all_results.extend(find_zdna(str(record.seq), 10))  # Example threshold value
        return pd.DataFrame(all_results)

    results_df = pd.DataFrame()

    if uploaded_files:
        results_df = process_uploaded_files(uploaded_files)
    elif pasted_sequence:
        results_df = pd.DataFrame(find_zdna(pasted_sequence, 10))

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
        motif_counts = results_df["len"].value_counts().reset_index()
        motif_counts.columns = ["Length", "Count"]
        
        st.subheader("Motif Frequency Bar Chart")
        fig_bar = px.bar(motif_counts, x="Length", y="Count", title="Frequency of Z-DNA Motifs", color="Length")
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
    st.write("This tool identifies Non-B DNA motifs such as Z-DNA in uploaded or pasted sequences.")

# ----------- CONTACT PAGE -----------
elif page == "Contact":
    st.title("Contact")
    st.write("Dr. Y V Rajesh: yvrajesh_bt@kluniversity.in")
    st.write("G. Aruna Sesha Chandrika: chandrikagummadi1@gmail.com")
