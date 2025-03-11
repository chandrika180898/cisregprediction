import streamlit as st
import pandas as pd
import plotly.express as px
from Bio import SeqIO
from io import StringIO
import re
from Bio.Seq import Seq

# Sidebar Navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "Visualization", "Download Report", "About", "Contact"])

# ----------- HOME PAGE -----------
if page == "Home":
    st.title("Welcome to NON-B DNA Motif Analysis Tool")
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
        "I-Motif": re.compile(r'((C[A,T]C){3,})'),
        "R-Loop": re.compile(r'(A{4,}[CG]{2,}A{4,})'),
        "Cruciform": re.compile(r'([ATGC]{4,})\1{2,}'),
        "G-Quadruplex": re.compile(r'(G{3,7})([ATCG]{1,7})(G{3,7})\1{2,}'),
        "Bipartite G-Quadruplex": re.compile(r'(G{3}N{1,3}G{3}N{1,3}G{3})N{1,7}(G{3}N{1,3}G{3}N{1,3}G{3})'),
        "G-Triplex DNA (G3-DNA)": re.compile(r'(G{3}N{1,7}){2}G{3}'),
        "G-Hairpin": re.compile(r'(G{3,})N{1,7}(G{3,})'),
        "G-Guanine Slip-Strand DNA": re.compile(r'(GGG){3,}')
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

# ----------- RESULTS PAGE -----------
elif page == "Results":
    st.title("Analysis Results")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        st.dataframe(results_df)
        motif_occurrence = results_df["Motif"].value_counts().reset_index()
        motif_occurrence.columns = ["Motif", "Total Count"]
        st.subheader("Motif Occurrence Summary")
        st.dataframe(motif_occurrence)
    else:
        st.warning("No results available. Please upload or paste a sequence first.")

# ----------- VISUALIZATION PAGE -----------
elif page == "Visualization":
    st.title("Visualization of Motif Analysis")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        motif_counts = results_df["Motif"].value_counts().reset_index()
        motif_counts.columns = ["Motif", "Count"]
        
        st.subheader("Motif Frequency Bar Chart")
        fig_bar = px.bar(motif_counts, x="Motif", y="Count", title="Frequency of Each Motif", color="Motif")
        st.plotly_chart(fig_bar)
        
        st.subheader("Motif Distribution Pie Chart")
        fig_pie = px.pie(motif_counts, names="Motif", values="Count", title="Distribution of Motifs")
        st.plotly_chart(fig_pie)
        
        st.subheader("Motif Positions in Sequences")
        fig_scatter = px.scatter(results_df, x="Start", y="End", color="Motif", title="Start vs. End Positions of Motifs")
        st.plotly_chart(fig_scatter)
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
    st.write("This tool helps analyze various non-B DNA motifs in uploaded DNA sequences.")

# ----------- CONTACT PAGE -----------
elif page == "Contact":
    st.title("Contact")
    st.write("Dr. Y V Rajesh: yvrajesh_bt@kluniversity.in")
    st.write("G. Aruna Sesha Chandrika: chandrikagummadi1@gmail.com")
