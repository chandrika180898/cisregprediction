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

# Results Page
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

# Visualization Page
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

# Download Report Page
elif page == "Download Report":
    st.title("Download Report")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        csv = results_df.to_csv(index=False)
        st.download_button("Download CSV", csv, file_name="motif_analysis_results.csv", mime="text/csv")
    else:
        st.warning("No data available. Please analyze sequences first.")
