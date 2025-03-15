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

    def analyze_sequence(dna_seq):
        results = []
        motifs = {
            "APR": r"([ATGC]{3,})\1{2,}",
            "Direct Repeat": r"(\w{3,})\1",
            "Mirror Repeat": r"(.{3,})\1[::-1]",
            "Short Tandem Repeat": r"([ATGC]{2,6})\1{2,}",
            "G-Quadruplex": r"(GGG\w{1,7}){3}GGG",
        }
        for motif_name, pattern in motifs.items():
            for match in re.finditer(pattern, dna_seq):
                results.append({"Motif": motif_name, "Start": match.start(), "End": match.end(), "Length": len(match.group(0))})
        return results

    def process_uploaded_files(uploaded_files):
        all_results = pd.DataFrame()
        for uploaded_file in uploaded_files:
            fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.getvalue().decode('utf-8')), 'fasta'))
            for record in fasta_sequences:
                results = analyze_sequence(str(record.seq))
                df = pd.DataFrame(results)
                df["Sequence ID"] = record.id
                all_results = pd.concat([all_results, df], ignore_index=True)
        return all_results

    def process_pasted_sequence(sequence):
        results = analyze_sequence(sequence)
        df = pd.DataFrame(results)
        df["Sequence ID"] = "Pasted_Sequence"
        return df

    if uploaded_files or pasted_sequence:
        try:
            results_df = pd.DataFrame()
            if uploaded_files:
                results_df = process_uploaded_files(uploaded_files)
            if pasted_sequence:
                results_df = pd.concat([results_df, process_pasted_sequence(pasted_sequence)], ignore_index=True)
            
            st.session_state["results_df"] = results_df
            st.success("Analysis completed! Go to 'Results' to view.")
        except Exception as e:
            st.error(f"An error occurred: {e}")

# Results Page
elif page == "Results":
    st.title("Analysis Results")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        st.dataframe(results_df)
        st.subheader("Motif Occurrence Summary")
        st.dataframe(results_df["Motif"].value_counts().reset_index().rename(columns={"index": "Motif", "Motif": "Total Count"}))
    else:
        st.warning("No results available. Please upload or paste a sequence first.")

# Visualization Page
elif page == "Visualization":
    st.title("Visualization of Motif Analysis")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        st.subheader("Motif Frequency Bar Chart")
        st.plotly_chart(px.bar(results_df, x="Motif", title="Frequency of Each Motif", color="Motif"))
        st.subheader("Motif Distribution Pie Chart")
        st.plotly_chart(px.pie(results_df, names="Motif", title="Distribution of Motifs"))
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
