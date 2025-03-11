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
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "Download Report", "About", "Contact"])

# Home Page
if page == "Home":
    st.title("Welcome to DNA Motif Analysis Tool")
    st.write("""
        This tool helps analyze DNA sequences to identify various **Non-B DNA motifs**.
         st.image("https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/New%20Microsoft%20PowerPoint%20Presentation.jpg")
    """)

# About Page
elif page == "About":
    st.title("About DNA Motif Analysis")
    st.write("""
        This tool analyzes **Non-B DNA motifs** in DNA sequences, including:
        - **Slipped DNA**
        - **Z-DNA**
        - **Short Tandem Repeats**
        - **I-Motif**
        - **R-Loop**
        - **Cruciform**
        - **G-Quadruplex**
        - **Hairpin**
        - **Triplex**
        - **H-DNA**
        - **Triplex-forming oligonucleotide (TFO)**
        
        Upload a FASTA file or paste a sequence and get detailed motif predictions!
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
    
    def find_motifs(sequence):
        results = []
        for motif_name, motif_pattern in motifs.items():
            for match in motif_pattern.finditer(str(sequence)):
                results.append({
                    "Motif": motif_name,
                    "Start": match.start() + 1,
                    "End": match.end(),
                    "Matched Sequence": str(sequence[match.start():match.end()])
                })
        return results
    
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
    
    def process_uploaded_files(uploaded_files):
        all_results = pd.DataFrame()
        for uploaded_file in uploaded_files:
            fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.getvalue().decode('utf-8')), 'fasta'))
            results_df = analyze_sequences_parallel(fasta_sequences)
            all_results = pd.concat([all_results, results_df], ignore_index=True)
        return all_results
    
    def process_pasted_sequence(sequence):
        fake_fasta_record = [SeqIO.SeqRecord(Seq(sequence), id="Pasted_Sequence", description="Pasted Sequence Analysis")]
        return analyze_sequences_parallel(fake_fasta_record)
    
    if uploaded_files or pasted_sequence:
        try:
            results_df = pd.DataFrame()
            if uploaded_files:
                results_df = process_uploaded_files(uploaded_files)
            if pasted_sequence:
                results_df = pd.concat([results_df, process_pasted_sequence(pasted_sequence)], ignore_index=True)
            
            if 'Matched Sequence' in results_df.columns:
                results_df['Matched Sequence'] = results_df['Matched Sequence'].astype(str)
            else:
                st.error("No motifs found or the 'Matched Sequence' column is missing!")
            
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
    else:
        st.warning("No results available. Please upload or paste sequences first.")

# Download Report Page
elif page == "Download Report":
    st.title("Download Report")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        csv = results_df.to_csv(index=False)
        st.download_button("Download CSV", csv, file_name="motif_analysis_results.csv", mime="text/csv")
    else:
        st.warning("No data available. Please analyze sequences first.")
