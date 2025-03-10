import streamlit as st
import pandas as pd
from Bio import SeqIO
from io import StringIO
import re
import plotly.express as px
from concurrent.futures import ProcessPoolExecutor
from reportlab.pdfgen import canvas
from Bio.Seq import Seq

st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "Download Report", "About", "Contact"])

if page == "About":
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
        
        Upload a FASTA file and get detailed motif predictions!
    """)

elif page == "Contact":
    st.title("Contact")
    st.write("""
        **Dr. Y V Rajesh**  
        📧 Email: yvrajesh_bt@kluniversity.in 
        
        **G. Aruna Sesha Chandrika**  
        📧 Email: chandrikagummadi1@gmail.com  
    """)

elif page == "Upload & Analyze":
    st.title('Upload and Analyze DNA Sequences')
    uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)

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

    def find_inverted_repeats(sequence):
        inverted_repeat_results = []
        pattern = r'([ATGC]{3,})[ATGC]{0,10}([ATGC]{3,})'
        for match in re.finditer(pattern, str(sequence)):
            part1 = match.group(1)
            part2 = match.group(2)[::-1]
            if part1 == part2:
                inverted_repeat_results.append({
                    "Motif": "Inverted Repeat",
                    "Start": match.start() + 1,
                    "End": match.end(),
                    "Matched Sequence": str(sequence[match.start():match.end()])
                })
        return inverted_repeat_results

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
        results.extend(find_inverted_repeats(sequence))
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

    if uploaded_files:
        try:
            results_df = process_uploaded_files(uploaded_files)

            if 'Matched Sequence' in results_df.columns:
                results_df['Matched Sequence'] = results_df['Matched Sequence'].apply(lambda x: str(x) if isinstance(x, Seq) else x)
            else:
                st.error("No motifs found or the 'Matched Sequence' column is missing!")

            st.session_state["results_df"] = results_df
            st.success("Analysis completed! Go to 'Results' to view.")

        except Exception as e:
            st.error(f"An error occurred: {e}")

elif page == "Results":
    st.title("Analysis Results")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        st.dataframe(results_df)
    else:
        st.warning("No results available. Please upload and analyze sequences first.")

elif page == "Download Report":
    st.title("Download Report")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]

        def generate_pdf(df):
            c = canvas.Canvas("motif_report.pdf")
            c.drawString(100, 800, "DNA Motif Analysis Report")
            y = 780
            for i, row in df.iterrows():
                c.drawString(100, y, f"{row['Sequence ID']} | {row['Motif']} | Start: {row['Start']} | End: {row['End']}")
                y -= 20
            c.save()

        if st.button("Generate PDF Report"):
            generate_pdf(results_df)
            with open("motif_report.pdf", "rb") as pdf:
                st.download_button("Download PDF Report", pdf, file_name="motif_analysis_report.pdf")

        csv = results_df.to_csv(index=False)
        st.download_button("Download CSV", csv, file_name="motif_analysis_results.csv", mime="text/csv")
    else:
        st.warning("No data available. Please analyze sequences first.")
