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
    
    # Fixed: Corrected GitHub Image Path (using the raw URL)
    st.image("https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/k.png")
    
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
            fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.read().decode('utf-8')), 'fasta'))
            for record in fasta_sequences:
                all_results.extend(find_motifs(record.seq, record.id))
        return pd.DataFrame(all_results)

    results_df = pd.DataFrame()

    if uploaded_files:
        results_df = process_uploaded_files(uploaded_files)
    elif pasted_sequence:
        results_df = pd.DataFrame(find_motifs(Seq(pasted_sequence)))

    if not results_df.empty:
        st.session_state["results_df"] = results_df
        st.success("Analysis completed! Go to 'Results' or 'Visualization'.")

# ----------- RESULTS PAGE -----------
elif page == "Results":
    st.title("Analysis Results")
    st.session_state.setdefault("results_df", None)
    if st.session_state["results_df"] is not None:
        st.dataframe(st.session_state["results_df"])
    else:
        st.warning("No results available. Please upload or paste a sequence first.")

# ----------- VISUALIZATION PAGE -----------
elif page == "Visualization":
    st.title("Visualization of Motif Analysis")
    st.session_state.setdefault("results_df", None)
    if st.session_state["results_df"] is not None:
        results_df = st.session_state["results_df"]

        motif_counts = results_df["Motif"].value_counts().reset_index()
        motif_counts.columns = ["Motif", "Count"]
        
        # Bar Chart
        st.subheader("Motif Frequency Bar Chart")
        fig_bar = px.bar(motif_counts, x="Motif", y="Count", title="Frequency of Each Motif", color="Motif")
        st.plotly_chart(fig_bar)
        
        # Pie Chart
        st.subheader("Motif Distribution Pie Chart")
        fig_pie = px.pie(motif_counts, names="Motif", values="Count", title="Distribution of Motifs")
        st.plotly_chart(fig_pie)
    else:
        st.warning("No data available for visualization.")

# ----------- DOWNLOAD REPORT PAGE -----------
elif page == "Download Report":
    st.title("Download Report")
    st.session_state.setdefault("results_df", None)
    if st.session_state["results_df"] is not None:
        results_df = st.session_state["results_df"]
        csv = results_df.to_csv(index=False)
        st.download_button("Download CSV", csv, file_name="motif_analysis_results.csv", mime="text/csv")
    else:
        st.warning("No data available for download.")

# ----------- ABOUT PAGE -----------
elif page == "About":
    st.title("About DNA Motif Analysis")
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

# ----------- CONTACT PAGE -----------
elif page == "Contact":
    st.title("Contact")
    st.write("Dr. Y V Rajesh: yvrajesh_bt@kluniversity.in")
    st.write("G. Aruna Sesha Chandrika: chandrikagummadi1@gmail.com")
