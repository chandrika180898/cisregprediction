import streamlit as st
import pandas as pd
from Bio import SeqIO
from io import StringIO
import re

# Sidebar Navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "Download Report", "About", "Contact"])

# Home Page
if page == "Home":
    st.title("Welcome to DNA Motif Analysis Tool")
    st.write("This tool helps analyze DNA sequences to identify various **Non-B DNA motifs**.")
    st.image("https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/New%20Microsoft%20PowerPoint%20Presentation.jpg")

# About Page
elif page == "About":
    st.title("About DNA Motif Analysis")
    st.write("""
    This tool detects various DNA motifs including:
    - **Direct Repeats (DR)**
    - **Inverted Repeats (IR)**
    - **Short Tandem Repeats (STR)**
    - **Z-DNA**
    - **G-Quadruplex (GQ)**
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

    class Motif:
        def __init__(self, start, end, motif_type):
            self.start = start
            self.end = end
            self.type = motif_type

    motifs = []

    def store_motif(start, end, motif_type):
        motifs.append(Motif(start, end, motif_type))

    def find_direct_repeats(dna, min_size=6, max_gap=5):
        seq_len = len(dna)
        for i in range(seq_len - min_size):
            for j in range(i + min_size + max_gap, seq_len):
                if dna[i:i + min_size] == dna[j:j + min_size]:
                    store_motif(i, j + min_size - 1, "Direct Repeat")

    def find_inverted_repeats(dna, min_size=6, max_gap=5):
        seq_len = len(dna)
        for i in range(seq_len - min_size):
            for j in range(i + min_size + max_gap, seq_len):
                if dna[i:i + min_size] == dna[j - min_size + 1:j + 1][::-1]:
                    store_motif(i, j, "Inverted Repeat")

    def find_short_tandem_repeats(dna, min_repeats=3):
        seq_len = len(dna)
        for i in range(seq_len - 2):
            for j in range(1, seq_len // 2):
                repeat_count = 0
                while dna[i:i + j] == dna[i + j * repeat_count:i + j * (repeat_count + 1)]:
                    repeat_count += 1
                    if repeat_count >= min_repeats:
                        store_motif(i, i + j * repeat_count - 1, "Short Tandem Repeat")
                        break

    def find_z_dna(dna, min_size=6):
        seq_len = len(dna)
        for i in range(seq_len - min_size):
            if re.fullmatch(r"[cg]+", dna[i:i + min_size]):
                store_motif(i, i + min_size - 1, "Z-DNA")

    def find_g_quadruplex(dna, minGQ=4, maxGQspacer=10):
        for match in re.finditer(r"(g{4,})(.{0," + str(maxGQspacer) + r"})(g{4,})", dna, re.IGNORECASE):
            store_motif(match.start(), match.end(), "G-Quadruplex")

    def read_fasta(file):
        sequences = []
        file_content = StringIO(file.getvalue().decode("utf-8"))
        for record in SeqIO.parse(file_content, "fasta"):
            sequences.append((record.id, str(record.seq).lower()))
        return sequences

    def analyze_sequences(sequences):
        motifs.clear()
        for seq_id, dna in sequences:
            find_direct_repeats(dna)
            find_inverted_repeats(dna)
            find_short_tandem_repeats(dna)
            find_z_dna(dna)
            find_g_quadruplex(dna)

    if uploaded_files or pasted_sequence:
        sequences = []
        if uploaded_files:
            for uploaded_file in uploaded_files:
                sequences.extend(read_fasta(uploaded_file))
        if pasted_sequence:
            sequences.append(("Pasted_Sequence", pasted_sequence.lower()))

        analyze_sequences(sequences)
        
        results_df = pd.DataFrame([{
            "Motif": motif.type,
            "Start": motif.start + 1,
            "End": motif.end + 1
        } for motif in motifs])

        if not results_df.empty:
            st.session_state["results_df"] = results_df
            st.success("Analysis completed! Go to 'Results' to view.")

# Results Page
elif page == "Results":
    st.title("Analysis Results")

    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]

        st.subheader("Motif Analysis Results")
        st.table(results_df)  # Display as table

        motif_counts = results_df["Motif"].value_counts().reset_index()
        motif_counts.columns = ["Motif", "Count"]
        st.subheader("Motif Occurrence Summary")
        st.table(motif_counts)

    else:
        st.warning("No results available. Please upload or paste a sequence first.")

# Download Report Page
elif page == "Download Report":
    st.title("Download Analysis Report")

    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]

        csv_data = results_df.to_csv(index=False).encode("utf-8")

        st.download_button(
            label="Download CSV Report",
            data=csv_data,
            file_name="motif_analysis_results.csv",
            mime="text/csv"
        )
    else:
        st.warning("No results available. Please upload or paste a sequence first.")
