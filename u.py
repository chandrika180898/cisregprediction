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

elif page == "Upload & Analyze":
    st.title('Upload and Analyze DNA Sequences')
    
    uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)
    pasted_sequence = st.text_area("Or Paste a DNA Sequence Here:")

    # Motif class to store motifs with start, end, and type
    class Motif:
        def __init__(self, start, end, motif_type, seq_id):
            self.start = start
            self.end = end
            self.type = motif_type
            self.seq_id = seq_id

    motifs = []

    def store_motif(start, end, motif_type, seq_id):
        motifs.append(Motif(start, end, motif_type, seq_id))

    # Motif Identification Functions
    def find_direct_repeats(dna, seq_id, min_size=6, max_gap=5):
        seq_len = len(dna)
        for i in range(seq_len - min_size):
            for j in range(i + min_size + max_gap, seq_len):
                if dna[i:i + min_size] == dna[j:j + min_size]:
                    store_motif(i, j + min_size - 1, "Direct Repeat", seq_id)

    def find_inverted_repeats(dna, seq_id, min_size=6, max_gap=5):
        seq_len = len(dna)
        for i in range(seq_len - min_size):
            for j in range(i + min_size + max_gap, seq_len):
                if dna[i:i + min_size] == dna[j - min_size + 1:j + 1][::-1]:
                    store_motif(i, j, "Inverted Repeat", seq_id)

    def find_short_tandem_repeats(dna, seq_id, min_repeats=3):
        seq_len = len(dna)
        for i in range(seq_len - 2):
            for j in range(1, seq_len // 2):
                repeat_count = 0
                while dna[i:i + j] == dna[i + j * repeat_count:i + j * (repeat_count + 1)]:
                    repeat_count += 1
                    if repeat_count >= min_repeats:
                        store_motif(i, i + j * repeat_count - 1, "Short Tandem Repeat", seq_id)
                        break

    def find_z_dna(dna, seq_id, min_size=6):
        seq_len = len(dna)
        for i in range(seq_len - min_size):
            if re.fullmatch(r"[cg]+", dna[i:i + min_size]):
                store_motif(i, i + min_size - 1, "Z-DNA", seq_id)

    def find_g_quadruplex(dna, seq_id, minGQ=4, maxGQspacer=10):
        for match in re.finditer(r"(g{4,})(.{0," + str(maxGQspacer) + r"})(g{4,})", dna, re.IGNORECASE):
            store_motif(match.start(), match.end(), "G-Quadruplex", seq_id)

    # Function to read FASTA
    def read_fasta(file):
        sequences = []
        file_content = StringIO(file.getvalue().decode("utf-8"))
        for record in SeqIO.parse(file_content, "fasta"):
            sequences.append((record.id, str(record.seq).lower()))
        return sequences

    def analyze_sequences(sequences):
        motifs.clear()
        for seq_id, dna in sequences:
            find_direct_repeats(dna, seq_id)
            find_inverted_repeats(dna, seq_id)
            find_short_tandem_repeats(dna, seq_id)
            find_z_dna(dna, seq_id)
            find_g_quadruplex(dna, seq_id)

    # Processing input
    if uploaded_files or pasted_sequence:
        sequences = []
        if uploaded_files:
            for uploaded_file in uploaded_files:
                sequences.extend(read_fasta(uploaded_file))
        if pasted_sequence:
            sequences.append(("Pasted_Sequence", pasted_sequence.lower()))

        analyze_sequences(sequences)

        if motifs:
            results_df = pd.DataFrame([{
                "Sequence ID": motif.seq_id,
                "Motif": motif.type,
                "Start": motif.start + 1,
                "End": motif.end + 1
            } for motif in motifs])

            st.session_state["results_df"] = results_df
            st.success("Analysis completed! Go to 'Results' or 'Visualization' to view.")
        else:
            st.warning("No motifs found in the provided sequences.")

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
elif page == "Visualization":
    st.title("Visualization of Motif Analysis")
    
    if "results_df" in st.session_state:
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
        
        # Scatter Plot
        st.subheader("Motif Positions in Sequences")
        fig_scatter = px.scatter(results_df, x="Start", y="End", color="Motif", title="Start vs. End Positions of Motifs")
        st.plotly_chart(fig_scatter)
        
        # Horizontal Thick Lines for Motif Positions
        st.subheader("Motif Start and End Positions")
        import plotly.graph_objects as go

        fig_lines = go.Figure()

        for _, row in results_df.iterrows():
            fig_lines.add_trace(go.Scatter(
                x=[row["Start"], row["End"]],
                y=[row["Motif"], row["Motif"]],
                mode="lines",
                line=dict(width=6),  # Thick lines for clarity
                name=row["Motif"]
            ))

        fig_lines.update_layout(
            title="Motif Prediction Start and End Positions",
            xaxis_title="Position in Sequence",
            yaxis_title="Motif",
            showlegend=False
        )

        st.plotly_chart(fig_lines)

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
