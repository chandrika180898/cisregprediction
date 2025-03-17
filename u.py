import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from Bio import SeqIO
from io import StringIO
import re
from Bio.Seq import Seq

# Initialize session state
if "results_df" not in st.session_state:
    st.session_state["results_df"] = None

# Sidebar Navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "Download Report", "Visualization", "About", "Contact"])

# Home Page
if page == "Home":
    st.title("Welcome to DNA Motif Analysis Tool")
    st.write("This tool helps analyze DNA sequences to identify various **Non-B DNA motifs**.")
    st.image("https://raw.githubusercontent.com/chandrikagummadi1/cisregprediction/main/images/New%20Microsoft%20PowerPoint%20Presentation.jpg")

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

    # **Motif Identification Functions**
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

    # **Processing FASTA Sequences**
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

        if motifs:
            results_df = pd.DataFrame([{
                "Motif": motif.type,
                "Start": motif.start + 1,
                "End": motif.end + 1
            } for motif in motifs])

            st.session_state["results_df"] = results_df
            st.success("Analysis completed! Go to 'Results' or 'Visualization' to view.")

# **Results Page**
elif page == "Results":
    st.title("Analysis Results")

    if st.session_state["results_df"] is not None:
        results_df = st.session_state["results_df"]
        st.subheader("Motif Analysis Results")
        st.dataframe(results_df)

        motif_counts = results_df["Motif"].value_counts().reset_index()
        motif_counts.columns = ["Motif", "Count"]
        st.subheader("Motif Occurrence Summary")
        st.dataframe(motif_counts)

    else:
        st.warning("No results available. Please upload or paste a sequence first.")

# **Visualization Page**
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

        # Horizontal Thick Lines for Motif Positions
        st.subheader("Motif Start and End Positions")

        fig_lines = go.Figure()
        motif_dict = {motif: i for i, motif in enumerate(results_df["Motif"].unique())}

        for _, row in results_df.iterrows():
            fig_lines.add_trace(go.Scatter(
                x=[row["Start"], row["End"]],
                y=[motif_dict[row["Motif"]], motif_dict[row["Motif"]]],
                mode="lines",
                line=dict(width=6),
                name=row["Motif"]
            ))

        fig_lines.update_layout(
            title="Motif Prediction Start and End Positions",
            xaxis_title="Position in Sequence",
            yaxis_title="Motif",
            yaxis=dict(
                tickmode="array",
                tickvals=list(motif_dict.values()),
                ticktext=list(motif_dict.keys())
            ),
            showlegend=True
        )

        st.plotly_chart(fig_lines)

    else:
        st.warning("No data available for visualization.")

# **Download Report Page**
elif page == "Download Report":
    st.title("Download Report")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        csv = results_df.to_csv(index=False)
        st.download_button("Download CSV", csv, file_name="motif_analysis_results.csv", mime="text/csv")
    else:
        st.warning("No data available. Please analyze sequences first.")
