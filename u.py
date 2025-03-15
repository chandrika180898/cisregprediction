import streamlit as st
import pandas as pd
import plotly.express as px
from Bio import SeqIO
from io import StringIO
import re
from reportlab.pdfgen import canvas
from Bio.Seq import Seq

def find_apr(dna):
    pattern = r"([ATGC]{3,})\1{2,}"
    matches = [(m.start(), len(m.group(0))) for m in re.finditer(pattern, dna)]
    return [{'start': m[0], 'len': m[1], 'motif': 'APR'} for m in matches]

def find_direct_repeats(dna):
    pattern = r"(\w{3,})\1"
    matches = [(m.start(), len(m.group(0))) for m in re.finditer(pattern, dna)]
    return [{'start': m[0], 'len': m[1], 'motif': 'Direct Repeat'} for m in matches]

def find_inverted_repeats(dna):
    results = []
    for i in range(len(dna)):
        for j in range(i + 3, len(dna)):
            if dna[i:j] == str(Seq(dna[i:j]).reverse_complement()):
                results.append({'start': i, 'len': j - i, 'motif': 'Inverted Repeat'})
    return results

def find_mirror_repeats(dna):
    results = []
    for i in range(len(dna)):
        for j in range(i + 3, len(dna)):
            if dna[i:j] == dna[i:j][::-1]:
                results.append({'start': i, 'len': j - i, 'motif': 'Mirror Repeat'})
    return results

def find_short_tandem_repeats(dna):
    pattern = r"([ATGC]{2,6})\1{2,}"
    matches = [(m.start(), len(m.group(0))) for m in re.finditer(pattern, dna)]
    return [{'start': m[0], 'len': m[1], 'motif': 'Short Tandem Repeat'} for m in matches]

def find_zdna(dna, min_z=10):
    total_bases = len(dna)
    npy = 1
    kvsum = 0
    zrep = []
    
    for i in range(total_bases - min_z):
        if dna[i:i+2] in ['AC', 'TG', 'CG', 'GC', 'CA', 'GT']:
            npy += 1
            kvsum += 3
        else:
            if npy >= min_z:
                zrep.append({'start': i - npy + 2, 'len': npy, 'motif': 'Z-DNA'})
            npy = 1
            kvsum = 0
    return zrep

def find_g_quadruplex(dna):
    pattern = r"(GGG\w{1,7}){3}GGG"
    matches = [(m.start(), len(m.group(0))) for m in re.finditer(pattern, dna)]
    return [{'start': m[0], 'len': m[1], 'motif': 'G-Quadruplex'} for m in matches]

st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "Visualization", "Download Report", "About", "Contact"])

if page == "Home":
    st.title("Welcome to NON-B DNA Motif Analysis Tool")
    st.image("https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/New%20Microsoft%20PowerPoint%20Presentation.jpg")
    st.write("Upload or paste DNA sequences to analyze Non-B DNA motifs.")

if page == "Upload & Analyze":
    st.title("Upload and Analyze NON-B DNA Sequences")
    uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)
    pasted_sequence = st.text_area("Or paste your DNA sequence here:")

    def process_uploaded_files(uploaded_files):
        all_results = []
        for uploaded_file in uploaded_files:
            fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.getvalue().decode('utf-8')), 'fasta'))
            for record in fasta_sequences:
                dna_seq = str(record.seq)
                all_results.extend(find_apr(dna_seq))
                all_results.extend(find_direct_repeats(dna_seq))
                all_results.extend(find_inverted_repeats(dna_seq))
                all_results.extend(find_mirror_repeats(dna_seq))
                all_results.extend(find_short_tandem_repeats(dna_seq))
                all_results.extend(find_zdna(dna_seq))
                all_results.extend(find_g_quadruplex(dna_seq))
        return pd.DataFrame(all_results)

    results_df = pd.DataFrame()

    if uploaded_files:
        results_df = process_uploaded_files(uploaded_files)
    elif pasted_sequence:
        results_df = pd.DataFrame(
            find_apr(pasted_sequence) +
            find_direct_repeats(pasted_sequence) +
            find_inverted_repeats(pasted_sequence) +
            find_mirror_repeats(pasted_sequence) +
            find_short_tandem_repeats(pasted_sequence) +
            find_zdna(pasted_sequence) +
            find_g_quadruplex(pasted_sequence)
        )

    if not results_df.empty:
        st.session_state["results_df"] = results_df
        st.success("Analysis completed! Go to 'Results' or 'Visualization'.")

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

if page == "About":
    st.title("About DNA Motif Analysis")
    st.write("This tool identifies various Non-B DNA motifs including Z-DNA, Direct Repeats, Inverted Repeats, Mirror Repeats, Short Tandem Repeats, G-Quadruplex, and APR.")

if page == "Contact":
    st.title("Contact")
    st.write("Dr. Y V Rajesh: yvrajesh_bt@kluniversity.in")
    st.write("G. Aruna Sesha Chandrika: chandrikagummadi1@gmail.com")
