import streamlit as st
import pandas as pd
from Bio import SeqIO
from io import StringIO
import re
import plotly.express as px
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from io import BytesIO
from Bio.Seq import Seq

# Title and Description
st.title('Advanced DNA Promoter Prediction and Non-B DNA Motif Analysis')
st.write('Upload multiple FASTA files to analyze DNA motifs, predict promoter regions, and visualize results.')

# Display Images of Motifs
st.header("Non-B DNA Motif Structures")
col1, col2, col3 = st.columns(3)

with col1:
    st.image("images/utr_image.jpg", caption="Untranslated Region (UTR)", use_column_width=True)
with col2:
    st.image("images/z_dna.jpg", caption="Z-DNA Structure", use_column_width=True)
with col3:
    st.image("images/g_quadruplex.jpg", caption="G-Quadruplex", use_column_width=True)

st.write("The images above represent different DNA motifs analyzed during the sequence processing.")

# Upload Section
uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)

# Motif Patterns
motifs = {
    "Slipped DNA": re.compile(r'([ATGC]{2,6})\1{1,}'),
    "Z-DNA": re.compile(r'(CG){6,}'),
    "Short Tandem Repeat": re.compile(r'([ATGC]{2,6})\1{2,}'),
    "I-Motif": re.compile(r'((C[A,T]C){3,})'),
    "R-Motif": re.compile(r'(A{4,}[CG]{2,}A{4,})'),
    "Cruciform": re.compile(r'([ATGC]{4,})\1{2,}'),
    "G-Quadruplex": re.compile(r'(G{3,}[ATGC]{1,5}G{3,}[ATGC]{1,5}G{3,}[ATGC]{1,5}G{3,})'),
    "Hairpin": re.compile(r'([ATGC]{4,})\1{1,}'),
    "Triplex": re.compile(r'(A{3,}[ATGC]{1,}A{3,})')
}

# Find Inverted Repeats
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
                "Matched Sequence": sequence[match.start():match.end()]
            })
    return inverted_repeat_results

# Find Motifs
def find_motifs(sequence):
    results = []
    for motif_name, motif_pattern in motifs.items():
        for match in motif_pattern.finditer(str(sequence)):
            results.append({
                "Motif": motif_name,
                "Start": match.start() + 1,
                "End": match.end(),
                "Matched Sequence": sequence[match.start():match.end()]
            })
    results.extend(find_inverted_repeats(sequence))
    return results

# Analyze Sequences
def analyze_sequences(sequences):
    data = []
    for record in sequences:
        motif_results = find_motifs(record.seq)
        for motif in motif_results:
            data.append({
                "Sequence ID": record.id,
                **motif,
                "Length": len(record.seq)
            })
    return pd.DataFrame(data)

# Visualization
def visualize_motifs(df):
    if not df.empty:
        fig = px.scatter(df, x='Start', y='Sequence ID', color='Motif',
                         hover_data=['Matched Sequence'],
                         title="Motif Distribution Across Sequences")
        st.plotly_chart(fig)
    else:
        st.warning("No motifs to visualize.")

# Generate PDF Report
def generate_pdf(df):
    buffer = BytesIO()
    c = canvas.Canvas(buffer, pagesize=letter)
    c.drawString(100, 800, "DNA Motif Analysis Report")
    y = 780
    for i, row in df.iterrows():
        c.drawString(100, y, f"{row['Sequence ID']} | {row['Motif']} | Start: {row['Start']} | End: {row['End']}")
        y -= 20
        if y < 50:  # Start a new page if space runs out
            c.showPage()
            y = 780
    c.save()
    buffer.seek(0)
    return buffer

# Process Uploaded Files
def process_uploaded_files(uploaded_files):
    all_results = pd.DataFrame()
    for uploaded_file in uploaded_files:
        try:
            fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.getvalue().decode('utf-8')), 'fasta'))
            if fasta_sequences:
                results_df = analyze_sequences(fasta_sequences)
                all_results = pd.concat([all_results, results_df], ignore_index=True)
            else:
                st.error(f"No valid sequences in {uploaded_file.name}")
        except Exception as e:
            st.error(f"Error processing {uploaded_file.name}: {str(e)}")
    return all_results

# Main Analysis Workflow
if uploaded_files:
    results_df = process_uploaded_files(uploaded_files)
    
    if not results_df.empty and 'Matched Sequence' in results_df.columns:
        results_df['Matched Sequence'] = results_df['Matched Sequence'].apply(lambda x: str(x) if isinstance(x, Seq) else x)
        st.write("### Motif Analysis Results")
        st.dataframe(results_df)
        visualize_motifs(results_df)
        
        if st.button("Generate PDF Report"):
            pdf_buffer = generate_pdf(results_df)
            st.download_button(
                "Download PDF Report", pdf_buffer, file_name="motif_analysis_report.pdf")

        csv = results_df.to_csv(index=False)
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name="motif_analysis_results.csv",
            mime="text/csv"
        )
    else:
        st.error("No motifs found or the 'Matched Sequence' column is missing!")
