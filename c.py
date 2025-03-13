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
        "G-Guanine Slip-Strand DNA": re.compile(r'(GGG){3,}'),
        "A-Tract": re.compile(r'(A{7,})'),

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
