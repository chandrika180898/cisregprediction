import streamlit as st
import pandas as pd
import re
from Bio.Seq import Seq
import plotly.express as px

# Sidebar Navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "Visualization", "Download Report", "About", "Contact"])

# Home Page
if page == "Home":
    st.title("Welcome to Non-B DNA Motif Analysis Tool")
    st.write("Analyze various non-B DNA motifs in your sequences with ease!")
    st.image("https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/New%20Microsoft%20PowerPoint%20Presentation.jpg")

# About Page
elif page == "About":
    st.title("About DNA Motif Analysis")
    st.write("""
    - **A-phased Repeats (APR)**  
    - **Direct Repeats (DR)**  
    - **G-quadruplexes (GQ)**  
    - **Inverted Repeats (IR)**  
    - **Mirror Repeats (MR)**  
    - **Short Tandem Repeats (STR)**  
    - **Z-DNA (Z)**  
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

# Function to analyze motifs
def find_motifs(dna):
    motifs = []

    # A-phased repeats (APR)
    for m in re.finditer(r"([ATGC]{3,})\1{2,}", dna):
        motifs.append({'start': m.start(), 'length': len(m.group(0)), 'motif': 'APR'})

    # Direct repeats (DR)
    for m in re.finditer(r"(\w{3,})\1", dna):
        motifs.append({'start': m.start(), 'length': len(m.group(0)), 'motif': 'DR'})

    # Inverted repeats (IR)
    for i in range(len(dna)):
        for j in range(i + 3, len(dna)):
            if dna[i:j] == str(Seq(dna[i:j]).reverse_complement()):
                motifs.append({'start': i, 'length': j - i, 'motif': 'IR'})

    # Mirror repeats (MR)
    for i in range(len(dna)):
        for j in range(i + 3, len(dna)):
            if dna[i:j] == dna[i:j][::-1]:
                motifs.append({'start': i, 'length': j - i, 'motif': 'MR'})

    # Short tandem repeats (STR)
    for m in re.finditer(r"([ATGC]{2,6})\1{2,}", dna):
        motifs.append({'start': m.start(), 'length': len(m.group(0)), 'motif': 'STR'})

    # Z-DNA (Z)
    for m in re.finditer(r"(GC){6,}", dna):
        motifs.append({'start': m.start(), 'length': len(m.group(0)), 'motif': 'Z'})

    # G-quadruplexes (GQ)
    for m in re.finditer(r"G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}", dna):
        motifs.append({'start': m.start(), 'length': len(m.group(0)), 'motif': 'GQ'})

    return motifs

# Upload & Analyze Page
elif page == "Upload & Analyze":
    st.title('Upload and Analyze DNA Sequences')
    uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)
    pasted_sequence = st.text_area("Or Paste a DNA Sequence Here:")

    if uploaded_files or pasted_sequence:
        try:
            if pasted_sequence:
                sequences = [{"name": "Pasted Sequence", "seq": pasted_sequence}]
            else:
                sequences = []
                for uploaded_file in uploaded_files:
                    for record in SeqIO.parse(uploaded_file, "fasta"):
                        sequences.append({"name": record.id, "seq": str(record.seq)})

            all_results = []
            for seq_data in sequences:
                motifs = find_motifs(seq_data["seq"])
                for motif in motifs:
                    all_results.append({"Sequence": seq_data["name"], **motif})

            results_df = pd.DataFrame(all_results)
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
        motif_occurrence = results_df["motif"].value_counts().reset_index()
        motif_occurrence.columns = ["Motif", "Total Count"]
        st.subheader("Motif Occurrence Summary")
        st.dataframe(motif_occurrence)
    else:
        st.warning("No results available. Please upload or paste a sequence first.")

# Visualization Page
elif page == "Visualization":
    st.title("Visualization of Motif Analysis")

    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        motif_counts = results_df["motif"].value_counts().reset_index()
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

# Download Report Page
elif page == "Download Report":
    st.title("Download Report")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        csv = results_df.to_csv(index=False)
        st.download_button("Download CSV", csv, file_name="motif_analysis_results.csv", mime="text/csv")
    else:
        st.warning("No data available. Please analyze sequences first.")
