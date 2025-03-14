import streamlit as st
import pandas as pd
import plotly.express as px
from Bio import SeqIO
from io import StringIO
import re
from concurrent.futures import ProcessPoolExecutor
from reportlab.pdfgen import canvas
from Bio.Seq import Seq

# Sidebar Navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "Visualization", "Download Report", "About", "Contact"])

# ----------- HOME PAGE -----------
if page == "Home":
    st.title("Welcome to NON-B DNA Motif Analysis Tool")
    
    # Fixed: Corrected GitHub Image Path (using the raw URL)
    st.image("https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/New%20Microsoft%20PowerPoint%20Presentation.jpg")

    st.write("Upload or paste DNA sequences to analyze Non-B DNA motifs.")

# ----------- UPLOAD & ANALYZE PAGE -----------
elif page == "Upload & Analyze":
    st.title("Upload and Analyze NON-B DNA Sequences")

    uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)
    pasted_sequence = st.text_area("Or paste your DNA sequence here:")

    # Define Motif Patterns
    motifs = {
       def get_a_tracts(dna, dna_rc, min_at, max_at):
    total_bases = len(dna)
    n_pat = 0
    p_aprs = []
    
    i = 0
    n_as = 0
    
    while i < total_bases:
        if dna[i] in ('a', 't'):
            n_as += 1
        else:
            if min_at <= n_as <= max_at:
                strt = i - n_as + 1
                at_end = strt + n_as
                
                max_at_len = max_t_len = 0
                max_at_len_rc = max_t_len_rc = 0
                max_at_end = max_at_end_rc = 0
                
                n_rc = total_bases - at_end
                at_len = alen = tlen = talen = 0
                at_len_rc = alen_rc = tlen_rc = talen_rc = 0
                
                for n in range(strt - 1, at_end - 1):
                    n_rc += 1
                    
                    if dna[n] == 'a':
                        tlen = talen = 0
                        alen = 0 if dna[n - 1] == 't' else alen + 1
                        at_len += 1
                    if dna_rc[n_rc] == 'a':
                        tlen_rc = talen_rc = 0
                        alen_rc = 0 if dna_rc[n_rc - 1] == 't' else alen_rc + 1
                        at_len_rc += 1
                    
                    if dna[n] == 't':
                        if talen < alen:
                            talen += 1
                            at_len += 1
                        else:
                            tlen += 1
                            talen = at_len = alen = 0
                    if dna_rc[n_rc] == 't':
                        if talen_rc < alen_rc:
                            talen_rc += 1
                            at_len_rc += 1
                        else:
                            tlen_rc += 1
                            talen_rc = at_len_rc = alen_rc = 0
                    
                    if max_at_len < at_len:
                        max_at_len = at_len
                        max_at_end = n
                    if max_t_len < tlen:
                        max_t_len = tlen
                    if max_at_len_rc < at_len_rc:
                        max_at_len_rc = at_len_rc
                        max_at_end_rc = n_rc
                    if max_t_len_rc < tlen_rc:
                        max_t_len_rc = tlen_rc
                
                if (max_at_len - max_t_len) >= min_at or (max_at_len_rc - max_t_len_rc) >= min_at:
                    p_aprs.append({
                        'start': strt,
                        'end': strt + n_as,
                        'a_center': ((max_at_end - ((max_at_len - 1) / 2)) + 1)
                        if (max_at_len - max_t_len) >= (max_at_len_rc - max_t_len_rc)
                        else total_bases - (max_at_end_rc - ((max_at_len_rc - 1) / 2))
                    })
                    n_pat += 1
            n_as = 0
        i += 1
    
    print(f"n potential a tracts = {n_pat}")
    return p_aprs

def find_apr(dna, dna_rc, min_apr, max_apr, min_atracts):
    p_aprs = get_a_tracts(dna, dna_rc, min_apr, max_apr)
    
    n_processed_ats = len(p_aprs)
    arep = []
    
    tracts = 1
    ndx = 0
    for i in range(n_processed_ats - (min_atracts + 1)):
        dist_to_next = p_aprs[i + 1]['a_center'] - p_aprs[i]['a_center']
        if 9.9 <= dist_to_next <= 11.1:
            tracts += 1
        else:
            if tracts >= min_atracts:
                arep.append({
                    'start': p_aprs[i - tracts + 1]['start'],
                    'loop': 0,
                    'num': tracts,
                    'strand': 0,
                    'len': tracts,
                    'end': p_aprs[i]['end'] - 1
                })
                ndx += 1
            tracts = 1
    return arep

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
            fasta_sequences = list(SeqIO.parse(StringIO(uploaded_file.getvalue().decode('utf-8')), 'fasta'))
            for record in fasta_sequences:
                all_results.extend(find_motifs(record.seq, record.id))
        return pd.DataFrame(all_results)

    results_df = pd.DataFrame()

    if uploaded_files:
        results_df = process_uploaded_files(uploaded_files)
    elif pasted_sequence:
        results_df = pd.DataFrame(find_motifs(pasted_sequence))

    if not results_df.empty:
        st.session_state["results_df"] = results_df
        st.success("Analysis completed! Go to 'Results' or 'Visualization'.")

# ----------- RESULTS PAGE -----------
elif page == "Results":
    st.title("Analysis Results")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        st.dataframe(results_df)
        
        # Motif Occurrence Table
        motif_occurrence = results_df["Motif"].value_counts().reset_index()
        motif_occurrence.columns = ["Motif", "Total Count"]
        st.subheader("Motif Occurrence Summary")
        st.dataframe(motif_occurrence)
    else:
        st.warning("No results available. Please upload or paste a sequence first.")

# ----------- VISUALIZATION PAGE -----------
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
    else:
        st.warning("No data available for visualization.")

# ----------- DOWNLOAD REPORT PAGE -----------
elif page == "Download Report":
    st.title("Download Report")
    if "results_df" in st.session_state:
        results_df = st.session_state["results_df"]
        
        csv = results_df.to_csv(index=False)
        st.download_button("Download CSV", csv, file_name="motif_analysis_results.csv", mime="text/csv")
    else:
        st.warning("No data available for download.")

# ----------- ABOUT PAGE -----------
elif page == "About":
    st.title("About DNA Motif Analysis")
    st.write("(Detailed description about motifs)")

# ----------- CONTACT PAGE -----------
elif page == "Contact":
    st.title("Contact")
    st.write("Dr. Y V Rajesh: yvrajesh_bt@kluniversity.in")
    st.write("G. Aruna Sesha Chandrika: chandrikagummadi1@gmail.com")
