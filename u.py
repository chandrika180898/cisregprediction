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
    st.title("Contact Information")
    st.write("📧 Email: yvrajesh_bt@kluniversity.in")

# Function to analyze motifs
def find_motifs(dna):
    motifs = []
    patterns = {
        "APR": r"([ATGC]{3,})\1{2,}",
        "DR": r"(\w{3,})\1",
        "STR": r"([ATGC]{2,6})\1{2,}",
        "Z": r"(GC){6,}",
        "GQ": r"G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}"
    }

    for motif, pattern in patterns.items():
        for match in re.finditer(pattern, dna):
            motifs.append({'start': match.start(), 'length': len(match.group(0)), 'motif': motif})

    return motifs

# Upload & Analyze Page
elif page == "Upload & Analyze":
    st.title('Upload & Analyze DNA Sequences')
    dna_sequence = st.text_area("Paste a DNA Sequence Here:")

    if dna_sequence:
        results = find_motifs(dna_sequence)
        results_df = pd.DataFrame(results)
        st.session_state["results_df"] = results_df
        st.success("Analysis completed! Go to 'Results' to view.")

# Results Page
elif page == "Results":
    st.title("Analysis Results")
    if "results_df" in st.session_state:
        st.dataframe(st.session_state["results_df"])
    else:
        st.warning("No results available. Please upload or paste a sequence first.")

# Download Report Page
elif page == "Download Report":
    st.title("Download Report")
    if "results_df" in st.session_state:
        csv = st.session_state["results_df"].to_csv(index=False)
        st.download_button("Download CSV", csv, file_name="motif_results.csv", mime="text/csv")
    else:
        st.warning("No data available.")
import streamlit as st
import pandas as pd
import re

# Sidebar Navigation
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "About", "Contact"])

# Home Page
if page == "Home":
    st.title("Welcome to Non-B DNA Motif Analysis Tool")
    st.write("Analyze various non-B DNA motifs in your sequences!")

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
    - **Z-DNA (ZDNA)**  
    """)

# Contact Page
elif page == "Contact":
    st.title("Contact Information")
    st.write("📧 Email: yvrajesh_bt@kluniversity.in")

# Function to detect motifs
def find_motifs(dna):
    motifs = []
    patterns = {
        "APR": r"([ATGC]{3,})\1{2,}",
        "DR": r"(\w{3,})\1",
        "STR": r"([ATGC]{2,6})\1{2,}",
        "ZDNA": r"(GC){6,}",
        "GQ": r"G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}"
    }

    for motif, pattern in patterns.items():
        for match in re.finditer(pattern, dna):
            motifs.append({'start': match.start(), 'length': len(match.group(0)), 'motif': motif})

    return motifs

# Upload & Analyze Page
elif page == "Upload & Analyze":  
    st.title('Upload & Analyze DNA Sequences')
    dna_sequence = st.text_area("Paste a DNA Sequence Here:")

    if dna_sequence:
        results = find_motifs(dna_sequence)
        results_df = pd.DataFrame(results)
        st.session_state["results_df"] = results_df
        st.success("Analysis completed! Go to 'Results' to view.")

# Results Page
elif page == "Results":
    st.title("Analysis Results")
    if "results_df" in st.session_state:
        st.dataframe(st.session_state["results_df"])
    else:
        st.warning("No results available. Please upload or paste a sequence first.")
