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
page = st.sidebar.radio("Go to", ["Home", "About", "Contact"])

if page == "About":
    st.title("About DNA Motif Analysis")
    st.write("""
        This tool analyzes **Non-B DNA motifs** in DNA sequences, including:
        - **Slipped DNA**
        - **Z-DNA**
        - **Short Tandem Repeats**
        - **I-Motif**
        - **R-Loop**
        - **Cruciform**
        - **G-Quadruplex**
        - **Hairpin**
        - **Triplex**
        - **H-DNA**
        - **Triplex-forming oligonucleotide (TFO)**
    
        Upload a FASTA file and get detailed motif predictions!
    """)
elif page == "Contact":
    st.title("Contact")
    st.write("""
        **Dr. Y V Rajesh**  
        📧 Email: yvrajesh_bt@kluniversity.in 
        
        **G. Aruna Sesha Chandrika**  
        📧 Email: chandrikagummadi1@gmail.com  
    """)
else:
    st.title('Advanced DNA Promoter Prediction and Non-B DNA Motif Analysis')
    st.write('Upload multiple FASTA files to analyze DNA motifs, predict promoter regions, and visualize results.')
    st
