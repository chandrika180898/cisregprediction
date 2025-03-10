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
page = st.sidebar.radio("Go to", ["Home", "Upload & Analyze", "Results", "Download Report", "About", "Contact"])

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

elif page == "Home":
    st.title('Advanced DNA Promoter Prediction and Non-B DNA Motif Analysis')
    st.write('Upload multiple FASTA files to analyze DNA motifs, predict promoter regions, and visualize results.')
    st.image('https://raw.githubusercontent.com/chandrika180898/cisregprediction/main/images/utr%20image.jpg', caption='DNA Structure')
    
    uploaded_files = st.file_uploader("Upload FASTA Files", type=['fasta'], accept_multiple_files=True)
    
    # Processing logic remains the same...

    if uploaded_files:
        try:
            results_df = process_uploaded_files(uploaded_files)
            
            if 'Matched Sequence' in results_df.columns:
                results_df['Matched Sequence'] = results_df['Matched Sequence'].apply(lambda x: str(x) if isinstance(x, Seq) else x)
            else:
                st.error("No motifs found or the 'Matched Sequence' column is missing!")
    
            st.write("### Motif Analysis Results")
            st.dataframe(results_df)
            visualize_motifs(results_df)
    
            if st.button("Generate PDF Report"):
                generate_pdf(results_df)
                with open("motif_report.pdf", "rb") as pdf:
                    st.download_button("Download PDF Report", pdf, file_name="motif_analysis_report.pdf")
    
            csv = results_df.to_csv(index=False)
            st.download_button(
                label="Download CSV",
                data=csv,
                file_name="motif_analysis_results.csv",
                mime="text/csv"
            )
        except Exception as e:
            st.error(f"An error occurred: {e}")
