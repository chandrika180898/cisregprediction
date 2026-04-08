import streamlit as st
import pandas as pd
from Bio import SeqIO
from io import StringIO
import re
import math

# Sidebar
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", ["Home","Upload & Analyze","Results"])

# ---------------- HOME ----------------
if page == "Home":
    st.title("Low Perplexity Non-B DNA Detector")

    st.write("""
    This tool detects **Non-B DNA motifs located in low-perplexity regions**.

    Pipeline:
    1. Sliding window (100 bp)
    2. Calculate perplexity
    3. Select bottom 5% windows
    4. Merge overlapping regions
    5. Detect Non-B DNA motifs
    """)

# ---------------- FUNCTIONS ----------------

def calc_perplexity_from_counts(counts):

    total = sum(counts.values())

    probs = [c/total for c in counts.values() if c>0]

    entropy = -sum(p*math.log2(p) for p in probs)

    return 2**entropy


def sliding_windows(seq, window=100):

    windows=[]
    perplexities=[]

    counts={"A":0,"C":0,"G":0,"T":0}

    for n in seq[:window]:
        if n in counts:
            counts[n]+=1

    p = calc_perplexity_from_counts(counts)

    windows.append((0,window,p))
    perplexities.append(p)

    for i in range(1,len(seq)-window+1):

        left=seq[i-1]
        right=seq[i+window-1]

        if left in counts:
            counts[left]-=1

        if right in counts:
            counts[right]+=1

        p = calc_perplexity_from_counts(counts)

        windows.append((i,i+window,p))
        perplexities.append(p)

    return windows,perplexities


def percentile(values,percent):

    values=sorted(values)

    index=int(len(values)*percent/100)

    return values[index]


def bottom_percentile_windows(windows,perplexities,percent=5):

    threshold=percentile(perplexities,percent)

    regions=[]

    for start,end,p in windows:
        if p<=threshold:
            regions.append((start,end))

    return regions


def merge_regions(regions):

    if not regions:
        return []

    regions=sorted(regions)

    merged=[list(regions[0])]

    for s,e in regions[1:]:

        last=merged[-1]

        if s<=last[1]:
            last[1]=max(last[1],e)

        else:
            merged.append([s,e])

    return merged


def build_nonb_regex():

    motifs={

        "PolyA_T":r"A{7,}|T{7,}",

        "STR_repeat":r"([ACGT]{1,6})\1{4,}",

        "G_quadruplex":r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}",

        "i_motif":r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}",

        "Z_DNA":r"(CG){4,}|(GC){4,}"
    }

    return {k:re.compile(v) for k,v in motifs.items()}


def scan_motifs(seq,regex_dict):

    hits=[]

    for name,regex in regex_dict.items():

        for m in regex.finditer(seq):

            hits.append((name,m.start(),m.end(),m.group()))

    return hits


def overlap(a_start,a_end,b_start,b_end):

    return max(a_start,b_start) < min(a_end,b_end)


# ---------------- ANALYSIS ----------------

def analyze_sequences(sequences):

    regex_dict = build_nonb_regex()

    results=[]

    for record in sequences:

        seq=str(record.seq).upper()

        if len(seq) < 100:
            continue

        if len(seq) > 500000:
            st.warning(f"{record.id} skipped (sequence too large)")
            continue

        windows,perplexities = sliding_windows(seq,100)

        low_regions = bottom_percentile_windows(windows,perplexities,5)

        merged_regions = merge_regions(low_regions)

        motifs = scan_motifs(seq,regex_dict)

        for name,ms,me,motif_seq in motifs:

            for rs,re in merged_regions:

                if overlap(ms,me,rs,re):

                    results.append({

                        "Sequence ID":record.id,
                        "Motif":name,
                        "Motif Start":ms,
                        "Motif End":me,
                        "Motif Sequence":motif_seq,
                        "LowP Start":rs,
                        "LowP End":re

                    })

    return pd.DataFrame(results)


# ---------------- UPLOAD PAGE ----------------

elif page == "Upload & Analyze":

    st.title("Upload DNA Sequences")

    uploaded_files = st.file_uploader(
        "Upload FASTA files",
        type=["fasta","fa","txt"],
        accept_multiple_files=True
    )

    if uploaded_files:

        if st.button("Run Analysis"):

            try:

                all_sequences=[]

                for uploaded_file in uploaded_files:

                    fasta_sequences=list(
                        SeqIO.parse(
                            StringIO(uploaded_file.getvalue().decode("utf-8")),
                            "fasta"
                        )
                    )

                    all_sequences.extend(fasta_sequences)

                results_df = analyze_sequences(all_sequences)

                st.session_state["results_df"] = results_df

                st.success("Analysis completed!")

            except Exception as e:

                st.error(f"Error: {e}")


# ---------------- RESULTS PAGE ----------------

elif page == "Results":

    st.title("Results")

    if "results_df" in st.session_state:

        df=st.session_state["results_df"]

        st.dataframe(df)

        csv=df.to_csv(index=False)

        st.download_button(
            "Download CSV",
            csv,
            "low_perplexity_nonB_results.csv",
            "text/csv"
        )

    else:

        st.warning("No results available. Run analysis first.")
