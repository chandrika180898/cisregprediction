# file: code.py

import math
import re
import streamlit as st
import pandas as pd
from io import StringIO


# ------------------------------
# Read uploaded sequence
# ------------------------------

def read_sequence(uploaded_file):

    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))

    seq = []

    for line in stringio:

        if line.startswith(">"):
            continue

        seq.append(line.strip().upper())

    sequence = "".join(seq)

    sequence = re.sub("[^ACGT]", "", sequence)

    return sequence


# ------------------------------
# Perplexity
# ------------------------------

def calculate_perplexity(seq):

    if len(seq) == 0:
        return 0

    counts = {n: seq.count(n) for n in "ACGT"}

    total = sum(counts.values())

    probs = [c/total for c in counts.values() if c > 0]

    entropy = -sum(p * math.log2(p) for p in probs)

    return 2 ** entropy


# ------------------------------
# Sliding window
# ------------------------------

def sliding_windows(seq, window=100):

    windows = []
    perplexities = []

    if len(seq) < window:
        return windows, perplexities

    for i in range(len(seq) - window + 1):

        sub = seq[i:i+window]

        p = calculate_perplexity(sub)

        windows.append((i, i+window, p))
        perplexities.append(p)

    return windows, perplexities


# ------------------------------
# Percentile
# ------------------------------

def percentile(values, percent):

    if not values:
        return None

    values = sorted(values)

    index = int(len(values)*percent/100)

    index = min(index, len(values)-1)

    return values[index]


# ------------------------------
# Low perplexity regions
# ------------------------------

def bottom_percentile_windows(windows, perplexities, percent=5):

    threshold = percentile(perplexities, percent)

    if threshold is None:
        return [], None

    regions = []

    for s,e,p in windows:

        if p <= threshold:
            regions.append((s,e))

    return regions, threshold


# ------------------------------
# Merge regions
# ------------------------------

def merge_regions(regions):

    if not regions:
        return []

    regions = sorted(regions)

    merged = [list(regions[0])]

    for s,e in regions[1:]:

        last = merged[-1]

        if s <= last[1]:

            last[1] = max(last[1],e)

        else:

            merged.append([s,e])

    return merged


# ------------------------------
# Non-B DNA regex
# ------------------------------

def build_nonb_regex():

    motifs = {

        "polyA_polyT": r"A{7,}|T{7,}",

        "STR_repeat": r"([ACGT]{1,6})\1{4,}",

        "G_quadruplex":
        r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}",

        "i_motif":
        r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}",

        "Z_DNA":
        r"(CG){4,}|(GC){4,}"
    }

    return {k: re.compile(v) for k,v in motifs.items()}


# ------------------------------
# Scan motifs
# ------------------------------

def scan_motifs(seq, regex_dict):

    hits = []

    for name,regex in regex_dict.items():

        for m in regex.finditer(seq):

            hits.append((name,m.start(),m.end(),m.group()))

    return hits


# ------------------------------
# Overlap
# ------------------------------

def overlap(a1,a2,b1,b2):

    return max(a1,b1) < min(a2,b2)


def intersect_motifs_lowP(motifs,regions):

    results=[]

    for name,ms,me,mseq in motifs:

        for rs,re in regions:

            if overlap(ms,me,rs,re):

                results.append({
                    "Motif":name,
                    "Motif_start":ms,
                    "Motif_end":me,
                    "Motif_sequence":mseq,
                    "LowP_start":rs,
                    "LowP_end":re
                })

    return results


# ------------------------------
# STREAMLIT UI
# ------------------------------

st.title("Low-Perplexity Non-B DNA Detector")

uploaded_file = st.file_uploader(
    "Upload FASTA / TXT file",
    type=["txt","fa","fasta"]
)

if uploaded_file is not None:

    seq = read_sequence(uploaded_file)

    st.write("Sequence length:",len(seq))

    windows,perplexities = sliding_windows(seq,100)

    low_regions,threshold = bottom_percentile_windows(
        windows,
        perplexities,
        5
    )

    merged = merge_regions(low_regions)

    regex_dict = build_nonb_regex()

    motifs = scan_motifs(seq,regex_dict)

    overlaps = intersect_motifs_lowP(motifs,merged)

    st.write("Perplexity threshold:",threshold)

    if overlaps:

        df = pd.DataFrame(overlaps)

        st.dataframe(df)

        st.download_button(
            "Download CSV",
            df.to_csv(index=False),
            "nonB_low_perplexity_results.csv",
            "text/csv"
        )

    else:

        st.warning("No motifs found in low-perplexity regions")
