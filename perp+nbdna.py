import streamlit as st 
import pandas as pd
from Bio import SeqIO
from io import StringIO
import re
import plotly.express as px
from concurrent.futures import ProcessPoolExecutor
from reportlab.pdfgen import canvas
from Bio.Seq import Seq

st.title("Low Perplexity Non-B DNA Detector")

uploaded_files = st.file_uploader(
    "Upload FASTA/TXT sequence",
    type=["fasta","fa","txt"],
    accept_multiple_files=True
)

# ----------------------------
# READ SEQUENCE
# ----------------------------
def read_sequences(uploaded_file):

    sequences = []

    text = uploaded_file.getvalue().decode("utf-8")

    if text.startswith(">"):
        fasta = list(SeqIO.parse(StringIO(text), "fasta"))
        for record in fasta:
            sequences.append((record.id, str(record.seq).upper()))
    else:
        seq = re.sub("[^ACGT]", "", text.upper())
        sequences.append(("sequence", seq))

    return sequences


# ----------------------------
# PERPLEXITY
# ----------------------------
def calculate_perplexity(seq):

    counts = {n: seq.count(n) for n in "ACGT"}

    total = sum(counts.values())

    probs = [c/total for c in counts.values() if c>0]

    entropy = -sum(p*math.log2(p) for p in probs)

    return 2**entropy


# ----------------------------
# SLIDING WINDOW
# ----------------------------
def sliding_windows(seq, window=100):

    windows = []
    perplexities = []

    for i in range(len(seq)-window+1):

        sub = seq[i:i+window]

        p = calculate_perplexity(sub)

        windows.append((i,i+window,p))
        perplexities.append(p)

    return windows, perplexities


# ----------------------------
# LOW PERPLEXITY
# ----------------------------
def percentile(values, percent):

    values = sorted(values)

    index = int(len(values)*percent/100)

    return values[index]


def bottom_percentile_windows(windows, perplexities, percent=5):

    threshold = percentile(perplexities, percent)

    regions = []

    for start,end,p in windows:

        if p <= threshold:
            regions.append((start,end))

    return regions, threshold


# ----------------------------
# MERGE REGIONS
# ----------------------------
def merge_regions(regions):

    if not regions:
        return []

    regions = sorted(regions)

    merged = [list(regions[0])]

    for s,e in regions[1:]:

        last = merged[-1]

        if s <= last[1]:
            last[1] = max(last[1], e)

        else:
            merged.append([s,e])

    return merged


# ----------------------------
# NON-B DNA REGEX
# ----------------------------
def build_nonb_regex():

    motifs = {

        "PolyA_T": r"A{7,}|T{7,}",

        "STR_repeat": r"([ACGT]{1,6})\1{4,}",

        "G_quadruplex":
        r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}",

        "i_motif":
        r"C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}",

        "Z_DNA":
        r"(CG){4,}|(GC){4,}"
    }

    return {k: re.compile(v) for k,v in motifs.items()}


# ----------------------------
# SCAN MOTIFS
# ----------------------------
def scan_motifs(seq, regex_dict):

    hits = []

    for name,regex in regex_dict.items():

        for m in regex.finditer(seq):

            hits.append((name,m.start(),m.end(),m.group()))

    return hits


# ----------------------------
# OVERLAP
# ----------------------------
def overlap(a_start,a_end,b_start,b_end):

    return max(a_start,b_start) < min(a_end,b_end)


def intersect_motifs_lowP(motifs, regions):

    results = []

    for name,ms,me,seq in motifs:

        for rs,re in regions:

            if overlap(ms,me,rs,re):

                results.append({
                    "Motif":name,
                    "Motif_start":ms,
                    "Motif_end":me,
                    "Sequence":seq,
                    "LowP_start":rs,
                    "LowP_end":re
                })

    return results


# ----------------------------
# MAIN ANALYSIS
# ----------------------------
if uploaded_files:

    all_results = []

    regex_dict = build_nonb_regex()

    for file in uploaded_files:

        sequences = read_sequences(file)

        for seq_id,seq in sequences:

            windows,perplexities = sliding_windows(seq,100)

            low_regions,threshold = bottom_percentile_windows(
                windows,
                perplexities,
                5
            )

            merged = merge_regions(low_regions)

            motifs = scan_motifs(seq,regex_dict)

            overlaps = intersect_motifs_lowP(motifs,merged)

            for r in overlaps:
                r["Sequence_ID"] = seq_id
                all_results.append(r)

    df = pd.DataFrame(all_results)

    st.session_state["results"] = df

    st.success("Analysis complete")


# ----------------------------
# RESULTS
# ----------------------------
if "results" in st.session_state:

    df = st.session_state["results"]

    st.subheader("Detected Motifs in Low Perplexity Regions")

    st.dataframe(df)

    csv = df.to_csv(index=False)

    st.download_button(
        "Download CSV",
        csv,
        "low_perplexity_nonB_results.csv",
        "text/csv"
    )
