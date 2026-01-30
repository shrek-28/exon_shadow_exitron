import re
from Bio import SeqIO
import pandas as pd
import argparse

# ----------------------------
# Regex to extract desired fields
# ----------------------------
HEADER_RE = re.compile(
    r"""
    ^.*_                # everything before last underscore
    (?P<frame>-?\d+)    # frame, allow negative
    ::
    (?P<chromosome>[^:]+)       # chromosome
    :
    (?P<shadow_start>\d+)-(?P<shadow_end>\d+)  # shadow coordinates
    \((?P<strand>[+-])\)        # strand
    """,
    re.VERBOSE
)

def parse_header(header):
    header = header.lstrip(">").strip()
    match = HEADER_RE.search(header)
    if not match:
        raise ValueError(f"Malformed header:\n{header}")
    return match.groupdict()

# ----------------------------
# GC and AT content calculation
# ----------------------------
def gc_at_percent(seq):
    seq = seq.upper()
    from collections import Counter
    counts = Counter(seq)
    valid = sum(counts[b] for b in "ATGC")
    if valid == 0:
        return 0.0, 0.0
    gc = (counts["G"] + counts["C"]) / valid * 100
    at = (counts["A"] + counts["T"]) / valid * 100
    return round(gc,3), round(at,3)

# ----------------------------
# FASTA parsing
# ----------------------------
def parse_fasta_to_csv(fasta_path, out_csv):
    rows = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        meta = parse_header(record.id)
        gc, at = gc_at_percent(str(record.seq))
        meta.update({
            "sequence_length": len(record.seq),
            "GC_percent": gc,
            "AT_percent": at
        })
        rows.append(meta)
    pd.DataFrame(rows).to_csv(out_csv, index=False)

# ----------------------------
# CLI
# ----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract frame, chromosome, shadow coordinates and strand from FASTA headers")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--out", required=True, help="Output CSV file")
    args = parser.parse_args()
    parse_fasta_to_csv(args.fasta, args.out)
