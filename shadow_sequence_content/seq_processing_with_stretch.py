import re
from Bio import SeqIO
import pandas as pd
from collections import Counter
import argparse
import subprocess

# ----------------------------
# Regex to extract header info
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
# GC and AT content
# ----------------------------
def gc_at_percent(seq):
    seq = seq.upper()
    counts = Counter(seq)
    valid = sum(counts[b] for b in "ATGC")
    if valid == 0:
        return 0.0, 0.0
    gc = (counts["G"] + counts["C"]) / valid * 100
    at = (counts["A"] + counts["T"]) / valid * 100
    return round(gc,3), round(at,3)

# ----------------------------
# Run Perl GC-stretch and capture output as list of lines
# ----------------------------
def run_gc_stretch_perl(fasta_file, perl_script_path):
    """
    Runs Perl GC-stretch script and returns a list of dictionaries per line.
    """
    try:
        output = subprocess.check_output(
            ["perl", perl_script_path, fasta_file],
            stderr=subprocess.STDOUT
        ).decode("utf-8").strip()
        lines = output.splitlines()
    except subprocess.CalledProcessError as e:
        print(f"Error running Perl script:\n{e.output.decode()}")
        return []

    dicts = []
    for line in lines:
        line = line.strip()
        if line.startswith("Species:") and "Average Length:" in line:
            try:
                parts = line.split(",")
                species_name = parts[0].split("Species:")[1].strip()
                avg_len = float(parts[1].split("Average Length:")[1].strip())
                dicts.append({"Species": species_name, "Average_Length": avg_len})
            except Exception:
                dicts.append({})
        else:
            dicts.append({})
    return dicts

# ----------------------------
# Parse FASTA, add GC/AT, and append Perl output
# ----------------------------
def parse_fasta_to_csv(fasta_path, out_csv, perl_script_path):
    perl_output_dicts = run_gc_stretch_perl(fasta_path, perl_script_path)

    rows = []
    for i, record in enumerate(SeqIO.parse(fasta_path, "fasta")):
        meta = parse_header(record.id)
        seq = str(record.seq)
        gc, at = gc_at_percent(seq)

        # Take corresponding Perl dictionary (if available)
        perl_dict = perl_output_dicts[i] if i < len(perl_output_dicts) else {}

        avg_len = perl_dict.get("Average_Length", None)  # extract value of Average_Length

        meta.update({
            "sequence_length": len(seq),
            "GC_percent": gc,
            "AT_percent": at,
            "GC_stretch": avg_len if avg_len is not None else 0
        })
        rows.append(meta)

    pd.DataFrame(rows).to_csv(out_csv, index=False)

# ----------------------------
# CLI
# ----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse FASTA, compute GC/AT, and append Perl GC-stretch output directly"
    )
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--out", required=True, help="Output CSV file")
    parser.add_argument("--perl", required=True, help="Path to Perl GC-stretch script")
    args = parser.parse_args()

    parse_fasta_to_csv(args.fasta, args.out, args.perl)
