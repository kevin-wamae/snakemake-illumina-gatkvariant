import os
import re
import pandas as pd
from subprocess import Popen, PIPE
import yaml


def get_stats(fastq_file):
    process = Popen(
        f"seqkit stats {fastq_file}", stdout=PIPE, stderr=PIPE, shell=True, text=True
    )
    header, data = process.communicate()[0].strip().split("\n")
    stats = {
        k: v.replace(",", "")
        for k, v in zip(header.split(), data.split())
        if k in ["num_seqs", "min_len", "avg_len", "max_len"]
    }
    return stats


# Read input directory from the config file
with open("config/config.yaml", "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

input_dir = config["input"]["fastq"]

output_file = "stats_1_raw_fastq.tsv"

stats_list = [
    {
        **get_stats(os.path.join(input_dir, f)),
        "id": re.sub(r"_R?[12]\.fastq\.gz", "", f),
    }
    for f in os.listdir(input_dir)
    if f.endswith(("_R1.fastq.gz", "_R2.fastq.gz", "_1.fastq.gz", "_2.fastq.gz"))
]

df = pd.DataFrame(stats_list).astype(
    {"num_seqs": int, "min_len": int, "avg_len": float, "max_len": int}
)
output_df = (
    df.groupby("id")
    .agg({"num_seqs": "sum", "min_len": "min", "avg_len": "mean", "max_len": "max"})
    .reset_index()
)
output_df.columns = [
    "id",
    "reads_raw_total",
    "reads_raw_len_min",
    "reads_raw_len_avg",
    "reads_raw_len_max",
]

output_df.to_csv(output_file, sep="\t", index=False)
