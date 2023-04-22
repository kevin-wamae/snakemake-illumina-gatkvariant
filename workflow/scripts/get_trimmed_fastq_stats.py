import os
import re
import pandas as pd
from subprocess import Popen, PIPE
import yaml


def get_stats(fastq_file, prefix):
    process = Popen(
        f"seqkit stats {fastq_file}", stdout=PIPE, stderr=PIPE, shell=True, text=True
    )
    header, data = process.communicate()[0].strip().split("\n")
    stats = {
        f"{prefix}_{k}": v.replace(",", "")
        for k, v in zip(header.split(), data.split())
        if k in ["num_seqs", "min_len", "avg_len", "max_len"]
    }
    return stats


def main():
    config_file = "config/config.yaml"
    with open(config_file, "r") as yaml_file:
        config = yaml.safe_load(yaml_file)

    output_dir = config["trim_reads"]["dir"]

    paired_dir = os.path.join(output_dir, "paired")
    unpaired_dir = os.path.join(output_dir, "unpaired")
    output_file = "stats_2_trimmed_fastq.tsv"

    paired_stats_list = []
    unpaired_stats_list = []

    for f in os.listdir(paired_dir):
        if f.endswith(
            (
                "_R1.trimmed.fastq.gz",
                "_R2.trimmed.fastq.gz",
                "_1.trimmed.fastq.gz",
                "_2.trimmed.fastq.gz",
            )
        ):
            paired_stats_list.append(
                {
                    **get_stats(os.path.join(paired_dir, f), "trimmed"),
                    "id": re.sub(r"_R?[12].trimmed.fastq.gz", "", f),
                }
            )

    for f in os.listdir(unpaired_dir):
        if f.endswith(
            (
                "_R1.unpaired.fastq.gz",
                "_R2.unpaired.fastq.gz",
                "_1.unpaired.fastq.gz",
                "_2.unpaired.fastq.gz",
            )
        ):
            unpaired_stats_list.append(
                {
                    **get_stats(os.path.join(unpaired_dir, f), "unpaired"),
                    "id": re.sub(r"_R?[12].unpaired.fastq.gz", "", f),
                }
            )

    df_paired = pd.DataFrame(paired_stats_list).astype(
        {
            "trimmed_num_seqs": int,
            "trimmed_min_len": int,
            "trimmed_avg_len": float,
            "trimmed_max_len": int,
        }
    )
    df_unpaired = pd.DataFrame(unpaired_stats_list).astype(
        {
            "unpaired_num_seqs": int,
            "unpaired_min_len": int,
            "unpaired_avg_len": float,
            "unpaired_max_len": int,
        }
    )

    output_df_paired = (
        df_paired.groupby("id")
        .agg(
            {
                "trimmed_num_seqs": "sum",
                "trimmed_min_len": "min",
                "trimmed_avg_len": "mean",
                "trimmed_max_len": "max",
            }
        )
        .reset_index()
    )
    output_df_unpaired = (
        df_unpaired.groupby("id")
        .agg(
            {
                "unpaired_num_seqs": "sum",
                "unpaired_min_len": "min",
                "unpaired_avg_len": "mean",
                "unpaired_max_len": "max",
            }
        )
        .reset_index()
    )

    output_df = output_df_paired.merge(output_df_unpaired, on="id")

    column_mapping = {
        "trimmed_num_seqs": "reads_trim_total",
        "trimmed_min_len": "reads_trim_len_min",
        "trimmed_avg_len": "reads_trim_len_avg",
        "trimmed_max_len": "reads_trim_len_max",
        "unpaired_num_seqs": "reads_unpaired_total",
        "unpaired_min_len": "reads_unpaired_len_min",
        "unpaired_avg_len": "reads_unpaired_len_avg",
        "unpaired_max_len": "reads_unpaired_len_max",
    }

    output_df = output_df.rename(columns=column_mapping)

    output_df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
