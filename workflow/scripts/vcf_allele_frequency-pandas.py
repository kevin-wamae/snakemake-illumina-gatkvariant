# ***************************************************************
# This code reads in a TSV file using Dask dataframe, performs
# data transformation and saves the result as another TSV file.

# The TSV file is assumed to have columns representing various
# features of a genomic variant, such as chromosome, position,
# reference allele, alternate allele, and so on. The remaining
# columns represent allele depth (AD) for each sample. AD
# represents the number of reads supporting the reference and
# alternate alleles for a given variant in each sample.

# The code then melts the allele depth columns into a long format
# using Dask's melt method. In the long format, each row represents
# a sample and a variant with its AD value. The calculate_proportion
# function is applied to the AD values for each sample to calculate
# the proportion of reads supporting the reference allele. The
# resulting proportions are rounded to 3 decimal places and saved
# to a TSV file using Dask's to_csv method.

# Overall, this code is performing allele frequency calculation for
# each variant across multiple samples.
# ***************************************************************


# ---------------------------------------------------------------
# load the necessary modules
# ---------------------------------------------------------------


import argparse
import pandas as pd


# ---------------------------------------------------------------
# argparse block to parse command-line arguments
# ---------------------------------------------------------------


parser = argparse.ArgumentParser(
    description="Calculate allele frequency for genomic variants."
)
parser.add_argument("input_file", help="Input TSV file path")
parser.add_argument("output_file", help="Output TSV file path")
args = parser.parse_args()


# ---------------------------------------------------------------
# load the data as a Pandas dataframe
# ---------------------------------------------------------------


allele_list = pd.read_csv(
    args.input_file,
    sep="\t",
    dtype={
        "AA.pos": "object",
        "Allele": "object",
        "Annotation": "object",
        "Annotation_Impact": "object",
        "CDS.pos": "object",
        "Distance": "float64",
        "Errors": "object",
        "Feature_ID": "object",
        "Feature_Type": "object",
        "Gene_ID": "object",
        "Gene_Name": "object",
        "HGVS.c": "object",
        "HGVS.p": "object",
        "Rank": "object",
        "Transcript_BioType": "object",
        "cDNA.pos": "object",
    },
)


# ---------------------------------------------------------------
# define a function to calculate the proportion from an allele depth
# (AD) value
# ---------------------------------------------------------------


def calculate_proportion(value):
    """Calculate proportion from allele depth (AD) value."""
    try:
        nums = list(map(int, value.split(",")))
    except ValueError:
        return None
    if sum(nums) == 0 or nums[1] == 0:
        return nums[0] / sum(nums) if sum(nums) != 0 else 0
    return nums[0] / (nums[0] + nums[1])


# ---------------------------------------------------------------
# reshape the dataframe from wide to long format, keeping the first
# 19 columns as id_vars
# ---------------------------------------------------------------


allele_freq = allele_list.melt(
    id_vars=allele_list.columns[:21],
    value_vars=allele_list.columns[21:],
    var_name="sample",
    value_name="AD",
).query("AD != '0,0'")

# ---------------------------------------------------------------
# apply the calculate_proportion function to the AD column and add the
# result as a new column to the dataframe
# ---------------------------------------------------------------


allele_freq["prop"] = allele_freq["AD"].apply(calculate_proportion).round(3)


# ---------------------------------------------------------------
# save the resulting dataframe(s) as a TSV file with tab-separated values
# ---------------------------------------------------------------


columns_to_keep = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "TYPE",
    "Annotation",
    "Gene_ID",
    "Feature_ID",
    "Transcript_BioType",
    "cDNA.pos",
    "CDS.pos",
    "HGVS.c",
    "AA.pos",
    "HGVS.p",
    "sample",
    "AD",
    "prop",
]

# cDNA.pos refers to the position of the variant within the coding DNA (cDNA) sequence of a gene.
# cDNA sequences include both exons and UTR (untranslated) regions of the gene.
# CDS.pos refers to the position of the variant within the coding sequence (CDS) of the gene.
# The CDS contains only the exons that code for the protein, excluding the UTR regions.

# ---------------------------------------------------------------
# Save the resulting dataframe as a TSV file with tab-separated values and no index
# ---------------------------------------------------------------


allele_freq.loc[:, columns_to_keep].to_csv(args.output_file, sep="\t", index=False)