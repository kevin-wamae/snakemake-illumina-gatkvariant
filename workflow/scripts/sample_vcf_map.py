#!/usr/bin/env python3
import os
import sys

# Get the directory path from the command line arguments
vcf_dir = sys.argv[1]
tsv_file = sys.argv[2]

# Get a list of all the VCF files in the directory
vcf_files = [f for f in os.listdir(vcf_dir) if f.endswith('.vcf.gz')]

# Create a list of tuples with the file paths and names
file_list = []
for vcf_file in vcf_files:
    file_name = os.path.splitext(vcf_file)[0]
    file_path = os.path.join(vcf_dir, vcf_file)
    file_list.append((file_path, file_name))

# Sort the list by the second element of each tuple (the file name)
file_list.sort(key=lambda x: x[1])

# Write the file names and paths to a TSV file
with open(tsv_file, 'w') as tsv_file:
    for file_path, file_name in file_list:
        tsv_file.write(os.path.splitext(file_name)[
                       0] + '\t' + file_path + '\n')
