#!/usr/bin/env python
import os
import glob

# change directory to where the fastq files are
os.chdir("input/fastq")

# loop through all files in the directory
for file in glob.glob("*.fastq.gz"):
    split_file = file.split('_')
    # check if the file name has at least 4 '_' characters, before renaming
    # if not, print a message and skip the file
    # the new file name should be the second 3 '_' separated fields
    if len(split_file) >= 4:
        new_file = split_file[0] + '_' + split_file[3] + '.fastq.gz'
        # test with print statement
        print(file, '\t', new_file)
        # rename
        os.rename(file, new_file)
    else:
        print(
            f"Skipping file {file} because it does not contain at least 4 '_' characters.")
