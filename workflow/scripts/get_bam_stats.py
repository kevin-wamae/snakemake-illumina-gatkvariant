import os
import re
import yaml


def main():
    # Read input directory from the config file
    with open("config/config.yaml", "r") as yaml_file:
        config = yaml.safe_load(yaml_file)

    dir_path = os.path.join(config["map_qual_stats"]["dir"], "samtools/flagstat")

    # Create an empty list to store the results
    results = []

    # Loop through all files in the directory
    for file in os.listdir(dir_path):
        if file.endswith(".bam.flagstat.txt"):
            with open(os.path.join(dir_path, file), "r") as f:
                content = f.read()

                # Extract sample ID
                sample_id = re.sub(r"\.bam\.flagstat\.txt$", "", file)

                # Extract total reads
                total_reads = int(
                    re.search(r"^(\d+) \+ \d+ in total", content, re.MULTILINE).group(1)
                )

                # Extract total mapped reads
                mapped_reads = int(
                    re.search(r"^(\d+) \+ \d+ mapped", content, re.MULTILINE).group(1)
                )

                # Calculate total unmapped reads
                unmapped_reads = total_reads - mapped_reads

                # Extract total duplicates
                duplicates = int(
                    re.search(r"^(\d+) \+ \d+ duplicates", content, re.MULTILINE).group(
                        1
                    )
                )

                # Append the results to the list
                results.append([sample_id, mapped_reads, duplicates, unmapped_reads])

    # Save the results to a file
    with open("stats_3_mapped_reads.tsv", "w") as f:
        # Write the header
        f.write(
            "sample_id\treads_mapped_total\treads_mapped_duplicates\treads_mapped_unmapped\n"
        )

        # Write the results in a tab-separated format
        for result in results:
            f.write("\t".join(map(str, result)) + "\n")


if __name__ == "__main__":
    main()
