#!/bin/bash

output_file=$1

awk -F '\t' 'BEGIN {
    OFS = FS
}
NR == 1 {
    # Print the header with new column names for ANN fields
    for (i = 1; i <= 5; ++i) {
        printf("%s%s", $i, OFS)
    }
    printf("Allele%sAnnotation%sAnnotation_Impact%sGene_Name%sGene_ID%sFeature_Type%sFeature_ID%sTranscript_BioType%sRank%sHGVS.c%sHGVS.p%scDNA.pos%sCDS.pos%sAA.pos%sDistance%sErrors", OFS, OFS, OFS, OFS, OFS, OFS, OFS, OFS, OFS, OFS, OFS, OFS, OFS, OFS, OFS)
    for (i = 7; i <= NF; ++i) {
        printf("%s%s", OFS, $i)
    }
    printf("\n")
    next
}
{
    # Print the data with ANN fields split and without the ANN field
    for (i = 1; i <= 5; ++i) {
        printf("%s%s", $i, OFS)
    }
    split($6, ann, "|")
    for (i = 1; i <= 16; ++i) {
        printf("%s%s", ann[i], OFS)
    }
    for (i = 7; i <= NF; ++i) {
        printf("%s%s", $i, (i == NF) ? "\n" : OFS)
    }
}' > /dev/stdout
