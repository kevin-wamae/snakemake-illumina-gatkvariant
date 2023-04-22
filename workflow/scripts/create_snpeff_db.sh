# *****************************************************************
# a bash script for creating a snpEff database in the output directory
# last update worked with SnpEff 5.1d (build 2022-04-19 15:49)
# *****************************************************************

# *****************************************************************
# we are going to create a P.falciparum database using data from PlasmoDB
# *****************************************************************

# variables for downloading genome files
GENOME_URL='https://plasmodb.org/common/downloads/release-51/Pfalciparum3D7/'
GENOME_FASTA='fasta/data/PlasmoDB-51_Pfalciparum3D7_Genome.fasta'
GENOME_CDS='fasta/data/PlasmoDB-51_Pfalciparum3D7_AnnotatedCDSs.fasta'
GENOME_PROT='fasta/data/PlasmoDB-51_Pfalciparum3D7_AnnotatedProteins.fasta'
GENOME_GFF='gff/data/PlasmoDB-51_Pfalciparum3D7.gff'

# *****************************************************************
# variables for snpEff database
#   - replace ORGANISM value  with the name of your organism
# *****************************************************************
DATADIR='output/1_snpeff_database'
ORGANISM='P.falciparum'

# *****************************************************************
# createsnpEff working directory inside DATADIR and download the
# genome files
# *****************************************************************
# create working directories
mkdir -p $DATADIR/genomes
mkdir -p $DATADIR/$ORGANISM

# download genome
# -----------------------------------------------------------------
wget --continue \
    $GENOME_URL$GENOME_FASTA \
    -O $DATADIR/genomes/$ORGANISM.fa

# download coding sequences
# -----------------------------------------------------------------
wget --continue \
    $GENOME_URL$GENOME_CDS \
    -O $DATADIR/$ORGANISM/cds.fa

# download protein sequences
# -----------------------------------------------------------------
wget --continue \
    $GENOME_URL$GENOME_PROT \
    -O $DATADIR/$ORGANISM/protein.fa

# ensure that the gene-names in cds.fa & protein.fa are identical
# -----------------------------------------------------------------
# 'snpEff build' will fail if the names don't match, e.g., we will
# - drop, the suffix '-p1' in protein.fa to match the names in cds.fa
if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' -e "s/1-p1/1/" $DATADIR/$ORGANISM/protein.fa
else
    sed -i -e "s/1-p1/1/" $DATADIR/$ORGANISM/protein.fa
fi

# download gff annotation files
# -----------------------------------------------------------------
wget --continue \
    $GENOME_URL$GENOME_GFF \
    -O $DATADIR/$ORGANISM/genes.gff

# create snpEff config file and populate with:
#   - DATADIR
#   - location of genome files
# -----------------------------------------------------------------
# echo "data.dir = $DATADIR" > $DATADIR/snpEff.config
# echo "$ORGANISM.genome : $ORGANISM" >> $DATADIR/snpEff.config
echo "data.dir = $DATADIR" >.snpEff.config
echo "$ORGANISM.genome : $ORGANISM" >>.snpEff.config

# finally, create snpEff database using the -gff3 option
# -----------------------------------------------------------------
snpEff build -gff3 -config .snpEff.config -v $ORGANISM
