# 02/03/2026
# Author: Iris Sammarco
# Purpose: annotate Festuca's assembled transcriptome with Trinotate using rice and wheat genomes

# Trinotate
# https://github.com/Trinotate/Trinotate/wiki/Software-installation-and-data-required

# Tutorials:
# https://neatseq-flow.readthedocs.io/projects/neatseq-flow-modules/en/latest/Workflow_docs/RNA_seq_Trinity.html
# https://docs.hpc.ufs.ac.za/training/transcriptomics_tutorial/transcriptomics_tutorial_1/

## Install all the tools in a mamba environment
module add mambaforge
export CONDA_ENVS_PATH=/path/to/conda/envs

mamba create -n trinotate_env -c bioconda -c conda-forge \
  trinotate transdecoder blast hmmer diamond sqlite trinity infernal \
  -y

mamba install -c bioconda eggnog-mapper -y # needed for improving the GO/KEGG annotations
mamba install seqkit -y # needed to split the fasta files in several chunks in order to run Pfam in parallel

# Download EggNOG databases
export TRINOTATE_DATA_DIR=/path/to/project/assembly/trinotate_data
EGGNOG_DATA="${TRINOTATE_DATA_DIR}/eggnog_data"

mkdir -p "${EGGNOG_DATA}"
cd "${EGGNOG_DATA}"

wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz  
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz

gunzip eggnog.db.gz eggnog_proteins.dmnd.gz
tar -zxf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz

## Install Pfam-A and create boilerplate for Trinotate (needed for Trinotate --create). This will create the output files:
# Trinotate.sqlite (boilerplate DB)
# uniprot_sprot.fa (for manual Diamond if needed)
cd /path/to/project/assembly/trinotate_data

# Pfam-A
wget -N "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
gunzip Pfam-A.hmm.gz
