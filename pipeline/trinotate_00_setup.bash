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

# Create boilerplate for Trinotate
export TRINITY_HOME=$(dirname $(dirname $(which Trinity))) || export TRINITY_HOME=${CONDA_PREFIX}
export TRINOTATE_HOME=$(dirname $(which Trinotate))

${TRINOTATE_HOME}/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate_new # In the pl script the eggnog path is broken: comment out the block that downloads and loads eggNOG, specifically the commands around wget http://eggnogdb.embl.de/.../NOG.annotations.tsv.gz and the downstream gunzip/import steps
# I've installed eggnog separately within the mamba environment (downloads ~2-5GB SwissProt/Pfam, takes 10-30min). The command should have included Trinotate as output rather than Trinotate.sqlite, Ive renamed the files to: Trinotate.sqlite  Trinotate.TaxonomyIndex  Trinotate.UniprotIndex
 
### Run the scripts (they're divided into several scripts as some steps are very long):
cd /storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_data
qsub trinotate_01_main_start.bash # gene→transcript map, TransDecoder LongOrfs, Direct Diamond SwissProt BLASTP
qsub trinotate_XX_blastx.bash # independent from the others, run at anytime; it took a few minutes to run
qsub trinotate_XX_rfam.bash # independent from the others, run at anytime; it took around 200 h to finish

# Pfam
cd /storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_output
bash trinotate_03_05_parallel_driver_pfam_eggnog.bash pfam 141    # split pfam input file into 141 chunks and runs pfam on each chunk; relies on trinotate_03a_pfam_array.pbs
bash trinotate_03_05_parallel_driver_pfam_eggnog.bash pfam 141 m # merges the output chunks, creates pfam.domtblout; relies on trinotate_03a_pfam_array.pbs

# transdecoder.predict
cd /storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_data
qsub trinotate_04_transdecoder.bash # run transdecoder.predict step (it relies on Pfam output and it's needed by eggnog)

# eggnog
cd /storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_output
bash trinotate_03_05_parallel_driver_pfam_eggnog.bash eggnog 72 # split eggnog input file and runs eggnog separately on each file; relies on trinotate_05a_eggnog_array.pbs
bash trinotate_03_05_parallel_driver_pfam_eggnog.bash eggnog 72 m  # merges the chunks, creates eggnog.emapper.annotations; relies on trinotate_05a_eggnog_array.pbs

# INIT, TMHMM 
cd /storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_data
qsub trinotate_06_init_tmhmm

# SignalP6
cd /storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_output
bash trinotate_07_parallel_driver_signalp6.bash 72 # Submit chunks
bash trinotate_07_parallel_driver_signalp6.bash 72 m # Merge and load after all jobs finish

# Rice/Wheat Poaceae Curation, final report
cd /storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_data
qsub trinotate_08_rice_curation.bash
qsub trinotate_09_wheat_curation.bash
qsub trinotate_10_load_final

# Run report summary from Trinotate:
cd /storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/final_annotation

cp /storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_output/Festuca_rubra.annotation_report_strict.xls.gz .
gunzip Festuca_rubra.annotation_report_strict.xls.gz

/storage/pruhonice1-ibot/home/irissammarco/.conda/envs/trinotate_env/bin/trinotate_report_summary.pl Festuca_rubra.annotation_report_strict.xls Festuca_rubra_report_strict

# Merge eggnog and rfam annotations to the final report, and create a list of GO and KEGG terms merging all the hits from several tools:
mamba activate python

cd /storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/final_annotation
cp -s /storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_output/eggnog.emapper.annotations .
cp -s /storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_output/rfam.tblout .

## Run in the queue without making an external pbs job
# Merge eggnog and rfam annotations to Trinotate report:
echo -e '#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=2:mem=8gb:scratch_local=20gb
#PBS -N Festuca_merge_annotations
#PBS -j oe

cd "$PBS_O_WORKDIR"
module add mambaforge
export CONDA_ENVS_PATH=/storage/pruhonice1-ibot/home/irissammarco/.conda/envs
mamba activate python

python3 Festuca_11_merge_annotations_filtered.py
' | qsub

# Combine GO terms from BLASTX, BLASTP, Pfam, and eggNOG:
echo -e '#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=2:mem=8gb:scratch_local=20gb
#PBS -N Festuca_GO_background
#PBS -j oe

cd "$PBS_O_WORKDIR"
module add mambaforge
export CONDA_ENVS_PATH=/storage/pruhonice1-ibot/home/irissammarco/.conda/envs
mamba activate python

python3 Festuca_12_GO_background.py
' | qsub

# Combine KEGG terms from BLASTX, BLASTP, Pfam, and eggNOG:
echo -e '#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=2:mem=8gb:scratch_local=20gb
#PBS -N Festuca_KEGG_background
#PBS -j oe

cd "$PBS_O_WORKDIR"
module add mambaforge
export CONDA_ENVS_PATH=/storage/pruhonice1-ibot/home/irissammarco/.conda/envs
mamba activate python

python3 Festuca_13_KEGG_background.py
' | qsub

## Remove empty columns from report (rnammer col 4, eggnog col 13, transcript col 18, peptide col 19) and create a final report:
awk -F'\t' 'BEGIN{OFS="\t"} {
    f=""; for(i=1;i<=NF;i++) if(i!=4 && i!=13 && i!=18 && i!=19) f = f (f==""?"":"\t") $i; print f
}' Festuca_rubra_annotation_report.extended.tsv > Festuca_rubra_annotation_report_final.tsv

## Check the report:
# ① Count total genes and empty genes
wc -l Festuca_rubra_annotation_report_final.tsv # 3549822
awk -F'\t' '$1!~"^#" && NF>0 {print}' Festuca_rubra_annotation_report_final.tsv | wc -l # non-header and non-empty rows: 3549821

awk -F'\t' '
  BEGIN{m=0}
  $1!~"^#" {
    empty=1
    for(i=3;i<=NF;i++) {
      if($i != "." && $i != "") {empty=0; break}
    }
    if(empty) m++
  }
  END{print "Genes with only dots from column 3 onward:", m}
' Festuca_rubra_annotation_report_final.tsv # genes with no BLAST hits at all 2701729, ~76%

# ② How many genes have Pfam, GO, KEGG, eggnog, SignalP, TmHMM?
awk -F'\t' '$1!~"^#" {print $9}' Festuca_rubra_annotation_report_final.tsv | grep -v "^.$" | wc -l    # Pfam 349266, ~10% of genes
awk -F'\t' '$1!~"^#" {print $13}' Festuca_rubra_annotation_report_final.tsv | grep -v "^.$" | wc -l   # KEGG 416087, ~12%
awk -F'\t' '$1!~"^#" {print $14}' Festuca_rubra_annotation_report_final.tsv | grep -v "^.$" | wc -l   # GO_BLASTX 260528, ~7%
awk -F'\t' '$1!~"^#" {print $15}' Festuca_rubra_annotation_report_final.tsv | grep -v "^.$" | wc -l   # GO_BLASTP 194487, ~7%
awk -F'\t' '$1!~"^#" {print $16}' Festuca_rubra_annotation_report_final.tsv | grep -v "^.$" | wc -l   # GO_Pfam 3549821
awk -F'\t' '$1!~"^#" {print $17}' Festuca_rubra_annotation_report_final.tsv | grep -v "^.$" | wc -l   # eggnog 3236549, ~91%
awk -F'\t' '$1!~"^#" {print $10}' Festuca_rubra_annotation_report_final.tsv | grep -c "sigP"        # SignalP 710507, ~20%
awk -F'\t' '$1!~"^#" {print $11}' Festuca_rubra_annotation_report_final.tsv | grep -c "ExpAA"       # TmHMM 96395, ~3% 

# 2) Sanity check GO / KEGG background files
wc -l Festuca_rubra_go_background.cleaned.tsv # 247029
wc -l Festuca_rubra_kegg_background.cleaned.tsv # 57265

# Count how many genes have GO vs KEGG
cut -f1 Festuca_rubra_go_background.cleaned.tsv | tail -n +2 | wc -l # 247028
cut -f1 Festuca_rubra_kegg_background.cleaned.tsv | tail -n +2 | wc -l # 57264

# Compress the file
gzip Festuca_rubra_annotation_report_final.tsv

