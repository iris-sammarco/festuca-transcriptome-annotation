# Date: 03/2026
# Author: Iris Sammarco
# Aim: Setup guide for the Festuca rubra transcriptome annotation pipeline using Trinotate. 
# Covers: conda environment creation, EggNOG/Pfam-A/Trinotate database downloads, SignalP v6 and TMHMM v2 manual installation, and the full pipeline run order.
# This script is NOT meant to be executed as a whole — run each section manually.
# Run: bash trinotate_00_setup.bash  (source sections manually as needed)
# Output: trinotate_env conda environment
#         ${TRINOTATE_DATA_DIR}/eggnog_data/ (EggNOG databases)
#         ${TRINOTATE_DATA_DIR}/Pfam-A.hmm (Pfam HMM database)
#         ${TRINOTATE_DATA_DIR}/Trinotate.sqlite (Trinotate boilerplate SQLite DB)

# =============================================================================
# USER CONFIGURATION — edit these before running anything
# =============================================================================
export CONDA_ENVS_PATH=/path/to/conda/envs
export TRINOTATE_DATA_DIR=/path/to/project/assembly/trinotate_data
OUTDIR=/path/to/project/assembly/trinotate_output
FINAL_DIR=/path/to/project/assembly/final_annotation
TMHMM_BIN=/path/to/project/tmhmm-2.0c/bin
# =============================================================================

# -----------------------------------------------------------------------------
# 1. CONDA ENVIRONMENT
# -----------------------------------------------------------------------------
module add mambaforge

mamba create -n trinotate_env -c bioconda -c conda-forge \
  trinotate transdecoder blast hmmer diamond sqlite trinity infernal \
  -y

mamba install -c bioconda eggnog-mapper -y # needed for improving the GO/KEGG annotations
mamba install seqkit -y # needed to split the fasta files in several chunks in order to run Pfam in parallel

# -----------------------------------------------------------------------------
# 2. EGGNOG DATABASES
# NOTE: download into ${TRINOTATE_DATA_DIR}/eggnog_data
# -----------------------------------------------------------------------------
EGGNOG_DATA="${TRINOTATE_DATA_DIR}/eggnog_data"
mkdir -p "${EGGNOG_DATA}"
cd "${EGGNOG_DATA}"

wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz  
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz

gunzip eggnog.db.gz eggnog_proteins.dmnd.gz
tar -zxf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz

# -----------------------------------------------------------------------------
# 3. PFAM-A + TRINOTATE BOILERPLATE
# -----------------------------------------------------------------------------
cd "${TRINOTATE_DATA_DIR}"

wget -N "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
gunzip Pfam-A.hmm.gz

# Locate Trinotate home and create the boilerplate SQLite DB.
# NOTE: in the Perl script, the eggnog download block is broken — comment out
# the block around: wget http://eggnogdb.embl.de/.../NOG.annotations.tsv.gz
# and the downstream gunzip/import steps. EggNOG is installed separately above.
export TRINITY_HOME=$(dirname $(dirname $(which Trinity))) || export TRINITY_HOME=${CONDA_PREFIX}
export TRINOTATE_HOME=$(dirname $(which Trinotate))

${TRINOTATE_HOME}/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
# Downloads ~2-5 GB (SwissProt/Pfam), takes 10-30 min.
# The output name will be Trinotate.sqlite

# -----------------------------------------------------------------------------
# 4. SIGNALP v6
# Download the fast model from: https://services.healthtech.dtu.dk/service.php?SignalP
# (fill the form, DTU will email a link valid for 4h; download the fast version,
# recommended by developers for de novo full-transcriptome annotations)
# Then install into the trinotate_env:
#
tar -xzvf signalp-6.0*.tar.gz
cd signalp-6-package
pip install --no-cache-dir .
#
# After install, find the package location with:
python3 -c "import signalp, os; print(os.path.dirname(signalp.__file__))"
# Then copy the fast model weights into the model_weights/ subdirectory there.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 5. TMHMM v2
# Download from: https://services.healthtech.dtu.dk/service.php?TMHMM-2.0
# (fill the form, DTU will email a link)
# Then:
tar -xzvf tmhmm-2.0c.Linux.tar.gz
# Edit the scripts tmhmm and tmhmmformat.pl following the Trinotate documentation:
# https://github.com/Trinotate/Trinotate/wiki/Software-installation-and-data-required
# The bin/ directory must be added to PATH — this is handled inside each run script.
# -----------------------------------------------------------------------------

# =============================================================================
# RUN ORDER
# All scripts live in pipeline/ except the two independent ones in independent/
# Run from the directory where the scripts are located, or provide full paths.
# =============================================================================

# --- INDEPENDENT (run at any time, do not block the main pipeline) ---
cd "${OUTDIR}"
qsub independent/trinotate_IND_blastx.bash  # BLASTx SwissProt; takes ~minutes
qsub independent/trinotate_IND_rfam.bash    # Rfam cmscan; takes ~200 h. It can be parallelized if needed

# --- STEP 01: gene→transcript map, TransDecoder LongOrfs, Diamond BLASTP ---
cd "${OUTDIR}"
qsub pipeline/trinotate_01_longorfs_blastp.bash

# --- STEP 02 (parallel): Pfam hmmscan ---
# Split input into N chunks, submit PBS array, then merge when all jobs finish.
# Adjust N_CHUNKS to your dataset size (used 141 here).
cd "${OUTDIR}"
bash pipeline/trinotate_02_parallel_pfam_eggnog.bash pfam 141    # split + submit
bash pipeline/trinotate_02_parallel_pfam_eggnog.bash pfam 141 m  # merge → pfam.domtblout

# --- STEP 03: TransDecoder Predict (requires Pfam output from step above) ---
cd "${OUTDIR}"
qsub pipeline/trinotate_03_transdecoder.bash

# --- EggNOG annotation (requires Trinity.fasta.transdecoder.pep from step 03) ---
# Adjust N_CHUNKS to your dataset size (used 72 here: less than Pfam because Pfam rns on longest_orfs.pep which is larger; EggNOG runs on the final transdecoder.pep which is smaller after TransDecoder.Predict filtering).
cd "${OUTDIR}"
bash pipeline/trinotate_03_parallel_pfam_eggnog.bash eggnog 72    # split + submit
bash pipeline/trinotate_03_parallel_pfam_eggnog.bash eggnog 72 m  # merge → eggnog.emapper.annotations

# --- STEP 05: Trinotate INIT + TMHMM ---
cd "${OUTDIR}"
qsub pipeline/trinotate_04_tmhmm.bash

# --- STEP 06 (parallel): SignalP6 ---
cd "${OUTDIR}"
bash pipeline/trinotate_05_parallel_signalp6.bash 72    # split + submit
bash pipeline/trinotate_05_parallel_signalp6.bash 72 m  # merge after all jobs finish

# --- STEPS 07-08: Poaceae curation (Rice + Wheat BLASTP) ---
cd "${OUTDIR}"
qsub pipeline/trinotate_06_blastp_rice.bash
qsub pipeline/trinotate_07_blastp_wheat.bash

# --- STEP 09: Load all results + generate Trinotate report ---
cd "${OUTDIR}"
qsub pipeline/trinotate_08_load_report.bash

# --- STEPS 10-12: Merge annotations, GO and KEGG background ---
cd "${FINAL_DIR}"
cp -s "${OUTDIR}/Festuca_rubra.annotation_report_strict.xls.gz" .
gunzip Festuca_rubra.annotation_report_strict.xls.gz
trinotate_report_summary.pl Festuca_rubra.annotation_report_strict.xls Festuca_rubra_report_strict

cp -s "${OUTDIR}/eggnog.emapper.annotations" .
cp -s "${OUTDIR}/rfam.tblout" .

cd "${FINAL_DIR}"
# Run in PBS if needed:
python3 pipeline/trinotate_09_merge_annotations.py
python3 pipeline/trinotate_10_go_background.py
python3 pipeline/trinotate_11_kegg_background.py

# =============================================================================
# POST-RUN: CLEAN UP REPORT + SANITY CHECKS
# Run manually in ${FINAL_DIR} after all steps complete.
# =============================================================================

## Remove empty columns from report (rnammer col 4, eggnog col 13, transcript col 18, peptide col 19; these columns are empty as: RNAmmer was not run; EggNOG was loaded via Python merge, not Trinotate; transcript/peptide sequences are excluded to keep file size manageable) and create a final report
# Run all checks above BEFORE compressing
awk -F'\t' 'BEGIN{OFS="\t"} {
    f=""; for(i=1;i<=NF;i++) if(i!=4 && i!=13 && i!=18 && i!=19) f = f (f==""?"":"\t") $i; print f
}' Festuca_rubra_annotation_report.extended.tsv > Festuca_rubra_annotation_report_final.tsv

## Check the report:
# Count total rows and unannotated genes:
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

# Per-tool annotation coverage: How many genes have Pfam, GO, KEGG, eggnog, SignalP, TmHMM?
awk -F'\t' '$1!~"^#" {print $9}' Festuca_rubra_annotation_report_final.tsv | grep -v "^.$" | wc -l    # Pfam 349266, ~10% of genes
awk -F'\t' '$1!~"^#" {print $13}' Festuca_rubra_annotation_report_final.tsv | grep -v "^.$" | wc -l   # KEGG 416087, ~12%
awk -F'\t' '$1!~"^#" {print $14}' Festuca_rubra_annotation_report_final.tsv | grep -v "^.$" | wc -l   # GO_BLASTX 260528, ~7%
awk -F'\t' '$1!~"^#" {print $15}' Festuca_rubra_annotation_report_final.tsv | grep -v "^.$" | wc -l   # GO_BLASTP 194487, ~7%
awk -F'\t' '$1!~"^#" {print $16}' Festuca_rubra_annotation_report_final.tsv | grep -v "^.$" | wc -l   # GO_Pfam 3549821
awk -F'\t' '$1!~"^#" {print $17}' Festuca_rubra_annotation_report_final.tsv | grep -v "^.$" | wc -l   # eggnog 3236549, ~91%
awk -F'\t' '$1!~"^#" {print $10}' Festuca_rubra_annotation_report_final.tsv | grep -c "sigP"        # SignalP 710507, ~20%
awk -F'\t' '$1!~"^#" {print $11}' Festuca_rubra_annotation_report_final.tsv | grep -c "ExpAA"       # TmHMM 96395, ~3% 

# Sanity check GO / KEGG background files
wc -l Festuca_rubra_go_background.cleaned.tsv # 247029
wc -l Festuca_rubra_kegg_background.cleaned.tsv # 57265

# Count how many genes have GO vs KEGG
cut -f1 Festuca_rubra_go_background.cleaned.tsv | tail -n +2 | wc -l # 247028
cut -f1 Festuca_rubra_kegg_background.cleaned.tsv | tail -n +2 | wc -l # 57264

# Compress the file
gzip Festuca_rubra_annotation_report_final.tsv

