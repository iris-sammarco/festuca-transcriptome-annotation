#!/bin/bash

#PBS -l walltime=200:00:0
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=28:mem=40gb:scratch_local=200gb
#PBS -N Festuca_Rfam
#PBS -j oe

# Author: Iris Sammarco
# Date: 06/03/2026
# Aim: Run the ncRNA step separately from the main script as it takes long. This step is independent from the other steps, so it can be run at any time before the final load
# Run: qsub Trinotate_rfam.bash

set -euo pipefail  # Exit on error, undefined vars, pipe failures
trap 'echo "[ERROR] Line $LINENO: ${BASH_COMMAND} failed on $(date)" >&2; exit 1' ERR # debug

# Environment & Paths
export TMPDIR=$SCRATCHDIR # This will force the application to place the temporary files into scratch directory instead of a /tmp directory
module add mambaforge
export CONDA_ENVS_PATH=/storage/pruhonice1-ibot/home/irissammarco/.conda/envs
mamba activate trinotate_env
export TRINITY_HOME="/storage/pruhonice1-ibot/home/irissammarco/.conda/envs/trinotate_env/opt/trinity-2.8.5"
export PATH="/storage/plzen1/home/irissammarco/Festuca_RNA_assembly/tmhmm-2.0c/bin:$PATH"
export TRINOTATE_DATA_DIR="/storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_data" # expects the files: Trinotate.sqlite and uniprot_sprot.fa
export EGGNOG_DATA_DIR="${TRINOTATE_DATA_DIR}/eggnog_data" 

# Config
ASSEMBLY="/storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/Trinity_output.Trinity.fasta"
PREFIX="Festuca_rubra"
OUTDIR="/storage/plzen1/home/irissammarco/Festuca_RNA_assembly/assembly/trinotate_output"
LOGDIR="${OUTDIR}/logs"
THREADS=28
EVAL="1e-5"  # Strict Diamond E-value
MAX_TARGETS=5  # Maximum number of target sequences to report alignments for (default=25)

mkdir -p "${OUTDIR}" "${LOGDIR}"
cd "${OUTDIR}"

# Per script logs:
exec 1> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.log")
exec 2> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.err")
echo "[INFO] Job ${PBS_JOBID} started in ${PWD} | Threads: ${THREADS}"

# Create symlinks of the db files in the OUTDIR:
ln -sf "${TRINOTATE_DATA_DIR}"/Trinotate.sqlite* .
ln -sf "${TRINOTATE_DATA_DIR}"/Pfam-A.hmm* .
ln -sf "${TRINOTATE_DATA_DIR}"/uniprot_sprot.pep .
ln -sf "${TRINOTATE_DATA_DIR}"/{uniprot_sprot.diamond,Rfam.cm,Rfam.clanin} . 2>/dev/null || true

## Sanity checks (fail-fast)
[[ -s "${ASSEMBLY}" ]] || { echo "[FATAL] Missing assembly: ${ASSEMBLY}"; exit 1; }
ln -sf "${ASSEMBLY}" Trinity.fasta

echo "[INFO] Rfam-only job started $(date)"

## Step 4C: Infernal cmscan Rfam for ncRNAs
RFAM_CM="${TRINOTATE_DATA_DIR}/Rfam.cm"
if [[ ! -s "rfam.tblout" || ! -s ".rfam.done" ]]; then
	# Calculate Z:
    echo "[INFO] Calculating Z for accurate E-values..." # For cmscan, Z is the length of the current query sequence [nucleotides] multiplied by 2 (because both strands of the sequence are searched [in each model]) and multiplied again by the number of CMs [covariance models] in the target CM database.
	# --Z = (total nucleotides in ALL Trinity.fasta × 2 strands) / 1,000,000 Mb
    TOTAL_NT=$(grep -v "^>" Trinity.fasta | awk '{total+=length($0)} END{print total+0}') # Remove headers and count number of nucleotides
	TOTAL_MB=$(echo "${TOTAL_NT} * 2 / 1000000" | bc -l | awk '{printf "%.6f", $1}')
	echo "[INFO] --Z ${TOTAL_MB}" # --Z 3182.349886
	
	TOTAL_MB=$(awk '/^>/ {if(seq) {total+=length(seq)}; seq=""} !/^>/ {seq=seq$0} END{if(seq) total+=length(seq); print (total*2)/1000000}' Trinity.fasta | awk '{printf "%.6f", $1}') # --Z 3182.349886
    echo "[INFO] --Z ${TOTAL_MB} (total Mb both strands)"

    echo "[INFO] Downloading + scanning Rfam (Infernal cmscan)..."
    rm -f rfam.tblout *.cmscan.log .rfam.done
    # Download/index if needed
    [[ -s "${RFAM_CM}" ]] || {
        wget -N "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz" -O "${TRINOTATE_DATA_DIR}/Rfam.cm.gz"
        gunzip "${TRINOTATE_DATA_DIR}/Rfam.cm.gz"
    }
    [[ -s "${TRINOTATE_DATA_DIR}/Rfam.clanin" ]] || {
        wget -N "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin" -O "${TRINOTATE_DATA_DIR}/Rfam.clanin"
    }
	# Index CM
    [[ -s "${RFAM_CM}.i1f" ]] || cmpress "${RFAM_CM}"
    # Scan (add -Z <total_Mb> for accurate E-values if known; est. total size)
    cmscan --cut_ga --rfam --nohmmonly --clanin "${TRINOTATE_DATA_DIR}/Rfam.clanin" \
           --oskip --tblout rfam.tblout --fmt 2 --notextw --cpu ${THREADS} -Z ${TOTAL_MB} \
           "${RFAM_CM}" Trinity.fasta \
           > "${LOGDIR}/rfam.cmscan.log" 2>&1
    #tail -1 rfam.cmscan.log  # Verify Z used
	touch .rfam.done
	echo "[SUCCESS] Rfam complete $(date)"
fi