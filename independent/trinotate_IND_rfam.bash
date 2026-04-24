#!/bin/bash

#PBS -l walltime=200:00:0
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=28:mem=40gb:scratch_local=200gb
#PBS -N Festuca_Rfam
#PBS -j oe

# Author: Iris Sammarco
# Date: 03/2026
# Aim: Scan the Trinity transcriptome assembly for ncRNAs using Infernal cmscan against the Rfam database. Uses --cut_ga (Rfam's gathering threshold) and dynamically computes --Z (total Mb × 2 strands) for accurate E-values.
# Position-independent — can run at any time before step 08. Note: slow (~200 h walltime), consider parallelizing this step.
# Rfam.cm and Rfam.clanin are downloaded automatically if absent.
# Run: qsub trinotate_IND_rfam.bash
# Input: Trinity_output.Trinity.fasta, Rfam.cm, Rfam.clanin
# Output: rfam.tblout

set -euo pipefail
trap 'echo "[ERROR] Line $LINENO: ${BASH_COMMAND} failed on $(date)" >&2; exit 1' ERR

export TMPDIR=$SCRATCHDIR
module add mambaforge
export CONDA_ENVS_PATH=/path/to/conda/envs
mamba activate trinotate_env
export TRINITY_HOME="/path/to/conda/envs/trinotate_env/opt/trinity-2.8.5"
export PATH="/path/to/project/tmhmm-2.0c/bin:$PATH"
export TRINOTATE_DATA_DIR="/path/to/project/assembly/trinotate_data"
export EGGNOG_DATA_DIR="${TRINOTATE_DATA_DIR}/eggnog_data"

ASSEMBLY="/path/to/project/assembly/Trinity_output.Trinity.fasta"
PREFIX="Festuca_rubra"
OUTDIR="/path/to/project/assembly/trinotate_output"
LOGDIR="${OUTDIR}/logs"
THREADS=28
EVAL="1e-5"
MAX_TARGETS=5

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

## Infernal cmscan Rfam for ncRNAs
RFAM_CM="${TRINOTATE_DATA_DIR}/Rfam.cm"
if [[ ! -s "rfam.tblout" || ! -s ".rfam.done" ]]; then
	# Calculate Z:
    echo "[INFO] Calculating Z for accurate E-values..." # For cmscan, Z is the length of the current query sequence [nucleotides] multiplied by 2 (because both strands of the sequence are searched [in each model]) and multiplied again by the number of CMs [covariance models] in the target CM database.
	# --Z = (total nucleotides in ALL Trinity.fasta × 2 strands) / 1,000,000 Mb
    TOTAL_NT=$(grep -v "^>" Trinity.fasta | awk '{total+=length($0)} END{print total+0}') # Remove headers and count number of nucleotides
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
