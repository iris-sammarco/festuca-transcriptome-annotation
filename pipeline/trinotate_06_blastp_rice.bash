#!/bin/bash

#PBS -l walltime=25:00:0
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=12:mem=45gb:scratch_local=200gb
#PBS -N Festuca_rice_blastp
#PBS -j oe

# Author: Iris Sammarco
# Date: 03/2026
# Aim: Download the Oryza sativa proteome from Ensembl Plants (release 62) and run BLAST+ blastp against it for Poaceae-specific annotation curation.
# Note: uses BLAST+ (not Diamond) as Trinotate 3.2.2 requires BLAST+ outfmt6 for LOAD_custom_blast. Can be run in parallel with step 07 (wheat).
# Run: qsub trinotate_06_blastp_rice.bash
# Input: Trinity.fasta.transdecoder.pep (step 03)
# Output: blastp.rice.outfmt6
#         rice.p{hr,in,sq} (BLAST+ database files)

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
THREADS=12
EVAL="1e-5"
MAX_TARGETS=5

mkdir -p "${OUTDIR}" "${LOGDIR}"
cd "${OUTDIR}"

# Per script logs:
exec 1> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.log")
exec 2> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.err")
echo "[INFO] Job ${PBS_JOBID} started in ${PWD} | Threads: ${THREADS}"

## Sanity checks (fail-fast)
[[ -s "${ASSEMBLY}" ]] || { echo "[FATAL] Missing assembly: ${ASSEMBLY}"; exit 1; }
ln -sf "${ASSEMBLY}" Trinity.fasta

## Rice Curation based on its proteome fasta file
PEP_FINAL="Trinity.fasta.transdecoder.pep"
if [[ ! -s "blastp.rice.outfmt6" || ! -s ".rice.done" ]]; then
    echo "[INFO] Rice proteome..."
    rm -f *.rice* .rice.done
    wget -N --no-check-certificate \ "https://ftp.ensemblgenomes.org/pub/plants/release-62/fasta/oryza_sativa/pep/Oryza_sativa.IRGSP-1.0.pep.all.fa.gz"
    gunzip -f Oryza_sativa.IRGSP-1.0.pep.all.fa.gz
    #diamond makedb --in Oryza_sativa.IRGSP-1.0.pep.all.fa -d rice --threads ${THREADS}
    #diamond blastp --db rice.dmnd --query "${PEP_FINAL}" \
    #    --out blastp.rice.outfmt6 --evalue ${EVAL} --max-target-seqs ${MAX_TARGETS} --threads ${THREADS} --outfmt 6 \
	makeblastdb -in Oryza_sativa.IRGSP-1.0.pep.all.fa -dbtype prot -out rice -blastdb_version 4 # run with BLAST+ as Trinotate takes this input
    blastp -db rice -query "${PEP_FINAL}" -out blastp.rice.outfmt6 \
    -evalue ${EVAL} -max_target_seqs ${MAX_TARGETS} -max_hsps 1 \
    -num_threads ${THREADS} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" > "${LOGDIR}/rice.log" 2>&1
    touch .rice.done
fi
