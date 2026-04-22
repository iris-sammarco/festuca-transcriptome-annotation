#!/bin/bash

#PBS -l walltime=25:00:0
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=12:mem=45gb:scratch_local=200gb
#PBS -N Festuca_rice_blastp
#PBS -j oe

# Author: Iris Sammarco
# Date: 06/03/2026
# Aim: Annotate Festuca rubra Trinity assembly with Trinotate + rice/wheat curation as no reference genome is available for F. rubra. When each step terminates successfully, it creates an output file with the "done" flag, so if the script needs to be resumed it will skip all the steps that have the files complete.
# Separate the longest steps in separate scripts.
# Outputs: annotation_report.xls with Pfam domains, TM/signal peptides, homology
# Run: qsub Trinotate_main.bash

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
THREADS=12
EVAL="1e-5"  # Strict Diamond E-value
MAX_TARGETS=5  # Maximum number of target sequences to report alignments for (default=25)

mkdir -p "${OUTDIR}" "${LOGDIR}"
cd "${OUTDIR}"

# Per script logs:
exec 1> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.log")
exec 2> >(tee -a "${LOGDIR}/${PBS_JOBNAME}.err")
echo "[INFO] Job ${PBS_JOBID} started in ${PWD} | Threads: ${THREADS}"

# Create symlinks of the db files in the OUTDIR:
#ln -sf "${TRINOTATE_DATA_DIR}"/Trinotate.* .
#ln -sf "${TRINOTATE_DATA_DIR}"/Pfam-A.hmm* .
#ln -sf "${TRINOTATE_DATA_DIR}"/uniprot_sprot.pep .
#ln -sf "${TRINOTATE_DATA_DIR}"/{uniprot_sprot.diamond,Rfam.cm,Rfam.clanin} . 2>/dev/null || true

## Sanity checks (fail-fast)
[[ -s "${ASSEMBLY}" ]] || { echo "[FATAL] Missing assembly: ${ASSEMBLY}"; exit 1; }
ln -sf "${ASSEMBLY}" Trinity.fasta

## STEP 11: Rice Curation based on their proteomes fasta files
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