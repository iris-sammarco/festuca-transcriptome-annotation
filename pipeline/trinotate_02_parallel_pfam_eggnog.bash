#!/bin/bash

# Author: Iris Sammarco
# Date: 03/2026
# Aim: Master (non-qsub) script to split a large TransDecoder peptide FASTA into chunks and submit PBS array jobs for parallel Pfam domain search (hmmscan) or EggNOG annotation (emapper.py). Includes a merge mode to concatenate chunked outputs.
# Calls: trinotate_02_pfam_array.pbs or trinotate_04_eggnog_array.pbs
# Run: bash trinotate_02_parallel_pfam_eggnog.bash pfam 141      # split + submit pfam (step 2a)
#      bash trinotate_02_parallel_pfam_eggnog.bash pfam 141 m    # merge chunks of pfam output (step 2b)
#      bash trinotate_02_parallel_pfam_eggnog.bash eggnog 72     # split + submit eggnog (step 4a)
#      bash trinotate_02_parallel_pfam_eggnog.bash eggnog 72 m   # merge chunks of eggnog output (step 4b)
# Input: Trinity.fasta.transdecoder_dir/longest_orfs.pep (pfam)
#        Trinity.fasta.transdecoder.pep (eggnog)
# Output: pfam.domtblout (merged Pfam domain table)
#         eggnog.emapper.annotations (merged EggNOG annotations)

set -euo pipefail
trap 'echo "[ERROR] Line $LINENO: ${BASH_COMMAND} failed on $(date)" >&2; exit 1' ERR

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 {pfam|eggnog} N_CHUNKS [m]" >&2
    exit 1
fi

STEP=$1
N_CHUNKS=$2
MERGE_ONLY=${3:-}
OUTDIR="/path/to/project/assembly/trinotate_output"
PEP_LONG="${OUTDIR}/Trinity.fasta.transdecoder_dir/longest_orfs.pep"

cd "${OUTDIR}"

case "${STEP}" in
    pfam)
        OUTPUT="pfam.domtblout"
        SPLIT_DIR="${OUTDIR}/Trinity.fasta.transdecoder_dir/longest_orfs.pep.split"
        PBS_SCRIPT="trinotate_02_pfam_array.pbs"
        INPUT_FILE="${PEP_LONG}"
		CHUNK_EXT="pep"
        ;;
    eggnog)
        OUTPUT="eggnog.emapper.annotations"
        SPLIT_DIR="${OUTDIR}/Trinity.fasta.transdecoder.pep.split"
        PBS_SCRIPT="trinotate_04_eggnog_array.pbs"
        INPUT_FILE="${OUTDIR}/Trinity.fasta.transdecoder.pep"
		CHUNK_EXT="pep"
        ;;
    *)
        echo "ERROR: Step must be 'pfam' or 'eggnog'" >&2
        exit 1
        ;;
esac

[[ -s "${INPUT_FILE}" ]] || { echo "[FATAL] Missing input: ${INPUT_FILE}" >&2; exit 1; }

if [[ -z "${MERGE_ONLY}" ]]; then
    echo "[INFO] === SPLIT + SUBMIT MODE: ${STEP^^} (${N_CHUNKS} chunks) ==="
    
    # 1. Use existing chunks OR create new ones
	if [[ ! -d "${SPLIT_DIR}" || ! -f "${SPLIT_DIR}/chunk_${STEP}_001.${CHUNK_EXT}" ]]; then
        echo "[INFO] Creating new chunks in ${SPLIT_DIR}"
        rm -rf "${SPLIT_DIR}"
        mkdir -p "${SPLIT_DIR}"
        seqkit split2 -s 10000 "${INPUT_FILE}" -o "${SPLIT_DIR}" --by-size-prefix "chunk_${STEP}_" # 10,000 sequences per chunk; for ~1.4M longest_orfs.pep this gives ~141 chunks; for ~720k transdecoder.pep ~72 chunks
    fi
	
	# 2. Count & verify chunks 
    CHUNKS=$(ls "${SPLIT_DIR}"/chunk_${STEP}_*.${CHUNK_EXT} 2>/dev/null | wc -l || echo 0)
    echo "[INFO] Found ${CHUNKS} chunks in ${SPLIT_DIR}"
    [[ ${CHUNKS} -ge 1 ]] || { echo "[FATAL] No chunks found!" >&2; exit 1; }

	# 3. Submit PBS array
	mkdir -p logs
	qsub -J 1-${CHUNKS}:1 \
		-N "${STEP}_parallel_${CHUNKS}" \
		-o "logs/${STEP}_%A_%a.out" \
		-e "logs/${STEP}_%A_%a.err" \
		"${PBS_SCRIPT}"
	echo "[INFO] Submitted ${CHUNKS} ${STEP^^} jobs. Monitor: qstat -u \$USER | grep ${STEP}"
	   
elif [[ "${MERGE_ONLY}" == "m" ]]; then
    echo "[INFO] === MERGE MODE: ${STEP^^} ==="
    
    if [[ -s "${OUTPUT}" && -s ".${STEP}.merged.done" ]]; then
        echo "[INFO] ${OUTPUT} already merged, skipping"
        exit 0
    fi
    # Pfam uses longest_orfs.pep (pre-TransDecoder.Predict); EggNOG uses transdecoder.pep (post-TransDecoder.Predict, smaller)
    case "${STEP}" in
        pfam)
            # Merge domtblout (HMMER3: 22-line header from first file only)
            rm -f "${OUTPUT}"
            FIRST=1
            for f in "${SPLIT_DIR}"/chunk_${STEP}_*.domtblout; do
                [[ -s "$f" ]] || continue
                if [[ $FIRST -eq 1 ]]; then
                    cat "$f" >> "${OUTPUT}"  # Header + data
                    FIRST=0
                else
                    tail -n +23 "$f" >> "${OUTPUT}"  # HMMER3 domtblout: 3 header lines + statistical summary block = skip first 22 lines from chunks 2+
                fi
            done
            [[ -s "${OUTPUT}" ]] && touch ".${STEP}.merged.done"
            echo "[INFO] Merged ${OUTPUT} ($(wc -l < "${OUTPUT}") lines)"
            ;;
        eggnog)
            # Merge eggnog (5-line header from first, data from all)
            rm -f "${OUTPUT}"
            FIRST=1
            for f in "${SPLIT_DIR}"/chunk_${STEP}_*.emapper.annotations; do
                [[ -s "$f" ]] || continue
                if [[ $FIRST -eq 1 ]]; then
                    head -5 "$f" > "${OUTPUT}"  # EggNOG emapper output: 4 comment lines + 1 column header = keep only from first chunk
                    tail -n +6 "$f" >> "${OUTPUT}"  # Data only
                    FIRST=0
                else
                    tail -n +6 "$f" >> "${OUTPUT}"  # Data only
                fi
            done
            [[ -s "${OUTPUT}" ]] && touch ".${STEP}.merged.done"
            echo "[INFO] Merged ${OUTPUT} ($(wc -l < "${OUTPUT}") lines)"
            ;;
    esac
fi
