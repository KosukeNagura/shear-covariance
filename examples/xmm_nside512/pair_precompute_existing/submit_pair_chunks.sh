#!/usr/bin/env bash
# Submit the established XMM/Nside=512 pair-worker over disjoint active-pixel ranges.
# The chunking, worker, and PBS resource defaults follow the original workflow.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WORK_ROOT="${SHEAR_PAIR_WORK_ROOT:-$SCRIPT_DIR}"
CHUNK_SIZE="${SHEAR_PAIR_CHUNK_SIZE:-2000}"
QUEUE="${SHEAR_PAIR_QUEUE:-tiny}"
MEMORY="${SHEAR_PAIR_MEMORY:-8gb}"
WALLTIME="${SHEAR_PAIR_WALLTIME:-48:00:00}"
CONDA_INIT="${SHEAR_CONDA_INIT:-/home/anaconda3/bin/activate}"
CONDA_ENV="${SHEAR_CONDA_ENV:-cellcov_3}"

for name in SHEAR_PAIR_CATALOG SHEAR_PAIR_THETA_NPZ; do
    if [[ -z "${!name:-}" || ! -r "${!name}" ]]; then
        echo "Set $name to a readable input file before submitting." >&2
        exit 2
    fi
done
if [[ ! -d "$WORK_ROOT" ]]; then
    echo "SHEAR_PAIR_WORK_ROOT must already exist: $WORK_ROOT" >&2
    exit 2
fi
if ! command -v qsub >/dev/null 2>&1; then
    echo "qsub was not found; this submission wrapper targets PBS/Torque." >&2
    exit 2
fi
if ! [[ "$CHUNK_SIZE" =~ ^[1-9][0-9]*$ ]]; then
    echo "SHEAR_PAIR_CHUNK_SIZE must be a positive integer." >&2
    exit 2
fi

export SHEAR_PAIR_WORK_ROOT="$WORK_ROOT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p "$WORK_ROOT/jobs_pair_chunks" "$WORK_ROOT/log"

n_active="$(PYTHONPATH="$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}" "${PYTHON:-python}" - <<'PY'
import numpy as np
from config_custom_10bins import pair_meta_path
with np.load(pair_meta_path()) as d:
    print(len(d['pix_idx']))
PY
)"
if ! [[ "$n_active" =~ ^[1-9][0-9]*$ ]]; then
    echo "Invalid active-pixel count from pair_meta.npz: $n_active" >&2
    exit 2
fi

echo "n_active=$n_active chunk_size=$CHUNK_SIZE queue=$QUEUE"
for ((start=0, chunk=0; start<n_active; start+=CHUNK_SIZE, chunk++)); do
    end=$((start + CHUNK_SIZE))
    jobfile="$WORK_ROOT/jobs_pair_chunks/pair_chunk_${chunk}.pbs"
    cat > "$jobfile" <<EOF_JOB
#!/usr/bin/env bash
#PBS -V
#PBS -N pc_${chunk}
#PBS -o ${WORK_ROOT}/log/pair_chunk_${chunk}.out
#PBS -e ${WORK_ROOT}/log/pair_chunk_${chunk}.err
#PBS -l select=1:ncpus=1:mem=${MEMORY}
#PBS -l walltime=${WALLTIME}
#PBS -q ${QUEUE}
set -euo pipefail
source ${CONDA_INIT}
conda activate ${CONDA_ENV}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd ${SCRIPT_DIR}
python -u pair_chunk_worker.py ${chunk} ${start} ${end} > ${WORK_ROOT}/log/log_pair_chunk_${chunk}.txt 2>&1
EOF_JOB
    chmod +x "$jobfile"
    qsub "$jobfile"
    echo "submitted chunk $chunk: i=[$start, $end)"
done
