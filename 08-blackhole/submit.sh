#!/bin/bash
set -euo pipefail

CASE="${1:-}"
MODE="${2:-dev}"

if [[ ! "$CASE" =~ ^(small|medium|large|huge)$ ]]; then
  echo "Usage: ./submit.sh <small|medium|large|huge> [dev|official]"
  echo "  dev      : 8-rank kp_interact run (default, for iteration)"
  echo "  official : official rank count on kp_run"
  exit 2
fi

if [[ ! "$MODE" =~ ^(dev|official)$ ]]; then
  echo "MODE must be dev or official"
  exit 2
fi

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p runs .slurm

case "$CASE" in
  small)
    OFFICIAL_RANKS=32; POINTS=10; FULL=0.15; ZERO=200 ;;
  medium)
    OFFICIAL_RANKS=128; POINTS=20; FULL=10; ZERO=500 ;;
  large)
    OFFICIAL_RANKS=32; POINTS=30; FULL=20; ZERO=600 ;;
  huge)
    OFFICIAL_RANKS=128; POINTS=40; FULL=30; ZERO=600 ;;
esac

if [[ "$MODE" == "dev" ]]; then
  PARTITION="kp_interact"
  QOS="kp_interact"
  NTASKS=1
  CPUS_PER_TASK=8
  RUN_RANKS=8
  WALLTIME="01:00:00"
else
  PARTITION="kp_run"
  QOS="kp_run"
  NTASKS="$OFFICIAL_RANKS"
  CPUS_PER_TASK=1
  RUN_RANKS="$OFFICIAL_RANKS"
  WALLTIME="00:20:00"
fi

JOBFILE=".slurm/blackhole-${CASE}-${MODE}.slurm"

cat > "$JOBFILE" <<EOF
#!/bin/bash
#SBATCH -J bh-${CASE}-${MODE}
#SBATCH -p ${PARTITION}
#SBATCH -q ${QOS}
#SBATCH -N 1
#SBATCH -n ${NTASKS}
#SBATCH -c ${CPUS_PER_TASK}
#SBATCH -t ${WALLTIME}
#SBATCH -o ${ROOT}/runs/${CASE}-${MODE}-%j.out
#SBATCH -e ${ROOT}/runs/${CASE}-${MODE}-%j.err

set -euo pipefail

cd "${ROOT}"
source ./env.sh

echo "===== BLACKHOLE JOB ====="
echo "case=${CASE}"
echo "mode=${MODE}"
echo "host=\$(hostname)"
echo "job_id=\$SLURM_JOB_ID"
echo "nodes=\$SLURM_JOB_NODELIST"
echo "run_ranks=${RUN_RANKS}"
echo "started=\$(date -Is)"

echo "===== SOURCE ====="
git rev-parse --short HEAD || true
git status --short -- src/surface_integral.C src/Makefile.inc || true
sha256sum src/surface_integral.C src/Makefile.inc

echo "===== BUILD ====="
if [[ "\${CLEAN_BUILD:-1}" == "1" ]]; then
  make -C src clean
fi
make -C src -j8 blackhole

echo "===== RUN ====="
rm -f benchmark-result.txt
RUNLOG="${ROOT}/runs/${CASE}-${MODE}-\$SLURM_JOB_ID.run.log"

if [[ "${MODE}" == "dev" ]]; then
  mpirun \
    --host "\$(hostname):8" \
    --oversubscribe \
    --bind-to core:overload-allowed \
    --map-by core \
    -np 8 \
    src/blackhole "configs/${CASE}.par" 2>&1 | tee "\$RUNLOG"
else
  mpirun \
    -np "${RUN_RANKS}" \
    --bind-to core \
    --map-by core \
    src/blackhole "configs/${CASE}.par" 2>&1 | tee "\$RUNLOG"
fi

echo "===== RESULT ====="
cat benchmark-result.txt
echo "===== BLACKHOLE ====="
grep '^BLACKHOLE_' "\$RUNLOG"

echo "===== CORRECTNESS ====="
python3 - "${CASE}" "${MODE}" "${POINTS}" "${FULL}" "${ZERO}" <<'PY'
import math, pathlib, re, sys

case, mode = sys.argv[1], sys.argv[2]
P, F, Z = map(float, sys.argv[3:6])

actual = list(map(float, pathlib.Path("benchmark-result.txt").read_text().split()))
expected = list(map(float, pathlib.Path(f"reference/{case}.expected").read_text().split()))

names = ["checksum", "wave_r", "wave_i", "mass", "energy"]
ok = len(actual) == len(expected)
max_abs = 0.0

for i, (x, y) in enumerate(zip(actual, expected)):
    diff = abs(x-y)
    max_abs = max(max_abs, diff)
    tol = 1e-8 + 1e-14*abs(y)
    passed = diff <= tol
    ok &= passed
    name = names[i] if i < len(names) else f"value_{i}"
    print(f"{name:10s} diff={diff:.3e} tol={tol:.3e} {'PASS' if passed else 'FAIL'}")

print(f"OVERALL: {'PASS' if ok else 'FAIL'}")
print(f"MAX_ABS_ERROR: {max_abs:.17g}")
PY

echo "===== PERFORMANCE ====="
python3 - "\$RUNLOG" "${MODE}" "${POINTS}" "${FULL}" "${ZERO}" <<'PY'
import pathlib, re, sys

runlog, mode = sys.argv[1], sys.argv[2]
P, F, Z = map(float, sys.argv[3:6])
text = pathlib.Path(runlog).read_text(errors="replace")
m = re.search(r"^BLACKHOLE_TIME max=([0-9eE+\-.]+)", text, re.M)
if not m:
    print("PERFORMANCE: BLACKHOLE_TIME not found")
    raise SystemExit(0)

t = float(m.group(1))
print(f"PERFORMANCE_SECONDS: {t:.17g}")
if mode == "official":
    raw = F * (Z - t) / (t * (Z - F))
    score = P * min(1.0, max(0.0, raw))
    print(f"CASE_SCORE: {score:.6f}/{P:g}")
else:
    print("CASE_SCORE: dev 8-rank run; not an official score")
PY
echo "finished=\$(date -Is)"
EOF

# Preserve BLACKHOLE lines in the normal Slurm stdout; after completion the
# helper below can extract time and (for official mode) compute the case score.

jid="$(sbatch --parsable "$JOBFILE")"
echo "Submitted job: $jid"
echo "Case: $CASE"
echo "Mode: $MODE"
echo "Watch:"
echo "  squeue -j $jid"
echo "After completion:"
echo "  grep -E '^BLACKHOLE_|OVERALL|MAX_ABS_ERROR|PERFORMANCE_SECONDS|CASE_SCORE' runs/${CASE}-${MODE}-${jid}.out"
echo "  cat runs/${CASE}-${MODE}-${jid}.err"

if [[ "$MODE" == "official" ]]; then
  echo "Official case: ranks=$OFFICIAL_RANKS, max_points=$POINTS, full_time=$FULL s, zero_time=$ZERO s"
fi
