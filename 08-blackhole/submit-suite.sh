#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p runs .slurm
JOBFILE=".slurm/blackhole-suite-dev.slurm"
cat > "$JOBFILE" <<'EOF'
#!/bin/bash
#SBATCH -J bh-suite-dev
#SBATCH -p kp_interact
#SBATCH -q kp_interact
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 1:00:00
#SBATCH -o runs/suite-dev-%j.out
#SBATCH -e runs/suite-dev-%j.err

set -euo pipefail
ROOT="/vault/home/acct-stu/stu1671/hellohpc-2nd/08-blackhole"
cd "$ROOT"
source ./env.sh
echo "===== BUILD ONCE ====="
make -C src clean
make -C src -j8 blackhole
for CASE in small medium large huge; do
  echo "===== CASE=$CASE ====="
  rm -f benchmark-result.txt
  RUNLOG="$ROOT/runs/${CASE}-suite-$SLURM_JOB_ID.run.log"
  mpirun --host "$(hostname):8" --oversubscribe --bind-to core:overload-allowed --map-by core -np 8 src/blackhole "configs/${CASE}.par" 2>&1 | tee "$RUNLOG"
  python3 - "$CASE" "$RUNLOG" <<'PY'
import pathlib,re,sys
case,runlog=sys.argv[1],sys.argv[2]
meta={"small":(10,0.15,200),"medium":(20,10,500),"large":(30,20,600),"huge":(40,30,600)}
P,F,Z=meta[case]
actual=list(map(float,pathlib.Path("benchmark-result.txt").read_text().split()))
expected=list(map(float,pathlib.Path(f"reference/{case}.expected").read_text().split()))
ok=len(actual)==len(expected); ma=0.0
for x,y in zip(actual,expected):
    d=abs(x-y); ma=max(ma,d); ok &= d <= 1e-8 + 1e-14*abs(y)
text=pathlib.Path(runlog).read_text(errors="replace")
m=re.search(r"^BLACKHOLE_TIME max=([0-9eE+\-.]+)",text,re.M)
t=float(m.group(1)) if m else float("nan")
print(f"SUITE_CASE {case} correctness={'PASS' if ok else 'FAIL'} max_abs={ma:.17g} time8={t:.17g} official_full={F:g}s official_points={P:g}")
PY
done
echo "===== SUITE COMPLETE ====="
EOF
jid="$(sbatch --parsable "$JOBFILE")"
echo "Submitted suite job: $jid"
echo "Watch: squeue -j $jid"
echo "After completion: grep '^SUITE_CASE ' runs/suite-dev-$jid.out"
