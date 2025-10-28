#!/usr/bin/env bash
set -euo pipefail

SYS="${1:-}"
if [[ -z "$SYS" ]]; then
  echo "Usage: bash code/run_one.sh <system-name>   # e.g., apo"
  exit 1
fi

BASE="systems/${SYS}"
PDB="${BASE}/inputs_min/start.pdb"
TOP="${BASE}/inputs_min/topol.top"
NDX="${BASE}/inputs_min/index.ndx"

MDP_MIN="mdp/minimization.mdp"
MDP_EQ="mdp/equilibration.mdp"   # NVT
MDP_MD="mdp/production.mdp"

# outputs will live inside systems/<sys>/
DEFFNM_MIN="${BASE}/step4.0_minimization"
DEFFNM_EQ="${BASE}/step4.1_equilibration"
DEFFNM_MD="${BASE}/step5_1"

# only add -n if index.ndx exists
NFLAG=()
if [[ -f "$NDX" ]]; then NFLAG=(-n "$NDX"); fi

echo "== [$SYS] 1) Energy minimization =="
gmx grompp -f "$MDP_MIN" -c "$PDB" -p "$TOP" -o "${DEFFNM_MIN}.tpr" -r "$PDB" "${NFLAG[@]}"
gmx mdrun -v -deffnm "${DEFFNM_MIN}"

echo "== [$SYS] 2) Equilibration (NVT) =="
gmx grompp -f "$MDP_EQ" -c "${DEFFNM_MIN}.gro" -p "$TOP" -o "${DEFFNM_EQ}.tpr" -r "${DEFFNM_MIN}.gro" "${NFLAG[@]}"
gmx mdrun -v -deffnm "${DEFFNM_EQ}"

echo "== [$SYS] 3) Production MD =="
gmx grompp -f "$MDP_MD" -c "${DEFFNM_EQ}.gro" -p "$TOP" -o "${DEFFNM_MD}.tpr" -r "${DEFFNM_EQ}.gro" "${NFLAG[@]}"
# remove -gpu_id 0 if no GPU
gmx mdrun -v -deffnm "${DEFFNM_MD}" -gpu_id 0

echo "== [$SYS] Done. Outputs in ${BASE}/ =="
