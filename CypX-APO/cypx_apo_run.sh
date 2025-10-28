#!/usr/bin/env bash
# CypX-APO one-click MD (GROMACS) — minimization → NVT → 1-μs production


set -euo pipefail

# ---- configurable paths (minimal inputs) ----
BASE_OUT="apo"                       # all outputs will live here
PDB="inputs_min/start.pdb"
TOP="inputs_min/topol.top"
NDX="inputs_min/index.ndx"           # optional; used only if exists

# Allow using gmx_mpi by: GMX=gmx_mpi bash cypx_apo_run.sh
GMX="${GMX:-gmx}"

# ---- create folders & sanity checks ----
mkdir -p "$BASE_OUT" mdp
[[ -f "$PDB" ]] || { echo "ERROR: $PDB not found"; exit 1; }
[[ -f "$TOP" ]] || { echo "ERROR: $TOP not found"; exit 1; }

# Use -n only if index exists
NFLAG=()
[[ -f "$NDX" ]] && NFLAG=(-n "$NDX")

# ---- write MDPs exactly as in logs ----
# 1) minimization (steep)
cat > mdp/minimization.mdp <<'EOF'
integrator              = steep
nsteps                  = 5000
emtol                   = 1000
emstep                  = 0.01

cutoff-scheme           = Verlet
nstlist                 = 10
rlist                   = 1.2
coulombtype             = PME
coulomb-modifier        = Potential-shift
rcoulomb                = 1.2

vdw-type                = Cut-off
vdw-modifier            = Force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2

constraints             = none
constraint-algorithm    = lincs
lincs-order             = 4
pbc                     = xyz
nstxout-compressed      = 0
nstenergy               = 1000
nstlog                  = 1000
EOF

# 2) equilibration (NVT, 125 ps)
cat > mdp/equilibration.mdp <<'EOF'
integrator              = md
dt                      = 0.001
nsteps                  = 125000    ; 125 ps

tcoupl                  = V-rescale
tc-grps                 = System
tau-t                   = 1.0
ref-t                   = 310.15

pcoupl                  = no

cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
coulombtype             = PME
coulomb-modifier        = Potential-shift
rcoulomb                = 1.2
fourierspacing          = 0.12
pme-order               = 4

vdw-type                = Cut-off
vdw-modifier            = Force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2

constraints             = h-bonds
constraint-algorithm    = lincs
lincs-order             = 4

nstxout-compressed      = 5000
nstenergy               = 1000
nstlog                  = 1000
pbc                     = xyz
EOF

# 3) production (1 μs, 2 fs)
cat > mdp/production.mdp <<'EOF'
integrator              = md
dt                      = 0.002
nsteps                  = 500000000   ; 1 microsecond

tcoupl                  = V-rescale
tc-grps                 = System
tau-t                   = 1.0
ref-t                   = 310.15

pcoupl                  = C-rescale
pcoupltype              = isotropic
tau-p                   = 5.0
ref-p                   = 1.0
compressibility         = 4.5e-5

cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
coulombtype             = PME
coulomb-modifier        = Potential-shift
rcoulomb                = 1.2
fourierspacing          = 0.12
pme-order               = 4
ewald-rtol-lj           = 0.001

vdw-type                = Cut-off
vdw-modifier            = Force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2

constraints             = h-bonds
constraint-algorithm    = lincs
lincs-order             = 4

nstxout-compressed      = 50000
nstenergy               = 50000
nstlog                  = 50000

pbc                     = xyz
EOF

# ---- run steps ----
MIN="${BASE_OUT}/step4.0_minimization"
EQ="${BASE_OUT}/step4.1_equilibration"
MD="${BASE_OUT}/step5_1"

echo "== [APO] 1) Energy minimization =="
$GMX grompp -f mdp/minimization.mdp -c "$PDB" -p "$TOP" -o "${MIN}.tpr" -r "$PDB" "${NFLAG[@]}"
$GMX mdrun  -v -deffnm "${MIN}"

echo "== [APO] 2) Equilibration (NVT) =="
$GMX grompp -f mdp/equilibration.mdp -c "${MIN}.gro" -p "$TOP" -o "${EQ}.tpr" -r "${MIN}.gro" "${NFLAG[@]}"
$GMX mdrun  -v -deffnm "${EQ}"

echo "== [APO] 3) Production MD (1 μs) =="
$GMX grompp -f mdp/production.mdp -c "${EQ}.gro" -p "$TOP" -o "${MD}.tpr" -r "${EQ}.gro" "${NFLAG[@]}"
# Remove -gpu_id 0 if no GPU is available
$GMX mdrun  -v -deffnm "${MD}" -gpu_id 0

echo "== DONE. Outputs in ./${BASE_OUT}/ =="
