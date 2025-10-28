#!/usr/bin/env bash
# CypX-EST one-click MD (GROMACS) — EM -> NVT -> 1-μs MD

set -euo pipefail

# ===== 用户需准备的最小输入 =====
OUTDIR="est"                          # 所有输出写入 ./est/
PDB="inputs_est/start.pdb"            # CypX + EST + 水/离子 坐标，需与拓扑一致
TOP="inputs_est/topol.top"            # 需 #include heme 与 EST 的 *.itp
NDX="inputs_est/index.ndx"            # 可选；存在则使用

# 可通过环境变量切换 gmx_mpi：GMX=gmx_mpi bash cypx_est_run.sh
GMX="${GMX:-gmx}"

mkdir -p "$OUTDIR" mdp
[[ -f "$PDB" ]] || { echo "ERROR: $PDB not found"; exit 1; }
[[ -f "$TOP" ]] || { echo "ERROR: $TOP not found"; exit 1; }

# 仅在存在 index.ndx 时添加 -n
NFLAG=()
[[ -f "$NDX" ]] && NFLAG=(-n "$NDX")

# ===== 写出与日志一致的 MDP =====

# 1) 能量最小化（steep；5000 步；PME + Force-switch 1.0→1.2 nm）
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

# 2) 平衡（NVT；125 ps；310.15 K）
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

# 3) 生产（1 μs；dt=2 fs；V-rescale 310.15K + C-rescale 1 bar）
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

# ===== 执行三步（与日志 deffnm/命令一致）=====
MIN="${OUTDIR}/step4.0_minimization"
EQ="${OUTDIR}/step4.1_equilibration"
MD="${OUTDIR}/step5_1"

echo "== [EST] 1) Energy minimization =="
$GMX grompp -f mdp/minimization.mdp -c "$PDB" -p "$TOP" -o "${MIN}.tpr" -r "$PDB" "${NFLAG[@]}"
$GMX mdrun  -v -deffnm "${MIN}"          # 最小化日志未使用 -gpu_id 0

echo "== [EST] 2) Equilibration (NVT) =="
$GMX grompp -f mdp/equilibration.mdp -c "${MIN}.gro" -p "$TOP" -o "${EQ}.tpr" -r "${MIN}.gro" "${NFLAG[@]}"
$GMX mdrun  -v -deffnm "${EQ}" -gpu_id 0 # 与平衡日志一致

echo "== [EST] 3) Production MD (1 μs) =="
$GMX grompp -f mdp/production.mdp -c "${EQ}.gro" -p "$TOP" -o "${MD}.tpr" -r "${EQ}.gro" "${NFLAG[@]}"
$GMX mdrun  -v -deffnm "${MD}" -gpu_id 0 # 与生产日志一致

echo "== DONE. Outputs in ./${OUTDIR}/ =="
