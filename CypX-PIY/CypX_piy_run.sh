#!/usr/bin/env bash
# CypX-PIY one-click MD (GROMACS) — minimization → NVT → 1-μs production
# T=310.15 K (V-rescale), P=1 bar (C-rescale), nsteps=5e8, dt=0.002, outputs every 50000 steps.

set -euo pipefail

# ===== 用户需准备的最小输入（PIY 体系）=====
OUTDIR="piy"                            # 所有输出写到 ./piy/
PDB="inputs_piy/start.pdb"              # 含 CypX + PIY + 溶剂/离子、与拓扑匹配
TOP="inputs_piy/topol.top"              # topol.top 要 #include heme 与 PIY 的 *.itp
NDX="inputs_piy/index.ndx"              # 可选；存在则使用

# 如需用 gmx_mpi：GMX=gmx_mpi bash cypx_piy_run.sh
GMX="${GMX:-gmx}"

mkdir -p "$OUTDIR" mdp
[[ -f "$PDB" ]] || { echo "ERROR: $PDB not found"; exit 1; }
[[ -f "$TOP" ]] || { echo "ERROR: $TOP not found"; exit 1; }

# 仅当 index.ndx 存在时才加 -n
NFLAG=()
[[ -f "$NDX" ]] && NFLAG=( -n "$NDX" )

# ===== 写出与日志一致的 MDP =====

# 1) Energy minimization（steep; 5k steps）
# 依据：mdrun -deffnm step4.0_minimization；PME + Force-switch 1.0→1.2 nm。 
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

# 2) Equilibration（NVT, 125 ps）
# 依据：mdrun -deffnm step4.1_equilibration；V-rescale 310.15 K；PME + Force-switch。 
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

# 3) Production（1 μs, dt=2 fs, V-rescale + C-rescale）
# 依据：mdrun -deffnm step5_1 -gpu_id 0；每 50000 步写坐标/能量/日志；PME/Force-switch。
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

# ===== 逐步执行（输出名与日志一致）=====
MIN="${OUTDIR}/step4.0_minimization"
EQ="${OUTDIR}/step4.1_equilibration"
MD="${OUTDIR}/step5_1"

echo "== [PIY] 1) Energy minimization =="
$GMX grompp -f mdp/minimization.mdp -c "$PDB" -p "$TOP" -o "${MIN}.tpr" -r "$PDB" "${NFLAG[@]}"
$GMX mdrun  -v -deffnm "${MIN}"              # 与日志命令一致：-deffnm step4.0_minimization

echo "== [PIY] 2) Equilibration (NVT) =="
$GMX grompp -f mdp/equilibration.mdp -c "${MIN}.gro" -p "$TOP" -o "${EQ}.tpr" -r "${MIN}.gro" "${NFLAG[@]}"
$GMX mdrun  -v -deffnm "${EQ}"               # 与日志命令一致：-deffnm step4.1_equilibration

echo "== [PIY] 3) Production MD (1 μs) =="
$GMX grompp -f mdp/production.mdp -c "${EQ}.gro" -p "$TOP" -o "${MD}.tpr" -r "${EQ}.gro" "${NFLAG[@]}"
# 如无 GPU，请把下一行的 -gpu_id 0 删除
$GMX mdrun  -v -deffnm "${MD}" -gpu_id 0     # 与日志命令一致：-deffnm step5_1 -gpu_id 0

echo "== DONE. All outputs in ./${OUTDIR}/ =="
