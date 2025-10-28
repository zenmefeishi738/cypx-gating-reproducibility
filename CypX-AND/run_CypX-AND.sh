#!/usr/bin/env bash
set -euo pipefail

# ========= 用户可改的变量 =========
GMX="${GMX:-gmx}"           # GROMACS 可执行名；有 MPI 版就改成 gmx_mpi
SYS="CypX-AND"
SYSTEM_GRO="system.gro"      # 初始整体结构（已溶剂化加离子）
TOP="topol.top"              # 主拓扑
POSRES="posre.itp"           # 位置约束（蛋白）
LIG_PREFIX="and"             # 配体文件前缀：and.itp / and.prm
T=300                        # 温度(K)
P=1.0                        # 压强(bar)

# GPU 自动探测（有 GPU 则附加 -nb gpu -pme gpu）
MDGPU_OPTS=()
if command -v nvidia-smi >/dev/null 2>&1; then
  MDGPU_OPTS=(-nb gpu -pme gpu)
fi

mkdir -p mdp logs

# ========== 1) 生成 .mdp 文件 ==========
# 1.1 能量最小化（与日志一致：integrator=steep, nsteps=5000, Verlet 截断）
cat > mdp/em.mdp <<'EOF'
integrator      = steep
emtol           = 1000
emstep          = 0.01
nsteps          = 5000

cutoff-scheme   = Verlet
ns_type         = grid
rlist           = 1.2
coulombtype     = PME
rcoulomb        = 1.2

vdwtype         = Cut-off
vdw-modifier    = Force-switch
rvdw-switch     = 1.0
rvdw            = 1.2

constraints     = none
pbc             = xyz
EOF

# 1.2 平衡（125 ps；dt=0.001；NPT：V-rescale恒温，Parrinello-Rahman恒压）
cat > mdp/equilibration.mdp <<EOF
integrator      = md
dt              = 0.001
nsteps          = 125000       ; 125 ps
continuation    = no
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4

cutoff-scheme   = Verlet
rlist           = 1.2
coulombtype     = PME
rcoulomb        = 1.2
vdwtype         = Cut-off
vdw-modifier    = Force-switch
rvdw-switch     = 1.0
rvdw            = 1.2

tcoupl          = V-rescale
tc-grps         = Protein_LIG Water_and_ions
tau_t           = 0.1   0.1
ref_t           = ${T}  ${T}

pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = ${P}
compressibility = 4.5e-5

gen_vel         = yes
gen_temp        = ${T}
gen_seed        = -1

pbc             = xyz
nstxout-compressed = 5000
nstenergy       = 1000
nstlog          = 1000
EOF

# 1.3 生产相（例如 1 μs：dt=0.002 -> nsteps=500000000）
cat > mdp/md.mdp <<'EOF'
integrator      = md
dt              = 0.002
nsteps          = 500000000     ; 1 us

constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4

cutoff-scheme   = Verlet
rlist           = 1.2
coulombtype     = PME
rcoulomb        = 1.2

vdwtype         = Cut-off
vdw-modifier    = Force-switch
rvdw-switch     = 1.0
rvdw            = 1.2

tcoupl          = V-rescale
tc-grps         = Protein_LIG Water_and_ions
tau_t           = 0.1   0.1
ref_t           = 300   300

pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5

pbc             = xyz
nstxout-compressed = 10000
nstenergy       = 5000
nstlog          = 5000
EOF

# ========== 2) 拓扑检查：确保包含配体参数 ==========
if ! grep -qiE "#include\\s+\"${LIG_PREFIX}\\.itp\"" "${TOP}"; then
  echo "警告：${TOP} 未包含 #include \"${LIG_PREFIX}.itp\"，请在 [ moleculetype ] 之前加入相应 include。"
fi

# ========== 3) grompp & mdrun 三步 ==========
# 3.1 EM
${GMX} grompp -f mdp/em.mdp -c "${SYSTEM_GRO}" -p "${TOP}" -o step4.0_minimization.tpr -r "${SYSTEM_GRO}" -maxwarn 1 2>&1 | tee logs/step4.0_grompp.log
${GMX} mdrun -v -deffnm step4.0_minimization "${MDGPU_OPTS[@]}" 2>&1 | tee logs/step4.0_mdrun.log

# 3.2 平衡
${GMX} grompp -f mdp/equilibration.mdp -c step4.0_minimization.gro -p "${TOP}" -o step4.1_equilibration.tpr -r step4.0_minimization.gro -maxwarn 1 2>&1 | tee logs/step4.1_grompp.log
${GMX} mdrun -v -deffnm step4.1_equilibration "${MDGPU_OPTS[@]}" 2>&1 | tee logs/step4.1_mdrun.log

# 3.3 生产
${GMX} grompp -f mdp/md.mdp -c step4.1_equilibration.gro -p "${TOP}" -o step5_1.tpr -r step4.1_equilibration.gro -maxwarn 1 2>&1 | tee logs/step5_1_grompp.log
${GMX} mdrun -v -deffnm step5_1 "${MDGPU_OPTS[@]}" 2>&1 | tee logs/step5_1_mdrun.log

echo "✅ ${SYS} | 全流程完成：step4.0_minimization → step4.1_equilibration → step5_1"
