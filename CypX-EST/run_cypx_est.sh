#!/usr/bin/env bash
set -euo pipefail

# ========= 用户可改的变量 =========
GMX="${GMX:-gmx}"           # GROMACS 可执行名；MPI 版可设为 gmx_mpi
SYS="CypX-EST"
SYSTEM_GRO="system.gro"      # 初始整体结构（已溶剂化加离子，EST 已结合）
TOP="topol.top"              # 主拓扑（含 EST 的 *.itp 引用）
LIG_PREFIX="est"             # 配体文件前缀：est.itp / est.prm（若有）
T=310.15                     # 温度(K)（日志为 310.15 K）
P=1.0                        # 压强(bar)

# 温压耦合分组（需在 index 或默认组中存在；二选一，按你体系实际情况修改）
TC_GRPS="Protein Non-Protein"
# TC_GRPS="Protein_LIG Water_and_ions"

# GPU 自动探测（有 GPU 则附加 -nb gpu -pme gpu）
MDGPU_OPTS=()
if command -v nvidia-smi >/dev/null 2>&1; then
  MDGPU_OPTS=(-nb gpu -pme gpu)
fi

mkdir -p mdp logs

# ========== 1) 生成 .mdp 文件（参数与日志匹配） ==========
# 1.1 能量最小化（steep；nsteps=5000；位能收敛到 emtol=1000）
cat > mdp/em.mdp <<'EOF'
integrator      = steep
emtol           = 1000
emstep          = 0.01
nsteps          = 5000

cutoff-scheme   = Verlet
rlist           = 1.2
coulombtype     = PME
coulomb-modifier= Potential-shift
rcoulomb        = 1.2

vdwtype         = Cut-off
vdw-modifier    = Force-switch
rvdw-switch     = 1.0
rvdw            = 1.2

constraints     = none
pbc             = xyz
nstenergy       = 1000
nstlog          = 1000
EOF

# 1.2 NPT 平衡（125 ps；dt=0.001；V-rescale + C-rescale；-DPOSRES）
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
coulomb-modifier= Potential-shift
rcoulomb        = 1.2
vdwtype         = Cut-off
vdw-modifier    = Force-switch
rvdw-switch     = 1.0
rvdw            = 1.2

tcoupl          = V-rescale
tc-grps         = ${TC_GRPS}
tau_t           = 1.0   1.0
ref_t           = ${T}  ${T}

pcoupl          = C-rescale
pcoupltype      = isotropic
tau_p           = 5.0
ref_p           = ${P}
compressibility = 4.5e-5

define          = -DPOSRES
gen_vel         = yes
gen_temp        = ${T}
gen_seed        = -1

pbc             = xyz
nstxout-compressed = 5000
nstenergy       = 1000
nstlog          = 1000
EOF

# 1.3 生产相（1 µs：dt=0.002 → nsteps=500000000；输出频率与日志一致）
cat > mdp/md.mdp <<EOF
integrator      = md
dt              = 0.002
nsteps          = 500000000     ; 1 µs

constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4

cutoff-scheme   = Verlet
rlist           = 1.2
coulombtype     = PME
coulomb-modifier= Potential-shift
rcoulomb        = 1.2

vdwtype         = Cut-off
vdw-modifier    = Force-switch
rvdw-switch     = 1.0
rvdw            = 1.2

tcoupl          = V-rescale
tc-grps         = ${TC_GRPS}
tau_t           = 1.0   1.0
ref_t           = ${T}  ${T}

pcoupl          = C-rescale
pcoupltype      = isotropic
tau_p           = 5.0
ref_p           = ${P}
compressibility = 4.5e-5

pbc             = xyz
nstxout         = 50000
nstvout         = 50000
nstfout         = 50000
nstenergy       = 50000
nstlog          = 50000
nstxout-compressed = 50000
compressed-x-precision = 1000

gen_vel         = no
continuation    = yes
EOF

# ========== 2) 拓扑检查：确保包含配体与 POSRES ==========
warn=false
if ! grep -qiE "#include\\s+\"${LIG_PREFIX}\\.itp\"" "${TOP}"; then
  echo "⚠️  警告：${TOP} 未包含 #include \"${LIG_PREFIX}.itp\"；请在适当位置加入以含 EST 力场参数。"
  warn=true
fi
if ! grep -qiE "POSRES" "${TOP}"; then
  echo "⚠️  警告：${TOP} 中未检测到 POSRES 宏；若需平衡阶段位置约束，请在拓扑中加入：
#ifdef POSRES
#include \"posre.itp\"
#endif"
  warn=true
fi
if [ "$warn" = true ]; then
  echo "—— 请根据警告修正拓扑后再继续 ——"
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
