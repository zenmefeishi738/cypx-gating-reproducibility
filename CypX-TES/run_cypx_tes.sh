#!/usr/bin/env bash
set -euo pipefail

# ========= User-tunable =========
GMX="${GMX:-gmx}"           # or gmx_mpi
SYS="CypX-TES"
SYSTEM_GRO="system.gro"      # solvated+ions with TES bound
TOP="topol.top"              # topology including TES parameters (tes.itp etc.)
LIG_PREFIX="tes"
T=310.15                     # K
P=1.0                        # bar
USE_POSRES=${USE_POSRES:-false}  # set true to define -DPOSRES during equilibration

# Temperature-coupling groups (adjust to your index groups)
TC_GRPS="${TC_GRPS:-Protein Non-Protein}"
# alt example: TC_GRPS="Protein_LIG Water_and_ions"

# GPU autodetect
MDGPU_OPTS=()
if command -v nvidia-smi >/dev/null 2>&1; then
  MDGPU_OPTS=(-nb gpu -pme gpu)
fi

mkdir -p mdp logs

# ----- 1) Write MDPs -----
# 1.1 Energy minimization (steep; nsteps=5000; rlist=1.26)
cat > mdp/em.mdp <<'EOF'
integrator      = steep
emtol           = 1000
emstep          = 0.01
nsteps          = 5000

cutoff-scheme   = Verlet
rlist           = 1.26
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

# 1.2 Equilibration (NVT, dt=0.001, nsteps=125000)
cat > mdp/equilibration.mdp <<EOF
integrator      = md
dt              = 0.001
nsteps          = 125000

constraints     = h-bonds
constraint_algorithm = lincs
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

pcoupl          = no

pbc             = xyz
EOF

# 1.3 Production (NPT, C-rescale, 1 µs, output = 50000)
cat > mdp/md.mdp <<EOF
integrator      = md
dt              = 0.002
nsteps          = 500000000     ; 1 µs

constraints     = h-bonds
constraint_algorithm = lincs
lincs_iter      = 1
lincs_order     = 4

cutoff-scheme   = Verlet
rlist           = 1.22
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

pme-order       = 4
fourierspacing  = 0.12

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

# ----- 2) Topology sanity checks -----
warn=false
if ! grep -qiE "#include\\s+\"${LIG_PREFIX}\\.itp\"" "${TOP}"; then
  echo "⚠️  ${TOP} missing #include \"${LIG_PREFIX}.itp\" (TES parameters)."
  warn=true
fi

if [ "${USE_POSRES}" = true ]; then
  if ! grep -qiE "POSRES" "${TOP}"; then
    cat <<'MSG'
⚠️  POSRES macro not detected in topology. To enable position restraints during equilibration, add:
#ifdef POSRES
#include "posre.itp"
#endif
MSG
    warn=true
  fi
fi

$warn && echo "Please fix warnings above if needed."

# ----- 3) Run -----
# 3.1 EM
${GMX} grompp -f mdp/em.mdp -c "${SYSTEM_GRO}" -p "${TOP}" -o step4.0_minimization.tpr -r "${SYSTEM_GRO}" -maxwarn 1 2>&1 | tee logs/step4.0_grompp.log
${GMX} mdrun -v -deffnm step4.0_minimization "${MDGPU_OPTS[@]}" 2>&1 | tee logs/step4.0_mdrun.log

# 3.2 Equilibration (NVT; optional POSRES via env USE_POSRES=true)
EQ_DEFINE=()
[ "${USE_POSRES}" = true ] && EQ_DEFINE=(-D "POSRES")
${GMX} grompp -f mdp/equilibration.mdp -c step4.0_minimization.gro -p "${TOP}" -o step4.1_equilibration.tpr -r step4.0_minimization.gro "${EQ_DEFINE[@]}" -maxwarn 1 2>&1 | tee logs/step4.1_grompp.log
${GMX} mdrun -v -deffnm step4.1_equilibration "${MDGPU_OPTS[@]}" 2>&1 | tee logs/step4.1_mdrun.log

# 3.3 Production (NPT C-rescale)
${GMX} grompp -f mdp/md.mdp -c step4.1_equilibration.gro -p "${TOP}" -o step5_1.tpr -r step4.1_equilibration.gro -maxwarn 1 2>&1 | tee logs/step5_1_grompp.log
${GMX} mdrun -v -deffnm step5_1 "${MDGPU_OPTS[@]}" 2>&1 | tee logs/step5_1_mdrun.log

echo "✅ ${SYS} finished: EM → NVT equil → 1 µs production (C-rescale)"
