# CypX Gating MD Reproducibility (GROMACS)

This repo contains exact run scripts and MDPs for CypX (CYP134A1) systems.
Current systems: `apo`. You can add more under `systems/<name>/`.

## Minimal inputs
Put minimal inputs in `systems/<name>/inputs_min/`:
- `start.pdb`, `topol.top`, required `*.itp` (heme), and optional `index.ndx`.

## How to run
Linux/WSL:
```bash
chmod +x code/run_one.sh code/run_all.sh
# run a single system:
bash code/run_one.sh apo
# or run all existing systems under systems/:
bash code/run_all.sh
