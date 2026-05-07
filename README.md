# RePISoSi — ReBCO Partially Insulated Solenoid Simulation

Transient electro-thermal network model of a partially-insulated, **Re**BCO solenoid coupled with a **Bi**ot–**Sa**vart magnetic-field solver. The code simulates quench propagation in a high-temperature superconducting (HTS) coil with axial shorts, helium gas convection, and copper current leads.

The model represents the conductor as a 3-D helix of inductively-coupled line elements. Each element carries a time-dependent current solved from a stiff ODE (`ode15s` / `ode23s`), a local temperature obtained from a finite-volume thermal substep loop, and material properties (resistivity, thermal conductivity, heat capacity) read from physically-motivated lookup tables (Bloch–Grüneisen, NIST OFHC copper, Al-5083, Fujikura tape parametrisation).

---

## Features

- 3-D Biot–Savart self- and mutual-inductance computation for an arbitrary helix mesh.
- Finite-volume axial heat flow with optional helium gas convection and copper current-lead conduction.
- Field- and angle-dependent critical current via the Fujikura FESC parametrisation (`parametrisation_fujikura.m`).
- Bloch–Grüneisen-based resistivity & thermal conductivity lookup tables, plus NIST fits for OFHC copper and Al-5083.
- Live 3-D field-map visualisation, optional MP4 video recording, and timestamped CSV-style data logs.

## Repository layout

The current flat layout is preserved for backward compatibility with the published thesis. The recommended target structure for future work is:

```
RePISoSi/
├── src/                # Core solver helpers (calc_*, generate_*, gather_*, myODE*, create_MLPInv_matrix)
├── physics/            # Material property models (BlochGrun, rhoCu_nist, CpCu_nist,
│                       # thermal_conductivity_*, calc_normal_*_resistivity,
│                       # heat_capacity_al_alloy_5083, parametrisation_fujikura,
│                       # conductivity_helium, n_value)
├── viz/                # Plotting & video helpers (drawFieldMap, quiverC3D[New],
│                       # makeplots, plotBET, plotIcTc, VIDEOplot, analyseaxialcurrent)
├── examples/           # Entry-point scripts (main_5turnMM.m and legacy variants)
├── data/               # Input datasets (HeConductivity_0.01MPa.csv, Alres.csv)
├── legacy/             # Older entry-point variants (main.m, main2.m, main3_lowerresolution.m)
├── tests/              # (Not yet populated) verification cases
├── README.md
├── LICENSE
└── CITATION.cff
```

## Requirements

- **MATLAB R2021b or newer** (uses `vecnorm`, string concatenation with `append`, name=value plot syntax such as `Color='b'`).
- Toolboxes:
  - **MATLAB** (base) — required.
  - **Symbolic Math Toolbox** — *not* required.
  - **Curve Fitting Toolbox** — only for opening `*.sfit` session files; not required for the simulation.
  - **Image Processing Toolbox** — not required for the core simulation; some legacy plotting helpers reference video writers (built-in `VideoWriter`).

## Installation

```matlab
% Clone or download the repository, then in MATLAB:
cd('path/to/RePISoSi')
addpath(genpath(pwd))    % puts every helper on the MATLAB path
```

## Usage

The canonical entry point is `main_5turnMM.m`.

```matlab
% 1. Open MATLAB and add the repo to the path (see Installation).
% 2. Edit the CONFIGURATION block at the top of main_5turnMM.m to set:
%       I0                       (driving current, A)
%       initial_temperature      (K)
%       N_shorts, numWindings    (geometry)
%       VHeater, RHeater         (heater pulse profile)
%       writedata, writevideo    (output toggles)
% 3. Run:
run('main_5turnMM.m')
```

Outputs:
- A 3-D field-map figure updated every `drawFigureAtIteration` iterations.
- Optional `<Description>.mp4` video of the field map.
- Optional `<Description><timestamp>.txt` log with `[time, centerB, Pheater, T_topRing, T_bottomRing, Efield, Tmax]` per row.
- A final `<Description><timestamp>.mat` workspace dump for post-processing with `makeplots.m` / `plotBET.m`.

### Minimal smoke test

To verify your install without running the full quench, reduce `maxIteration` to e.g. 20 and `numWindings`/`N_shorts` to small integers in `main_5turnMM.m`.

## Known issues / things to fix before extending

- `calc_HTS_critical_current_refined.m` references undefined symbols (`meas_Jc`, `Ic0`) and is currently dead code — see the `% TODO` banner at the top of the file. Use `parametrisation_fujikura.m` instead.
- `n_value.m` originally used a broken element-wise string comparison; fixed to `contains(HTStapeName, 'Fujikura')`. Behaviour for non-Fujikura tapes is undefined and returns no output — caller must handle that.
- `Conductivity.m` is a research scratchpad with a syntax-error line preserved verbatim; it is not on any execution path.
- Several hardcoded paths/numbers in the entry-point scripts are flagged with `% TODO: parameterize` (notably the hardcoded radius `0.12` inside the E-field metric).
- All four `main*.m` scripts perform `clear; close all; clc;` at start-up. The simulation state is therefore lost when the script ends — wrap calls inside a function if you need a reproducible API.

## Citing

If you use this code in academic work, please cite via [`CITATION.cff`](CITATION.cff). The model and validation data are described in J.L. Van den Eijnden's MSc thesis (TU Delft / CERN, 2023).

## License

[MIT](LICENSE) © 2026 J.L. Van den Eijnden, with original contributions by M. Mentink and A. Vaskuri.
