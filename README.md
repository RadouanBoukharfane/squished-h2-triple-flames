# Squished hydrogen/air triple flames from 1 to 20 bar — reproducibility package

This repository contains the inputs, code modifications, scripts and derived data needed to reproduce the
simulations and figures of

> *<Squished hydrogen--air triple flames: a tip-curvature law and its pressure dependence, R. Boukharfane (2026))>*

All simulations were performed with **PeleLMeX** (2D, low-Mach, AMR) with two transport formulations: mixture-averaged (MA) and constant Schmidt numbers (SC).

## Contents

| Directory | Content |
| --- | --- |
| `pelelmex/` | Code versions (`VERSIONS.md`), the PelePhysics patch adding the constant-Schmidt transport, the San Diego H2 mechanism for PelePhysics |
| `pele/C1/` | Template case (Chaplet et al. case C1): `GNUmakefile`, problem files `pelelmex_prob.{H,cpp}`, input file, job script. Executable with MA transport |
| `pele/TripleFlame-SC-build/` | Same problem compiled against the patched PelePhysics (SC transport) |
| `pele/FlameSheetH2*/` | 1D flames at φ = 1.6 (MA and SC) used for s_L, δ_L and s_T,th |
| `pele/tools/` | Case generator and post-processing (see below) |
| `pele/build_sc_campaign.py`, `pele/build_sc_1d.py` | Generation of the constant-Schmidt cases and 1D flames |
| `cases/` | Input file and job script of **every run** used in the paper (see table) |
| `chem/` | Cantera mechanisms (San Diego, Burke et al. 2012) and constant-Schmidt transport data per condition (`transport/Sc*.json`) |
| `theory/` | Cantera evaluation of Markstein numbers, Lewis numbers, flame thickness; tabulated s_L(φ); Cantera flame-speed sweeps (`cantera/`) and s_T,th of the mixture-averaged series (`sTth_cantera_ma.csv`) |
| `paper/tools/` | Figures 1–3 of the paper, tip-radius measurement, ε* measured on the cold runs (`eps_measured.py`), leading-edge analysis (`tip_analysis.py`), illustrative figures |
| `paper/analysis/` | Robustness analyses quoted in the paper (intercept, measured ε*, tip Markstein number) with their outputs |
| `paper/data/` | Derived data behind the figures and Table 1 |


## Requirements

* PeleLMeX, PelePhysics and AMReX at the commits of `pelelmex/VERSIONS.md`; a C++17 compiler and MPI
  (runs were made with GCC + OpenMPI, 112 MPI ranks per case, 224 for the 960 × 576 base grids).
* Python ≥ 3.10 with Cantera 3.2, NumPy, SciPy, Matplotlib, yt 4.4, contourpy (versions used: Python 3.12.13,
  Cantera 3.2.0, NumPy 2.5.4, SciPy 1.18.1, Matplotlib 3.11.2, yt 4.4.3, contourpy 1.4.0).
* `fextract` from AMReX (`Tools/Plotfile`) for the 1D flame and cold-layer post-processing.

## 1. Build PeleLMeX

```bash
git clone --recursive https://github.com/AMReX-Combustion/PeleLMeX.git && cd PeleLMeX
git checkout 111026b1eafbffd69d77f0e692c3ea27751e0f79 && git submodule update --init --recursive
export PELE_HOME=$PWD
# PelePhysics: San Diego mechanism + constant-Schmidt transport
cp -r <this repo>/pelelmex/Mechanisms/SanDiegoH2 Submodules/PelePhysics/Mechanisms/
```

Two executables are needed:

```bash
cd <this repo>/pele/C1 && make -j TPL && make -j                 # mixture-averaged transport
cd $PELE_HOME/Submodules/PelePhysics && git apply <this repo>/pelelmex/pelephysics_constant_schmidt.patch
cd <this repo>/pele/TripleFlame-SC-build && make -j TPL && make -j # constant Schmidt numbers
```

The patch only acts when `transport.viscosity_power` and `transport.schmidt_numbers` are given in the input
file: μ = μ_ref (T/T_ref)^n, λ = μ c_p / Pr, and ρD_k = μ / Sc_k · W_k / W̄ (the Schmidt numbers of
Chaplet et al., Table 1). Without these keys the Sutherland/simple transport is unchanged.

## 2. Run a case

`cases/<ID>/input.case` and `run.sbatch` are the files of each run. Link the executable into the case
directory, generate the equilibrium table the input refers to, and submit:

```bash
cd cases/C4
ln -s ../../pele/C1/PeleLMeX2d.gnu.MPI.SanDiegoH2.ex .            # *-SC cases: ../../pele/TripleFlame-SC-build/...
python ../../pele/tools/make_eq_table.py eq_z_p1_T300.dat 1 300 1001   # name, p (bar), Tu (K): see pmf.datafile
sbatch run.sbatch                                                  # fill in account/partition first
```

New cases are built with `pele/tools/make_pele_case.py` (lengths and times scaled with δ_L(p, T_u)/δ_L(1 atm,
300 K); see its docstring), e.g. `python pele/tools/make_pele_case.py C4 3.2 18.0 7.6`. The mechanism path is
read from `$PELE_HOME`. The 1D flames need an initial profile from `pele/tools/make_pmf.py` (see its usage line).

The tip speed s_T is the inlet velocity of the active control at the end of the run (`AC_History.dat`); the
companion cold runs (`COLD*`, no ignition, same inlet and velocity) give the mixing-layer thickness Δ at the
tip height (`pele/tools/delta_profile.py`).

### Runs

| Group | Cases | Purpose |
| --- | --- | --- |
| 1 atm, 300 K (Chaplet et al. series) | `C1`, `C2`, `C4`, `C5`, `C6` | MA reference series |
| 10 bar, 579 K and 20 bar, 706 K (isentropic compression) | `B1-p10-*`, `B1-p20-*` | pressure series, MA |
| constant Schmidt, all conditions | `*-SC` (`C*`, `A6*`, `B1-*`, `B2-*` = 10 bar/300 K, `B3-*` = 1 bar/579 K) | second transport formulation |
| 5 bar, 475 K; 10 bar, 300 K; 1 bar, 579 K (mixture-averaged) | `B5-MA-*`, `B2-MA-*`, `B3-MA-*` | isentropic path at 5 bar, pressure only, temperature only |
| mixing-layer shape and inlet distance | `A6P-erf`, `A6P-tanh`, `A6P-tanh10`, `S6-*`, `S45-*`, `D6-tanh-far`, `D1-erf-near` | tanh vs erf layers, tip close to/far from the inlet |
| mesh convergence | `MESH-C4-L3`, `MESH-C4-L5`, `MESH-p10D5b-L3`, `MESH-p10D5b-L5` | finest spacing δ_L/8 and δ_L/32 (production: δ_L/16) |
| Soret effect | `C1-SORET`, `C4-SORET` | thermal diffusion on/off |
| domain sensitivity | `C1-WIDE`, `C1-TALL` | |
| cold companions | `COLD-*`, `COLDF-*` | Δ⁺ at the tip height |
| refined cold companions | `COLDR-*` | layer refined to the finest level over the whole height, for the measured ε* |
| 1D flames | `pele/FlameSheetH2*`, `FlameSheetH2-SC-*` | s_L, δ_L, s_T,th with SC transport |

## 3. Theory (Cantera)

```bash
python theory/theory_groups.py chem/sandiego_h2.yaml theory/groups.json   # Le, α_u, sL, δL, σ, T_ad
python theory/markstein_defs.py                                           # Markstein numbers (consumption, unburnt-gas), San Diego and Burke 2012
python theory/markstein_by_transport.py                                   # Markstein numbers with MA and SC transport (uses sL_sc_pelelmex.csv)
```

`theory/sL_cantera_phi.csv`: Cantera mixture-averaged s_L(φ), δ_L at the conditions of the paper.
`theory/sTth_cantera_ma.csv`: s_T,th = (ρ_u/ρ_b)^{1/2}_{φ=1} s_L(φ = 1.6) of the mixture-averaged series, both factors from
the Cantera sweeps of `theory/cantera/` (San Diego, mixture-averaged). For the constant-Schmidt series s_T,th uses the
PeleLMeX-SC 1D flames (`sL_sc_pelelmex.csv`).
`theory/markstein_sweep.py`: consumption Markstein number along the isentropic compression path (continuous
prediction of figure 3c).
`theory/sL_sc_pelelmex.csv`: PeleLMeX 1D flames with the SC transport (`pele/tools/post_flamesheet.py`).

## 4. Figures

```bash
RUNS=<directory holding the run directories> python paper/tools/make_paper_figures.py   # fig1-3 -> paper/figures/, fits -> paper/data/fits.json
RUNS=... python paper/tools/tip_radius.py out.json C4 B1-p10-D5b ...                      # tip radius of an isotherm (T_TIP, default 1000 K)
```

Figure 1(a) needs the last plotfile of case C4; figures 2 and 3 only use `paper/data/`.

## 5. Robustness analyses (`paper/analysis/`)

```bash
cd paper/analysis
python fit_fixed_A.py                 # intercept A free, = 1, = widest layer, common per transport -> fits_fixed_A.json
python fit_eps_measured.py            # slopes with eps* measured on the refined cold runs -> fits_eps_measured.json
RUNS=... python tip_post.py           # tip Markstein number and eps* seen by the tip -> tip_post.json
```

Their inputs are produced by `paper/tools/eps_measured.py` (ε* measured on a cold run: z(x) at the tip height with
`fextract`, s_L(z(x)), osculating parabola at the maximum; outputs `eps_refined_COLDR.json`, `eps_refined_COLDF.json`)
and `paper/tools/tip_analysis.py` (consumption speed, stretch and curvature at the leading edge of the reacting runs;
output `tip_analysis.json`); job scripts `paper/tools/run_*.sbatch`. `tip_radius_meshC4_*.json`: tip radius of the
mesh-convergence runs.

## Data

* `paper/data/results_current.json`: per condition and transport (`ma_*`, `sc_*`), inlet Δ₀⁺, Δ⁺ at the tip,
  ε*, s_T and s_T/s_T,th; s_T,th per condition; predictions of the slope ratios.
* `paper/data/chaplet_eps.json`: values of Chaplet, Selle & Poinsot (2026) used for comparison (from their
  tables/figures; cite the original paper).
* `paper/data/tip_radius*.json`: tip radii on the 1000/1500/2000 K isotherms.
* `paper/data/fits.json`: lines fitted at each condition and transport (written by `make_paper_figures.py`).

Comments in the input files that mention dates refer to the authors' run log.

## License

BSD 3-Clause (see `LICENSE`), as PeleLMeX. The mechanism files keep the terms of their original sources
(San Diego mechanism, UCSD; Burke et al., *Int. J. Chem. Kinet.* 44 (2012) 444).

## Citation

TBD later.