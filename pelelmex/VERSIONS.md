# Code versions

| Code | Commit | Date |
| --- | --- | --- |
| PeleLMeX | `111026b1eafbffd69d77f0e692c3ea27751e0f79` | 2026-10-21 |
| PelePhysics | `10a59e5e18bd728f66f1f87ae1588318f8618dc7` | 2026-09-27 |
| AMReX | `5a08d94fc6ab8acc24ad62eedbe0f87b832ca9d6` | 2026-09-15 |

PeleLMeX: https://github.com/AMReX-Combustion/PeleLMeX (BSD-3-Clause), PelePhysics and AMReX as its submodules.
`pelephysics_constant_schmidt.patch` (apply in Submodules/PelePhysics: `git apply`) extends the Sutherland transport
model with a power-law viscosity (`transport.viscosity_power`) and species Schmidt numbers
(`transport.schmidt_numbers`); without these keys the model is unchanged. `Mechanisms/SanDiegoH2` is the San Diego
mechanism converted for PelePhysics (copy into Submodules/PelePhysics/Mechanisms/).
