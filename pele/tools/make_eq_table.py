"""HP-equilibrium table of fresh H2/air mixtures versus mixture fraction z, written in the PelePhysics PMF
format with z in place of X (columns z, T [K], u = 0, rho [g/cm3], mole fractions in PelePhysics order).
Used by PELE/C1 (pelelmex_prob.H) to initialise the burnt band with the local equilibrium state.
Usage: python make_eq_table.py OUT.dat P_BAR TU [NZ=401]"""
import os
import re, sys
import numpy as np
import cantera as ct

MECH = os.path.join(os.environ["PELE_HOME"], "Submodules/PelePhysics/Mechanisms/SanDiegoH2/")
out, p_bar, tu = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
nz = int(sys.argv[4]) if len(sys.argv) > 4 else 401
P = ct.one_atm if p_bar == 1.0 else p_bar * 1e5
order = [s for _, s in sorted((int(i), s) for s, i in re.findall(r"#define (\w+)_ID (\d+)", open(MECH + "mechanism.H").read()))]
g = ct.Solution(MECH + "mechanism.yaml")
rows = []
for z in np.concatenate([np.linspace(0.0, 0.2, 4 * nz // 5), np.linspace(0.2, 1.0, nz // 5 + 1)[1:]]):
    g.TPY = tu, P, {"H2": z, "O2": 0.233 * (1 - z), "N2": 0.767 * (1 - z)} if 0 < z < 1 else \
        ({"H2": 1.0} if z >= 1 else {"O2": 0.233, "N2": 0.767})
    g.equilibrate("HP")
    rows.append([z, g.T, 0.0, g.density * 1e-3] + [g.X[g.species_index(s)] for s in order])
with open(out, "w") as h:
    h.write('VARIABLES = "X" "temp" "u" "rho" ' + " ".join('"%s"' % s for s in order) + "\n")
    h.write(" ZONE I=%d FORMAT=POINT\n" % len(rows))
    np.savetxt(h, np.array(rows), fmt="%.16e")
print("%s: %d points, Tmax = %.0f K" % (out, nz, max(r[1] for r in rows)))
