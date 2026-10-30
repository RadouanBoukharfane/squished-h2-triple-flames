"""PMF file (PelePhysics format) of a 1D premixed H2/air flame computed with Cantera (San Diego, the same
YAML as PelePhysics/Mechanisms/SanDiegoH2): columns X [cm], T [K], u [cm/s], rho [g/cm3], then mole
fractions in the PelePhysics species order (read from mechanism.H).
Usage: python make_pmf.py OUT.dat P_BAR TU PHI [mixture-averaged|multicomponent] [soret]"""
import os
import re, sys
import numpy as np
import cantera as ct

MECH = os.path.join(os.environ["PELE_HOME"], "Submodules/PelePhysics/Mechanisms/SanDiegoH2/")
out, p_bar, tu, phi = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
model = sys.argv[5] if len(sys.argv) > 5 else "mixture-averaged"
soret = len(sys.argv) > 6 and sys.argv[6] == "soret"
P = ct.one_atm if p_bar == 1.0 else p_bar * 1e5
order = [s for _, s in sorted((int(i), s) for s, i in re.findall(r"#define (\w+)_ID (\d+)", open(MECH + "mechanism.H").read()))]
g = ct.Solution(MECH + "mechanism.yaml")
g.set_equivalence_ratio(phi, "H2", "O2:1, N2:3.76")
g.TP = tu, P
f = ct.FreeFlame(g, width=0.03 * ct.one_atm / P)
f.transport_model = model
f.soret_enabled = soret
f.set_refine_criteria(ratio=2, slope=0.03, curve=0.06)
f.solve(loglevel=0, auto=True)
dL = (f.T[-1] - f.T[0]) / np.max(np.gradient(f.T, f.grid))
X = np.array([f.X[g.species_index(s)] for s in order]).T
data = np.column_stack([f.grid * 100, f.T, f.velocity * 100, f.density * 1e-3, X])
with open(out, "w") as h:
    h.write('VARIABLES = "X" "temp" "u" "rho" ' + " ".join('"%s"' % s for s in order) + "\n")
    h.write(" ZONE I=%d FORMAT=POINT\n" % len(f.grid))
    np.savetxt(h, data, fmt="%.16e")
print("%s: p=%g bar Tu=%g phi=%g %s soret=%s sL=%.4f m/s dL=%.1f um L=%.2f cm npts=%d order=%s"
      % (out, p_bar, tu, phi, model, soret, f.velocity[0], dL * 1e6, f.grid[-1] * 100, len(f.grid), order))
