"""Build a PeleLMeX triple-flame case PELE/<ID>/ from PELE/C1 (same executable, problem files and run
script): input file scaled with the flame thickness (lengths given as 1-atm equivalents, multiplied
by ls = deltaL(p, Tu)/deltaL(1 atm, 300 K)), HP-equilibrium table vs z for the burnt band.
Usage: python make_pele_case.py ID DELTA0_PLUS Y_TARGET_MM V0 [--p BAR] [--tu K] [--soret] [--hours H]
  C4        : python make_pele_case.py C4 3.2 18.7 7.6
  B1-p10-D20: python make_pele_case.py B1-p10-D20 19 18.7 16.0 --p 10 --tu 579"""
import argparse, os, re, shutil, subprocess, sys
import numpy as np
import cantera as ct

HERE = os.path.dirname(os.path.abspath(__file__))
PELE = os.path.dirname(HERE)
MECH = os.path.join(os.environ["PELE_HOME"], "Submodules/PelePhysics/Mechanisms/SanDiegoH2/mechanism.yaml")
ap = argparse.ArgumentParser()
ap.add_argument("id"); ap.add_argument("d0p", type=float); ap.add_argument("yt_mm", type=float); ap.add_argument("v0", type=float)
ap.add_argument("--p", type=float, default=1.0); ap.add_argument("--tu", type=float, default=300.0)
ap.add_argument("--soret", action="store_true"); ap.add_argument("--hours", type=float, default=24.0)
ap.add_argument("--cold", action="store_true", help="no reaction, no hot band, no active control: Delta(y) at V0")
ap.add_argument("--exe-dir", default="C1",
                help="PELE/<dir> holding the executable to link (e.g. TripleFlame-SC-build for the constant-Schmidt transport); "
                     "the case gets a symbolic link: never copy over it (cp writes through the link)")
a = ap.parse_args()


def flame(p, tu):
    g = ct.Solution(MECH)
    g.set_equivalence_ratio(1.6, "H2", "O2:1, N2:3.76")
    P = ct.one_atm if p == 1.0 else p * 1e5
    g.TP = tu, P
    f = ct.FreeFlame(g, width=0.03 * ct.one_atm / P)
    f.transport_model = "mixture-averaged"
    f.set_refine_criteria(ratio=3, slope=0.05, curve=0.1)
    f.solve(loglevel=0, auto=True)
    return (f.T[-1] - f.T[0]) / np.max(np.gradient(f.T, f.grid)), f.velocity[0]


ref = a.p == 1.0 and a.tu == 300.0
dl1, sl1 = flame(1.0, 300.0)
dlp, slp = (dl1, sl1) if ref else flame(a.p, a.tu)
ls, ts = dlp / dl1, (dlp / slp) / (dl1 / sl1)
dL = 317.0e-6 * ls
case = os.path.join(PELE, a.id)
os.makedirs(case)
for f in ("PeleLMeX2d.gnu.MPI.SanDiegoH2.ex",):
    os.symlink(os.path.join(PELE, a.exe_dir, f), os.path.join(case, f))
eq = "eq_z_p%g_T%g.dat" % (a.p, a.tu)
subprocess.run([sys.executable, os.path.join(HERE, "make_eq_table.py"), os.path.join(case, eq), repr(a.p), repr(a.tu), "1001"], check=True)
s = open(os.path.join(PELE, "C1", "input.c1")).read()
P = 101325.0 if a.p == 1.0 else a.p * 1e5
L = lambda m: "%.6e" % (m * ls)
rep = {r"geometry.prob_hi\s*=.*": "geometry.prob_hi     = %.15e %.15e" % (0.075 * ls, 0.075 * ls * 144 / 240),   # exact 0.6 ratio: PeleLMeX asserts dx = dy to 10 ulp
       r"prob.P_mean = .*": "prob.P_mean = %.1f" % P, r"prob.xc = .*": "prob.xc = %s" % L(0.0375),
       r"prob.delta0 = .*": "prob.delta0 = %.6e" % (a.d0p * dL), r"prob.T_in = .*": "prob.T_in = %g" % a.tu,
       r"prob.V_in = .*": "prob.V_in = %g" % a.v0, r"prob.y_hot = .*": "prob.y_hot = %s" % L((a.yt_mm + 1.0) * 1e-3),
       r"prob.w_hot = .*": "prob.w_hot = %s" % L(0.001), r"amr.stop_time = .*": "amr.stop_time = %.4e" % (0.025 * ts),
       r"amr.max_wall_time = .*": "amr.max_wall_time = %.1f" % (a.hours - 0.5),
       r"amr.plot_per = .*": "amr.plot_per = %.4e" % (5e-4 * ts), r"amr.check_per = .*": "amr.check_per = %.4e" % (2.5e-3 * ts),
       r"active_control.tau = .*": "active_control.tau = %.4e" % (1e-3 * ts),
       r"active_control.height = .*": "active_control.height = %s" % L(a.yt_mm * 1e-3),
       r"active_control.velMax = .*": "active_control.velMax = %.1f" % (2.0 * a.v0),
       r"active_control.changeMax = .*": "active_control.changeMax = %.4f" % (0.05 * a.v0 / 8.6),
       r'pmf.datafile = .*': 'pmf.datafile = "%s"' % eq}
for k, v in rep.items():
    s, n = re.subn("^" + k + "$", v, s, flags=re.M)
    assert n == 1, k
s = s.replace("# Case C1 of Chaplet", "# Case %s (built by PELE/tools/make_pele_case.py from C1): p = %g bar, Tu = %g K, Delta0+ = %g,\n# y_target = %g mm (1-atm equivalent), soret = %s.\n# Template: case C1 of Chaplet" % (a.id, a.p, a.tu, a.d0p, a.yt_mm, a.soret), 1)
# resolve the inlet mixing layer when it is thin (Delta0 < 3 base cells)
if a.d0p * dL < 10 * 0.075 * ls / 240:      # Delta0 below 10 base cells
    s += "amr.refinement_indicators = hot gradT HR layer\namr.layer.max_level = %d\namr.layer.adjacent_difference_greater = 0.05\namr.layer.field_name = Y(H2)\n" % (4 if a.d0p * dL < 10 * 0.075 * ls / 240 / 4 else 2)   # thin layer: refine to the finest level
    s = re.sub(r"^amr.refinement_indicators = hot gradT HR$", "", s, count=1, flags=re.M)
# thin inlet layer: the base grid itself must resolve it (>= 3 cells over Delta0), otherwise the run aborts
# at t ~ 0.15 ms (C5/C6, B1-p10-D3, JOURNAL 25/10 23h); refine level 0 and drop levels to keep the finest dx
k = 0
while a.d0p * dL < 3 * 0.075 * ls / (240 * 2 ** k) and k < 2:
    k += 1
if k:
    s = re.sub(r"^amr.n_cell\s*=.*$", "amr.n_cell          = %d %d" % (240 * 2 ** k, 144 * 2 ** k), s, flags=re.M)
    for key in ("amr.max_level", "amr.gradT.max_level", "amr.HR.max_level", "amr.layer.max_level"):
        s = re.sub(r"^(%s\s*=\s*)(\d+)$" % re.escape(key), lambda m: m.group(1) + str(max(int(m.group(2)) - k, 1)), s, flags=re.M)
    s = re.sub(r"^amr.ref_ratio\s*=.*$", "amr.ref_ratio       = " + " ".join(["2"] * (4 - k)), s, flags=re.M)
    s = re.sub(r"^amr.n_error_buf\s*=.*$", "amr.n_error_buf     = " + " ".join(["2"] * (4 - k)), s, flags=re.M)
if a.soret:
    s += "transport.use_soret = 1\n"
if a.cold:      # cold companion (A3-like): 2 flow-throughs of the domain at V0, fresh layer everywhere
    s = re.sub(r"^active_control.on = .*$", "active_control.on = 0", s, flags=re.M)
    s = re.sub(r"^prob.y_hot = .*$", "prob.y_hot = 1.0e3", s, flags=re.M)
    s = re.sub(r"^amr.stop_time = .*$", "amr.stop_time = %.4e" % (2 * 0.045 * ls / a.v0), s, flags=re.M)
    s = re.sub(r"^amr.plot_per = .*$", "amr.plot_per = %.4e" % (0.25 * 0.045 * ls / a.v0), s, flags=re.M)
    # chemistry left on (do_react = 0 crashes in this PeleLMeX version): no hot gas, so nothing reacts
open(os.path.join(case, "input.case"), "w").write(s)
sb = open(os.path.join(PELE, "C1", "run.sbatch")).read().replace("pele_C1", "pele_" + a.id).replace("input.c1", "input.case")
sb = re.sub(r"--time=\S+", "--time=%02d:00:00" % int(a.hours), sb)
open(os.path.join(case, "run.sbatch"), "w").write(sb)
print("%s: p=%g Tu=%g ls=%.4f ts=%.4f dL=%.2f um Delta0=%.4f mm y_target=%.3f mm V0=%g soret=%s" %
      (case, a.p, a.tu, ls, ts, dL * 1e6, a.d0p * dL * 1e3, a.yt_mm * ls, a.v0, a.soret))
