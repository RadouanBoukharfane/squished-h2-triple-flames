"""1D flames phi = 1.6 with the constant-Schmidt transport (constant Schmidt numbers) at the four other conditions, for
sT,th of the constant-Schmidt campaign (28/10/2026). From FlameSheetH2-SC (1 atm): lengths x ls, times x ts (same
values as make_pele_case.py), pressure, PMF from tools/make_pmf.py (mixture-averaged Cantera, initial state only),
standoff placing the flame at the same fraction of the domain as at 1 atm, transport of the condition (chem/transport/Sc*.json).
Usage: python build_sc_1d.py"""
import json, os, re, shutil, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CHEM = os.path.join(os.path.dirname(HERE), "chem", "transport")
REF = os.path.join(HERE, "FlameSheetH2-SC")
ORDER = ["N2", "H", "O2", "OH", "O", "H2", "H2O", "HO2", "H2O2"]
COND = {"p10_T579": (10, 579, 0.08021700269911883, 0.03195295426140144, "Sc_p10_T579.json"),
        "p20_T706": (20, 706, 0.0383466, 0.0118, "Sc_p20_T706.json"),
        "p10_T300": (10, 300, 0.0674, 0.0850, "Sc_p10_T300.json"),
        "p1_T579": (1, 579, 1.3735, 0.4978, "Sc_p1_T579.json")}


def x_mid(pmf):
    d = np.loadtxt(pmf, skiprows=2)
    x, T = d[:, 0] * 1e-2, d[:, 1]
    Tm = 0.5 * (T[0] + T[-1])
    k = np.where(T >= Tm)[0][0]
    return x[k - 1] + (Tm - T[k - 1]) / (T[k] - T[k - 1]) * (x[k] - x[k - 1])


ref = open(os.path.join(REF, "input.1d_phi1.6")).read()
g = lambda key, s=ref: re.search(r"^%s\s*=\s*(.*)$" % re.escape(key), s, re.M).group(1).strip()
Lx, Ly = [float(v) for v in g("geometry.prob_hi").split()]
frac = (float(g("prob.standoff")) + x_mid(os.path.join(REF, g("pmf.datafile").strip('"')))) / Ly
for tag, (p, tu, ls, ts, scf) in COND.items():
    case = os.path.join(HERE, "FlameSheetH2-SC-" + tag)
    if os.path.exists(case):
        print(tag, "exists"); continue
    os.makedirs(case)
    pmf = "pmf_%s_phi1.6_MA.dat" % tag
    subprocess.run([sys.executable, os.path.join(HERE, "tools", "make_pmf.py"), os.path.join(case, pmf), repr(p), repr(tu), "1.6"],
                   check=True, stdout=subprocess.DEVNULL)
    s = json.load(open(os.path.join(CHEM, scf)))
    t = ref
    rep = {"geometry.prob_hi": "%.6e %.6e" % (Lx * ls, Ly * ls),
           "prob.P_mean": "%.1f" % (101325.0 if p == 1 else p * 1e5),
           "prob.standoff": "%.6e" % (frac * Ly * ls - x_mid(os.path.join(case, pmf))),
           "pmf.datafile": '"%s"' % pmf,
           "amr.stop_time": "%.4e" % (float(g("amr.stop_time")) * ts),
           "amr.plot_per": "%.4e" % (float(g("amr.plot_per")) * ts),
           "amr.check_per": "%.4e" % (float(g("amr.check_per")) * ts),
           "transport.viscosity_mu_ref": "%.12e" % (10.0 * s["mu_ref"]),
           "transport.viscosity_T_ref": "%.1f" % s["T_ref"],
           "transport.viscosity_power": "%.12f" % s["power"],
           "transport.Prandtl_number": "%.12f" % s["Pr"],
           "transport.schmidt_numbers": " ".join("%.2f" % s["Sc"][k] for k in ORDER)}
    for k, v in rep.items():
        t, n = re.subn(r"^(%s\s*=\s*).*$" % re.escape(k), lambda m: m.group(1) + v, t, flags=re.M)
        assert n == 1, k
    open(os.path.join(case, "input.1d_phi1.6"), "w").write(t)
    os.symlink(os.path.join(REF, "PeleLMeX2d.gnu.MPI.SanDiegoH2.ex"), os.path.join(case, "PeleLMeX2d.gnu.MPI.SanDiegoH2.ex"))
    sb = open(os.path.join(HERE, "FlameSheetH2", "run.sbatch")).read()
    sb = sb.replace("--job-name=pele_1d", "--job-name=pele_1d_SC_" + tag).replace("--ntasks-per-node=16", "--ntasks-per-node=32")
    sb = re.sub(r"^# PeleLMeX verification.*$", "# 1D flame phi = 1.6, %s, constant-Schmidt transport (constant Schmidt numbers)" % tag, sb, flags=re.M)
    open(os.path.join(case, "run.sbatch"), "w").write(sb)
    print(tag, "built: Ly = %.3e m, standoff = %s, stop = %s" % (Ly * ls, rep["prob.standoff"], rep["amr.stop_time"]), flush=True)
