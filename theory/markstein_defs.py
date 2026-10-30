"""Mechanism sensitivity and Markstein definitions (27/10/2026, JOURNAL), phi = 1.6, Cantera mixture-averaged.
For each mechanism (San Diego 2016, Burke et al. 2012) and condition (p, Tu):
  free flame       : sL, deltaL, sigma = rho_u/rho_b (phi = 1.6) and sigma(phi = 1) (sT,th), sT,th = sigma1^0.5 sL;
  twin counterflow : fresh against fresh, strain K = max |du/dx| ahead of the flame (as R0/cantera_ma_ext.py);
     L_c  : consumption speed S_c = -int w_D dx / (rho_u Y_D,u);
     L_d  : density-weighted displacement speed of the T = 1000 K isotherm (the leading-edge definition of the
            triple-flame runs), S_d = rho(x_1000) |u(x_1000)| / rho_u;
     L_d0 : same at the T = Tu + 10 K isotherm (unburnt-side edge);
  Markstein lengths = -dS/dK (linear fit, Ka < 0.3); Markstein numbers in deltaL units; predicted slope
  b = lambda Ma_delta, lambda = 1/(1 + alpha), alpha = 1 - 1/sigma(phi = 1.6).
Usage (compute node): python markstein_defs.py OUT.json NPROC"""
import json, os, sys, warnings
import numpy as np
import cantera as ct
from multiprocessing import Pool
warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MECHS = {"sandiego": ROOT + "/chem/sandiego_h2.yaml", "burke2012": ROOT + "/chem/burke2012_h2.yaml"}
CASES = [(1, 300), (10, 579), (20, 706), (1, 579), (10, 300)]
KA = (0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.14, 0.18)
PHI = 1.6
Pa = lambda p: 101325.0 if p == 1 else p * 1e5


def free(mech, p, Tu, phi):
    g = ct.Solution(mech)
    g.set_equivalence_ratio(phi, "H2", "O2:1, N2:3.76")
    g.TP = Tu, Pa(p)
    rho_u = g.density
    f = ct.FreeFlame(g, width=0.03 * 101325.0 / Pa(p))
    f.transport_model = "mixture-averaged"
    f.set_refine_criteria(ratio=3, slope=0.04, curve=0.08)
    f.solve(loglevel=0, auto=True)
    return f.velocity[0], (f.T[-1] - f.T[0]) / np.max(np.gradient(f.T, f.grid)), rho_u / f.density[-1]


def iso(x, T, val):
    k = np.where(T >= val)[0][0]
    return k - 1 + (val - T[k - 1]) / (T[k] - T[k - 1])


def run(args):
    name, p, Tu = args
    mech = MECHS[name]
    try:
        sL, dL, sig = free(mech, p, Tu, PHI)
        _, _, sig1 = free(mech, p, Tu, 1.0)
        g = ct.Solution(mech)
        g.set_equivalence_ratio(PHI, "H2", "O2:1, N2:3.76")
        g.TP = Tu, Pa(p)
        rho_u = g.density
        kD = g.species_index("O2")
        YD, WD = g.Y[kD], g.molecular_weights[kD]
        rows = []
        for ka in KA:
            a = ka * sL / dL
            W = (2.5 / ka + 15.0) * dL
            g.set_equivalence_ratio(PHI, "H2", "O2:1, N2:3.76")
            g.TP = Tu, Pa(p)
            f = ct.CounterflowTwinPremixedFlame(g, width=W)
            f.transport_model = "mixture-averaged"
            f.reactants.mdot = rho_u * a * W
            f.set_refine_criteria(ratio=2.5, slope=0.02, curve=0.04, prune=0.005)
            try:
                f.solve(loglevel=0, auto=True)
            except ct.CanteraError:
                continue
            x, u, T, rho = f.grid, f.velocity, f.T, f.density
            if T.max() < Tu + 500:
                continue
            kf = np.argmax(np.gradient(T, x))
            kmin = np.argmin(u[:kf + 1])
            K = np.max(np.abs(np.gradient(u, x))[:max(kmin, 2)])
            Sc = -np.trapezoid(f.net_production_rates[kD] * WD, x) / (rho_u * YD)
            out = [K, Sc]
            for Tiso in (1000.0, Tu + 10.0):
                s = iso(x, T, Tiso)
                k0 = int(np.floor(s)); w = s - k0
                ui = (1 - w) * u[k0] + w * u[k0 + 1]; ri = (1 - w) * rho[k0] + w * rho[k0 + 1]
                out.append(ri * abs(ui) / rho_u)
            rows.append(out)
        A = np.array(rows)
        sel = A[:, 0] * dL / sL < 0.3
        res = dict(mech=name, p=p, Tu=Tu, sL=sL, dL=dL, sigma=sig, sigma1=sig1, sTth=np.sqrt(sig1) * sL,
                   lam=1.0 / (1.0 + (1.0 - 1.0 / sig)), K=list(A[:, 0]))
        for j, key in ((1, "c"), (2, "d1000"), (3, "d0")):
            c = np.polyfit(A[sel, 0], A[sel, j], 1)
            res["L_" + key] = -c[0]
            res["Ma_" + key] = -c[0] / dL
            res["b_" + key] = res["lam"] * (-c[0] / dL)
            res["S_" + key] = list(A[:, j])
        return res
    except Exception as e:                                    # noqa: BLE001
        return dict(mech=name, p=p, Tu=Tu, error=str(e))


if __name__ == "__main__":
    jobs = [(m, p, T) for m in MECHS for p, T in CASES]
    with Pool(int(sys.argv[2])) as pool:
        res = pool.map(run, jobs)
    json.dump(res, open(sys.argv[1], "w"), indent=1)
    for r in res:
        if "error" in r:
            print(r["mech"], r["p"], r["Tu"], "ERROR", r["error"]); continue
        print("%-9s p=%2d Tu=%3d sL=%.3f dL=%.1fum sTth=%.3f lam=%.3f | Ma_c=%.3f b=%.3f | Ma_d1000=%.3f b=%.3f | Ma_d0=%.3f b=%.3f"
              % (r["mech"], r["p"], r["Tu"], r["sL"], r["dL"] * 1e6, r["sTth"], r["lam"], r["Ma_c"], r["b_c"],
                 r["Ma_d1000"], r["b_d1000"], r["Ma_d0"], r["b_d0"]))
