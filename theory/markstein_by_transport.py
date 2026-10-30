"""Point 1 (27/10/2026): does the transport model explain the different slopes b of the two transport formulations?
Asymptotic unburnt-gas Markstein number (Clavin & Williams 1982, effective Lewis number of Matalon, Cui &
Bechtold 2003), phi = 1.6, evaluated with the two transport models:
  MA   : Cantera mixture-averaged (= PeleLMeX), Lewis numbers and alpha_u of the fresh gas (groups.json);
  SC   : constant Schmidt numbers (chem/transport/Sc*.json, Chaplet et al. Table 1), Le_k = Sc_k / Pr,
         alpha_u = mu(Tu) / (rho_u Pr), mu = power law; sL and deltaL from the PeleLMeX constant-Schmidt 1D flames
         (sL_sc_pelelmex.csv).
Ze and sigma are taken from Cantera for both (same chemistry and thermodynamics).
Ma_delta = Ma_d * (alpha_u / sL) / deltaL. Usage: python markstein_by_transport.py"""
import csv, json, os
import numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
G = json.load(open(os.path.join(HERE, "groups.json")))
SC1D = list(csv.DictReader(open(os.path.join(HERE, "sL_sc_pelelmex.csv"))))
SCF = {(1, 300): "Sc.json", (10, 579): "Sc_p10_T579.json", (20, 706): "Sc_p20_T706.json",
       (1, 579): "Sc_p1_T579.json", (10, 300): "Sc_p10_T300.json"}
I = lambda s: quad(lambda x: np.log(x) / (x - 1), 1, s)[0]


def ma_d(s, Ze, LeH2, LeO2, phi=1.6):
    LeD, LeE, Phi = (LeO2, LeH2, phi) if phi >= 1 else (LeH2, LeO2, 1 / phi)
    A = 1 + Ze * (Phi - 1)
    Leff = 1 + ((LeE - 1) + (LeD - 1) * A) / (1 + A)
    return s / (s - 1) * np.log(s) + Ze * (Leff - 1) / (2 * (s - 1)) * I(s), Leff


out = {}
for (p, Tu), scf in SCF.items():
    a, b = [r for r in G if r["p"] == p and r["Tu"] == Tu and r["phi"] == 1.6]
    m0, m1 = a["rho_u"] * a["sL"], b["rho_u"] * b["sL"]
    Ze = -2 * (np.log(m1) - np.log(m0)) / (1 / b["Tad"] - 1 / a["Tad"]) * (a["Tad"] - Tu) / a["Tad"] ** 2
    s = a["rr"]
    # mixture-averaged
    Ma_MA, Le_MA = ma_d(s, Ze, a["Le_H2"], a["Le_O2"])
    d_MA = a["alpha_u"] / a["sL"]
    MAd_MA = Ma_MA * d_MA / a["dL"]
    # SC
    S = json.load(open(os.path.join(ROOT, "chem", "transport", scf)))
    Pr = S["Pr"]
    LeH2, LeO2 = S["Sc"]["H2"] / Pr, S["Sc"]["O2"] / Pr
    mu_u = S["mu_ref"] * (Tu / S["T_ref"]) ** S["power"]
    alpha_A = mu_u / (a["rho_u"] * Pr)
    r = [r for r in SC1D if float(r["p_bar"]) == p and float(r["Tu_K"]) == Tu][0]
    sL_A, dL_A = float(r["sL_m_s"]), float(r["deltaL_um"]) * 1e-6
    Ma_A, Le_A = ma_d(s, Ze, LeH2, LeO2)
    MAd_A = Ma_A * (alpha_A / sL_A) / dL_A
    out["p%d_T%d" % (p, Tu)] = dict(sigma=s, Ze=Ze, MA=dict(LeH2=a["Le_H2"], LeO2=a["Le_O2"], Leff=Le_MA, Ma_d=Ma_MA,
                                    d_over_dL=d_MA / a["dL"], Ma_delta=MAd_MA),
                                    SC=dict(LeH2=LeH2, LeO2=LeO2, Leff=Le_A, Ma_d=Ma_A, d_over_dL=alpha_A / sL_A / dL_A,
                                              Ma_delta=MAd_A), ratio_SC_MA=MAd_A / MAd_MA)
    print("p=%2d Tu=%3d  sig=%.2f Ze=%.2f | MA: LeH2 %.2f LeO2 %.2f Leff %.2f Ma_d %.2f d/dL %.4f Ma_del %.4f | "
          "SC: LeH2 %.2f LeO2 %.2f Leff %.2f Ma_d %.2f d/dL %.4f Ma_del %.4f | SC/MA %.2f"
          % (p, Tu, s, Ze, a["Le_H2"], a["Le_O2"], Le_MA, Ma_MA, d_MA / a["dL"], MAd_MA,
             LeH2, LeO2, Le_A, Ma_A, alpha_A / sL_A / dL_A, MAd_A, MAd_A / MAd_MA))
json.dump(out, open(os.path.join(HERE, "markstein_by_transport.json"), "w"), indent=1)
