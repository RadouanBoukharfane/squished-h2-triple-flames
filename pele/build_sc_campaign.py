"""PeleLMeX constant-Schmidt campaign (28/10/2026): one case per point of the paper, same inlet Delta0+,
leading edge held at the target tip height (1-atm-equivalent mm), V0 = initial guess of sT; constant-Schmidt transport
(constant Schmidt numbers, power-law viscosity, Pr; chem/transport/Sc*.json, CGS) with the executable of
TripleFlame-SC-build (linked, --exe-dir); positive regrid interpolation and species floor. 4 nodes when the
base grid is refined (thin inlet layer). Builds only; submission separate. Usage: python build_sc_campaign.py"""
import json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CHEM = os.path.join(os.path.dirname(HERE), "chem", "transport")
SC = {(1, 300): "Sc.json", (10, 579): "Sc_p10_T579.json", (20, 706): "Sc_p20_T706.json",
      (10, 300): "Sc_p10_T300.json", (1, 579): "Sc_p1_T579.json"}
ORDER = ["N2", "H", "O2", "OH", "O", "H2", "H2O", "HO2", "H2O2"]          # SanDiegoH2 species order
CASES = [  # id, Delta0+, y_tip (1-atm eq. mm), V0 = initial guess of sT, p, Tu
    ("C1-SC", 32, 22.82, 8.850, 1, 300), ("C3-SC", 5.3, 18.57, 8.2475, 1, 300), ("C5-SC", 0.8, 18.77, 7.6225, 1, 300),
    ("C6-SC", 0.8, 2.08, 5.530, 1, 300), ("A6-tanh-SC", 4.4, 9.67, 7.815, 1, 300), ("A6-erf-SC", 1.5, 23.55, 7.7175, 1, 300),
    ("B1-p10-D20-SC", 19, 20.70, 16.900, 10, 579), ("B1-p10-D7-SC", 5, 9.10, 14.712, 10, 579),
    ("B1-p10-D4-SC", 2, 6.77, 13.312, 10, 579), ("B1-p10-D3-SC", 1.5, 4.22, 12.406, 10, 579),
    ("B1-p20-D20-SC", 18.5, 20.60, 19.700, 20, 706), ("B1-p20-D7-SC", 5, 7.15, 16.957, 20, 706),
    ("B1-p20-D5-SC", 3, 5.50, 15.938, 20, 706), ("B1-p20-D4-SC", 2, 3.95, 14.900, 20, 706),
    ("B1-p20-D3-SC", 1.2, 2.62, 13.600, 20, 706),
    ("B2-D32-SC", 32, 18.60, 6.8885, 10, 300), ("B2-D7-SC", 5, 3.12, 5.330, 10, 300), ("B2-D5-SC", 3, 1.77, 4.550, 10, 300),
    ("B3-D32-SC", 32, 22.05, 18.2965, 1, 579), ("B3-D7-SC", 5, 14.97, 16.6602, 1, 579), ("B3-D3-SC", 1.5, 6.75, 14.5171, 1, 579)]


def transport(p, tu):
    s = json.load(open(os.path.join(CHEM, SC[(p, tu)])))
    return ("# constant-Schmidt transport at %g bar / %g K (constant Schmidt numbers, chem/transport/%s; CGS)\n" % (p, tu, SC[(p, tu)]) +
            "transport.viscosity_mu_ref = %.12e\n" % (10.0 * s["mu_ref"]) +
            "transport.viscosity_T_ref = %.1f\n" % s["T_ref"] +
            "transport.viscosity_power = %.12f\n" % s["power"] +
            "transport.Prandtl_number = %.12f\n" % s["Pr"] +
            "transport.schmidt_numbers = %s\n" % " ".join("%.2f" % s["Sc"][k] for k in ORDER) +
            "transport.const_bulk_viscosity = 0.0\n"
            "peleLM.regrid_interp_method = 0\npeleLM.floor_species = 1\n")


if __name__ == "__main__":
    for cid, d0, yt, v0, p, tu in CASES:
        if os.path.exists(os.path.join(HERE, cid)):
            print(cid, "exists, skipped"); continue
        cmd = [sys.executable, os.path.join(HERE, "tools", "make_pele_case.py"), cid, repr(d0), repr(yt), repr(v0),
               "--exe-dir", "TripleFlame-SC-build", "--hours", "24"]
        if (p, tu) != (1, 300):
            cmd += ["--p", repr(p), "--tu", repr(tu)]
        out = subprocess.run(cmd, capture_output=True, text=True)
        print(out.stdout.strip().splitlines()[-1] if out.stdout.strip() else out.stderr[-500:], flush=True)
        inp = os.path.join(HERE, cid, "input.case")
        open(inp, "a").write(transport(p, tu))
        txt = open(inp).read()
        if re.search(r"^amr.n_cell\s*=\s*960", txt, re.M):
            sb = os.path.join(HERE, cid, "run.sbatch")
            s = open(sb).read()                       # read first: open(sb, "w") would empty the file before the read
            open(sb, "w").write(s.replace("--nodes=2", "--nodes=4"))
