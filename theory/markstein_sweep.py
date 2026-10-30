"""Pressure sweep along the isentropic compression path (29/10/2026): Tu = 300 p^0.2856 (300, 579, 706 K at 1, 10,
20 bar), phi = 1.6, San Diego, Cantera mixture-averaged; consumption Markstein number of twin counterflow flames and
predicted slope b = lambda Ma_c,delta exactly as markstein_defs.py (same functions). Gives the continuous prediction
b(p)/b(1 atm) of figure 3(c). Usage (compute node): python markstein_sweep.py OUT.json NPROC"""
import json, sys
from multiprocessing import Pool
import markstein_defs as M

P = [1, 2, 3, 5, 7, 10, 13, 16, 20]
TU = {1: 300.0, 10: 579.0, 20: 706.0}

if __name__ == "__main__":
    jobs = [("sandiego", p, TU.get(p, round(300.0 * p ** 0.2856, 1))) for p in P]
    with Pool(int(sys.argv[2])) as pool:
        res = pool.map(M.run, jobs)
    json.dump(res, open(sys.argv[1], "w"), indent=1)
    b1 = [r for r in res if r["p"] == 1][0]["b_c"]
    for r in res:
        if "error" in r:
            print(r["p"], r["Tu"], "ERROR", r["error"]); continue
        print("p=%5.1f Tu=%5.1f sL=%.3f dL=%.1fum lam=%.3f Ma_c=%.3f b=%.4f b/b1=%.3f"
              % (r["p"], r["Tu"], r["sL"], r["dL"] * 1e6, r["lam"], r["Ma_c"], r["b_c"], r["b_c"] / b1))
