"""Ze (N2-dilution method), Lewis numbers, Le_eff and asymptotic Markstein number
for the six campaign conditions, phi = 0.5, 1.0, 1.6 (San Diego, mixture-averaged)."""
import cantera as ct, numpy as np, json, sys, warnings
from multiprocessing import Pool
from scipy.integrate import quad
warnings.filterwarnings('ignore')
MECH = sys.argv[1]
CASES = [(1, 300), (5, 475), (10, 579), (20, 706), (10, 300), (1, 579)]
PHIS = [0.5, 1.0, 1.6]
DIL = [1.0, 1.06]          # N2/O2 multiplier on 3.76 -> changes Tad at fixed phi


def flame(args):
    p, Tu, phi, d = args
    g = ct.Solution(MECH)
    g.set_equivalence_ratio(phi, 'H2', {'O2': 1, 'N2': 3.76 * d})
    g.TP = Tu, p * 1e5
    rho_u = g.density
    ff = ct.FreeFlame(g, width=0.03 / p)
    ff.transport_model = 'mixture-averaged'
    ff.set_refine_criteria(ratio=3, slope=0.04, curve=0.08)
    ff.solve(loglevel=0, auto=True)
    T, x = ff.T, ff.grid
    out = dict(p=p, Tu=Tu, phi=phi, d=d, sL=ff.velocity[0], rho_u=rho_u,
               rr=rho_u / ff.density[-1], Tad=T[-1],
               dL=(T[-1] - T[0]) / np.max(np.gradient(T, x)))
    if d == 1.0:   # Lewis numbers in the unburnt mixture
        g.TP = Tu, p * 1e5
        g.transport_model = 'mixture-averaged'
        a = g.thermal_conductivity / (g.density * g.cp_mass)
        D = g.mix_diff_coeffs
        out['Le_H2'] = a / D[g.species_index('H2')]
        out['Le_O2'] = a / D[g.species_index('O2')]
        out['alpha_u'] = a
        out['D_H2'] = D[g.species_index('H2')]
    return out


if __name__ == '__main__':
    jobs = [(p, T, phi, d) for p, T in CASES for phi in PHIS for d in DIL]
    with Pool(6) as pool:
        res = pool.map(flame, jobs)
    json.dump(res, open(sys.argv[2], 'w'), indent=1)
    I = lambda s: quad(lambda x: np.log(x) / (x - 1), 1, s)[0]
    print(f"{'p':>3} {'Tu':>4} {'phi':>4} {'sL':>6} {'dL_um':>6} {'sig':>5} {'Ze':>5} "
          f"{'LeH2':>5} {'LeO2':>5} {'Leff':>5} {'Ma':>6}")
    for p, T in CASES:
        for phi in PHIS:
            a, b = [r for r in res if r['p'] == p and r['Tu'] == T and r['phi'] == phi]
            m0, m1 = a['rho_u'] * a['sL'], b['rho_u'] * b['sL']
            EaR = -2 * (np.log(m1) - np.log(m0)) / (1 / b['Tad'] - 1 / a['Tad'])
            Ze = EaR * (a['Tad'] - T) / a['Tad'] ** 2
            # deficient (D) / excess (E) reactant, Matalon-Cui-Bechtold effective Le
            if phi < 1: LeD, LeE, Phi = a['Le_H2'], a['Le_O2'], 1 / phi
            else:       LeD, LeE, Phi = a['Le_O2'], a['Le_H2'], phi
            A = 1 + Ze * (Phi - 1)
            Leff = 1 + ((LeE - 1) + (LeD - 1) * A) / (1 + A)
            s = a['rr']
            Ma = s / (s - 1) * np.log(s) + Ze * (Leff - 1) / (2 * (s - 1)) * I(s)
            print(f"{p:3d} {T:4d} {phi:4.1f} {a['sL']:6.2f} {a['dL']*1e6:6.1f} {s:5.2f} {Ze:5.1f} "
                  f"{a['Le_H2']:5.2f} {a['Le_O2']:5.2f} {Leff:5.2f} {Ma:6.2f}")
