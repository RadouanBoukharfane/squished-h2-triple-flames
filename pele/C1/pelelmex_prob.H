#ifndef PELELM_PROB_H
#define PELELM_PROB_H
// H2/air triple flame of Chaplet, Selle & Poinsot (CNF 283, 2026), case C1 and family, for PeleLMeX
// Adapted from Exec/RegTests/TripleFlame.
// - inlet (y = 0): tanh mixing layer z(x) = 0.5 (1 + tanh(2 (x - xc) / Delta0)), pure H2 for x > xc,
//   air (O2 0.233 / N2 0.767 by mass) for x < xc, T_in, uniform V (active control: FCData.ctrl_V_in);
// - initial state: the same layer everywhere, plus a hot band y > y_hot (T = T_hot, fresh composition,
//   tanh in y of width w_hot) that autoignites and leaves the flame at the target height (no descent).
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>
#include "PMFData.H"
#include <PMF.H>
#include "PelePhysics.H"
#include <PeleLMeX_ProblemSpecificFunctions.H>
#include "PeleLMeX_Index.H"
#include <PeleLMeX_FlowControllerData.H>

struct MyProbParm : public ProbParmDefault
{
  amrex::Real xc = 0.0375;      // centre of the mixing layer [m]
  amrex::Real delta0 = 0.010144; // Delta0 = 1 / max dz/dx at the inlet [m]
  amrex::Real T_in = 300.0;
  amrex::Real V_in = 8.6;
  amrex::Real y_hot = 0.021;
  amrex::Real w_hot = 0.0005;
  amrex::Real T_hot = 1500.0;
  FlowControllerData FCData;
};

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
set_Y_layer(MyProbParm const& pp, const amrex::Real x, amrex::Real Y[])
{
  const amrex::Real z = 0.5 * (1.0 + std::tanh(2.0 * (x - pp.xc) / pp.delta0));
  for (int n = 0; n < NUM_SPECIES; ++n) {
    Y[n] = 0.0;
  }
  Y[H2_ID] = z;
  Y[O2_ID] = 0.233 * (1.0 - z);
  Y[N2_ID] = 0.767 * (1.0 - z);
}

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void initdata(
    int i, int j, int k, int /*is_incompressible*/,
    amrex::Array4<amrex::Real> const& state,
    amrex::Array4<amrex::Real> const& /*aux*/,
    amrex::GeometryData const& geomdata,
    MyProbParm const& pp,
    pele::physics::PMF::PmfData::DataContainer const* pmf_data)
  {
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real Y[NUM_SPECIES] = {0.0};
    set_Y_layer(pp, x, Y);
    // hot band: HP-equilibrium state of the local mixture (Cantera table vs z read through the PMF
    // reader, pmf.datafile = eq_z_*.dat, PELE/tools/make_eq_table.py), blended with the fresh gas by a
    // tanh in y: no autoignition and no dissociation transient
    const amrex::Real w = 0.5 * (1.0 + std::tanh((y - pp.y_hot) / pp.w_hot));
    const amrex::Real z = 0.5 * (1.0 + std::tanh(2.0 * (x - pp.xc) / pp.delta0));
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> eq = {0.0};
    pele::physics::PMF::pmf(pmf_data, z, z, eq);
    amrex::Real Xb[NUM_SPECIES] = {0.0}, Yb[NUM_SPECIES] = {0.0};
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Xb[n] = eq[3 + n];
    }
    eos.X2Y(Xb, Yb);
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y[n] = (1.0 - w) * Y[n] + w * Yb[n];
    }
    state(i, j, k, TEMP) = (1.0 - w) * pp.T_in + w * eq[0];
    state(i, j, k, VELX) = 0.0;
    state(i, j, k, VELY) = pp.V_in;
    amrex::Real rho_cgs, h_cgs;
    const amrex::Real P_cgs = pp.P_mean * 10.0;
    eos.PYT2R(P_cgs, Y, state(i, j, k, TEMP), rho_cgs);
    state(i, j, k, DENSITY) = rho_cgs * 1.0e3;
    eos.TY2H(state(i, j, k, TEMP), Y, h_cgs);
    state(i, j, k, RHOH) = h_cgs * 1.0e-4 * state(i, j, k, DENSITY);
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, FIRSTSPEC + n) = Y[n] * state(i, j, k, DENSITY);
    }
  }

  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void bcnormal(
    const amrex::Real x[AMREX_SPACEDIM],
    const amrex::CellData<amrex::Real>& /*s_in*/,
    amrex::Real s_ext[NVAR],
    const int /*idir*/,
    const int sgn,
    const amrex::Real time,
    amrex::GeometryData const& /*geomdata*/,
    MyProbParm const& pp,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    if (sgn != 1) {
      return;
    }
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real Y[NUM_SPECIES] = {0.0};
    s_ext[VELX] = 0.0;
    if ((pp.FCData.ctrl_active == 0U) || time <= 0.0) {
      s_ext[VELY] = pp.V_in;
    } else {
      s_ext[VELY] = pp.FCData.ctrl_V_in + (time - pp.FCData.ctrl_tBase) * pp.FCData.ctrl_dV;
    }
    s_ext[TEMP] = pp.T_in;
    set_Y_layer(pp, x[0], Y);
    amrex::Real rho_cgs, h_cgs;
    const amrex::Real P_cgs = pp.P_mean * 10.0;
    eos.PYT2R(P_cgs, Y, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;
    eos.TY2H(s_ext[TEMP], Y, h_cgs);
    s_ext[RHOH] = h_cgs * 1.0e-4 * s_ext[DENSITY];
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC + n] = Y[n] * s_ext[DENSITY];
    }
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;
#endif
