#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");
  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("xc", PeleLM::prob_parm->xc);
  pp.query("delta0", PeleLM::prob_parm->delta0);
  pp.query("T_in", PeleLM::prob_parm->T_in);
  pp.query("V_in", PeleLM::prob_parm->V_in);
  pp.query("y_hot", PeleLM::prob_parm->y_hot);
  pp.query("w_hot", PeleLM::prob_parm->w_hot);
  pp.query("T_hot", PeleLM::prob_parm->T_hot);
  PeleLM::pmf_data.initialize();   // equilibrium table vs z (pmf.datafile), used for the burnt band
}

void
PeleLM::freeProbParm()
{
}
