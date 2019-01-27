///////////////////////////////////////////////////////////////////////////////
// job config for RMC rejection analysis
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name rmc_0010 ("rmc_ana");

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void  rmc_ana(double KMax = 90, int DebugBit = -1) {

  m_tcm = (TTrackCompModule*) g.x->AddModule("TTrackCompModule",0);  
  m_tcm->SetPdgCode      (-11);
  m_tcm->SetGeneratorCode(41);

  m_tcm->SetKMaxRMC(KMax);

  if (DebugBit >= 0) m_tcm->SetDebugBit(DebugBit,1);
}

