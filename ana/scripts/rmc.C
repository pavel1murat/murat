///////////////////////////////////////////////////////////////////////////////
// job config for RMC rejection analysis
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name rmc_0010 ("rmc_ana_000");
def_name rmc_0020 ("rmc_ana");

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void  rmc_ana_000(double KMax = 90, int DebugBit = -1) {

  m_tcm = (TTrackCompModule*) g.x->AddModule("TTrackCompModule",0);  
  m_tcm->SetPdgCode      (-11);
  m_tcm->SetGeneratorCode(41);

  m_tcm->SetKMaxRMC(KMax);

  if (DebugBit >= 0) m_tcm->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void  rmc_ana(double KMax = 90, int DebugBit = -1) {

  m_rmc = (TRMCAnaModule*) g.x->AddModule("TRMCAnaModule",0);  
  m_rmc->SetPdgCode      (-11);
  m_rmc->SetGeneratorCode(41);

  m_rmc->SetKMaxRMC(KMax);

  if (DebugBit >= 0) m_rmc->SetDebugBit(DebugBit,1);
}

