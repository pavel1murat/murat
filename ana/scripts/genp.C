///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name genp_001("genp_ana");
def_name genp_002("spmc_ana");
def_name genp_003("bflash_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void  smpc_ana(int DebugBit = -1) {
  m_spmc = (TStepPointMCAnaModule*) g.x->AddModule("TStepPointMCAnaModule",0);  
  if (DebugBit >= 0) m_spmc->SetDebugBit(DebugBit,1);
}

void  bflash_ana(int DebugBit = -1) {
  m_bfl = (TBeamFlashAnaModule*) g.x->AddModule("TBeamFlashAnaModule",0);  
  if (DebugBit >= 0) m_bfl->SetDebugBit(DebugBit,1);
}

