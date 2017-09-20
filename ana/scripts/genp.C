///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name genp_001("genp_ana");
def_name genp_002("spmc_ana");
def_name genp_003("mustop_ana");
def_name genp_004("bflash_ana_spmc");
def_name genp_005("bflash_ana_vdet");
def_name genp_006("g4val_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void  spmc_ana(int DebugBit = -1) {
  m_spmc = (TStepPointMCAnaModule*) g.x->AddModule("TStepPointMCAnaModule",0);  
  if (DebugBit >= 0) m_spmc->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
// output of stage 3 has only one StepPointMCCollection - ::virtualdetectors
//-----------------------------------------------------------------------------
void  bflash_ana_spmc(int Stage = 1, int DebugBit = -1) {
  m_bfl = (TBeamFlashAnaModule*) g.x->AddModule("TBeamFlashAnaModule",0);
  m_bfl->SetSpmcBlockName("SpmcBlock");
  if (DebugBit >= 0) m_bfl->SetDebugBit(DebugBit,1);
}

void  bflash_ana_vdet(int Stage = 1, int DebugBit = -1) {
  m_bfl = (TBeamFlashAnaModule*) g.x->AddModule("TBeamFlashAnaModule",0);
  m_bfl->SetSpmcBlockName("VdetBlock");
  if (DebugBit >= 0) m_bfl->SetDebugBit(DebugBit,1);
}

void  mustop_ana(int DebugBit = -1) {
  m_must = (TMuonStopAnaModule*) g.x->AddModule("TMuonStopAnaModule",0);  
  if (DebugBit >= 0) m_must->SetDebugBit(DebugBit,1);
}

void  g4val_ana(int DebugBit = -1) {
  m_g4val = (TG4ValidationModule*) g.x->AddModule("TG4ValidationModule",0);  
  if (DebugBit >= 0) m_g4val->SetDebugBit(DebugBit,1);
}

