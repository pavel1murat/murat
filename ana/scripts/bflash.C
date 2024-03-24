///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name murat_bflash_004   ("murat_bflash_ana_spmc");
def_name murat_bflash_005   ("murat_bflash_ana_vdet");
///////////////////////////////////////////////////////////////////////////////

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
  m_bfl->SetSpmcBlockName("VDetBlock");
  if (DebugBit >= 0) m_bfl->SetDebugBit(DebugBit,1);
}
