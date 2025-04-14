///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name murat_spmc_002("murat_spmc_ana");
def_name murat_spmc_003("murat_spmc_ana_old");

//-----------------------------------------------------------------------------
void  murat_spmc_ana(int DebugBit = -1, const char* VDetBlockName = "SpmcBlockVDet") {
  murat::m_spmc = (murat::TSpmcAnaModule*) g.x->AddModule("murat::TSpmcAnaModule",0);  
  murat::m_spmc->SetVDetBlockName(VDetBlockName);
  //  murat::m_spmc->SetSpmcBlockName("SpmcBlockVDet");
  if (DebugBit >= 0) murat::m_spmc->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
// before Dec'2018, the data blocks were called 'VdetBlock'
//-----------------------------------------------------------------------------
void  murat_spmc_ana_old(int Stage = 1, int DebugBit = -1) {
  murat::m_spmc = (murat::TSpmcAnaModule*) g.x->AddModule("murat::TSpmcAnaModule",0);  
  if (Stage == 3) murat::m_spmc->SetVDetBlockName("VdetBlock");
  if (Stage == 2) murat::m_spmc->SetVDetBlockName("VdetBlock");

  if (DebugBit >= 0) murat::m_spmc->SetDebugBit(DebugBit,1);
}

