///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name murat_genp_001   ("murat_genp_ana");
def_name murat_genp_002   ("murat_spmc_ana");
def_name murat_genp_002_01("murat_spmc_ana_old");
def_name murat_genp_003   ("murat_mustop_ana");
def_name murat_genp_004   ("murat_bflash_ana_spmc");
def_name murat_genp_005   ("murat_bflash_ana_vdet");
def_name murat_genp_006   ("murat_g4val_ana");
def_name murat_genp_007   ("murat_dose_ana");
def_name murat_genp_008   ("murat_pbarabs_dose_ana");
def_name murat_genp_009   ("murat_coll31_dose_ana");
def_name murat_genp_010   ("murat_coll32_dose_ana");
def_name murat_genp_011   ("murat_coll3_dose_ana");
def_name murat_genp_012   ("murat_coll1_dose_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void  murat_genp_ana(int DebugBit = -1) {
  m_gen = (TGenAnaModule*) g.x->AddModule("TGenAnaModule",0);  
  if (DebugBit >= 0) m_gen->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  murat_spmc_ana(int DebugBit = -1) {
  murat::m_spmc = (murat::TSpmcAnaModule*) g.x->AddModule("murat::TSpmcAnaModule",0);  
  // murat::m_spmc->SetSpmcBlockName("SpmcBlockVDet");
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

void  mustop_ana(const char* VDetBlockName = "VDetBlock", int DebugBit = -1) {
  m_must = (TMuonStopAnaModule*) g.x->AddModule("TMuonStopAnaModule",0);
  m_must->SetVDetBlockName(VDetBlockName);
  if (DebugBit >= 0) m_must->SetDebugBit(DebugBit,1);
}

void  g4val_ana(int DebugBit = -1) {
  m_g4val = (TG4ValidationModule*) g.x->AddModule("TG4ValidationModule",0);  
  if (DebugBit >= 0) m_g4val->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  dose_ana(const char* Name="SpmcBlock",int DebugBit = -1) {
  m_dose = (TDoseAnaModule*) g.x->AddModule("TDoseAnaModule",0);  
  m_dose->SetSpmcBlockName(Name);
  if (DebugBit >= 0) m_g4val->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  pbarabs_dose_ana(int DebugBit = -1) {
  m_dose = (TDoseAnaModule*) g.x->AddModule("TDoseAnaModule",0);  
  m_dose->SetSpmcBlockName("PbarAbsSpmcBlock");
  if (DebugBit >= 0) m_dose->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  coll31_dose_ana(int DebugBit = -1) {
  m_dose = (TDoseAnaModule*) g.x->AddModule("TDoseAnaModule",0);  
  m_dose->SetSpmcBlockName("Coll31SpmcBlock");
  if (DebugBit >= 0) m_dose->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  coll32_dose_ana(int DebugBit = -1) {
  m_dose = (TDoseAnaModule*) g.x->AddModule("TDoseAnaModule",0);  
  m_dose->SetSpmcBlockName("Coll32SpmcBlock");
  if (DebugBit >= 0) m_dose->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  coll3_dose_ana(int DebugBit = -1) {
  m_coll3 = (TColl3DoseAnaModule*) g.x->AddModule("TColl3DoseAnaModule",0);  
  if (DebugBit >= 0) m_coll3->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  coll1_dose_ana(int DebugBit = -1) {
  m_coll1 = (TColl1DoseAnaModule*) g.x->AddModule("TColl1DoseAnaModule",0);  
  //  m_coll1->SetSpmcBlockName("PbarAbsDiskSpmcBlock");
  if (DebugBit >= 0) m_coll1->SetDebugBit(DebugBit,1);
}

