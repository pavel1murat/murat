///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name genp_001("gen_ana");
def_name genp_002("spmc_ana");
def_name genp_002_01("spmc_ana_old");
def_name genp_003("mustop_ana");
def_name genp_004("bflash_ana_spmc");
def_name genp_005("bflash_ana_vdet");
def_name genp_006("g4val_ana");
def_name genp_007("dose_ana");
def_name genp_008("pbarabs_dose_ana");
def_name genp_009("coll31_dose_ana");
def_name genp_010("coll32_dose_ana");
def_name genp_011("coll3_dose_ana");
def_name genp_012("coll1_dose_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void  gen_ana(int DebugBit = -1) {
  m_gen = (TGenAnaModule*) g.x->AddModule("TGenAnaModule",0);  
  if (DebugBit >= 0) m_gen->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  spmc_ana(int DebugBit = -1) {
  m_spmc = (TStepPointMCAnaModule*) g.x->AddModule("TStepPointMCAnaModule",0);  
  if (DebugBit >= 0) m_spmc->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
// before Dec'2018, the data blocks were called 'VdetBlock'
//-----------------------------------------------------------------------------
void  spmc_ana_old(int Stage = 1, int DebugBit = -1) {
  m_spmc = (TStepPointMCAnaModule*) g.x->AddModule("TStepPointMCAnaModule",0);  
  if (Stage == 3) m_spmc->SetVDetBlockName("VdetBlock");
  if (Stage == 2) m_spmc->SetVDetBlockName("VdetBlock");

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

