///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name murat_dose_007   ("murat_dose_ana");
def_name murat_dose_008   ("murat_pbarabs_dose_ana");
def_name murat_dose_009   ("murat_coll31_dose_ana");
def_name murat_dose_010   ("murat_coll32_dose_ana");
def_name murat_dose_011   ("murat_coll3_dose_ana");
def_name murat_dose_012   ("murat_coll1_dose_ana");
///////////////////////////////////////////////////////////////////////////////

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

