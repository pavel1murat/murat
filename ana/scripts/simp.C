///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name simp_001("murat_simp_ana");

//-----------------------------------------------------------------------------
void  murat_simp_ana(int DebugBit = -1) {
  murat::m_sim = (murat::TSimpAnaModule*) g.x->AddModule("murat::TSimpAnaModule",0);  
  if (DebugBit >= 0) murat::m_sim->SetDebugBit(DebugBit,1);
}

