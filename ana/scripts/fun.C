///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name fun_001("murat_fun_ana");

//-----------------------------------------------------------------------------
void  murat_fun_ana(int DebugBit = -1) {
  murat::m_fun = (murat::TFunAnaModule*) g.x->AddModule("murat::TFunAnaModule",0);  
  if (DebugBit >= 0) murat::m_fun->SetDebugBit(DebugBit,1);
}

