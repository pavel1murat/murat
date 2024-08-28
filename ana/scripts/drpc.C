///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name murat_drpc_0010 ("murat_drpc_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void  murat_drpc_ana(int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  murat::m_drpc = (murat::TDegraderRpcAnaModule*) g.x->AddModule("murat::TDegraderRpcAnaModule",0);  

  printf("DebugBit:%i\n",DebugBit);
  if (DebugBit >= 0) murat::m_drpc->SetDebugBit(DebugBit,1);
}
