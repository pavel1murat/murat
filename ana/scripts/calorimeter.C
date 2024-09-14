///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name cal_001("murat_cal_ana");
def_name cal_002("murat_cluster_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// configure jobs
//-----------------------------------------------------------------------------
void  murat_cal_ana(int DebugBit = -1) {
  murat::m_cal = (murat::TCalAnaModule*) g.x->AddModule("murat::TCalAnaModule",0);  
  if (DebugBit >= 0) murat::m_cal->SetDebugBit(DebugBit,1);
}


//-----------------------------------------------------------------------------
void  murat_cluster_ana(int DebugBit = -1) {
  murat::m_cls = (murat::TClusterAnaModule*) g.x->AddModule("murat::TClusterAnaModule",0);  
  if (DebugBit >= 0) murat::m_cls->SetDebugBit(DebugBit,1);
}
