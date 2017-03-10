///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name cal_001("cal_ana");
def_name cal_002("cluster_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// configure jobs
//-----------------------------------------------------------------------------
void  cal_ana(int DebugBit = -1) {
  m_cal = (TCalAnaModule*) g.x->AddModule("TCalAnaModule",0);  
  if (DebugBit >= 0) m_cal->SetDebugBit(DebugBit,1);
}


//-----------------------------------------------------------------------------
void  cluster_ana(int DebugBit = -1) {
  m_cls = (TClusterAnaModule*) g.x->AddModule("TClusterAnaModule",0);  
  if (DebugBit >= 0) m_cls->SetDebugBit(DebugBit,1);
}





