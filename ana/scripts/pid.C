///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name pid_001("pid_ana");
def_name pid_002("pid_tpana");
def_name pid_003("pid_emuana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void  pid_ana(int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  murat::m_pid = (murat::TPidAnaModule*) g.x->AddModule("murat::TPidAnaModule",0);  
  if (DebugBit >= 0) murat::m_pid->SetDebugBit(DebugBit,1);
}


//-----------------------------------------------------------------------------
void  pid_tpana(int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  m_tpa = (TTrackPidAnaModule*) g.x->AddModule("TTrackPidAnaModule",0);  
  if (DebugBit >= 0) m_tpa->SetDebugBit(DebugBit,1);
}




//-----------------------------------------------------------------------------
void  pid_emuana(int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  murat::m_emu = (murat::TEmuModule*) g.x->AddModule("murat::TEmuModule",0);  
  if (DebugBit >= 0) murat::m_emu->SetDebugBit(DebugBit,1);
}





