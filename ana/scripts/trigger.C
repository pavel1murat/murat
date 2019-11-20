///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name trigger_0010("trigger_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void  trigger_ana(int PdgCode = 11, int ProcessCode = 7, int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  m_trig = (TTriggerAnaModule*) g.x->AddModule("TTriggerAnaModule",0);
  m_trig->SetPdgCode(PdgCode);
  m_trig->SetProcessCode(ProcessCode);
  if (DebugBit >= 0) {
    m_trig->SetDebugBit(DebugBit,1);
  }
}

