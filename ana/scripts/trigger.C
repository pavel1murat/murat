///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name trigger_0010("trigger_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void  trigger_ana(const char* TrackBlockName = "TrackBlockPar", int PdgCode = 11, int MCProcessCode = 7, int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  murat::m_trig = (murat::TTriggerAnaModule*) g.x->AddModule("murat::TTriggerAnaModule",0);

  murat::m_trig->SetTrackBlockName(TrackBlockName);
  murat::m_trig->SetPDGCode(PdgCode);
  murat::m_trig->SetMCProcessCode(MCProcessCode);

  if (DebugBit >= 0) {
    murat::m_trig->SetDebugBit(DebugBit,1);
  }
}

