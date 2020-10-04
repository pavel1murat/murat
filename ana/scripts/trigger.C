///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name trigger_0010("trigger_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void  trigger_ana(int PdgCode = 11, int GeneratorCode = 7, int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  murat::m_trig = (murat::TTriggerAnaModule*) g.x->AddModule("murat::TTriggerAnaModule",0);
  murat::m_trig->SetPdgCode(PdgCode);
  murat::m_trig->SetGeneratorCode(GeneratorCode);
  if (DebugBit >= 0) {
    murat::m_trig->SetDebugBit(DebugBit,1);
  }
}

