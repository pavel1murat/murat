///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name helix_001("helix_ana");
def_name helix_002("helix_ana_old");

//-----------------------------------------------------------------------------
void  helix_ana(int DebugBit = -1) {
  m_hel = (THelixAnaModule*) g.x->AddModule("THelixAnaModule",0);  

  if (DebugBit >= 0) m_hel->SetDebugBit(DebugBit,1);
}
//-----------------------------------------------------------------------------
void  helix_ana_old(int DebugBit = -1) {
  m_hel = (THelixAnaModule*) g.x->AddModule("THelixAnaModule",0);

  m_hel->SetHelixBlockName(0,"HelixBlockTpr");
  m_hel->SetHelixBlockName(1,"HelixBlockCpr");
  m_hel->SetHelixBlockName(2,"HelixBlock"   );

  m_hel->SetTrackSeedBlockName("TrackSeedBlockTpr");

  if (DebugBit >= 0) m_hel->SetDebugBit(DebugBit,1);
}
