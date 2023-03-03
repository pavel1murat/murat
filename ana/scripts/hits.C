///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name hit_001("straw_hit_ana");
def_name hit_002("track_straw_hit_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// by default, assume murat/test/read_reco_stn_tcn.fcl and TrkPatRec tracks
//-----------------------------------------------------------------------------
void  straw_hit_ana(int DebugBit = -1) {
  m_strh = (TStrawHitAnaModule*) g.x->AddModule("TStrawHitAnaModule",0);  

  if (DebugBit >= 0) m_strh->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  track_straw_hit_ana(int DebugBit = -1) {
  m_tsh = (TTrackStrawHitAnaModule*) g.x->AddModule("TTrackStrawHitAnaModule",0);  

  m_tsh->SetHitBlockName("TrackHitBlock");
  
  if (DebugBit >= 0) m_tsh->SetDebugBit(DebugBit,1);
}



//-----------------------------------------------------------------------------
// by default, assume murat/test/read_reco_stn_tcn.fcl and TrkPatRec tracks
//-----------------------------------------------------------------------------
// void  vst_ana(int DebugBit = -1) {
//   m_vst = (TVstAnaModule*) g.x->AddModule("TVstAnaModule",0);  

//   if (DebugBit >= 0) m_vst->SetDebugBit(DebugBit,1);
// }



