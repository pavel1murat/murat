///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name reco_eff_ana_001("reco_eff_ana");
def_name reco_eff_ana_002("reco_eff_ana_no_mc");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void  reco_eff_ana(int DebugBit = -1, int TrkPatRecOnly = 0) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  m_eff = (TTrackRecoEffAnaModule*) g.x->AddModule("TTrackRecoEffAnaModule",0);  
  
  m_eff->SetTrkPatRecOnly(TrkPatRecOnly);

  if (DebugBit >= 0) m_eff->SetDebugBit(DebugBit,1);

}

//-----------------------------------------------------------------------------
void  reco_eff_ana_no_mc(int DebugBit = -1, int TrkPatRecOnly = 0) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  m_eff = (TTrackRecoEffAnaModule*) g.x->AddModule("TTrackRecoEffAnaModule",0);  
  
  m_eff->SetTrkPatRecOnly(TrkPatRecOnly);

  if (DebugBit >= 0) m_eff->SetDebugBit(DebugBit,1);

  m_eff->SetMinNMCHits(-1);
  m_eff->SetMinMCMomentum(-1);
  m_eff->SetMinMCPitch(-1e12);
  m_eff->SetMaxMCPitch( 1e12);
}


