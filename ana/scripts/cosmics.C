///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name cosmics_001("cosmics_ana");
def_name cosmics_002("cosmics_ralf");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
void  cosmics_ana(const char* TrackBlockName, int DebugBit = -1) {

  murat::m_cos = (murat::TCosmicsAnaModule*) g.x->AddModule("murat::TCosmicsAnaModule",0);

  if (TrackBlockName) murat::m_cos->SetTrackBlockName(TrackBlockName);

  if (DebugBit >= 0 ) murat::m_cos->SetDebugBit(DebugBit,1);
}


void cosmics_ralf() {
//-----------------------------------------------------------------------------
// strip cosmic-like Z's 
//-----------------------------------------------------------------------------
  int run_event [] = {
    15792036,287384,    // looks OK, but may be cosmics
    15792043,373517,
    15792072,438616,
    15792072,321389,
    15794482,707397,
    15794484,242583,
    15794485,186161,
    15794485,801503,
    15794485,758874,
    15799260,961365,
    15799263,240006,
    15799263,899123,
    15799263,718730,
    15799263,146445,
    15799266, 11122,
    15799266,600477,
    15799279,319063,
    15800060,452944,
    15800077,448781,
    15800077,863692,
    16467078,672455,
    16467078,903744,
    16467079,318358,
    16467079,137813,
    16467079, 24884,
    16467097,418894,
    16467477,502784,
    16467477,342543,
    16467478,317445,
    16467479,625785,
    16467479,171935,
    16467479,912702,
    16467479,367288,
    16467480,113124,
    16467482,225144,
    16467483,616276,
    16467483,990933,
    16467483,756733,
    16467648,137165,
    16467648,208953,
    -1
  };
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  g.x->SetNEventsToReport(5000);
  g.x->SetEventList(run_event);

  murat::m_cos = (murat::TCosmicsAnaModule*) g.x->AddModule("murat::TCosmicsAnaModule",0);
  murat::m_cos->SetDebugBit(39,1);
}
