///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name pid_001("pid_ana");
def_name pid_002("pid_tpana");
def_name pid_003("pid_emuana");
def_name pid_005("pid_emuana_write_mva_tree");
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
// use MVA-based PID 
// - TrackBlockName and MVA training code need to be in sync - training depends on the track reco
// - if MinTrq < 0, the default stored in the training results is used
// defaults: TrackBlockName[0]: "TrackBlockDarDe"
//           TrackBlockName[1]: "TrackBlockDarDmu"
//-----------------------------------------------------------------------------
void  pid_emuana(const char* TrackBlockName=nullptr, int MVATrainingCode = 1070, float MinTrq=-1, int DebugBit = -1) {
  murat::m_emu = (murat::TEmuAnaModule*) g.x->AddModule("murat::TEmuAnaModule",0);  

  if (TrackBlockName) murat::m_emu->SetTrackBlockName(0,TrackBlockName);

  if (MVATrainingCode >= 0) {
    int imva = 0;
    murat::m_emu->SetTrqMVA(imva,"fele2s51b1",MVATrainingCode);

    float min_trq = MinTrq;
    if (min_trq < 0) min_trq = murat::m_emu->fTrqMVA[imva]->CutValue();

    murat::m_emu->fTrackID[0]->SetMinTrkQual(min_trq);
  }

  if (DebugBit >= 0) murat::m_emu->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  pid_emuana_write_mva_tree(int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  murat::m_emu = (murat::TEmuAnaModule*) g.x->AddModule("murat::TEmuAnaModule",0);  
  murat::m_emu->SetWriteMvaTree(1);
  if (DebugBit >= 0) murat::m_emu->SetDebugBit(DebugBit,1);
}


