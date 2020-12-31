///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name pid_001("pid_ana");
def_name pid_002("pid_tpana");
def_name pid_003("pid_emuana");
def_name pid_004("pid_emuana_use_trq_mva");
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
void  pid_emuana(const char* TrackBlockName=nullptr, int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  murat::m_emu = (murat::TEmuModule*) g.x->AddModule("murat::TEmuModule",0);  

  if (TrackBlockName) murat::m_emu->SetTrackBlockName(0,TrackBlockName);

  if (DebugBit >= 0) murat::m_emu->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  pid_emuana_use_trq_mva(const char* TrackBlockName=nullptr, int MVATrainingCode = -1, float MinTrq=0.2, int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  murat::m_emu = (murat::TEmuModule*) g.x->AddModule("murat::TEmuModule",0);  

  if (TrackBlockName) murat::m_emu->SetTrackBlockName(0,TrackBlockName);

  if (MVATrainingCode >= 0) {
    murat::m_emu->SetTrqMVA(0,"fele2s51b1",MVATrainingCode);
    murat::m_emu->fTrackID[0]->SetMinTrkQual(MinTrq);
  }

  if (DebugBit >= 0) murat::m_emu->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  pid_emuana_write_mva_tree(int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  murat::m_emu = (murat::TEmuModule*) g.x->AddModule("murat::TEmuModule",0);  
  murat::m_emu->SetWriteMvaTree(1);
  if (DebugBit >= 0) murat::m_emu->SetDebugBit(DebugBit,1);
}





