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
// use MVA-based PID , assume DAR tracks
// - TrackBlockName and MVA training code need to be in sync - training depends on the track reco
// - if MinTrq < 0, the default stored in the training results is used
// defaults: TrackBlockName[0]: "TrackBlockDarDe"
//           TrackBlockName[1]: "TrackBlockDarDmu"
//-----------------------------------------------------------------------------
void  pid_emuana(const char* TrackBlockName=nullptr, int RunningMode = 1010, int DebugBit = -1) {
  murat::m_emu = (murat::TEmuAnaModule*) g.x->AddModule("murat::TEmuAnaModule",0);  

  if (TrackBlockName) murat::m_emu->SetTrackBlockName(0,TrackBlockName);

  
  if (RunningMode >= 0) {
    int imva = 0;
    
    int track_type  = (RunningMode / 1000);
    int channel     = (RunningMode % 1000) / 100;
    int use_trq_mva = (RunningMode %  100) / 10;

    if (use_trq_mva) {
      if      (channel == 0) murat::m_emu->SetTrqMVA(imva,"fele2s51b1",1070);   // 105 MeV/c
      else if (channel == 1) murat::m_emu->SetTrqMVA(imva,"fpos2s51b1",1170);   //  92 MeV/c
    }

    float min_trq = murat::m_emu->fTrqMVA[imva]->CutValue();

    murat::m_emu->fTrackID[0]->SetMinTrkQual(min_trq);

    // PID MVA is always used

    if      (channel == 0) murat::m_emu->SetPidMVA("ele00s61b0",1000);
    else if (channel == 1) murat::m_emu->SetPidMVA("ele01s51b0",1100);
  }


  if (DebugBit >= 0) murat::m_emu->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  pid_emuana_write_mva_tree(int MVATrainingCode = 1070, int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  murat::m_emu = (murat::TEmuAnaModule*) g.x->AddModule("murat::TEmuAnaModule",0);  
  murat::m_emu->SetWriteMvaTree(1);

  if (MVATrainingCode >= 0) {
    int imva = 0;
    
    int channel_code = MVATrainingCode % 1000;
    
    int channel = channel_code / 100;
    
    TString training_dataset = "fele2s51b1";             // 105 MeV/c
    if (channel == 1) training_dataset = "fpos2s51b1";   //  92 MeV/c
  
    murat::m_emu->SetTrqMVA(imva,training_dataset.Data(),MVATrainingCode);

    float min_trq = murat::m_emu->fTrqMVA[imva]->CutValue();

    murat::m_emu->fTrackID[0]->SetMinTrkQual(min_trq);
  }

  if (DebugBit >= 0) murat::m_emu->SetDebugBit(DebugBit,1);
}


