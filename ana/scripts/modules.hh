#ifndef murat_ana_scripts_modules_hh
#define murat_ana_scripts_modules_hh

#include "murat/ana/TBeamFlashAnaModule.hh"
#include "murat/ana/TCalAnaModule.hh"
#include "murat/ana/TClusterAnaModule.hh"
#include "murat/ana/TColl1DoseAnaModule.hh"
#include "murat/ana/TColl3DoseAnaModule.hh"
#include "murat/ana/TCosmicsAnaModule.hh"
#include "murat/ana/TDioCalibModule.hh"
#include "murat/ana/TDoseAnaModule.hh"
#include "murat/ana/TG4ValidationModule.hh"
#include "murat/ana/TMuonStopAnaModule.hh"
#include "murat/ana/TPidAnaModule.hh"
#include "murat/ana/TStepPointMCAnaModule.hh"
#include "murat/ana/TStrawHitAnaModule.hh"
#include "murat/ana/TTrackAnaModule.hh"
#include "murat/ana/TTrackAnaModuleA.hh"
#include "murat/ana/TTrackCompModule.hh"
#include "murat/ana/TTrackPidAnaModule.hh"
#include "murat/ana/TTrackRecoEffAnaModule.hh"
#include "murat/ana/TTrackStrawHitAnaModule.hh"
#include "murat/ana/TTriggerAnaModule.hh"
#include "murat/ana/TValCalPatRecModule.hh"
#include "murat/ana/TVdetAnaModule.hh"
//-----------------------------------------------------------------------------
// hide CDF  names inside the namespace to preserve them
//-----------------------------------------------------------------------------
namespace cdf_stntuple {
  class TClcAnaModule;
  class TDFCModule;
  class TCalAnaModule;
  class TCesAnaModule;
  class TClusterAnaModule;
  class TConversionFilterModule;
  class TEmFilterModule;
  class TFwdDetAnaModule;
  class TJetAnaModule;
  class TMcAnaModule;
  class TMetAnaModule;
  class TMuoAnaModule;
  class TPesAnaModule;
  class TPhotonAnaModule;
  class TSvtAnaModule;
  class TTrackAnaModule;
  class TTrigAnaModule;
  class TWenuMonModule;
  class TJpsiMonModule;

  TClcAnaModule*           m_clc   = NULL;
  TCalAnaModule*           m_cal   = NULL;
  TCesAnaModule*           m_ces   = NULL;
  TClusterAnaModule*       m_clu   = NULL;
  TConversionFilterModule* m_cnv   = NULL;
  TStnDebugModule*         m_dbg   = NULL;
  TEmFilterModule*         m_emf   = NULL;
  TFwdDetAnaModule*        m_fwd   = NULL;
  TJetAnaModule*           m_jet   = NULL;
  TMcAnaModule*            m_mc    = NULL;
  TMetAnaModule*           m_met   = NULL;   
  TMuoAnaModule*           m_muo   = NULL;
  TPesAnaModule*           m_pes   = NULL;
  TPhotonAnaModule*        m_pho   = NULL;
  TSvtAnaModule*           m_svt   = NULL;
  TTrackAnaModule*         m_trk   = NULL;
  TTrigAnaModule*          m_trg   = NULL;
  
  TWenuMonModule*          m_wen   = NULL;
  TJpsiMonModule*          m_jps   = NULL;
};

//-----------------------------------------------------------------------------
// modules defined in 'murat' package (libmurat_ana.so)
//-----------------------------------------------------------------------------
TBeamFlashAnaModule*       m_bfl   = NULL;
TCalAnaModule*             m_cal   = NULL;
TClusterAnaModule*         m_cls   = NULL;
TColl1DoseAnaModule*       m_coll1 = NULL;
TColl3DoseAnaModule*       m_coll3 = NULL;
TCosmicsAnaModule*         m_cos   = NULL;
TDioCalibModule*           m_dio   = NULL;
TDoseAnaModule*            m_dose  = NULL;
TG4ValidationModule*       m_g4val = NULL;
TGenAnaModule*             m_gen   = NULL;
TMuonStopAnaModule*        m_must  = NULL;
TPidAnaModule*             m_pid   = NULL;
TStrawHitAnaModule*        m_strh  = NULL;
TStepPointMCAnaModule*     m_spmc  = NULL;
TTrackStrawHitAnaModule*   m_tsh   = NULL;
TTrackAnaModule*           m_trk   = NULL;
TTrackAnaModuleA*          m_trka  = NULL;
TTrackCompModule*          m_tcm   = NULL;
TTrackPidAnaModule*        m_tpa   = NULL;
TTrackRecoEffAnaModule*    m_eff   = NULL;
TTriggerAnaModule*         m_trig  = NULL;
TValCalPatRecModule*       m_vcpr  = NULL;
TVdetAnaModule*            m_vdt   = NULL;

// TStnTrackID*             trk_id  = NULL;
// TAnalysisDataset*        a_dset  = NULL;

#endif

