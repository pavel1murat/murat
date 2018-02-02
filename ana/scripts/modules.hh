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
// #include "murat/ana/TGenAnaModule.hh"
// #include "murat/ana/TPhotosAnaModule.hh"
#include "murat/ana/TPidAnaModule.hh"
#include "murat/ana/TStepPointMCAnaModule.hh"
#include "murat/ana/TStrawHitAnaModule.hh"
// #include "murat/ana/TStnGeneratorModule.hh"
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
// TPhotosAnaModule*          m_pho   = NULL;
TPidAnaModule*             m_pid   = NULL;
// TStnGeneratorModule*       m_stg   = NULL;
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

