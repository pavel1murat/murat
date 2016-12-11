#ifndef murat_ana_scripts_modules_hh
#define murat_ana_scripts_modules_hh

#include "murat/ana/TCalAnaModule.hh"
#include "murat/ana/TCosmicsAnaModule.hh"
#include "murat/ana/TDioCalibModule.hh"
// #include "murat/ana/TGenAnaModule.hh"
#include "murat/ana/TPidAnaModule.hh"
#include "murat/ana/TStrawHitAnaModule.hh"
#include "murat/ana/TTrackAnaModule.hh"
#include "murat/ana/TTrackAnaModuleA.hh"
#include "murat/ana/TTrackCompModule.hh"
#include "murat/ana/TTrackPidAnaModule.hh"
#include "murat/ana/TTrackRecoEffAnaModule.hh"
#include "murat/ana/TValCalPatRecModule.hh"
#include "murat/ana/TVdetAnaModule.hh"
#include "murat/ana/TEventDisplayModule.hh"


//-----------------------------------------------------------------------------
// modules defined in 'murat' package (libmurat_ana.so)
//-----------------------------------------------------------------------------
TCalAnaModule*           m_cal   = NULL;
TCosmicsAnaModule*       m_cos   = NULL;
TDioCalibModule*         m_dio   = NULL;
// TGenAnaModule*           m_gen   = NULL;
TEventDisplayModule*     m_evd   = NULL;
TPidAnaModule*           m_pid   = NULL;
TStrawHitAnaModule*      m_strh  = NULL;
TTrackAnaModule*         m_trk   = NULL;
TTrackAnaModuleA*        m_trka  = NULL;
TTrackCompModule*        m_tcm   = NULL;
TTrackPidAnaModule*      m_tpa   = NULL;
TTrackRecoEffAnaModule*  m_eff   = NULL;
TValCalPatRecModule*     m_vcpr  = NULL;
TVdetAnaModule*          m_vdt   = NULL;

// TStnTrackID*             trk_id  = NULL;
// TAnalysisDataset*        a_dset  = NULL;

#endif

