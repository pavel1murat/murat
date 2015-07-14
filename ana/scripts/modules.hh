#ifndef murat_ana_scripts_modules_hh
#define murat_ana_scripts_modules_hh
//-----------------------------------------------------------------------------
// modules defined in 'murat' package (libmurat_ana.so)
//-----------------------------------------------------------------------------
TCalAnaModule*           m_cal   = NULL;
TCosmicsAnaModule*       m_cos   = NULL;
TDioCalibModule*         m_dio   = NULL;
// TStrawHitAnaModule*      m_str   = NULL;
TTrackAnaModule*         m_trk   = NULL;
TTrackCompModule*        m_tcm   = NULL;
TTrackPidAnaModule*      m_pid   = NULL;
TTrackRecoEffAnaModule*  m_eff   = NULL;
TValCalPatRecModule*     m_vcpr  = NULL;
TVdetAnaModule*          m_vdt   = NULL;

TStnTrackID*             trk_id  = NULL;

TAnalysisDataset*        a_dset  = NULL;
#endif

