///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/scripts/global_vars.h"
#include "modules.hh"

def_name murat_det_001("murat_calo_time_ana");
def_name murat_det_002("murat_crv_tc_ana");
def_name murat_det_003("murat_crv_time_ana");
def_name murat_det_004("murat_det_time_ana");
def_name murat_det_005("murat_panel_dt_ana");
def_name murat_det_006("murat_event_dump"  );

//-----------------------------------------------------------------------------
void  murat_calo_time_ana(int DebugBit = -1) {
  murat::m_calt = (murat::TCaloTimeAnaModule*) g.x->AddModule("murat::TCaloTimeAnaModule",0);  

  if (DebugBit >= 0) murat::m_calt->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  murat_det_time_ana(int DebugBit = -1) {
  murat::m_det = (murat::TDetTimeAnaModule*) g.x->AddModule("murat::TDetTimeAnaModule",0);  

  if (DebugBit >= 0) murat::m_det->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  murat_crv_tc_ana(int DebugBit = -1) {
  murat::m_ctc = (murat::TCrvTcAnaModule*) g.x->AddModule("murat::TCrvTcAnaModule",0);  

  if (DebugBit >= 0) murat::m_ctc->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  murat_crv_time_ana(int DebugBit = -1) {
  murat::m_crvt = (murat::TCrvTimeAnaModule*) g.x->AddModule("murat::TCrvTimeAnaModule",0);  

  if (DebugBit >= 0) murat::m_crvt->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  murat_panel_dt_ana(int DebugBit = -1) {
  murat::m_pdt = (murat::TPanelDtAnaModule*) g.x->AddModule("murat::TPanelDtAnaModule",0);  

  if (DebugBit >= 0) murat::m_pdt->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  murat_event_dump(int DebugBit = -1) {
  murat::m_dump = (murat::TEventDumpModule*) g.x->AddModule("murat::TEventDumpModule",0);  

  if (DebugBit >= 0) murat::m_dump->SetDebugBit(DebugBit,1);
}
