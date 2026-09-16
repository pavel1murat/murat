///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/scripts/global_vars.h"
#include "modules.hh"

def_name murat_det_002("murat_det_time_ana");

//-----------------------------------------------------------------------------
void  murat_det_time_ana(int DebugBit = -1) {
  murat::m_det = (murat::TDetTimeAnaModule*) g.x->AddModule("murat::TDetTimeAnaModule",0);  

  if (DebugBit >= 0) murat::m_det->SetDebugBit(DebugBit,1);
}
