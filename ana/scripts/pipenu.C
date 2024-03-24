///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name murat_pipenu_0010 ("murat_pipenu_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void  murat_pipenu_ana(const char* TrackBlockName = "TrackBlock", int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  murat::m_pipenu = (murat::TPipenuAnaModule*) g.x->AddModule("murat::TPipenuAnaModule",0);  

  murat::m_pipenu->SetTrackBlockName(TrackBlockName);
  
  printf("DebugBit:%i\n",DebugBit);
  if (DebugBit >= 0) murat::m_pipenu->SetDebugBit(DebugBit,1);
}
