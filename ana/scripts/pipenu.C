///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name murat_pipenu_0010 ("murat_pipenu_ana");
def_name murat_pipenu_0020 ("murat_strip_tracks");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void  murat_pipenu_ana(int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  murat::m_pipenu = (murat::TPipenuAnaModule*) g.x->AddModule("murat::TPipenuAnaModule",0);  

  printf("DebugBit:%i\n",DebugBit);
  if (DebugBit >= 0) murat::m_pipenu->SetDebugBit(DebugBit,1);
}

//-----------------------------------------------------------------------------
void  murat_strip_tracks(const char* Project, const char* ODsid) {
//-----------------------------------------------------------------------------
// configure analysis module *g.x" is the TStnAna thing
//-----------------------------------------------------------------------------
  murat::m_filter = (murat::TFilterModule*) g.x->AddModule("murat::TFilterModule",0);  
  murat::m_filter->SetFilteringMode(1);   // 0: disabled (default), 1: filter, 2: veto

  TStnOutputModule* m = new TStnOutputModule(Form("nts.mu2e.%s.%s.000000_00000000.stn",ODsid,Project));
  g.x->SetOutputModule(m);
}
