///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/scripts/modules.hh"

def_name murat_momscale_0010 ("murat_momscale_ana");
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void  murat_momscale_ana(int DebugBit = -1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  murat::m_moms = (murat::TMomscaleAnaModule*) g.x->AddModule("murat::TMomscaleAnaModule",0);  

  //  murat::m_moms->SetSignalParticle(-11,14);

  printf("DebugBit:%i\n",DebugBit);
  if (DebugBit >= 0) murat::m_moms->SetDebugBit(DebugBit,1);
}

