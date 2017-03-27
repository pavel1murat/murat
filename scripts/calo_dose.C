///////////////////////////////////////////////////////////////////////////////
// estimate radiation dose due to different processes in the Mu2e calorimeter 
// ProcessName = "FLASH", "DIO" etc ...for example
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/TCaloDose.hh"

TStnAna*      x(NULL);
TCaloDose*    mcd(NULL); 

int calo_dose(const char* ProcessName, int NEvents = -1, int SaveHistograms = 0) {

  x = new TStnAna("none","gen");

  mcd = new TCaloDose(ProcessName);

  x->AddModule(mcd,0);

  x->GetFolder()->Add(mcd->GetFolder());

  mcd->Loop(NEvents);

  if (SaveHistograms) x->SaveHist(Form("calo_dose_%s.hist",ProcessName));

  return 0;
}
