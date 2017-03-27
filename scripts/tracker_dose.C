///////////////////////////////////////////////////////////////////////////////
// Backgrounds are named, see murat/ana/TTrackerDose.cc. 
// background name = "FLASH", for example
// example: tracker_dose("FLASH",-1,1)
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/TTrackerDose.hh"

TStnAna*      x(NULL);
TTrackerDose* mtd(NULL); 

int tracker_dose(const char* ProcessName, int NEvents = -1, int SaveHistograms = 0) {

  x = new TStnAna("none","gen");


  mtd = new TTrackerDose(ProcessName);

  x->AddModule(mtd,0);

  x->GetFolder()->Add(mtd->GetFolder());

  mtd->Loop(NEvents);

  if (SaveHistograms) {
    x->SaveHist(Form("tracker_dose_%s.hist",ProcessName));
  }

  return 0;
}
