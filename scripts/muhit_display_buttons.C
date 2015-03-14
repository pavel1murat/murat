//-----------------------------------------------------------------------------
// 2014-06-05 P.Murat: a poor-man's "gui" for MuHitDisplay
//-----------------------------------------------------------------------------


void PrintCalPatRec() {
  TAnaDump* d = TAnaDump::Instance();
  printf("\n");
  d->printKalRepCollection("CalPatRec");
}

void PrintTrkPatRec() {
  TAnaDump* d = TAnaDump::Instance();
  printf("\n");
  d->printKalRepCollection("TrkPatRec");
}


void PrintCalTimePeakCollection() {
  TAnaDump* d = TAnaDump::Instance();
  printf("\n");
  d->printCalTimePeakCollection("CalPatRec");
}

void PrintStrawHits() {
  TAnaDump* d = TAnaDump::Instance();
  printf("\n");
  d->printStrawHitCollection("makeSH"); 
}

//-----------------------------------------------------------------------------
void muhit_display_buttons() {
   // This macro shows a control bar o run some of the ROOT tutorials.
   // To execute an item, click with the left mouse button.

  //   gROOT->Reset();

   //Add the tutorials directory to the macro path
   //This is necessary in case this macro is executed from another user directory
   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("demos.C","");
   dir.ReplaceAll("/./","");
   const char *current = gROOT->GetMacroPath();
   gROOT->SetMacroPath(Form("%s:%s",current,dir.Data()));
   
   TControlBar *bar = new TControlBar("vertical", "Demos",10,10);
   bar->AddButton("browser",   "new TBrowser;",         "Start the ROOT Browser");
   bar->AddButton("Print CalPatRec", "PrintCalPatRec()","print CalPatRec tracks");
   bar->AddButton("Print TrkPatRec", "PrintTrkPatRec()","print TrkPatRec tracks");
   bar->AddButton("Print CalTimePeaks","PrintCalTimePeakCollection()","print Cal Time Peaks");
   bar->AddButton("Print Straw Hits","PrintStrawHits(); > aaa.txt", "print straw hits");
   bar->AddButton("Next event",".q", "next event");
   bar->SetButtonWidth(90);
   bar->Show();
   //   gROOT->SaveContext();
}
