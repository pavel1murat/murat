///////////////////////////////////////////////////////////////////////////////
// plots for pet-0001
// -------------------------
// directory with the histogram files: by default: $WORK_DIR/results
// can be redefined with "pet.HistDir" in .rootrc
// 
// figures:
// ---------
// Fig.   1: distributions in R and Z
// Fig.   2: point source : distribution in E_max and E_tot
// Fig.   3: point source : E2:E1, crystal and wedge
// Fig.   4: point source : W2:W1, crystal and wedge
//
// Fig. 302: dist source, vacuum: number of hit crystals and wedges in "trigger" events
// Fig. 304: dist source, vacuum: W2:W1, crystal and wedge for distributed source

// Fig. 402: dist source : number of hit crystals and wedges in "trigger" events
// Fig. 404: point source, water: W2:W1, crystal and wedge 
//
// Fig. 504: point source : W2:W1, crystal and wedge for distributed source, water
//
// Fig. 600-699: geometry: 32x5mm, distributed source, water , single event
// Fig. 602: distributed source, water: number of hit crystals and wedges in "trigger" events
// Fig. 605: 32x05 imager, water phantom, single event: dt_min for 2 trigger wedges
//
// Fig. 700-799: geometry: 32x5mm, distributed source, water, dose=1mCi
//
// Fig. 704: 32x05 imager, water phantom, 1 mCi: W2:W1, crystal and wedge
// Fig. 705: 32x05 imager, water phantom, 1 mCi: dt_min for 2 trigger wedges
// Fig. 706: 32x05 imager, water phantom, 1 mCi: evt_4/W2:W1, crystal and wedge
// Fig. 707: 32x05 imager, water phantom, 1 mCi: evt_4/dist (distribution for scatterers)
// Fig. 708: 32x05 imager, water phantom, 1 mCi: coincidence type
//
// Fig. 800-899: geometry: 32x5x5x15mm3 LYSO, distributed water phantom, dose 10mCi
// Fig. 807: 32x05 imager, water phantom, 10 mCi: evt_4/dist (distribution for scatterers)
// Fig. 808: 32x05 imager, water phantom, 10 mCi: coincidence type
//
// Fig. 805:     dt_min (10mCi vs 1mCi)
//
// Fig. 911: energies evt_2 BaF2 vs LYSO
////////////////////////////////////////////////////////////////////////////////
#include "plot/TPet0001.hh"

#include "TPaveLabel.h"
#include "TArrow.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TProfile.h"
#include <iostream>

ClassImp(TPet0001)

//_____________________________________________________________________________
TPet0001::TPet0001(int PlotMode, int BlessingMode): TPlotNote(PlotMode,BlessingMode) {

  const char* hist_dir;
  char         res_dir[500];
  
  fFiguresDir   = Form("%s/pet_0001",gEnv->GetValue("pet.Figures","."));
  fWorkDir      = gSystem->Getenv("WORK_DIR");

  sprintf(res_dir,"%s/results",fWorkDir.Data());
  hist_dir      = gEnv->GetValue("pet.HistDir",res_dir);

  //  fPet002   = Form("%s/v4_0_3/petana002.hist.2013-11-06",hist_dir); 
  fPet002   = Form("%s/v4_0_3/petana002.hist",hist_dir); 
  //  fPet003   = Form("%s/v4_0_3/petana003.hist.2013-11-06",hist_dir); 
  fPet003   = Form("%s/v4_0_3/petana003.hist",hist_dir); 
  //  fPet004   = Form("%s/v4_0_3/petana004.hist.2013-11-06",hist_dir); 
  fPet004   = Form("%s/v4_0_3/petana004.hist",hist_dir); 
  //  fPet005   = Form("%s/v4_0_3/petana005.hist.2013-11-06",hist_dir); 
  fPet005   = Form("%s/v4_0_3/petana005.hist",hist_dir); 
  //  fPet006   = Form("%s/v4_0_3/petana006.hist.2013-11-06",hist_dir); 
  fPet006   = Form("%s/v4_0_3/petana006.hist",hist_dir); 
  //  fPet007   = Form("%s/v4_0_3/petana007.hist.2013-11-06",hist_dir); 
  fPet007   = Form("%s/v4_0_3/petana007.hist",hist_dir); 
  //  fPet008   = Form("%s/v4_0_3/petana008.hist.2013-11-06",hist_dir); 
  fPet008   = Form("%s/v4_0_3/petana008.hist",hist_dir); 
  //  fPet009   = Form("%s/v4_0_3/petana009.hist.2013-11-06",hist_dir); 
  fPet009   = Form("%s/v4_0_3/petana009.hist",hist_dir); 
  //  fPet010   = Form("%s/v4_0_3/petana010.hist.2013-11-06",hist_dir); 
  fPet010   = Form("%s/v4_0_3/petana010.hist",hist_dir); 

  sprintf(res_dir,"%s",fWorkDir.Data());
//-----------------------------------------------------------------------------
// different settings (colors etc) for talk, note and paper modes
//-----------------------------------------------------------------------------
  // if (PlotMode == kNoteMode) {
  // }
  // else if (PlotMode == kTalkMode) {
  // }
  // else if (PlotMode == kPaperMode) {
  // }

};

//_____________________________________________________________________________
TPet0001::~TPet0001() {
};

//-----------------------------------------------------------------------------
const char* TPet0001::GetFiguresDir() {
  return fFiguresDir.Data();
}

void TPet0001::get_filename(int Figure, char* Filename) {

  if      (Figure == -1) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_dt");
  else {
//-----------------------------------------------------------------------------
// undefined, just figure
//-----------------------------------------------------------------------------
    sprintf(Filename,"fig_%i",Figure);
  }
}
//-----------------------------------------------------------------------------
// determine name of the .eps and .gif files in the ./figures directory
// today's default set is 'frr_03', all the plots should be coming from it
//-----------------------------------------------------------------------------
void TPet0001::remake_plots() {

  int fig [] = { 610, 

		 -1
  };

  int figi;

  for (int i=0; fig[i]>0; i++) {
    figi = fig[i];
    plot (figi); 
    print(figi); 
  }
}


//-----------------------------------------------------------------------------
void TPet0001::plot(Int_t Figure, const char* CanvasName) {
//-----------------------------------------------------------------------------
//  make sure all the needed scripts are loaded
//-----------------------------------------------------------------------------
  char        macro[200];

  const char* script[] = { 
//    "plot/make_fake_rate_hist.C",
//    "plot/make_ratio_hist.C",
//    "plot/plot_ntrk10.C",
    0 
  };

  const char* work_dir = gSystem->Getenv("PWD");

  for (int i=0; script[i] != 0; i++) {
    sprintf(macro,"%s/murat/%s",work_dir,script[i]);
    if (! gInterpreter->IsLoaded(macro)) gInterpreter->LoadMacro(macro);
  }
//-----------------------------------------------------------------------------
  char        name[200], title[200];

  //  int old_opt_stat = gStyle->GetOptStat();
  
  if ((! CanvasName) || (CanvasName[0] == '0')) {
    sprintf(name,"fig_%i",Figure);
    sprintf(title,"figure_%i",Figure);
  }
  else {
    sprintf(name,"%s",CanvasName);
    sprintf(title,"figure %i",Figure);
  }
//-----------------------------------------------------------------------------
// Fig.   1: distributions in R and Z
//-----------------------------------------------------------------------------
  if (Figure == 1) {

    TH1F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    gPad->SetLogy(1);
    h1 = gh1(fPet002.Data(),"PetAna001","pho_0/r");
    h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    h1->GetXaxis()->SetTitle("Radius, mm");
    h1->Draw();

    h2 = gh1(fPet003.Data(),"PetAna001","pho_0/r");
    h2->SetFillStyle(3002);
    h2->SetFillColor(kBlue);
    h2->Draw("same");

    fP1->cd(2);
    gPad->SetLogy(1);
    h1 = gh1(fPet002.Data(),"PetAna001","pho_0/z");
    h1->GetXaxis()->SetTitle("Z, mm");
    h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    h1->Draw();

    h2 = gh1(fPet003.Data(),"PetAna001","pho_0/z");
    h2->SetFillStyle(3002);
    h2->SetFillColor(kBlue);
    h2->Draw("same");
  }    
//-----------------------------------------------------------------------------
// Fig.   2: point source : distribution in E_max and E_tot
//-----------------------------------------------------------------------------
  if (Figure == 2) {

    TH1F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh1(fPet002.Data(),"PetAna001","evt_2/emax");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    h1->GetXaxis()->SetTitle("Crystal Energy, MeV");
    h1->Draw();

    fP1->cd(2);
    h1 = gh1(fPet002.Data(),"PetAna001","evt_2/etot");

    h1->GetXaxis()->SetTitle("Module Energy, MeV");
    h1->Draw();
  }    
//-----------------------------------------------------------------------------
// Fig.   3: point source : E2:E1, crystal and wedge
//-----------------------------------------------------------------------------
  if (Figure == 3) {

    TH2F  *h1, *h2;
    TBox  *b1, *b2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh2(fPet002.Data(),"PetAna001","evt_2/e2_vs_e1_cr");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Crystal Energy, MeV");
    h1->Draw();

    b1 = new TBox(0.42,0.42,0.6,0.6);
    b1->SetFillStyle(0);
    b1->SetLineColor(2);
    b1->SetLineWidth(2);
    b1->Draw("same");

    fP1->cd(2);
    h1 = gh2(fPet002.Data(),"PetAna001","evt_2/e2_vs_e1_wed");
    h1->GetXaxis()->SetTitle("Module Energy, MeV");
    h1->SetMarkerStyle(7);
    h1->Draw();

    b2 = new TBox(0.42,0.42,0.6,0.6);
    b2->SetFillStyle(0);
    b2->SetLineColor(2);
    b2->SetLineWidth(2);
    b2->Draw("same");
  }    
//-----------------------------------------------------------------------------
// Fig.   4: point source : W2:W1, crystal and wedge
//-----------------------------------------------------------------------------
  if (Figure == 4) {

    TH2F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh2(fPet002.Data(),"PetAna001","evt_2/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");

    fP1->cd(2);
    h1 = gh2(fPet002.Data(),"PetAna001","evt_3/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");
  }    
//-----------------------------------------------------------------------------
// Fig. 302: point source : number of hit crystals and wedges in "trigger" events
//-----------------------------------------------------------------------------
  if (Figure == 302) {

    TH1F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh1(fPet003.Data(),"PetAna001","evt_2/nhcr");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("number of hit crystals");
    h1->Draw("");

    fP1->cd(2);
    h1 = gh1(fPet003.Data(),"PetAna001","evt_2/nhmod");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("number of wedges with hit crystals");
    h1->Draw("");
  }    
//-----------------------------------------------------------------------------
// Fig. 304: point source : W2:W1, crystal and wedge
//-----------------------------------------------------------------------------
  if (Figure == 304) {

    TH2F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh2(fPet003.Data(),"PetAna001","evt_2/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");

    fP1->cd(2);
    h1 = gh2(fPet003.Data(),"PetAna001","evt_3/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");
  }    
//-----------------------------------------------------------------------------
// Fig. 402: point source, water : number of hit crystals and wedges in "trigger" events
//-----------------------------------------------------------------------------
  if (Figure == 402) {

    TH1F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh1(fPet004.Data(),"PetAna001","evt_2/nhcr");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("number of hit crystals");
    h1->Draw("");

    fP1->cd(2);
    h1 = gh1(fPet004.Data(),"PetAna001","evt_2/nhmod");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("number of wedges with hit crystals");
    h1->Draw("");
  }    
//-----------------------------------------------------------------------------
// Fig. 404: point source : W2:W1, crystal and wedge
//-----------------------------------------------------------------------------
  if (Figure == 404) {

    TH2F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh2(fPet004.Data(),"PetAna001","evt_2/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");

    fP1->cd(2);
    h1 = gh2(fPet004.Data(),"PetAna001","evt_3/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");
  }    
//-----------------------------------------------------------------------------
// Fig. 502: point source, water: number of hit crystals and wedges in "trigger" events
//-----------------------------------------------------------------------------
  if (Figure == 502) {

    TH1F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh1(fPet005.Data(),"PetAna001","evt_2/nhcr");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("number of hit crystals");
    h1->Draw("");

    fP1->cd(2);
    h1 = gh1(fPet005.Data(),"PetAna001","evt_2/nhmod");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("number of wedges with hit crystals");
    h1->Draw("");
  }    
//-----------------------------------------------------------------------------
// Fig. 504: point source : W2:W1, crystal and wedge
//-----------------------------------------------------------------------------
  if (Figure == 504) {

    TH2F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh2(fPet005.Data(),"PetAna001","evt_2/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");

    fP1->cd(2);
    h1 = gh2(fPet005.Data(),"PetAna001","evt_3/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");
  }    
//-----------------------------------------------------------------------------
// Fig. 602: distributed source, water: number of hit crystals and wedges in "trigger" events
//-----------------------------------------------------------------------------
  if (Figure == 602) {

    TH1F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh1(fPet006.Data(),"PetAna001","evt_3/nhcr");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("number of hit crystals");
    h1->Draw("");

    fP1->cd(2);
    h2 = gh1(fPet006.Data(),"PetAna001","evt_3/nhmod");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h2->GetXaxis()->SetTitle("number of wedges with hit crystals");
    h2->Draw("");
  }    
//-----------------------------------------------------------------------------
// Fig. 603: 32x05 imager, water phantom, 1 event: E2:E1, crystal and wedge
//-----------------------------------------------------------------------------
  if (Figure == 603) {

    TH2F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh2(fPet006.Data(),"PetAna001","evt_2/e2_vs_e1_cr");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Crystal Energy, MeV");
    h1->Draw();

    fP1->cd(2);
    h1 = gh2(fPet006.Data(),"PetAna001","evt_2/e2_vs_e1_wed");
    h1->GetXaxis()->SetTitle("Module Energy, MeV");
    h1->SetMarkerStyle(7);
    h1->Draw();
  }    
//-----------------------------------------------------------------------------
// Fig. 604: 32x05 imager, water phantom, 1 event : W2:W1, crystal and wedge
//-----------------------------------------------------------------------------
  if (Figure == 604) {

    TH2F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh2(fPet006.Data(),"PetAna001","evt_2/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");

    fP1->cd(2);
    h1 = gh2(fPet006.Data(),"PetAna001","evt_3/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");
  }    
//-----------------------------------------------------------------------------
// Fig. 605: 32x05 imager, water phantom, single event: dt_min for 2 trigger wedges
//-----------------------------------------------------------------------------
  if (Figure == 605) {

    TH2F  *h1;

    fCanvas  = new_slide(name,title,1,1,800,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh2(fPet006.Data(),"PetAna001","evt_3/dt_min");
    h1->GetXaxis()->SetRangeUser(-10,9.99);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("#Delta_{T}, ns");
    h1->Draw("");
  }    
//-----------------------------------------------------------------------------
// Fig. 702: 32_05_01_mci: distributed source, water: number of hit crystals and wedges in "trigger" events
//-----------------------------------------------------------------------------
  if (Figure == 702) {

    TH1F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh1(fPet007.Data(),"PetAna001","evt_3/nhcr");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("number of hit crystals");
    h1->Draw("");

    fP1->cd(2);
    h2 = gh1(fPet007.Data(),"PetAna001","evt_3/nhmod");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h2->GetXaxis()->SetTitle("number of wedges with hit crystals");
    h2->Draw("");
  }    
//-----------------------------------------------------------------------------
// Fig. 704: 32x05 imager, water phantom, 1 event : W2:W1, crystal and wedge
//-----------------------------------------------------------------------------
  if (Figure == 704) {

    TH2F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh2(fPet007.Data(),"PetAna001","evt_2/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");

    fP1->cd(2);
    h1 = gh2(fPet007.Data(),"PetAna001","evt_3/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");
  }    
//-----------------------------------------------------------------------------
// Fig. 705: 32x05 imager, water phantom, dose=1mCi: dt_min for 2 trigger wedges
//-----------------------------------------------------------------------------
  if (Figure == 705) {

    TH2F  *h1;

    fCanvas  = new_slide(name,title,1,1,800,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh2(fPet007.Data(),"PetAna001","evt_3/dt_min");
    h1->GetXaxis()->SetRangeUser(-10,9.99);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("#Delta_{T}, ns");
    h1->Draw("");
  }    
//-----------------------------------------------------------------------------
// Fig. 706: 32x05 imager, water phantom, 1 event : evt_4/W2:W1, crystal and wedge
//-----------------------------------------------------------------------------
  if (Figure == 706) {

    TH2F  *h1, *h2;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh2(fPet007.Data(),"PetAna001","evt_3/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");

    fP1->cd(2);
    h1 = gh2(fPet007.Data(),"PetAna001","evt_4/w2_vs_w1");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Wedge number");
    h1->Draw("box");
  }    
//-----------------------------------------------------------------------------
// Fig. 707: 32x05 imager, water phantom, 1 mCi: distribution for scatterers
//-----------------------------------------------------------------------------
  if (Figure == 707) {

    TH1F  *h1, *h2;

    fCanvas  = new_slide(name,title,1,1,800,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    gPad->SetLogy(1);
    h1 = gh1(fPet007.Data(),"PetAna001","evt_4/dist");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("scattered distance, 1 mCi");
    h1->Draw("");

    TArrow* arr = new TArrow(8,5,8,1,0.015);
    arr->Draw();
  }    
//-----------------------------------------------------------------------------
// Fig. 708: 32x05 imager, water phantom, 1 mCi: coincidence type
//-----------------------------------------------------------------------------
  if (Figure == 708) {

    TH1F  *h1, *h2;

    fCanvas  = new_slide(name,title,1,1,800,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    gPad->SetLogy(1);
    h1 = gh1(fPet007.Data(),"PetAna001","evt_4/ctype");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("coincidence type, 1 mCi");
    h1->Draw("");
  }    
//-----------------------------------------------------------------------------
// Fig. 805: 32x05 imager, 10mCi vs 1mCi
//-----------------------------------------------------------------------------
  if (Figure == 805) {

    TH1F  *h1, *h2;

    fCanvas  = new_slide(name,title,1,1,900,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh1(fPet007.Data(),"PetAna001","evt_3/dt_min");
    h1->GetXaxis()->SetRangeUser(-2.,1.999);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("#Delta T");
    h1->Draw("");

    h2 = gh1(fPet008.Data(),"PetAna001","evt_3/dt_min");
    h2->SetLineColor(2);
    h2->SetFillStyle(3002);
    h2->SetFillColor(2);
    h2->GetXaxis()->SetRangeUser(-2.,1.999);
    h2->Draw("sames");

    TLegend* leg;

    leg = new TLegend(0.55,0.75,0.75,0.85,"");
    leg->AddEntry(h1,"1  mCi","f");
    leg->AddEntry(h2,"10 mCi","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }    
//-----------------------------------------------------------------------------
// Fig. 807: 32x05 imager, water phantom, 10 mCi: distribution for scatterers
//-----------------------------------------------------------------------------
  if (Figure == 807) {

    TH1F  *h1, *h2;

    fCanvas  = new_slide(name,title,1,1,800,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    gPad->SetLogy(1);
    h1 = gh1(fPet008.Data(),"PetAna001","evt_4/dist");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("scattered distance, 10 mCi");
    h1->Draw("");

    TArrow* arr = new TArrow(8,5,8,1,0.015);
    arr->Draw();
  }    
//-----------------------------------------------------------------------------
// Fig. 808: 32x05 imager, water phantom, 10 mCi: coincidence type
//-----------------------------------------------------------------------------
  if (Figure == 808) {

    TH1F  *h1, *h2;

    fCanvas  = new_slide(name,title,1,1,800,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    gPad->SetLogy(1);
    h1 = gh1(fPet008.Data(),"PetAna001","evt_4/ctype");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("coincidence type, 1 mCi");
    h1->Draw("");
  }    
//-----------------------------------------------------------------------------
// Fig. 911: 32x05 imager, water phantom, evt_2/emax: BaF2 vs LYSO
//-----------------------------------------------------------------------------
  if (Figure == 911) {

    TH1F  *h1, *h2;

    fCanvas  = new_slide(name,title,1,1,800,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    h1 = gh1(fPet006.Data(),"PetAna001","evt_2/emax");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h1->GetXaxis()->SetTitle("Energy, MeV");
    h1->Draw("");

    h2 = gh1(fPet009.Data(),"PetAna001","evt_2/emax");
    // h1->GetYaxis()->SetRangeUser(0.1,5.e4);
    //    h1->SetMarkerStyle(7);
    h2->SetFillStyle(3002);
    h2->SetFillColor(2);
    h2->SetLineColor(2);
    h2->Draw("same");
  }    

}

