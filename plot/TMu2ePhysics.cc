///////////////////////////////////////////////////////////////////////////////
// plots for Mu2e-Physics: Mu2e sensitivity studies
// -------------------------
// directory with the histogram files: by default: $WORK_DIR/results
// can be redefined with "mu2e.HistDir" in .rootrc
// 
// figures:
// ---------
// Fig.  1: expected momentum spectrum at 105 MeV - smearing
// Fig.  2: momentum spectrum of the
// fig.  3: track reco efficiency vs track momentum
//
// Fig. 100: DIO plots
//
// Fig. 200: RCP plots
// Fig. 201: e- momentum from PCP's
//
// Fig. 300: cosmics    plots
// Fig. 400: antiproton plots
// 
////////////////////////////////////////////////////////////////////////////////
#include "plot/TMu2ePhysics.hh"

#include "TPaveLabel.h"
#include "TArrow.h"
#include "TObjArray.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TObjString.h"
#include <iostream>

ClassImp(TMu2ePhysics)

//_____________________________________________________________________________
TMu2ePhysics::TMu2ePhysics(int PlotMode, int BlessingMode): TPlotNote(PlotMode,BlessingMode) {

  const char* hist_dir;
  char         res_dir[500];
  
  fFiguresDir   = Form("%s/mu2e_physics",gEnv->GetValue("mu2e.Figures","."));
  fWorkDir      = gSystem->Getenv("WORK_DIR");

  sprintf(res_dir,"%s/results",fWorkDir.Data());
  hist_dir      = gEnv->GetValue("mu2e.HistDir",res_dir);

  sprintf(f_conv  ,"%s/v4_0_6/conversion_01_tcalm002.hist", hist_dir); 
  sprintf(f_rpc   ,"%s/v4_0_6/rpc00102_tcalm002.hist", hist_dir); 

  sprintf(f_data  ,"%s/limits/data01.mu2e_limits.hist", hist_dir); 
  sprintf(f_dio   ,"%s/limits/dio02.mu2e_limits.hist" , hist_dir); 

  sprintf(res_dir,"%s",fWorkDir.Data());
  hist_dir      = gEnv->GetValue("mu2e.HistDir",res_dir);

					// integral of the crystal ball function
  fNConvEle = 1.0;
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
TMu2ePhysics::~TMu2ePhysics() {
};

//-----------------------------------------------------------------------------
const char* TMu2ePhysics::GetFiguresDir() {
  return fFiguresDir.Data();
}

void TMu2ePhysics::get_filename(int Figure, char* Filename) {

  if      (Figure ==   1) sprintf(Filename,"fig_%03i_%s",Figure,"treff_vs_zv");
  else if (Figure ==   2) sprintf(Filename,"fig_%03i_%s",Figure,"treff_vs_rv");
  else if (Figure == 201) sprintf(Filename,"fig_%03i_%s",Figure,"emom_rpc");
  else if (Figure == -12) sprintf(Filename,"fig_%03i_%s",Figure,"prt_energy");
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
void TMu2ePhysics::remake_plots() {

  int fig [] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
		 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
		 21, 22, 23,
		 103, 113, 123,
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
void TMu2ePhysics::plot(Int_t Figure, const char* CanvasName) {
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
// Fig.  1: track reconstruction efficiencies (all tracks and "Set C" tracks) vs Z
//-----------------------------------------------------------------------------
  if (Figure == 1) {
    TH1F      *h_eff_1, *h_eff_6, *h_zv;
    TLegend   *leg;

    fCanvas = new TCanvas(name,"",1000,700);

    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    TPad *pad2 = new TPad("pad2","",0,0,1,1);

    pad1->Draw();
    pad1->cd();

    h_eff_1 = create_eff_hist(f_conv,"TCalm002","evt_0/zv","evt_1/zv","h_eff_1",
			      "Track reconstruction efficiency (all) vs Z(v)",20,1);
    h_eff_6 = create_eff_hist(f_conv,"TCalm002","evt_0/zv","evt_6/zv","h_eff_6",
			      "Track reconstruction efficiency Set C vs Z(v) vs Zv",24,1);

    h_eff_1->GetXaxis()->SetRangeUser(5400,6399);
    h_eff_1->GetXaxis()->SetTitle("Z, mm");
    h_eff_1->GetYaxis()->SetRangeUser(0.3,1.0);
    h_eff_1->SetStats(0);
    h_eff_1->Draw();

    h_eff_6->Draw("same");
    
   //compute the pad range with suitable margins

    h_zv  = gh1(f_conv,"TCalm002","evt_0/zv");

    Double_t ymin = 0;
    Double_t ymax = h_zv->GetMaximum()*1.1;
    Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
    Double_t xmin = 5400;
    Double_t xmax = 6400;
    Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
    pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
    
    pad2->SetFillStyle(4000);   //will be transparent

    pad2->Draw();
    pad2->cd();

    h_zv->SetLineColor(kRed);
    h_zv->SetStats(0);
    h_zv->Draw("][sames");
    
    TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
    axis->SetLabelColor(kRed);
    axis->Draw();

    leg = new TLegend(0.35,0.75,0.65,0.85,"");

    leg->AddEntry(h_eff_1,"all tracks","pe");
    leg->AddEntry(h_eff_6,"Set C tracks","pe");
    leg->AddEntry(h_zv   ,"Z(conversion)","f");

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();
    
    //    fCanvas->Draw();

  }
//-----------------------------------------------------------------------------
// Fig.  2: track reconstruction efficiencies (all tracks and "Set C" tracks) vs R
//-----------------------------------------------------------------------------
  if (Figure == 2) {
    TH1F      *h_eff_1, *h_eff_6, *h_rv;
    TLegend   *leg;

    fCanvas = new TCanvas(name,"",1000,700);

    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    TPad *pad2 = new TPad("pad2","",0,0,1,1);

    pad1->Draw();
    pad1->cd();

    h_eff_1 = create_eff_hist(f_conv,"TCalm002","evt_0/rv","evt_1/rv","h_eff_1",
			      "Track reconstruction efficiency (all) vs R(v)",20,1);
    h_eff_6 = create_eff_hist(f_conv,"TCalm002","evt_0/rv","evt_6/rv","h_eff_6",
			      "Track reconstruction efficiency Set C vs R(v) vs Zv",24,1);

    h_eff_1->GetXaxis()->SetRangeUser(0,99.9);
    h_eff_1->GetXaxis()->SetTitle("R, mm");
    h_eff_1->GetYaxis()->SetRangeUser(0.3,1.0);
    h_eff_1->SetStats(0);
    h_eff_1->Draw();

    h_eff_6->Draw("same");
    
   //compute the pad range with suitable margins

    h_rv  = gh1(f_conv,"TCalm002","evt_0/rv");

    Double_t ymin = 0;
    Double_t ymax = h_rv->GetMaximum()*1.1;
    Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
    Double_t xmin = 0;
    Double_t xmax = 99.9;
    Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
    pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
    
    pad2->SetFillStyle(4000);   //will be transparent

    pad2->Draw();
    pad2->cd();

    h_rv->SetLineColor(kRed);
    h_rv->SetStats(0);
    h_rv->Draw("][sames");
    
    TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
    axis->SetLabelColor(kRed);
    axis->Draw();

    leg = new TLegend(0.65,0.75,0.90,0.85,"");

    leg->AddEntry(h_eff_1,"all tracks","pe");
    leg->AddEntry(h_eff_6,"Set C tracks","pe");
    leg->AddEntry(h_rv   ,"R(conversion)","f");

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();
    
    //    fCanvas->Draw();

  }
//-----------------------------------------------------------------------------
// Fig.  201: RPC momentum spectrum
//-----------------------------------------------------------------------------
  if (Figure == 201) {
    TH1D   *h1;
    fCanvas  = new_slide(name,title,1,1,1100,800);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    h1 = (TH1D*) gh1(f_rpc ,"TCalm002","trk_1/p2");

    h1->SetTitle("electrons from RPC, v4_0_6"); 
    h1->GetXaxis()->SetTitle("p, MeV/c"); 
    h1->Draw();
  }
}

