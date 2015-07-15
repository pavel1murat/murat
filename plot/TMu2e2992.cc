///////////////////////////////////////////////////////////////////////////////
// plots for cdf-10008
// -------------------------
// directory with the histogram files: by default: $WORK_DIR/results
// can be redefined with "zzx.HistDir" in .rootrc
// 
// figures:
// ---------
// Fig.   1:  delta(T) for Set C tracks with clusters
// Fig.   2:  E/P      for Set C tracks with clusters 
// Fig.   3:  DY       for Set C tracks with clusters (vane geometry)
// Fig.   4:  DZ       for Set C tracks with clusters (vane geometry)
// Fig.   5:  D0 for conversion electrons and cosmic muons
// Fig.   6:  Z0 for conversion electrons and cosmic muons
// Fig.  10:  trk_1/Dt  for electrons and muons
// Fig.  11:  trk_1/Ep  for electrons and muons
//
// Fig.  30:  trk_21/d0 for cb000401
// Fig.  31  (a) trk_24/d0 for cb000401
//           (b) trk_22/Ep
//
// Fig.  34  (a) trk_25/p for cb000401
// Fig.  35  (a) trk_30/lhr_slope for cb000401  - no cluster
//           (b) trk_31/p                       - no cluster
//           (c) trk_27/lhr_slope for cb000401  - no cluster
//           (d) trk_28/p                       - no cluster
//
// Fig.  41  (a) trk_1/xs        for e0000001 and m0000001
// Fig.  42  (a) trk_1/llhr_dedx for e0000001 and m0000001
// Fig.  44: (a) trk_1/llhr_trk cut efficiency for e0000001 and m0000001
//
// Fig.  60: track reconstruction efficiencies vs luminosity
//
// Fig.  62: LogLH_CAL selection efficiency vs luminosity
// Fig.  63: LogLH_TRK selection efficiency vs luminosity
////////////////////////////////////////////////////////////////////////////////
#include "plot/TMu2e2992.hh"

#include "TPaveLabel.h"
#include "TArrow.h"
#include "TObjArray.h"
#include "TGraphErrors.h"
#include "TObjString.h"
#include <iostream>

ClassImp(TMu2e2992)

//_____________________________________________________________________________
TMu2e2992::TMu2e2992(int PlotMode, int BlessingMode): TPlotNote(PlotMode,BlessingMode) {

  const char* hist_dir;
  char         res_dir[500];
  
  fFiguresDir   = Form("%s/mu2e_2992",gEnv->GetValue("mu2e.Figures","."));
  fWorkDir      = gSystem->Getenv("WORK_DIR");

  sprintf(res_dir,"%s/results",fWorkDir.Data());
  hist_dir      = gEnv->GetValue("mu2e.HistDir",res_dir);

  e0000001_tcalm002 = Form("%s/v3_0_0/e0000001_tcalm002_vane.hist", hist_dir); 
  cb000101_tcalm002 = Form("%s/v2_1_0/cb000101.tcalm002.hist", hist_dir); 
  m0000001_tcalm002 = Form("%s/v3_0_0/m0000001_tcalm002_vane.hist", hist_dir); 
  me000001_tcalm002 = Form("%s/v2_1_0/me000001.tcalm002.hist", hist_dir); 
  cb000401_tcalm004 = Form("%s/v3_0_0/cb000401_tcalm004.hist", hist_dir); 

  // electron and muon guns with different (0,x1,x2,x4) background mixes

  egun0001_tcalm002 = Form("%s/v3_0_0/electronGun_600_1700_tcalm002_vane.hist", hist_dir); 
  egun0101_tcalm002 = Form("%s/v3_0_0/electronGun_600_1700_mixing_01_vane_x1_tcalm002.hist", hist_dir); 
  egun0201_tcalm002 = Form("%s/v3_0_0/electronGun_600_1700_mixing_01_vane_x2_tcalm002.hist", hist_dir); 
  egun0401_tcalm002 = Form("%s/v3_0_0/electronGun_600_1700_mixing_01_vane_x4_tcalm002.hist", hist_dir); 

  mgun0001_tcalm002 = Form("%s/v3_0_0/muonGun_600_1700_tcalm002_vane.hist", hist_dir); 
  mgun0101_tcalm002 = Form("%s/v3_0_0/muonGun_600_1700_mixing_01_vane_x1_tcalm002.hist", hist_dir); 
  mgun0201_tcalm002 = Form("%s/v3_0_0/muonGun_600_1700_mixing_01_vane_x2_tcalm002.hist", hist_dir); 
  mgun0401_tcalm002 = Form("%s/v3_0_0/muonGun_600_1700_mixing_01_vane_x4_tcalm002.hist", hist_dir); 

  sprintf(res_dir,"%s",fWorkDir.Data());
  hist_dir      = gEnv->GetValue("mu2e.HistDir",res_dir);
  cb000301_cosmicAnalyzer = Form("%s/cb000301.cosmicAnalyzer.hist", hist_dir); 
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
TMu2e2992::~TMu2e2992() {
};

//-----------------------------------------------------------------------------
const char* TMu2e2992::GetFiguresDir() {
  return fFiguresDir.Data();
}

void TMu2e2992::get_filename(int Figure, char* Filename) {

  if      (Figure == 1) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_dt");
  else if (Figure == 2) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_ep");
  else if (Figure == 3) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_dy");
  else if (Figure == 4) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_dz");
  else if (Figure == 5) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_d0");
  else if (Figure == 6) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_z0");
  else if (Figure == 100) sprintf(Filename,"fig_%03i_%s",Figure,"EoverP");
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
void TMu2e2992::remake_plots() {

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
void TMu2e2992::plot(Int_t Figure, const char* CanvasName) {
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
//  Fig.1: delta(T) distributions for electrons and muons, vane geometry
//-----------------------------------------------------------------------------
  if (Figure == 1) {
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(e0000001_tcalm002,"TCalm002", "trk_1/dt");
    fH1[1] = gh1(m0000001_tcalm002,"TCalm002", "trk_1/dt");
    
    fP1 ->cd(1);

    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("Tracks passing Set C cuts");
    fH1[0]->GetXaxis()->SetTitle("T_{track}-T_{cluster}, ns");
    fH1[0]->Draw();

    fH1[1]->SetFillColor(kBlue+3);
    fH1[1]->SetFillStyle(3001);
    fH1[1]->Draw("same");

    TLegend* leg = new TLegend(0.2,0.7,0.4,0.8,"");
    leg->AddEntry(fH1[0],"electrons","f");
    leg->AddEntry(fH1[1],"muons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }
//-----------------------------------------------------------------------------
//  Fig.2: E/P distributions for MC electrons and muons, vane geometry may enter
//         through clustering
//-----------------------------------------------------------------------------
  if (Figure == 2) {
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(e0000001_tcalm002,"TCalm002", "trk_1/ep");
    fH1[1] = gh1(m0000001_tcalm002,"TCalm002", "trk_1/ep");
    
    fP1 ->cd(1);

    //    fH1[0]->SetOptStat(0);

    double h_m0 = fH1[0]->GetMaximum();
    double h_m1 = fH1[1]->GetMaximum();

    if (h_m0 < h_m1) fH1[0]->SetMaximum(h_m1*1.1);

    fH1[0]->SetTitle("Tracks passing Set C cuts");
    fH1[0]->GetXaxis()->SetTitle("E_{cluster}/P_{track}");
    fH1[0]->GetXaxis()->SetRangeUser(0.,1.2);
    fH1[0]->Draw();

    fH1[1]->SetFillColor(kBlue+3);
    fH1[1]->SetFillStyle(3001);
    fH1[1]->Draw("same");

    fH1[0]->Draw("same");

    TLegend* leg = new TLegend(0.15,0.7,0.25,0.8,"");
    leg->AddEntry(fH1[0],"electrons","f");
    leg->AddEntry(fH1[1],"muons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }
//-----------------------------------------------------------------------------
//  Fig.3: DY = Y_track - Y_cluster for e0000001 trk_1/dy
//-----------------------------------------------------------------------------
  if (Figure == 3) {
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(e0000001_tcalm002,"TCalm002", "trk_1/dy");
    
    fP1 ->cd(1);

    //    fH1[0]->SetOptStat(0);

    fH1[0]->SetTitle("Tracks passing Set C cuts");
    fH1[0]->GetXaxis()->SetTitle("Y_{track}-Y_{cluster}, mm");
    fH1[0]->GetXaxis()->SetRangeUser(-200.,200);
    fH1[0]->Fit("gaus");
  }
//-----------------------------------------------------------------------------
//  Fig.4: DZ = Z_track - Z_cluster, vane geometry
//-----------------------------------------------------------------------------
  if (Figure == 4) {
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(e0000001_tcalm002,"TCalm002", "trk_1/dz");
    
    fP1 ->cd(1);

    //    fH1[0]->SetOptStat(0);

    fH1[0]->SetTitle("Tracks passing Set C cuts");
    fH1[0]->GetXaxis()->SetTitle("Z_{track}-Z_{cluster}, mm");
    fH1[0]->GetXaxis()->SetRangeUser(-200.,200);
    fH1[0]->Fit("gaus");
  }
//-----------------------------------------------------------------------------
//  Fig.5: distributions in |D0| cosmic MC and conversion electrons
//-----------------------------------------------------------------------------
  if (Figure == 5) {
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(me000001_tcalm002,"TCalm002", "trk_1/d0");
    fH1[1] = gh1(cb000101_tcalm002,"TCalm002", "trk_1/d0");
    
    fP1->cd(1);
    gPad->SetLogy(1);

    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("Tracks passing Set C cuts");
    fH1[0]->GetXaxis()->SetTitle("Track D0, mm");
    fH1[0]->Draw();

    fH1[1]->SetLineColor(kBlue-7);
    fH1[1]->SetFillColor(kBlue-7);
    fH1[1]->SetFillStyle(3003);
    fH1[1]->Draw("same");

    fH1[0]->Draw("same");

    TArrow* arr;

    arr = new TArrow(-120.,4,-120.,0.7,0.015);
    arr->SetLineWidth(3);
    arr->Draw();
    arr = new TArrow( 120.,4, 120.,0.7,0.015);
    arr->SetLineWidth(3);
    arr->Draw();

    TLegend* leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"conversion electrons","f");
    leg->AddEntry(fH1[1],"tracks from cosmics" ,"f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }
//-----------------------------------------------------------------------------
//  Fig.6: distributions in |Z0| cosmic MC and conversion electrons
//-----------------------------------------------------------------------------
  if (Figure == 6) {
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(e0000001_tcalm002,"TCalm002", "trk_1/z0");
    fH1[1] = gh1(cb000101_tcalm002,"TCalm002", "trk_1/z0");
    
    fP1->cd(1);
    gPad->SetLogy(1);

    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("Tracks passing Set C cuts");
    fH1[0]->GetXaxis()->SetTitle("Track Z0, mm");
    fH1[0]->GetXaxis()->SetRangeUser(-2500,2500);
    fH1[0]->SetMinimum(0.5);
    fH1[0]->Draw();

    fH1[1]->SetLineColor(kBlue-7);
    fH1[1]->SetFillColor(kBlue-7);
    fH1[1]->SetFillStyle(3003);
    fH1[1]->Draw("same");

    fH1[0]->Draw("same");

    TArrow* arr;

    arr = new TArrow(-1000.,4,-1000.,0.7,0.015);
    arr->SetLineWidth(3);
    arr->Draw();
    arr = new TArrow( 1000.,4, 1000.,0.7,0.015);
    arr->SetLineWidth(3);
    arr->Draw();

    TLegend* leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"conversion electrons","f");
    leg->AddEntry(fH1[1],"tracks from cosmics" ,"f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }
//-----------------------------------------------------------------------------
//  Fig. 10: distributions in Delta(T) for 100 MeV electrons and muons
//-----------------------------------------------------------------------------
  if (Figure == 10) {
    TPaveLabel*  label;
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    fH1[0] = gh1(e0000001_tcalm002,"TCalm002", "trk_1/dt");
    fH1[1] = gh1(m0000001_tcalm002,"TCalm002", "trk_1/dt");
    
    //    fH1[0]->SetOptStat(0);
    //    fH1[0]->SetTitle("Tracks passing Set C cuts");
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("Track #Delta_{T}, ns");
    fH1[0]->Draw();

    fH1[1]->SetLineColor(kBlue-7);
    fH1[1]->SetFillColor(kBlue-7);
    fH1[1]->SetFillStyle(3003);
    fH1[1]->Draw("same");

    fH1[0]->Draw("same");

    // TArrow* arr;

    // arr = new TArrow(-1000.,4,-1000.,0.7,0.015);
    // arr->SetLineWidth(3);
    // arr->Draw();
    // arr = new TArrow( 1000.,4, 1000.,0.7,0.015);
    // arr->SetLineWidth(3);
    // arr->Draw();

    TLegend* leg = new TLegend(0.15,0.73,0.40,0.83,"");
    leg->AddEntry(fH1[0],"100 MeV electrons","f");
    leg->AddEntry(fH1[1],"100 MeV muons" ,"f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);
  }

//-----------------------------------------------------------------------------
//  Fig. 11: distributions in E/P for 100 MeV electrons and muons
//-----------------------------------------------------------------------------
  if (Figure == 11) {
    TPaveLabel*  label;

    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    fH1[0] = gh1(e0000001_tcalm002,"TCalm002", "trk_1/ep");
    fH1[1] = gh1(m0000001_tcalm002,"TCalm002", "trk_1/ep");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("E_{cluster}/P_{track}");
    fH1[0]->Draw();

    fH1[1]->SetLineColor(kBlue-7);
    fH1[1]->SetFillColor(kBlue-7);
    fH1[1]->SetFillStyle(3003);
    fH1[1]->Draw("same");

    fH1[0]->Draw("same");

    // TArrow* arr;

    // arr = new TArrow(-1000.,4,-1000.,0.7,0.015);
    // arr->SetLineWidth(3);
    // arr->Draw();
    // arr = new TArrow( 1000.,4, 1000.,0.7,0.015);
    // arr->SetLineWidth(3);
    // arr->Draw();

    TLegend* leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"100 MeV electrons","f");
    leg->AddEntry(fH1[1],"100 MeV muons" ,"f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);
  }
//-----------------------------------------------------------------------------
//  Fig. 12: log LHR distributions for electrons and muons
//-----------------------------------------------------------------------------
  if (Figure == 12) {
    TPaveLabel*  label;

    fCanvas  = new_slide(name,title,2,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    fH1[0] = gh1(e0000001_tcalm002,"TCalm002", "trk_1/lhr");
    fH1[1] = gh1(m0000001_tcalm002,"TCalm002", "trk_1/lhr");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("log(LHR)");
    fH1[0]->GetXaxis()->SetRangeUser(-100.,70.);
    fH1[0]->Draw();

    fH1[1]->SetLineColor(kBlue-7);
    fH1[1]->SetFillColor(kBlue-7);
    fH1[1]->SetFillStyle(3003);
    fH1[1]->Draw("same");

    fH1[0]->Draw("same");

    // TArrow* arr;

    // arr = new TArrow(-1000.,4,-1000.,0.7,0.015);
    // arr->SetLineWidth(3);
    // arr->Draw();
    // arr = new TArrow( 1000.,4, 1000.,0.7,0.015);
    // arr->SetLineWidth(3);
    // arr->Draw();

    TLegend* leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"100 MeV electrons","f");
    leg->AddEntry(fH1[1],"100 MeV muons" ,"f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);
  }
//-----------------------------------------------------------------------------
// Fig.  30: (a) trk_21/d0 for cb000401
//           (b) trk_22/Ep
//-----------------------------------------------------------------------------
  if (Figure == 30) {
    TPaveLabel*  label;
    TLegend*     leg;

    fCanvas  = new_slide(name,title,2,1,1200,500);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    fH1[0] = gh1(cb000401_tcalm004,"TCalm004", "trk_21/d0");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("D0, mm");
    //    fH1[0]->GetXaxis()->SetRangeUser(-15.,10.);
    fH1[0]->Draw();

    // TArrow* arr;

    // arr = new TArrow(-1000.,4,-1000.,0.7,0.015);
    // arr->SetLineWidth(3);
    // arr->Draw();
    // arr = new TArrow( 1000.,4, 1000.,0.7,0.015);
    // arr->SetLineWidth(3);
    // arr->Draw();

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"tracks from cosmics","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);


    fP1->cd(2);
    fH1[0] = gh1(cb000401_tcalm004,"TCalm004", "trk_22/ep");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("E/P");
    //    fH1[0]->GetXaxis()->SetRangeUser(-15.,10.);
    fH1[0]->Draw();

    // TArrow* arr;

    // arr = new TArrow(-1000.,4,-1000.,0.7,0.015);
    // arr->SetLineWidth(3);
    // arr->Draw();
    // arr = new TArrow( 1000.,4, 1000.,0.7,0.015);
    // arr->SetLineWidth(3);
    // arr->Draw();

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"tracks from cosmics","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);
  }
//-----------------------------------------------------------------------------
// Fig.  31  (a) trk_24/d0 for cb000401
//           (b) trk_22/Ep
//-----------------------------------------------------------------------------
  if (Figure == 31) {
    TPaveLabel*  label;
    TLegend*     leg;

    fCanvas  = new_slide(name,title,2,1,1200,500);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    fH1[0] = gh1(cb000401_tcalm004,"TCalm004", "trk_23/dy");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("DY, mm");
    //    fH1[0]->GetXaxis()->SetRangeUser(-15.,10.);
    fH1[0]->Draw();

    TArrow* arr;

    arr = new TArrow(-150.,50,-150.,5,0.015);
    arr->SetLineWidth(2);
    arr->Draw();
    arr = new TArrow( 150.,50, 150.,5,0.015);
    arr->SetLineWidth(2);
    arr->Draw();

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"tracks from cosmics","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);


    fP1->cd(2);
    fH1[0] = gh1(cb000401_tcalm004,"TCalm004", "trk_23/dz");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("#Delta Z, mm");
    //    fH1[0]->GetXaxis()->SetRangeUser(-15.,10.);
    fH1[0]->Draw();

    arr = new TArrow(-150.,50,-150.,5.,0.015);
    arr->SetLineWidth(2);
    arr->Draw();
    arr = new TArrow( 150.,50, 150.,5.,0.015);
    arr->SetLineWidth(2);
    arr->Draw();

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"tracks from cosmics","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);
  }
//-----------------------------------------------------------------------------
// Fig.  32  (a) trk_24/lhr for cb000401
//-----------------------------------------------------------------------------
  if (Figure == 32) {
    TPaveLabel*  label;
    TLegend*     leg;

    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    fH1[0] = gh1(cb000401_tcalm004,"TCalm004", "trk_24/lhr");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("Log LH_{e} - Log LH_{#mu}");
    fH1[0]->GetXaxis()->SetRangeUser(-100.,100.);
    fH1[0]->Draw();

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"cosmics background, MC","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);
  }
//-----------------------------------------------------------------------------
// Fig.  33  (a) trk_25/pdg for cb000401
//           (b) trk_26/pdg
//-----------------------------------------------------------------------------
  if (Figure == 33) {
    TPaveLabel*  label;
    TLegend*     leg;

    fCanvas  = new_slide(name,title,2,1,1200,500);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    gPad->SetLogy(1);
    fH1[0] = gh1(cb000401_tcalm004,"TCalm004", "trk_25/pdg");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("PDG code");
    fH1[0]->GetXaxis()->SetRangeUser(-25.,25.);
    fH1[0]->Draw();

    //    TArrow* arr;

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"identified electrons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);

    fP1->cd(2);
    gPad->SetLogy(1);
    fH1[0] = gh1(cb000401_tcalm004,"TCalm004", "trk_26/pdg");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("PDG code");
    fH1[0]->GetXaxis()->SetRangeUser(-25.,25.);
    fH1[0]->Draw();

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"identified muons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);
  }
//-----------------------------------------------------------------------------
// Fig.  34  (a) trk_25/p for cb000401
//-----------------------------------------------------------------------------
  if (Figure == 34) {
    TPaveLabel*  label;
    TLegend*     leg;

    fCanvas  = new_slide(name,title,2,1,1200,500);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    gPad->SetLogy(1);
    fH1[0] = gh1(cb000401_tcalm004,"TCalm004", "trk_25/p");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("track momentum,  MeV/c");
    //    fH1[0]->GetXaxis()->SetRangeUser(-25.,25.);
    fH1[0]->Draw();

    //    TArrow* arr;

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"identified electrons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);
  }
//-----------------------------------------------------------------------------
// Fig.  35  (a) trk_30/lhr_slope for cb000401  - no cluster
//           (b) trk_31/p                       - no cluster
//           (c) trk_27/lhr_slope for cb000401  - no cluster
//           (d) trk_28/p                       - no cluster
//-----------------------------------------------------------------------------
  if (Figure == 35) {
    //    TArrow*      arr;
    TPaveLabel*  label;
    TLegend*     leg;

    fCanvas  = new_slide(name,title,2,2,1200,1000);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    fH1[0] = gh1(cb000401_tcalm004,"TCalm004", "trk_30/llhr_slope");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("LHR(slope)");
    //    fH1[0]->GetXaxis()->SetRangeUser(-25.,25.);
    fH1[0]->Draw();

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"identified electrons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);

    fP1->cd(2);
    fH1[0] = gh1(cb000401_tcalm004,"TCalm004", "trk_31/p");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("track momentum,  MeV/c");
    //    fH1[0]->GetXaxis()->SetRangeUser(-25.,25.);
    fH1[0]->Draw();

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"identified electrons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);

    fP1->cd(3);
    fH1[0] = gh1(cb000401_tcalm004,"TCalm004", "trk_27/llhr_slope");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("LHR(slope)");
    //    fH1[0]->GetXaxis()->SetRangeUser(-25.,25.);
    fH1[0]->Draw();

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"identified electrons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);

    fP1->cd(4);
    fH1[0] = gh1(cb000401_tcalm004,"TCalm004", "trk_28/p");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetTitle("track momentum,  MeV/c");
    //    fH1[0]->GetXaxis()->SetRangeUser(-25.,25.);
    fH1[0]->Draw();

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"identified electrons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);
  }
//-----------------------------------------------------------------------------
// Fig.  41  (a) trk_1/xs for e0000001 and m0000001
//-----------------------------------------------------------------------------
  if (Figure == 41) {
    //    TArrow*      arr;
    TPaveLabel*  label;
    TLegend*     leg;

    fCanvas  = new_slide(name,title,1,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    fH1[0] = (TH1F*) gh1(e0000001_tcalm002,"TCalm002", "trk_1/xslope")->Clone("h_fig_41");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetRangeUser(-10,10);
    fH1[0]->GetXaxis()->SetTitle("slope/#sigma(slope)");
    //    fH1[0]->GetXaxis()->SetRangeUser(-25.,25.);
    fH1[0]->Scale(1./fH1[0]->Integral());
    fH1[0]->Fit("gaus");
    //    fH1[0]->DrawNormalized();

    fH1[1] = gh1(m0000001_tcalm002,"TCalm002", "trk_1/xslope");
    fH1[1]->SetLineColor(kBlue-7);
    fH1[1]->SetFillColor(kBlue-7);
    fH1[1]->SetFillStyle(3003);
    fH1[1]->DrawNormalized("same");

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"electrons","f");
    leg->AddEntry(fH1[1],"muons"    ,"f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);
  }

//-----------------------------------------------------------------------------
// Fig.  42  (a) trk_1/llhr_dedx for e0000001 and m0000001
//-----------------------------------------------------------------------------
  if (Figure == 42) {
    //    TArrow*      arr;
    TPaveLabel*  label;
    TLegend*     leg;

    fCanvas  = new_slide(name,title,1,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    fH1[0] = gh1(e0000001_tcalm002,"TCalm002", "trk_1/llhr_dedx");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetRangeUser(-15,15);
    fH1[0]->GetXaxis()->SetTitle("log LHR(dE/dX)");
    //    fH1[0]->GetXaxis()->SetRangeUser(-25.,25.);
    fH1[0]->DrawNormalized();

    fH1[1] = gh1(m0000001_tcalm002,"TCalm002", "trk_1/llhr_dedx");
    fH1[1]->SetLineColor(kBlue-7);
    fH1[1]->SetFillColor(kBlue-7);
    fH1[1]->SetFillStyle(3003);
    fH1[1]->DrawNormalized("same");

    leg = new TLegend(0.15,0.7,0.30,0.8,"");
    leg->AddEntry(fH1[0],"electrons","f");
    leg->AddEntry(fH1[1],"muons"    ,"f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);
  }

//-----------------------------------------------------------------------------
// Fig.  43  (a) trk_1/llhr_trk for e0000001 and m0000001
//-----------------------------------------------------------------------------
  if (Figure == 43) {
    //    TArrow*      arr;
    TPaveLabel*  label;
    TLegend*     leg;

    fCanvas  = new_slide(name,title,1,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    fH1[0] = gh1(e0000001_tcalm002,"TCalm002", "trk_1/llhr_trk");
    
    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("");
    fH1[0]->GetXaxis()->SetRangeUser(-15,15);
    fH1[0]->GetXaxis()->SetTitle("track-only log LHR");
    //    fH1[0]->GetXaxis()->SetRangeUser(-25.,25.);
    fH1[0]->DrawNormalized();

    fH1[1] = gh1(m0000001_tcalm002,"TCalm002", "trk_1/llhr_trk");
    fH1[1]->SetLineColor(kBlue-7);
    fH1[1]->SetFillColor(kBlue-7);
    fH1[1]->SetFillStyle(3003);
    fH1[1]->DrawNormalized("same");

    leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"electrons","f");
    leg->AddEntry(fH1[1],"muons"    ,"f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.15,0.85,0.5,0.95,52);
    label->SetTextSize(0.4);
  }

//-----------------------------------------------------------------------------
// Fig.  44: log LHR (trk) rejection vs efficiency, electrons and muons
//-----------------------------------------------------------------------------
  if (Figure == 44) {
    TPaveLabel*  label;
    TLegend*     leg;

    fCanvas  = new TCanvas(name,title,1100,700);
    fCanvas->Divide(1,1);
    fP1      = (TPad*) fCanvas->GetPad(0);

    fP1->cd(1);
    fH1[0] = gh1(e0000001_tcalm002,"TCalm002", "trk_1/llhr_trk");
    fH1[1] = gh1(m0000001_tcalm002,"TCalm002", "trk_1/llhr_trk");
    
    TH1* h1 = (TH1*) fH1[0]->Clone("h_fig_20_ele_efficiency");
    TH1* h2 = (TH1*) fH1[1]->Clone("h_fig_20_muo_rejection");

    int nbins = h1->GetNbinsX();

    double h1_integral = fH1[0]->Integral();
    double h2_integral = fH1[1]->Integral();

    double eff1, ome2;

    for (int i=1; i<=nbins; i++) { 
      eff1 = 1 - fH1[0]->Integral(1,i)/h1_integral; 
      ome2 = 1 - fH1[1]->Integral(1,i)/h2_integral; 

      h1->SetBinContent(i,eff1);
      h2->SetBinContent(i,ome2);
    }
    
    h1->GetXaxis()->SetRangeUser(-15,10);
    h1->SetLineColor(2);
    h1->SetTitle("Probability to identity a particle");
    h1->GetXaxis()->SetTitle("log(L(e)/L(#mu))");

    if (fPlotMode == TPlotNote::kTalkMode) {
      h1->SetStats(0);
    }

    h1->SetStats(0);
    h1->Draw();
    h2->Draw("same");
    
    leg = new TLegend(0.15,0.45,0.40,0.55,"");
    leg->AddEntry(h1,"electrons","f");
    leg->AddEntry(h2,"muons"    ,"f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"tracks passing Set C cuts",0.55,0.85,0.90,0.95,52);
    label->SetTextSize(0.4);
  }

//-----------------------------------------------------------------------------
// Fig.  60: track reconstruction efficiencies vs luminosity
//-----------------------------------------------------------------------------
  if (Figure == 60) {
    //    TPaveLabel*  label;
    TLegend*     leg;

    fCanvas  = new_slide(name,title,2,2,1300,900);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    fH1[0] = create_eff_hist(egun0001_tcalm002,"TCalm002", "evt_0/ce_costh","evt_1/ce_costh","eff_t0","",20,1.0);
    fH1[1] = create_eff_hist(egun0001_tcalm002,"TCalm002", "evt_0/ce_costh","evt_6/ce_costh","eff_t1","",24,1.0);
    
    fH1[0]->GetXaxis()->SetRangeUser(-0.6,0.6);
    fH1[0]->SetTitle("Efficiency, background: x0");
    fH1[0]->Draw("ep");
    fH1[1]->SetLineColor(kBlue-7);
    fH1[1]->SetMarkerColor(kBlue-7);
    fH1[1]->Draw("same,ep");
    
    leg = new TLegend(0.15,0.75,0.40,0.85,"");
    leg->AddEntry(fH1[0],"all tracks","h");
    leg->AddEntry(fH1[1],"Set C tracks","h");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
					// x1 background
    fP1->cd(2);
    fH1[0] = create_eff_hist(egun0101_tcalm002,"TCalm002", "evt_0/ce_costh","evt_1/ce_costh","eff_t0","",20,1.0);
    fH1[1] = create_eff_hist(egun0101_tcalm002,"TCalm002", "evt_0/ce_costh","evt_6/ce_costh","eff_t1","",24,1.0);
    
    fH1[0]->GetXaxis()->SetRangeUser(-0.6,0.6);
    fH1[0]->SetTitle("Efficiency, background: x1");
    fH1[0]->Draw("ep");
    fH1[1]->SetLineColor(kBlue-7);
    fH1[1]->SetMarkerColor(kBlue-7);
    fH1[1]->Draw("same,ep");
    
    leg = new TLegend(0.15,0.75,0.40,0.85,"");
    leg->AddEntry(fH1[0],"all tracks","h");
    leg->AddEntry(fH1[1],"Set C tracks","h");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
					// x2 background
    fP1->cd(3);
    fH1[0] = create_eff_hist(egun0201_tcalm002,"TCalm002", "evt_0/ce_costh","evt_1/ce_costh","eff_t0","",20,1.0);
    fH1[1] = create_eff_hist(egun0201_tcalm002,"TCalm002", "evt_0/ce_costh","evt_6/ce_costh","eff_t1","",24,1.0);
    
    fH1[0]->GetXaxis()->SetRangeUser(-0.6,0.6);
    fH1[0]->SetTitle("Efficiency, background: x2");
    fH1[0]->Draw("ep");
    fH1[1]->SetLineColor(kBlue-7);
    fH1[1]->SetMarkerColor(kBlue-7);
    fH1[1]->Draw("same,ep");
    
    leg = new TLegend(0.15,0.75,0.40,0.85,"");
    leg->AddEntry(fH1[0],"all tracks","h");
    leg->AddEntry(fH1[1],"Set C tracks","h");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
					// x4 background
    fP1->cd(4);
    fH1[0] = create_eff_hist(egun0401_tcalm002,"TCalm002", "evt_0/ce_costh","evt_1/ce_costh","eff_t0","",20,1.0);
    fH1[1] = create_eff_hist(egun0401_tcalm002,"TCalm002", "evt_0/ce_costh","evt_6/ce_costh","eff_t1","",24,1.0);
    
    fH1[0]->GetXaxis()->SetRangeUser(-0.6,0.6);
    fH1[0]->SetTitle("Efficiency, background: x4");
    fH1[0]->Draw("ep");
    fH1[1]->SetLineColor(kBlue-7);
    fH1[1]->SetMarkerColor(kBlue-7);
    fH1[1]->Draw("same,ep");
    
    leg = new TLegend(0.15,0.75,0.40,0.85,"");
    leg->AddEntry(fH1[0],"all tracks","h");
    leg->AddEntry(fH1[1],"Set C tracks","h");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  61: track reconstruction efficiencies [-0.2,0.2] vs luminosity
//-----------------------------------------------------------------------------
  if (Figure == 61) {
    //    TPaveLabel*  label;
    TLegend*     leg;

    double   x[4] = {0., 1., 2., 4.};
    double  ex[4] = {0., 0., 0., 0.};
    double  eff1[4], eff6[4], err1[4], err6[4];
    double  q0, q1, q6;

    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    
    q0 = gh1(egun0001_tcalm002,"TCalm002","evt_0/ce_costh")->Integral(41,60);
    q1 = gh1(egun0001_tcalm002,"TCalm002","evt_1/ce_costh")->Integral(41,60);
    q6 = gh1(egun0001_tcalm002,"TCalm002","evt_6/ce_costh")->Integral(41,60);

    eff1[0] = q1/q0;
    err1[0] = sqrt(q0-q1)/q0;
    eff6[0] = q6/q0;
    err6[0] = sqrt(q0-q6)/q0;

    q0 = gh1(egun0101_tcalm002,"TCalm002","evt_0/ce_costh")->Integral(41,60);
    q1 = gh1(egun0101_tcalm002,"TCalm002","evt_1/ce_costh")->Integral(41,60);
    q6 = gh1(egun0101_tcalm002,"TCalm002","evt_6/ce_costh")->Integral(41,60);

    eff1[1] = q1/q0;
    err1[1] = sqrt(q0-q1)/q0;
    eff6[1] = q6/q0;
    err6[1] = sqrt(q0-q6)/q0;

    q0 = gh1(egun0201_tcalm002,"TCalm002","evt_0/ce_costh")->Integral(41,60);
    q1 = gh1(egun0201_tcalm002,"TCalm002","evt_1/ce_costh")->Integral(41,60);
    q6 = gh1(egun0201_tcalm002,"TCalm002","evt_6/ce_costh")->Integral(41,60);

    eff1[2] = q1/q0;
    err1[2] = sqrt(q0-q1)/q0;
    eff6[2] = q6/q0;
    err6[2] = sqrt(q0-q6)/q0;

    q0 = gh1(egun0401_tcalm002,"TCalm002","evt_0/ce_costh")->Integral(41,60);
    q1 = gh1(egun0401_tcalm002,"TCalm002","evt_1/ce_costh")->Integral(41,60);
    q6 = gh1(egun0401_tcalm002,"TCalm002","evt_6/ce_costh")->Integral(41,60);

    eff1[3] = q1/q0;
    err1[3] = sqrt(q0-q1)/q0;
    eff6[3] = q6/q0;
    err6[3] = sqrt(q0-q6)/q0;


    TGraphErrors* gr1 = new TGraphErrors(4,x,eff1,ex,err1);
    TGraphErrors* gr6 = new TGraphErrors(4,x,eff6,ex,err6);
    
    TH2F* h2 = new TH2F("h2","",1,-1,5,1,0,1.1);
    h2->SetTitle("track reconstruction efficiency vs background occupancy");
    h2->GetXaxis()->SetTitle("nominal background, scale factor");
    h2->SetStats(0);
    h2->Draw();

    gr1->SetMarkerStyle(20);
    gr1->Draw("ep");

    gr6->SetMarkerStyle(24);
    gr6->Draw("ep");

    leg = new TLegend(0.15,0.35,0.40,0.45,"");
    leg->AddEntry(gr1,"all tracks","ep");
    leg->AddEntry(gr6,"Set C tracks","ep");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  62: LogLH_CAL selection efficiency vs luminosity
//-----------------------------------------------------------------------------
  if (Figure == 62) {
    //    TPaveLabel*  label;
    TLegend*     leg;

    double   x[4] = {0., 1., 2., 4.};
    double  ex[4] = {0., 0., 0., 0.};
    double  effe[4], effm[4], erre[4], errm[4];
    double  q0, q1; //, q6;

    fCanvas  = new TCanvas(name,title,1000,700);
    fP1      = (TPad*) fCanvas->GetPad(0);

    //    fP1->cd(1);
    
    q1 = gh1(egun0001_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral(101,200);
    q0 = gh1(egun0001_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral();

    effe[0] = q1/q0;
    erre[0] = sqrt(q0-q1)/q0;

    q1 = gh1(egun0101_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral(101,200);
    q0 = gh1(egun0101_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral();

    effe[1] = q1/q0;
    erre[1] = sqrt(q0-q1)/q0;

    q1 = gh1(egun0201_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral(101,200);
    q0 = gh1(egun0201_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral();

    effe[2] = q1/q0;
    erre[2] = sqrt(q0-q1)/q0;

    q1 = gh1(egun0401_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral(101,200);
    q0 = gh1(egun0401_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral();

    effe[3] = q1/q0;
    erre[3] = sqrt(q0-q1)/q0;

    TGraphErrors* gre = new TGraphErrors(4,x,effe,ex,erre);

    q1 = gh1(mgun0001_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral(101,200);
    q0 = gh1(mgun0001_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral();

    effm[0] = q1/q0;
    errm[0] = sqrt(q1)/q0;

    q1 = gh1(mgun0101_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral(101,200);
    q0 = gh1(mgun0101_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral();

    effm[1] = q1/q0;
    errm[1] = sqrt(q1)/q0;

    q1 = gh1(mgun0201_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral(101,200);
    q0 = gh1(mgun0201_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral();

    effm[2] = q1/q0;
    errm[2] = sqrt(q1)/q0;

    q1 = gh1(mgun0401_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral(101,200);
    q0 = gh1(mgun0401_tcalm002,"TCalm002","trk_1/llhr_cal")->Integral();

    effm[3] = q1/q0;
    errm[3] = sqrt(q1)/q0;

    TGraphErrors* grm = new TGraphErrors(4,x,effm,ex,errm);
    
    TH2F* h2 = new TH2F("h2","",1,-1,5,100,0.0,1.05);
    h2->SetTitle("");
    h2->GetXaxis()->SetTitle("Background Scale Factor");
    h2->GetYaxis()->SetTitle("Identification Efficiency");
    h2->GetYaxis()->SetRangeUser(1.e-4,1.2);
    h2->SetStats(0);
    h2->Draw();

    if (fPlotMode == TPlotNote::kNoteMode) {
      gre->SetMarkerStyle(20);
      grm->SetMarkerStyle(24);
    }
    else if (fPlotMode == kTalkMode) {
      gre->SetMarkerStyle(20);
      gre->SetMarkerColor(2);
      grm->SetMarkerStyle(20);
      grm->SetMarkerColor(4);
    }

    gre->Draw("ep");
    grm->Draw("ep");

    leg = new TLegend(0.15,0.35,0.40,0.45,"");
    leg->AddEntry(gre,"#epsilon(e)"  ,"ep");
    leg->AddEntry(grm,"#epsilon(#mu)","ep");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  63: LogLH_TRK selection efficiency vs luminosity
//-----------------------------------------------------------------------------
  if (Figure == 63) {
    //    TPaveLabel*  label;
    TLegend*     leg;

    double   x[4] = {0., 1., 2., 4.};
    double  ex[4] = {0., 0., 0., 0.};
    double  effe[4], effm[4], erre[4], errm[4];
    double  q0, q1; //, q6;

    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);
    
    q1 = gh1(egun0001_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral(91,200);
    q0 = gh1(egun0001_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral();

    effe[0] = q1/q0;
    erre[0] = sqrt(q0-q1)/q0;

    q1 = gh1(egun0101_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral(91,200);
    q0 = gh1(egun0101_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral();

    effe[1] = q1/q0;
    erre[1] = sqrt(q0-q1)/q0;

    q1 = gh1(egun0201_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral(91,200);
    q0 = gh1(egun0201_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral();

    effe[2] = q1/q0;
    erre[2] = sqrt(q0-q1)/q0;

    q1 = gh1(egun0401_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral(91,200);
    q0 = gh1(egun0401_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral();

    effe[3] = q1/q0;
    erre[3] = sqrt(q0-q1)/q0;

    TGraphErrors* gre = new TGraphErrors(4,x,effe,ex,erre);

    q1 = gh1(mgun0001_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral(91,200);
    q0 = gh1(mgun0001_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral();

    effm[0] = q1/q0;
    errm[0] = sqrt(q1)/q0;

    q1 = gh1(mgun0101_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral(91,200);
    q0 = gh1(mgun0101_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral();

    effm[1] = q1/q0;
    errm[1] = sqrt(q1)/q0;

    q1 = gh1(mgun0201_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral(91,200);
    q0 = gh1(mgun0201_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral();

    effm[2] = q1/q0;
    errm[2] = sqrt(q1)/q0;

    q1 = gh1(mgun0401_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral(91,200);
    q0 = gh1(mgun0401_tcalm002,"TCalm002","trk_1/llhr_trk")->Integral();

    effm[3] = q1/q0;
    errm[3] = sqrt(q1)/q0;

    TGraphErrors* grm = new TGraphErrors(4,x,effm,ex,errm);
    
    TH2F* h2 = new TH2F("h2","",1,-1,5,100,0.0,1.05);
    h2->SetTitle("LLHR(TRK) : background occupancy dependence");
    h2->GetXaxis()->SetTitle("nominal background, scale factor");
    h2->SetStats(0);
    h2->Draw();

    gre->SetMarkerStyle(20);
    gre->Draw("ep");

    grm->SetMarkerStyle(24);
    grm->Draw("ep");

    leg = new TLegend(0.15,0.35,0.40,0.45,"");
    leg->AddEntry(gre,"#epsilon(e)","ep");
    leg->AddEntry(grm,"#epsilon(#mu)","ep");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }
//-----------------------------------------------------------------------------
//  Fig.100: distribution of the E/P ratio for downstream tracks
//           which are matched with an EMC with energy>10MeV and
//           which have a distance from the EMC (for each coordinate)
//           less than 15 mm (half crystal size)
//-----------------------------------------------------------------------------
  if (Figure == 100) {
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(cb000301_cosmicAnalyzer,"ex00", "EoverP_eminus");
    fH1[1] = gh1(cb000301_cosmicAnalyzer,"ex00", "EoverP_muminus");
    
    fP1->cd(1);
    //gPad->SetLogy(1);

    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("E/P ratio for downstream tracks ");
    fH1[0]->GetXaxis()->SetTitle("E/P");
    //fH1[0]->GetXaxis()->SetRangeUser(-2500,2500);
    //fH1[0]->SetMinimum(0.5);
    fH1[0]->Draw();

    fH1[1]->SetLineColor(kBlue-7);
    fH1[1]->SetFillColor(kBlue-7);
    fH1[1]->SetFillStyle(3003);
    fH1[1]->Draw("same");

  //   TArrow* arr;

//     arr = new TArrow(-1000.,4,-1000.,0.7,0.015);
//     arr->SetLineWidth(3);
//     arr->Draw();
//     arr = new TArrow( 1000.,4, 1000.,0.7,0.015);
//     arr->SetLineWidth(3);
//     arr->Draw();

    TLegend* leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"electrons","f");
    leg->AddEntry(fH1[1],"muons" ,"f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }

  if (Figure == 101) {
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(cb000301_cosmicAnalyzer,"ex00", "deltatime_eminus_0");
    fH1[1] = gh1(cb000301_cosmicAnalyzer,"ex00", "deltatime_muminus_0");
    
    fP1->cd(1);
    //gPad->SetLogy(1);

    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("#Deltat distribution for downstream tracks");
    fH1[0]->GetXaxis()->SetTitle("E/P");
    //fH1[0]->GetXaxis()->SetRangeUser(-2500,2500);
    //fH1[0]->SetMinimum(0.5);
    fH1[0]->Draw();

    fH1[1]->SetLineColor(kBlue-7);
    //fH1[1]->SetFillColor(kBlue-7);
    //fH1[1]->SetFillStyle(3003);
    fH1[1]->Draw("same");

    // fH1[0]->Draw("same");

  //   TArrow* arr;

//     arr = new TArrow(-1000.,4,-1000.,0.7,0.015);
//     arr->SetLineWidth(3);
//     arr->Draw();
//     arr = new TArrow( 1000.,4, 1000.,0.7,0.015);
//     arr->SetLineWidth(3);
//     arr->Draw();

    TLegend* leg = new TLegend(0.15,0.7,0.40,0.8,"");
    leg->AddEntry(fH1[0],"electrons","f");
    leg->AddEntry(fH1[1],"muons" ,"f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }

}

