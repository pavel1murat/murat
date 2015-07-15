///////////////////////////////////////////////////////////////////////////////
// plots for Mu2e-4256: Mu2e Calorimeter-based PID
// -------------------------
// directory with the histogram files: by default: $WORK_DIR/results
// can be redefined with "mu2e.HistDir" in .rootrc
// 
// figures:
// ---------
// Fig.   1: Du
// Fig.   2: Dv
// Fig.   3: Du vs Path
// Fig.   4: Dv vs Path
// Fig.   5: Du(corrected) vs Path - ProjectionY
// Fig.   6: Dv(corrected) vs Path - ProjectionY - electrons
// Fig.   7: Dt for electrons 
// Fig.   8: chi2 (tcm) selection for electrons 
// Fig.   9: E/P selection - CE  without background(e00s1212) and with the default background(e00s1412)
//
// Fig.  11: E/P for different occupancies and datasets : TRK_19
// Fig.  12: E/P for different occupancies and datasets : TRK_19
// Fig.  13: E/P  vs DT TRK_19
//
// Fig. 101: LLHR_CAL rejection vs luminosity
// Fig. 102: muon rejection vs electron efficiency, vary calorimeter performance 
//           parameters 
////////////////////////////////////////////////////////////////////////////////
#include "plot/TMu2e4256.hh"

#include "TPaveLabel.h"
#include "TArrow.h"
#include "TObjArray.h"
#include "TGraphErrors.h"
#include "TObjString.h"
#include <iostream>

ClassImp(TMu2e4256)

//_____________________________________________________________________________
TMu2e4256::TMu2e4256(int PlotMode, int BlessingMode): TPlotNote(PlotMode,BlessingMode) {

  const char* hist_dir;
  char         res_dir[500];
  
  fFiguresDir   = Form("%s/mu2e_4256",gEnv->GetValue("mu2e.Figures","."));
  fWorkDir      = gSystem->Getenv("WORK_DIR");

  sprintf(res_dir,"%s/results",fWorkDir.Data());
  hist_dir      = gEnv->GetValue("mu2e.HistDir",res_dir);

  sprintf(f_e00s1212,"%s/v4_2_1/e00s1212.track_ana.hist", hist_dir); 
  sprintf(f_m00s1212,"%s/v4_2_1/m00s1212.track_ana.hist", hist_dir); 
  sprintf(f_e00s1412,"%s/v4_2_1/e00s1412.track_ana.hist", hist_dir); 
  sprintf(f_m00s1412,"%s/v4_2_1/m00s1412.track_ana.hist", hist_dir); 
  sprintf(f_e00s1512,"%s/v4_2_1/e00s1512.track_ana.hist", hist_dir); 
  sprintf(f_m00s1512,"%s/v4_2_1/m00s1512.track_ana.hist", hist_dir); 
  sprintf(f_e00s1612,"%s/v4_2_1/e00s1612.track_ana.hist", hist_dir); 
  sprintf(f_m00s1612,"%s/v4_2_1/m00s1612.track_ana.hist", hist_dir); 

  sprintf(res_dir,"%s",fWorkDir.Data());
  hist_dir = gEnv->GetValue("mu2e.HistDir",res_dir);
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
TMu2e4256::~TMu2e4256() {
};

//-----------------------------------------------------------------------------
const char* TMu2e4256::GetFiguresDir() {
  return fFiguresDir.Data();
}

void TMu2e4256::get_filename(int Figure, char* Filename) {

  if      (Figure ==  -1) sprintf(Filename,"fig_%03i_%s",Figure,"momentum");
  else if (Figure == -12) sprintf(Filename,"fig_%03i_%s",Figure,"prt_energy");
  else {
//-----------------------------------------------------------------------------
// undefined, just figure
//-----------------------------------------------------------------------------
    sprintf(Filename,"fig_%03i",Figure);
  }
}
//-----------------------------------------------------------------------------
// determine name of the .eps and .gif files in the ./figures directory
// today's default set is 'frr_03', all the plots should be coming from it
//-----------------------------------------------------------------------------
void TMu2e4256::remake_plots() {

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


//---------------------------------------------------------------------------------------
// example: FnEle: "~/hist/mu2e/v4_2_1/e00s1212.track_ana.hist"
//          FnMuo: "~/hist/mu2e/v4_2_1/m00s1212.track_ana.hist"
// 'MuonRecoEff' - reconstruction efficiency at zero background occupancy
//---------------------------------------------------------------------------------------
void TMu2e4256::create_llhr_rejection_graph(const char* FnEle, 
					    const char* FnMuo, 
					    double      MuonRecoEff,
					    TGraph*&    Graph      ) {

  TH1* h_llhr_e = (TH1*) gh1(FnEle,"TrackAna","trk_25/llhr_cal")->Clone("h_llhr_e");
  TH1* h_llhr_m = (TH1*) gh1(FnMuo,"TrackAna","trk_25/llhr_cal")->Clone("h_llhr_m");

  double pe[1000], rm[1000];

  double qnm_13, qnm_25, qne_13, qne_25;

  int nbx = h_llhr_e->GetNbinsX();

  qne_25 = gh1(FnEle,"TrackAna","trk_25/llhr_cal")->GetEntries();
  qne_13 = gh1(FnEle,"TrackAna","trk_13/llhr_cal")->GetEntries();

  qnm_25 = gh1(FnMuo,"TrackAna","trk_25/llhr_cal")->GetEntries();
  qnm_13 = gh1(FnMuo,"TrackAna","trk_13/llhr_cal")->GetEntries();

//-----------------------------------------------------------------------------
// reconstruction efficiency above 100 MeV/c depends on occupancy, 
// Ralf's numbers correspond to zero background level, account for that
//-----------------------------------------------------------------------------
  double muon_reco_eff, nt, nev, reco_eff_sf(1);

  nt  = gh1(FnMuo,"TrackAna","trk_13/ep")->GetEntries();
  nev = gh1(FnMuo,"TrackAna","evt_0/ntrk")->GetEntries();

  muon_reco_eff = nt/nev;

  reco_eff_sf   = MuonRecoEff/muon_reco_eff;

  double se = h_llhr_e->Integral();
  double sm = h_llhr_m->Integral();

  double qe, qm;

  for (int i=0; i<nbx; i++) {
    qe = h_llhr_e->Integral(1,i+1);
    qm = h_llhr_m->Integral(1,i+1);

    pe[i] =(1-qe/se)*(qne_25/qne_13);
    rm[i] = 1./(1-qm/sm + 1.e-6)*(qnm_13/qnm_25)*reco_eff_sf;

    if (GetDebugBit(1) == 1) {
      printf(" i, llhr , qe, qm, prob(e) , rej(mu) : %3i %10.3f %10.3f %10.3f %10.5f %10.3f\n",
	    i,h_llhr_e->GetBinCenter(i+1),qe,qm,pe[i],rm[i]);
    }
  }

  if (Graph != NULL) delete Graph;

  Graph = new TGraph(nbx,pe,rm);
}


//-----------------------------------------------------------------------------
void TMu2e4256::plot(Int_t Figure, const char* CanvasName) {
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
// Fig.  1: single electrons , DU
//-----------------------------------------------------------------------------
  if (Figure == 1) {
    TH1F   *h1;

    fCanvas  = new TCanvas(name,title,1,1,1100,800);

    fCanvas->cd(1);

    h1 = gh1(f_e00s1212,"TrackAna","trk_1/du");

    h1->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  2: single electrons , DV
//-----------------------------------------------------------------------------
  if (Figure == 2) {
    TH1F   *h1;

    fCanvas  = new TCanvas(name,title,1,1,1100,800);

    fCanvas->cd(1);

    h1 = gh1(f_e00s1212,"TrackAna","trk_1/dv");

    h1->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  3: DU vs path
//-----------------------------------------------------------------------------
  if (Figure == 3) {
    TH1F   *h1;

    fCanvas  = new TCanvas(name,title,1,1,1100,800);

    fCanvas->cd(1);

    h1 = gh1(f_e00s1212,"TrackAna","trk_1/du_vs_path");

    h1->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  4: DV vs path
//-----------------------------------------------------------------------------
  if (Figure == 4) {
    TH2F   *h2;

    fCanvas  = new TCanvas(name,title,1,1,1100,800);

    fCanvas->cd(1);

    h2 = gh2(f_e00s1212,"TrackAna","trk_1/dv_vs_path");

    h2->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  5: DU vs path
//-----------------------------------------------------------------------------
  if (Figure == 5) {
    TH2F   *h2;

    fCanvas  = new TCanvas(name,title,1,1,1100,800);

    fCanvas->cd(1);

    h2 = gh2(f_e00s1212,"TrackAna","trk_1/duc_vs_path");

    h2->ProjectionY()->Fit("gaus");
  }
//-----------------------------------------------------------------------------
// Fig.  6: DV vs path
//-----------------------------------------------------------------------------
  if (Figure == 6) {
    TH2F   *h2;

    fCanvas  = new TCanvas(name,title,1,1,1100,800);

    fCanvas->cd(1);

    h2 = gh2(f_e00s1212,"TrackAna","trk_1/dvc_vs_path");

    h2->ProjectionY()->Fit("gaus");
  }
//-----------------------------------------------------------------------------
// Fig.  7: Dt selection
//-----------------------------------------------------------------------------
  if (Figure == 7) {
    TH1F   *h_1212, *h_1412;

    fCanvas  = new TCanvas(name,title,1,1,1100,800);

    fCanvas->cd(1);

    h_1212 = gh1(f_e00s1212,"TrackAna","trk_19/dt");
    h_1412 = gh1(f_e00s1412,"TrackAna","trk_19/dt");

    double qn_1412=h_1412->Integral();
    h_1212->SetLineColor(2);
    h_1212->SetFillColor(2);
    h_1212->SetFillStyle(3002);
    h_1212->DrawNormalized("",qn_1412);

    h_1412->Draw("same");

    TArrow  *a1, *a2;

    a1 = new TArrow(-3.,100,-3,10);
    a1->SetArrowSize(0.015);
    a1->Draw();

    a2 = new TArrow(3.,100,3,10);
    a2->SetArrowSize(0.015);
    a2->Draw();

    TLegend* leg = new TLegend(0.15,0.65,0.4,0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(h_1212,"CE only"    ,"ep");
    leg->AddEntry(h_1412,"CE+MIXP3-x1","ep");

    leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  8: Chi2(TCM) selection
//-----------------------------------------------------------------------------
  if (Figure == 8) {
    TH1F   *h_1212, *h_1412;

    fCanvas  = new TCanvas(name,title,1,1,1100,800);
    gPad->SetLogy(1);

    fCanvas->cd(1);

    h_1212 = gh1(f_e00s1212,"TrackAna","trk_13/chi2tcm");
    h_1412 = gh1(f_e00s1412,"TrackAna","trk_13/chi2tcm");

    double qn_1412=h_1412->Integral();
    h_1212->SetLineColor(2);
    h_1212->SetFillColor(2);
    h_1212->SetFillStyle(3002);
    h_1212->DrawNormalized("",qn_1412);

    h_1412->Draw("same");

    TArrow  *a1; //, *a2;

    a1 = new TArrow(100.,100,100.,10);
    a1->SetArrowSize(0.015);
    a1->Draw();

    // a2 = new TArrow(3.,100,3,10);
    // a2->SetArrowSize(0.015);
    // a2->Draw();

    TLegend* leg = new TLegend(0.15,0.65,0.4,0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(h_1212,"CE only"    ,"f");
    leg->AddEntry(h_1412,"CE+MIXP3-x1","f");

    leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  9: E/P selection
//-----------------------------------------------------------------------------
  if (Figure == 9) {
    TH1F   *h_1212, *h_1412;

    fCanvas  = new TCanvas(name,title,1,1,1200,800);
    gPad->SetLogy(1);

    fCanvas->cd(1);

    h_1212 = gh1(f_e00s1212,"TrackAna","trk_19/ep");
    h_1412 = gh1(f_e00s1412,"TrackAna","trk_19/ep");

    double qn_1412=h_1412->Integral();
    h_1212->SetLineColor(2);
    h_1212->SetFillColor(2);
    h_1212->SetFillStyle(3001);
    h_1212->DrawNormalized("",qn_1412);

    h_1412->Draw("same");

    TArrow  *a1; //, *a2;

    a1 = new TArrow(1.15,100,1.15,10);
    a1->SetArrowSize(0.015);
    a1->Draw();

    // a2 = new TArrow(3.,100,3,10);
    // a2->SetArrowSize(0.015);
    // a2->Draw();

    TLegend* leg = new TLegend(0.15,0.65,0.4,0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(h_1212,"CE only"    ,"f");
    leg->AddEntry(h_1412,"CE+MIXP3-x1","f");

    leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  11: delta(T) for different datasets, TRK_19 (track+cluster)
//-----------------------------------------------------------------------------
  if (Figure == 11) {
    TH1F   *he_1212, *hm_1212, *hm_1412, *hm_1512;

    fCanvas  = new TCanvas(name,title,1,1,1200,800);
    gPad->SetLogy(1);

    fCanvas->cd(1);

    he_1212 = gh1(f_e00s1212,"TrackAna","trk_19/dt");

    hm_1212 = gh1(f_m00s1212,"TrackAna","trk_19/dt");
    hm_1412 = gh1(f_m00s1412,"TrackAna","trk_19/dt");
    hm_1512 = gh1(f_m00s1512,"TrackAna","trk_19/dt");

				// drawing part
    hm_1512->SetMaximum(0.12);
    hm_1512->SetStats(0);
    hm_1512->DrawNormalized("",1);

    hm_1412->SetLineColor(2);
    hm_1412->SetFillColor(2);
    hm_1412->SetFillStyle(3002);
    hm_1412->DrawNormalized("same",1);

    hm_1212->SetFillColor(kBlue-7);
    hm_1212->SetFillStyle(3001);
    hm_1212->DrawNormalized("same",1);

    he_1212->SetMarkerSize(1);
    he_1212->SetMarkerStyle(20);
    he_1212->DrawNormalized("same,p",1);

    TLegend* leg = new TLegend(0.15,0.65,0.4,0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(he_1212,"Electrons 105 MeV/c","p");
    leg->AddEntry(hm_1212,"muons 105 MeV/c"    ,"f");
    leg->AddEntry(hm_1412,"M105+MIXP3-x1","f");
    leg->AddEntry(hm_1512,"M105+MIXP3-x2","f");

    leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  12: E/P for different datasets, TRK_25 (track+cluster)
//-----------------------------------------------------------------------------
  if (Figure == 12) {
    TH1F   *he_1212, *hm_1212, *hm_1412, *hm_1512;

    fCanvas  = new TCanvas(name,title,1,1,1200,800);
    gPad->SetLogy(1);

    fCanvas->cd(1);

    he_1212 = gh1(f_e00s1212,"TrackAna","trk_25/ep");
    hm_1212 = gh1(f_m00s1212,"TrackAna","trk_25/ep");

    hm_1412 = gh1(f_m00s1412,"TrackAna","trk_25/ep");


    hm_1512 = gh1(f_m00s1512,"TrackAna","trk_25/ep");

    hm_1512->SetMaximum(0.12);
    hm_1512->SetStats(0);
    hm_1512->DrawNormalized("",1);


    hm_1412->SetLineColor(2);
    hm_1412->SetFillColor(2);
    hm_1412->SetFillStyle(3002);
    hm_1412->DrawNormalized("same",1);

    hm_1212->SetFillColor(kBlue-7);
    hm_1212->SetFillStyle(3001);
    hm_1212->DrawNormalized("same",1);

    he_1212->SetMarkerSize(1);
    he_1212->SetMarkerStyle(20);
    he_1212->DrawNormalized("same,p",1);

    TLegend* leg = new TLegend(0.15,0.65,0.4,0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(he_1212,"Electrons 105 MeV/c","p");
    leg->AddEntry(hm_1212,"muons 105 MeV/c"    ,"f");
    leg->AddEntry(hm_1412,"M105+MIXP3-x1","f");
    leg->AddEntry(hm_1512,"M105+MIXP3-x2","f");

    leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  13: E/P vs DT
//-----------------------------------------------------------------------------
  if (Figure == 13) {
    TH1F   *he_1212, *hm_1212; // , *hm_1412, *hm_1512;

    fCanvas  = new TCanvas(name,title,1,1,1400,700);
    //    gPad->SetLogy(1);

    fCanvas->Divide(2,1);
    fCanvas->cd(1);

    he_1212 = gh1(f_e00s1212,"TrackAna","trk_19/ep_vs_dt");
    he_1212->SetMarkerColor(kRed+2);
    he_1212->SetStats(0);
    he_1212->Draw();

    TLegend* leg1 = new TLegend(0.15,0.75,0.4,0.85);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);

    leg1->AddEntry(he_1212,"Electrons, 105 MeV/c","f");

    leg1->Draw();

    fCanvas->cd(2);
    hm_1212 = gh1(f_m00s1212,"TrackAna","trk_19/ep_vs_dt");

    hm_1212->SetMarkerColor(kBlue+2);
    hm_1212->SetStats(0);
    hm_1212->Draw();

    TLegend* leg2 = new TLegend(0.15,0.75,0.4,0.85);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);

    leg2->AddEntry(hm_1212,"muons, 105 MeV/c"    ,"f");

    leg2->Draw();
  }
//-----------------------------------------------------------------------------
// Fig. 101: single electrons , DV
//-----------------------------------------------------------------------------
  if (Figure == 101) {
    TGraph   *gr_x0(0), *gr_x1(0), *gr_x2(0); // , *gr_x4(0);

    fCanvas  = new TCanvas(name,title,1,1,1100,800);

    fCanvas->cd(1);

    TH2F* h2 = new TH2F("h2","Muon Rejection Vs Electron Efficiency",500,0.8,1.0,1,10,2000);

    h2->GetXaxis()->SetTitle("Electron efficiency");
    h2->GetYaxis()->SetTitle("Muon Rejection");
    h2->SetStats(0);

    h2->Draw();

    gPad->SetLogy(1);
    gPad->SetGridy(1);

    double nt_1212  = gh1(f_m00s1212,"TrackAna","trk_13/ep")->GetEntries();

    double nev_1212 = gh1(f_m00s1212,"TrackAna","evt_0/ntrk")->GetEntries();

    double muon_reco_eff_1212 = nt_1212/nev_1212;

    create_llhr_rejection_graph(f_e00s1212,f_m00s1212,muon_reco_eff_1212,gr_x0);

    gr_x0->SetMarkerStyle(20);
    gr_x0->SetMarkerSize(1);
    gr_x0->Draw("LP");

    create_llhr_rejection_graph(f_e00s1412,f_m00s1412,muon_reco_eff_1212,gr_x1);

    gr_x1->SetMarkerStyle(24);
    gr_x1->SetMarkerSize(1);
    gr_x1->Draw("same,LP");

    create_llhr_rejection_graph(f_e00s1512,f_m00s1512,muon_reco_eff_1212,gr_x2);

    gr_x2->SetMarkerStyle(25);
    gr_x2->SetMarkerSize(1);
    gr_x2->Draw("same,LP");

    // create_llhr_rejection_graph(f_e00s1612,f_m00s1612,muon_reco_eff_1212,gr_1612);

    // gr_1612->SetMarkerStyle(26);
    // gr_1612->SetMarkerSize(1);
    // gr_1612->Draw("same,LP");

    TLegend* leg = new TLegend(0.15,0.25,0.4,0.4);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(gr_x0,"CE only"    ,"ep");
    leg->AddEntry(gr_x1,"CE+MIXP3-x1","ep");
    leg->AddEntry(gr_x2,"CE+MIXP3-x2","ep");
    //    leg->AddEntry(gr_1612,"CE+MIXP3-x4","ep");

    leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig. 102: results of the toy MC
//-----------------------------------------------------------------------------
  if (Figure == 102) {

    fCanvas  = new TCanvas(name,title,1,1,1200,800);

    fCanvas->cd(1);

    gInterpreter->LoadMacro("murat/plot/plot_lhr_rejection.C");
    gInterpreter->ProcessLine("plot_lhr_rejection(1412)");
  }
}

