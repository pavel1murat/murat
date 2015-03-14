///////////////////////////////////////////////////////////////////////////////
// plots for Mu2e-4375: Mu2e update on the Calorimeter-based PID
// used offline 4_2_4
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
// Fig.   5: E/P for electrons
// Fig.   6: 
// Fig.   7: Dt for electrons 
// Fig.   8: chi2 (tcm) selection for electrons 
// Fig.   9: E/P selection - CE  without background(e00s2012) and with the default background(e00s1412)
//
// Fig.  11: E/P for different occupancies and datasets : TRK_19
// Fig.  12: E/P for different occupancies and datasets : TRK_19
// Fig.  13: E/P  vs DT TRK_19
//
// Fig.  21: efficiency of the P preselection cut
// Fig.  22: efficiency of the E/P preselection cut 
// Fig.  23: efficiency of the Dt preselection cut
// Fig.  24: efficiency of the chi2(tcm) preselection cut
//
// Fig. 101: LLHR_CAL rejection vs luminosity
// Fig. 102: muon rejection vs electron efficiency, vary calorimeter performance 
//           parameters 
////////////////////////////////////////////////////////////////////////////////
#include "plot/TMu2e4375.hh"

#include "TPaveLabel.h"
#include "TArrow.h"
#include "TObjArray.h"
#include "TGraphErrors.h"
#include "TObjString.h"
#include "TPaveStats.h"

// #include <iostream>

ClassImp(TMu2e4375)

//_____________________________________________________________________________
TMu2e4375::TMu2e4375(int PlotMode, int BlessingMode): TPlotNote(PlotMode,BlessingMode) {

  const char* hist_dir;
  char         res_dir[500];
  
  fFiguresDir   = Form("%s/mu2e_4375",gEnv->GetValue("mu2e.Figures","."));
  fWorkDir      = gSystem->Getenv("WORK_DIR");

  sprintf(res_dir,"%s/results",fWorkDir.Data());
  hist_dir      = gEnv->GetValue("mu2e.HistDir",res_dir);

  sprintf(f_e00s0160,"%s/v4_2_6/e00s0160.track_ana.hist", hist_dir); 
  sprintf(f_e00s0161,"%s/v4_2_6/e00s0161.track_ana.hist", hist_dir); 
  sprintf(f_e00s0162,"%s/v4_2_6/e00s0162.track_ana.hist", hist_dir); 

  sprintf(f_m00s0160,"%s/v4_2_6/m00s0160.track_ana.hist", hist_dir); 
  sprintf(f_m00s0161,"%s/v4_2_6/m00s0161.track_ana.hist", hist_dir); 
  sprintf(f_m00s0162,"%s/v4_2_6/m00s0162.track_ana.hist", hist_dir); 

  sprintf(f_cnvs2012,"%s/v4_2_6/cnvs2012.track_ana.hist", hist_dir); 
  sprintf(f_cnvs2022,"%s/v4_2_6/cnvs2022.track_ana.hist", hist_dir); 
  sprintf(f_cnvs2032,"%s/v4_2_6/cnvs2032.track_ana.hist", hist_dir); 
  sprintf(f_cnvs2042,"%s/v4_2_6/cnvs2042.track_ana.hist", hist_dir); 

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
TMu2e4375::~TMu2e4375() {
};

//-----------------------------------------------------------------------------
const char* TMu2e4375::GetFiguresDir() {
  return fFiguresDir.Data();
}

void TMu2e4375::get_filename(int Figure, char* Filename) {

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
void TMu2e4375::remake_plots() {

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
// example: FnEle: "~/hist/mu2e/v4_2_1/e00s2012.track_ana.hist"
//          FnMuo: "~/hist/mu2e/v4_2_1/m00s2012.track_ana.hist"
// 'MuonRecoEff' - reconstruction efficiency at zero background occupancy
//---------------------------------------------------------------------------------------
void TMu2e4375::create_llhr_rejection_graph(const char* FnEle, 
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
// trk_13: Set C tracks in the window 100 < P < 110 MeV/c
//-----------------------------------------------------------------------------
  double muon_reco_eff, nt, nev, reco_eff_sf(1);

  nt  = gh1(FnMuo,"TrackAna","trk_13/ep")->GetEntries();
  nev = gh1(FnMuo,"TrackAna","evt_0/ntrk")->GetEntries();

  muon_reco_eff = nt/nev;

  reco_eff_sf   = MuonRecoEff/muon_reco_eff;

  printf(" muon reco_eff_sf : %10.5f\n",reco_eff_sf);

  double se = h_llhr_e->Integral();
  double sm = h_llhr_m->Integral();

  double qe, qm;

  for (int i=0; i<nbx; i++) {
    qe = h_llhr_e->Integral(1,i+1);
    qm = h_llhr_m->Integral(1,i+1);

    pe[i] =(1-qe/se)*(qne_25/qne_13);
    rm[i] = 1./(1-qm/sm + 1.e-6)*(qnm_13/qnm_25); // *reco_eff_sf;

    if (GetDebugBit(1) == 1) {
      printf(" i, llhr , qe, qm, prob(e) , rej(mu) : %3i %10.3f %10.3f %10.3f %10.5f %10.3f\n",
	    i,h_llhr_e->GetBinCenter(i+1),qe,qm,pe[i],rm[i]);
    }
  }

  if (Graph != NULL) delete Graph;

  Graph = new TGraph(nbx,pe,rm);
}


//-----------------------------------------------------------------------------
void TMu2e4375::plot(Int_t Figure, const char* CanvasName) {
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
    //    TH1F   *h1;

    // fCanvas  = new TCanvas(name,title,1,1,1100,800);

    // fCanvas->cd(1);

    // h1 = gh1(f_e00s2012,"TrackAna","trk_1/du");

    // h1->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  2: single electrons , DV
//-----------------------------------------------------------------------------
  if (Figure == 2) {
    //    TH1F   *h1;

    // fCanvas  = new TCanvas(name,title,1,1,1100,800);

    // fCanvas->cd(1);

    // h1 = gh1(f_e00s2012,"TrackAna","trk_1/dv");

    // h1->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  3: DU vs path
//-----------------------------------------------------------------------------
  if (Figure == 3) {
    //    TH1F   *h1;

    // fCanvas  = new TCanvas(name,title,1,1,1100,800);

    // fCanvas->cd(1);

    // h1 = gh1(f_e00s2012,"TrackAna","trk_1/du_vs_path");

    // h1->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  4: DV vs path
//-----------------------------------------------------------------------------
  if (Figure == 4) {
    //    TH2F   *h2;

    // fCanvas  = new TCanvas(name,title,1,1,1100,800);

    // fCanvas->cd(1);

    // h2 = gh2(f_e00s2012,"TrackAna","trk_1/dv_vs_path");

    // h2->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  5: DU vs path
//-----------------------------------------------------------------------------
  if (Figure == 5) {
    TH1F   *h1;

    fCanvas  = new TCanvas(name,title,1,1,1100,800);

    fCanvas->cd(1);

    h1 = gh1(f_e00s0160,"TrackAna","trk_19/ep");

    h1->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  6: DV vs path
//-----------------------------------------------------------------------------
  if (Figure == 6) {
    TH2F   *h2;

    fCanvas  = new TCanvas(name,title,1,1,1100,800);

    fCanvas->cd(1);

    h2 = gh2(f_e00s0160,"TrackAna","trk_1/dvc_vs_path");

    h2->ProjectionY()->Fit("gaus");
  }
//-----------------------------------------------------------------------------
// Fig.  7: Dt selection
//-----------------------------------------------------------------------------
  if (Figure == 7) {
    TH1F   *h_2012, *h_1442;

    fCanvas  = new TCanvas(name,title,1,1,1100,800);

    fCanvas->cd(1);

    h_2012 = gh1(f_cnvs2012,"TrackAna","trk_19/dt");
    h_1442 = gh1(f_e00s0160,"TrackAna","trk_19/dt");

    double qn_1442=h_1442->Integral();
    h_2012->SetLineColor(2);
    h_2012->SetFillColor(2);
    h_2012->SetFillStyle(3002);
    h_2012->DrawNormalized("",qn_1442);

    h_1442->Draw("same");

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

    leg->AddEntry(h_2012,"CE only"    ,"ep");
    leg->AddEntry(h_1442,"CE+MIXP3-x1","ep");

    leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  8: Chi2(TCM) selection
//-----------------------------------------------------------------------------
  if (Figure == 8) {
    TH1F   *h_1442, *h_2012;

    fCanvas  = new TCanvas(name,title,1,1,1100,800);
    gPad->SetLogy(1);

    fCanvas->cd(1);

    h_1442 = gh1(f_e00s0160,"TrackAna","trk_19/chi2tcm");
    h_2012 = gh1(f_cnvs2012,"TrackAna","trk_19/chi2tcm");

    double qn_1442=h_1442->Integral();
    h_2012->SetLineColor(2);
    h_2012->SetFillColor(2);
    h_2012->SetFillStyle(3002);
    h_2012->SetStats(0);

    h_2012->SetMinimum(0.01);
    h_2012->DrawNormalized("",qn_1442);

    h_1442->Draw("same");

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

    leg->AddEntry(h_1442,"CE only"    ,"f");
    leg->AddEntry(h_2012,"CE+MIXP3-x1","f");

    leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  9: E/P selection
//-----------------------------------------------------------------------------
  if (Figure == 9) {
    TH1F   *h_2012, *h_1442;

    fCanvas  = new TCanvas(name,title,1,1,1200,800);
    //    gPad->SetLogy(1);

    fCanvas->cd(1);

    h_2012 = gh1(f_cnvs2012,"TrackAna","trk_19/ep");
    h_1442 = gh1(f_e00s0160,"TrackAna","trk_19/ep");

    double qn_1442=h_1442->Integral();
    h_2012->SetLineColor(2);
    h_2012->SetFillColor(2);
    h_2012->SetFillStyle(3001);

    h_1442->Draw();
				// then plot low-statistics CE+MIXP plot

    h_2012->DrawNormalized("same",qn_1442);

				// plot the CE-only distribution again to have it 
				// on top of the filled histogram
    h_1442->Draw("same");

    TArrow  *a1; //, *a2;

    a1 = new TArrow(1.15,2000,1.15,100);
    a1->SetArrowSize(0.015);
    a1->Draw();

    // a2 = new TArrow(3.,100,3,10);
    // a2->SetArrowSize(0.015);
    // a2->Draw();

    TLegend* leg = new TLegend(0.15,0.65,0.3,0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(h_1442,"CE only"    ,"f");
    leg->AddEntry(h_2012,"CE+MIXP3-x1","f");

    leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  11: delta(T) for different datasets, TRK_19 (track+cluster)
//-----------------------------------------------------------------------------
  if (Figure == 11) {
    // TH1F   *he_2012, *hm_2012, *hm_1442, *hm_1512;

    // fCanvas  = new TCanvas(name,title,1,1,1200,800);
    // gPad->SetLogy(1);

    // fCanvas->cd(1);

    // he_2012 = gh1(f_e00s2012,"TrackAna","trk_19/dt");

    // hm_2012 = gh1(f_m00s2012,"TrackAna","trk_19/dt");
    // hm_1442 = gh1(f_m00s0160,"TrackAna","trk_19/dt");
    // hm_1512 = gh1(f_m00s1512,"TrackAna","trk_19/dt");

    // 				// drawing part
    // hm_1512->SetMaximum(0.12);
    // hm_1512->SetStats(0);
    // hm_1512->DrawNormalized("",1);

    // hm_1442->SetLineColor(2);
    // hm_1442->SetFillColor(2);
    // hm_1442->SetFillStyle(3002);
    // hm_1442->DrawNormalized("same",1);

    // hm_2012->SetFillColor(kBlue-7);
    // hm_2012->SetFillStyle(3001);
    // hm_2012->DrawNormalized("same",1);

    // he_2012->SetMarkerSize(1);
    // he_2012->SetMarkerStyle(20);
    // he_2012->DrawNormalized("same,p",1);

    // TLegend* leg = new TLegend(0.15,0.65,0.4,0.85);
    // leg->SetBorderSize(0);
    // leg->SetFillStyle(0);

    // leg->AddEntry(he_2012,"Electrons 105 MeV/c","p");
    // leg->AddEntry(hm_2012,"muons 105 MeV/c"    ,"f");
    // leg->AddEntry(hm_1442,"M105+MIXP3-x1","f");
    // leg->AddEntry(hm_1512,"M105+MIXP3-x2","f");

    // leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  12: E/P for different datasets, TRK_25 (track+cluster)
//-----------------------------------------------------------------------------
  if (Figure == 12) {
    // TH1F   *he_2012, *hm_2012, *hm_1442, *hm_1512;

    // fCanvas  = new TCanvas(name,title,1,1,1200,800);
    // gPad->SetLogy(1);

    // fCanvas->cd(1);

    // he_2012 = gh1(f_e00s2012,"TrackAna","trk_25/ep");
    // hm_2012 = gh1(f_m00s2012,"TrackAna","trk_25/ep");

    // hm_1442 = gh1(f_m00s0160,"TrackAna","trk_25/ep");


    // hm_1512 = gh1(f_m00s1512,"TrackAna","trk_25/ep");

    // hm_1512->SetMaximum(0.12);
    // hm_1512->SetStats(0);
    // hm_1512->DrawNormalized("",1);


    // hm_1442->SetLineColor(2);
    // hm_1442->SetFillColor(2);
    // hm_1442->SetFillStyle(3002);
    // hm_1442->DrawNormalized("same",1);

    // hm_2012->SetFillColor(kBlue-7);
    // hm_2012->SetFillStyle(3001);
    // hm_2012->DrawNormalized("same",1);

    // he_2012->SetMarkerSize(1);
    // he_2012->SetMarkerStyle(20);
    // he_2012->DrawNormalized("same,p",1);

    // TLegend* leg = new TLegend(0.15,0.65,0.4,0.85);
    // leg->SetBorderSize(0);
    // leg->SetFillStyle(0);

    // leg->AddEntry(he_2012,"Electrons 105 MeV/c","p");
    // leg->AddEntry(hm_2012,"muons 105 MeV/c"    ,"f");
    // leg->AddEntry(hm_1442,"M105+MIXP3-x1","f");
    // leg->AddEntry(hm_1512,"M105+MIXP3-x2","f");

    // leg->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  21: efficiency of P preselection cut
//-----------------------------------------------------------------------------
  if (Figure == 21) {
    TH1F  *h1, *h2;

    fCanvas  = new TCanvas(name,title,1,1,1200,700);

    h1 = gh1(f_e00s0161,"TrackAna","trk_1/p_2");
    h2 = gh1(f_e00s0161,"TrackAna","trk_13/p_2");

    fCanvas->cd();
    gPad->SetLogy(1);

    h2->SetLineColor(kBlue-6);
    h2->SetFillColor(kBlue-6);
    h2->SetFillStyle(3001);

    h1->GetXaxis()->SetRangeUser(60,119.999);
    h1->GetXaxis()->SetTitle("P, MeV/c");
    h1->Draw();
    
    h2->Draw("sames");

    fCanvas->Modified();
    fCanvas->Update();
  
    TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

    s1->SetY1NDC(0.65);
    s1->SetY2NDC(0.90);
    s1->Draw();

    TPaveStats* s2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");
    
    s2->SetY1NDC(0.40);
    s2->SetY2NDC(0.65);
    s2->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  22: efficiency of E/P preselection cut
//-----------------------------------------------------------------------------
  if (Figure == 22) {
    TH1F  *h1, *h2;

    fCanvas  = new TCanvas(name,title,1,1,1200,700);

    h1 = gh1(f_e00s0161,"TrackAna","trk_13/ep");
    h2 = gh1(f_e00s0161,"TrackAna","trk_29/ep");

    fCanvas->cd();
    gPad->SetLogy(1);

    h2->SetLineColor(kBlue-6);
    h2->SetFillColor(kBlue-6);
    h2->SetFillStyle(3001);

    //    h1->GetXaxis()->SetRangeUser(60,119.999);
    h1->GetXaxis()->SetTitle("E/P");
    h1->Draw();
    
    h2->Draw("sames");

    fCanvas->Modified();
    fCanvas->Update();
  
    TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

    s1->SetY1NDC(0.65);
    s1->SetY2NDC(0.90);
    s1->Draw();

    TPaveStats* s2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");
    
    s2->SetY1NDC(0.40);
    s2->SetY2NDC(0.65);
    s2->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  23: efficiency of DT preselection cut
//-----------------------------------------------------------------------------
  if (Figure == 23) {
    TH1F  *h1, *h2;

    fCanvas  = new TCanvas(name,title,1,1,1200,700);

    h1 = gh1(f_e00s0161,"TrackAna","trk_29/dt");
    h2 = gh1(f_e00s0161,"TrackAna","trk_32/dt");

    fCanvas->cd();
    gPad->SetLogy(1);

    h2->SetLineColor(kBlue-6);
    h2->SetFillColor(kBlue-6);
    h2->SetFillStyle(3001);

    //    h1->GetXaxis()->SetRangeUser(60,119.999);
    h1->GetXaxis()->SetTitle("Delta T, ns");
    h1->Draw();
    
    h2->Draw("sames");

    fCanvas->Modified();
    fCanvas->Update();
  
    TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

    s1->SetY1NDC(0.65);
    s1->SetY2NDC(0.90);
    s1->Draw();

    TPaveStats* s2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");
    
    s2->SetY1NDC(0.40);
    s2->SetY2NDC(0.65);
    s2->Draw();
  }
//-----------------------------------------------------------------------------
// Fig.  24: efficiency of chi2(Track-Cluster Match) preselection cut
//-----------------------------------------------------------------------------
  if (Figure == 24) {
    TH1F  *h1, *h2;

    fCanvas  = new TCanvas(name,title,1,1,1200,700);

    h1 = gh1(f_e00s0161,"TrackAna","trk_32/chi2tcm");
    h2 = gh1(f_e00s0161,"TrackAna","trk_25/chi2tcm");

    fCanvas->cd();
    gPad->SetLogy(1);

    h2->SetLineColor(kBlue-6);
    h2->SetFillColor(kBlue-6);
    h2->SetFillStyle(3001);

    //    h1->GetXaxis()->SetRangeUser(60,119.999);
    h1->GetXaxis()->SetTitle("#chi^{2}(track-cluster match)");
    h1->Draw();
    
    h2->Draw("sames");

    fCanvas->Modified();
    fCanvas->Update();
  
    TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

    s1->SetY1NDC(0.65);
    s1->SetY2NDC(0.90);
    s1->Draw();

    TPaveStats* s2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");
    
    s2->SetY1NDC(0.40);
    s2->SetY2NDC(0.65);
    s2->Draw();
  }
//-----------------------------------------------------------------------------
// Fig. 101: muon rejection vs electron efficiency, different background levels
//-----------------------------------------------------------------------------
  if (Figure == 101) {
    TGraph   *gr_x0(0), *gr_x1(0), *gr_x2(0); // *gr_1612(0);

    fCanvas  = new TCanvas(name,title,1,1,1100,800);

    fCanvas->cd(1);

    TH2F* h2 = new TH2F("h2","Muon Rejection Vs Electron Efficiency",500,0.8,1.0,1,10,2000);
    
    h2->GetXaxis()->SetTitle("Electron efficiency");
    h2->GetYaxis()->SetTitle("Muon Rejection");
    h2->SetStats(0);
    
    h2->Draw();
    
    gPad->SetLogy(1);
    gPad->SetGridy(1);
    
    double nt_x0  = gh1(f_m00s0160,"TrackAna","trk_13/ep")->GetEntries();
    
    double nev_x0 = gh1(f_m00s0160,"TrackAna","evt_0/ntrk")->GetEntries();
    
    double muon_reco_eff_x0 = nt_x0/nev_x0;
    
    create_llhr_rejection_graph(f_e00s0160,f_m00s0160,muon_reco_eff_x0,gr_x0);
    
    gr_x0->SetMarkerStyle(20);
    gr_x0->SetMarkerSize(1);
    gr_x0->Draw("LP");
    
    create_llhr_rejection_graph(f_e00s0161,f_m00s0161,muon_reco_eff_x0,gr_x1);
    
    gr_x1->SetMarkerStyle(24);
    gr_x1->SetMarkerSize(1);
    gr_x1->Draw("same,LP");
    
    create_llhr_rejection_graph(f_e00s0162,f_m00s0162,muon_reco_eff_x0,gr_x2);
    
    gr_x2->SetMarkerStyle(25);
    gr_x2->SetMarkerSize(1);
    gr_x2->Draw("same,LP");

    // create_llhr_rejection_graph(f_e00s1612,f_m00s1612,muon_reco_eff_2012,gr_1612);

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
    gInterpreter->ProcessLine("plot_lhr_rejection(0041)");
  }
}

