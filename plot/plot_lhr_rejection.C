///////////////////////////////////////////////////////////////////////////////
// 2014-06-11: plot result produced by murat/ana/lhr_rejection.cc
///////////////////////////////////////////////////////////////////////////////

#include "murat/scripts/datasets.hh"

//-----------------------------------------------------------------------------
// templates tell the dataset name: 1212 - use e00s1212 and m00s1212 (no background added)
//                                  1412 - use e00s1412 and m00s1412
//
// plot results by the toy simulation
//-----------------------------------------------------------------------------
void plot_lhr_rejection(int Templates = 1212) {

  float  sigt [5] = { 0.05,  0.10, 0.20, 0.50, 1.00};

  float  ex [5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  float  err[5] = {0.002, 0.002, 0.002, 0.002, 0.002};

  float  eff_002[5] = { 0.995, 0.995, 0.994 , 0.991, 0.984 };
  float  eff_005[5] = { 0.993, 0.992, 0.990 , 0.988, 0.978 };
  float  eff_010[5] = { 0.988, 0.989, 0.987 , 0.980, 0.963 };
  float  eff_015[5] = { 0.980, 0.980, 0.980 , 0.965, 0.937 };
  float  eff_020[5] = { 0.965, 0.964, 0.960 , 0.937, 0.886 };

  // float  eff_002_1412[5] = { 0.989, 0.988, 0.988 , 0.986, 0.976 };
  // float  eff_005_1412[5] = { 0.986, 0.985, 0.984 , 0.980, 0.970 };
  // float  eff_010_1412[5] = { 0.976, 0.975, 0.974 , 0.969, 0.954 };
  // float  eff_015_1412[5] = { 0.960, 0.960, 0.957 , 0.945, 0.917 };
  // float  eff_020_1412[5] = { 0.929, 0.928, 0.928 , 0.909, 0.855 };
//-----------------------------------------------------------------------------
// 2014-06-20: absolute numbers
//-----------------------------------------------------------------------------
  float  eff_002_1412[5] = { 0.970, 0.970, 0.969 , 0.966, 0.961 };
  float  eff_005_1412[5] = { 0.968, 0.967, 0.967 , 0.963, 0.956 };
  float  eff_010_1412[5] = { 0.963, 0.962, 0.961 , 0.954, 0.943 };
  float  eff_015_1412[5] = { 0.955, 0.955, 0.952 , 0.939, 0.917 };
  float  eff_020_1412[5] = { 0.942, 0.940, 0.935 , 0.914, 0.870 };
//-----------------------------------------------------------------------------
// 2014-07-31: absolute numbers for '0041' datasets - v4_2_4 105 MeV/c 
//             electron and muon datasets with the standard background mix
//             the numbers are by 0.5-1% lower than the previos version 
//             what does it mean? 
//-----------------------------------------------------------------------------
  float  eff_002_0041[5] = { 0.964, 0.964, 0.963 , 0.960, 0.950 };
  float  eff_005_0041[5] = { 0.962, 0.961, 0.960 , 0.955, 0.943 };
  float  eff_010_0041[5] = { 0.957, 0.956, 0.954 , 0.944, 0.923 };
  float  eff_015_0041[5] = { 0.947, 0.946, 0.942 , 0.923, 0.883 };
  float  eff_020_0041[5] = { 0.930, 0.930, 0.921 , 0.888, 0.823 };

  TGraphErrors *gr_002, *gr_005, *gr_010, *gr_015, *gr_020;

  int np = 5;

  TH2F* h2 = new TH2F("h2","Electron efficiency for muon rejection of 200",1,0,1.1,1,0.80,1.0 );

  h2->SetStats(0);
  h2->GetXaxis()->SetTitle("#sigma_{T}, ns");
  h2->Draw();

  if (Templates == 1212) {
    gr_002 = new TGraphErrors(np, sigt,eff_002,ex,err);
    gr_005 = new TGraphErrors(np, sigt,eff_005,ex,err);
    gr_010 = new TGraphErrors(np, sigt,eff_010,ex,err);
    gr_015 = new TGraphErrors(np, sigt,eff_015,ex,err);
    gr_020 = new TGraphErrors(np, sigt,eff_020,ex,err);
  }
  else if (Templates == 1412) {
    gr_002 = new TGraphErrors(np, sigt,eff_002_1412,ex,err);
    gr_005 = new TGraphErrors(np, sigt,eff_005_1412,ex,err);
    gr_010 = new TGraphErrors(np, sigt,eff_010_1412,ex,err);
    gr_015 = new TGraphErrors(np, sigt,eff_015_1412,ex,err);
    gr_020 = new TGraphErrors(np, sigt,eff_020_1412,ex,err);
  }
  else if (Templates == 0041) {
    gr_002 = new TGraphErrors(np, sigt,eff_002_0041,ex,err);
    gr_005 = new TGraphErrors(np, sigt,eff_005_0041,ex,err);
    gr_010 = new TGraphErrors(np, sigt,eff_010_0041,ex,err);
    gr_015 = new TGraphErrors(np, sigt,eff_015_0041,ex,err);
    gr_020 = new TGraphErrors(np, sigt,eff_020_0041,ex,err);
  }

  gr_002->SetMarkerSize (1);
  gr_002->SetMarkerStyle(20);
  gr_002->Draw("same,LP");

  gr_005->SetMarkerSize (1);
  gr_005->SetMarkerStyle(21);
  gr_005->Draw("same,LP");

  gr_010->SetMarkerSize (1);
  gr_010->SetMarkerStyle(22);
  gr_010->Draw("same,LP");

  gr_015->SetMarkerSize (1);
  gr_015->SetMarkerStyle(23);
  gr_015->Draw("same,LP");

  gr_020->SetMarkerSize (1);
  gr_020->SetMarkerStyle(24);
  gr_020->Draw("same,LP");

  TLegend* leg = new TLegend(0.15,0.15,0.4,0.4);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);


  leg->AddEntry(gr_002,"#sigma_{E}/E = 0.02","ep");
  leg->AddEntry(gr_005,"#sigma_{E}/E = 0.05","ep");
  leg->AddEntry(gr_010,"#sigma_{E}/E = 0.10","ep");
  leg->AddEntry(gr_015,"#sigma_{E}/E = 0.15","ep");
  leg->AddEntry(gr_020,"#sigma_{E}/E = 0.20","ep");

  leg->Draw();
}


//-----------------------------------------------------------------------------
// plot muon rejection vs electron efficiency for different datasets
// it looks that requiring P > 100 helps - yes, it rejects the DIO electrons !
// need to require uniformly 
// ultimately need to switch from TRK_13 to TRK_29 for normalization
//-----------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
// example: FnEle: "~/hist/mu2e/v4_2_1/e00s1212.track_ana.hist"
//          FnMuo: "~/hist/mu2e/v4_2_1/m00s1212.track_ana.hist"
// Print = 1: print few values around LLHR = 0 to determine the operational point
//---------------------------------------------------------------------------------------
void create_llhr_cal_rejection_graph(const dataset_t* DsEle,
				     const dataset_t* DsMuo,
				     TGraph*&         Graph,
				     int              Print = 0) {

  TH1* h_llhr_e = (TH1*) gh1(DsEle->fn_track_ana,"TrackAna","trk_25/llhr_cal")->Clone("h_llhr_e");
  TH1* h_llhr_m = (TH1*) gh1(DsMuo->fn_track_ana,"TrackAna","trk_25/llhr_cal")->Clone("h_llhr_m");

  double pe[1000], rm[1000];

  double qnm_13, qnm_25, qne_13, qne_25;

  int nbx = h_llhr_e->GetNbinsX();

  qne_25 = gh1(DsEle->fn_track_ana,"TrackAna","trk_25/llhr_cal")->GetEntries();
  qne_13 = gh1(DsEle->fn_track_ana,"TrackAna","trk_13/llhr_cal")->GetEntries();

  qnm_25 = gh1(DsMuo->fn_track_ana,"TrackAna","trk_25/llhr_cal")->GetEntries();
  qnm_13 = gh1(DsMuo->fn_track_ana,"TrackAna","trk_13/llhr_cal")->GetEntries();

  double se = h_llhr_e->Integral();
  double sm = h_llhr_m->Integral();

  double qe, qm;

  if (Print != 0) printf(" ------------------- print efficiency and rejection numbers for %s and %s \n",
			 DsEle->name,DsMuo->name);
  for (int i=0; i<nbx; i++) {
    qe = h_llhr_e->Integral(1,i+1);
    qm = h_llhr_m->Integral(1,i+1);

    pe[i] =(1-qe/se)*(qne_25/qne_13);
    rm[i] = 1./(1-qm/sm + 1.e-6)*(qnm_13/qnm_25);

    if ((Print != 0) && (pe[i] > 0.85) && (rm[i] > 10)) {
      printf(" i, llhr , qe, qm, prob(e) , rej(mu) : %3i %10.3f %10.3f %10.3f %10.5f %10.3f\n",
	     i,h_llhr_e->GetBinCenter(i+1),qe,qm,pe[i],rm[i]);
    }
  }

  if (Graph != NULL) delete Graph;

  Graph = new TGraph(nbx,pe,rm);
}


//-----------------------------------------------------------------------------
void plot_llhr_cal_rejection_2(int OffVer = 421, int Print = 0) {

  TGraph  *gr_x0(0), *gr_x1(0);

  dataset_t   *ele_x0, *muo_x0, *ele_x1, *muo_x1;

  char fn_ele_x0[200], fn_muo_x0[200], fn_ele_x1[200], fn_muo_x1[200];

  if (OffVer == 421) {
    strcpy(fn_ele_x0,"~/hist/mu2e/v4_2_1/e00s1212.track_ana.hist");
    strcpy(fn_muo_x0,"~/hist/mu2e/v4_2_1/m00s1212.track_ana.hist");

    strcpy(fn_ele_x1,"~/hist/mu2e/v4_2_1/e00s1412.track_ana.hist");
    strcpy(fn_muo_x1,"~/hist/mu2e/v4_2_1/m00s1412.track_ana.hist");
  }
  else if (OffVer == 572) {
    ele_x0 = &e40s5720;
    muo_x0 = &m40s5720;
    ele_x1 = &e42s5721;
    muo_x1 = &m40s5721;
  }

  printf(" OffVer = %3i fn_ele: %s fn_muo: %s\n",
	 OffVer,ele_x0->fn_track_ana,muo_x0->fn_track_ana);

  TH2F* h2 = new TH2F("h2","Muon Rejection Vs Electron Efficiency",1000,0.8,1.0,1,10,5000);

  h2->GetXaxis()->SetTitle("Electron efficiency");
  h2->GetYaxis()->SetTitle("Muon Rejection");
  h2->SetStats(0);

  h2->Draw();

  gPad->SetLogy(1);
  gPad->SetGridy(1);

  create_llhr_cal_rejection_graph(ele_x0,muo_x0,gr_x0,Print);

  gr_x0->SetMarkerStyle(20);
  gr_x0->SetMarkerSize(1);
  gr_x0->Draw("LP");

  create_llhr_cal_rejection_graph(ele_x1,muo_x1,gr_x1,Print);

  gr_x1->SetMarkerStyle(24);
  gr_x1->SetMarkerSize(1);
  gr_x1->Draw("same,LP");

  TLegend* leg = new TLegend(0.15,0.25,0.4,0.4);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  leg->AddEntry(gr_x0,"Signal only"           ,"ep");
  leg->AddEntry(gr_x1,"Signal + Overlays (x1)","ep");

  leg->Draw();
}
