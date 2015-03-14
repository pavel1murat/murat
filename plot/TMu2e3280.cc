///////////////////////////////////////////////////////////////////////////////
// plots for Mu2e-3280: calorimeter occupancy studies
// background indexing: [0]:DIO, [1]:EPR, [2]:ENU, [3]:EPH
// -------------------------
// directory with the histogram files: by default: $WORK_DIR/results
// can be redefined with "zzx.HistDir" in .rootrc
// 
// figures:
// ---------
// Fig.  1: DIO energy             (ibgr=0)
// Fig.  2: proton energy          (ibgr=1)
// Fig.  3: neutron energy         (ibgr=2)
// Fig.  4: photon energy          (ibgr=3)
// Fig.  5: DIO number of hits     (ibgr=0)
// Fig.  6: proton number of hits  (ibgr=1)
// Fig.  7: neutron number of hits (ibgr=2)
// Fig.  8: photon number of hits  (ibgr=3)
// ... 9-20
// Fig. 21: single hit energy distributions for all 4 background sources
// Fig. 22: single Hit energy distributions for all 4 background sources, E > 0.1 MeV 
//
// Fig. 23: total occupancies vs R : E>0 vs E>1 MeV
// Fig. 24: total occupancies: new PA vs old PA
// Fig. 25: total occupancies: new INA vs no INA
// Fig. 26: total occupancies: old PA vs new PA for hit E > 1 MeV : 
// Fig. 27: total occupancies, E > 1 MeV, new PA, NO_NABS vs NEW_NABS
//
// Fig. 31: DIO energy     (ibgr = 0), NEW_PA, NO_NABS vs NEW_NABS
// Fig. 32: proton  energy (ibgr = 1), NEW_PA, NO_NABS vs NEW_NABS
// Fig. 33: neutron energy (ibgr = 2), NEW_PA, NO_NABS vs NEW_NABS
// Fig. 34: photon  energy (ibgr = 3), NEW_PA, NO_NABS vs NEW_NABS
// Fig. 35: DIO     nhits  (ibgr = 0), NEW_PA, NO_NABS vs NEW_NABS
// Fig. 36: proton  nhits  (ibgr = 1), NEW_PA, NO_NABS vs NEW_NABS
// Fig. 37: neutron nhits  (ibgr = 2), NEW_PA, NO_NABS vs NEW_NABS
// Fig. 38: photon  nhits  (ibgr = 3), NEW_PA, NO_NABS vs NEW_NABS
// Fig. 39: DIO     nhits  (ibgr = 0), NEW_PA, NO_NABS vs NEW_NABS, T>700 ns 
// Fig. 40: proton  nhits  (ibgr = 1), NEW_PA, NO_NABS vs NEW_NABS, T>700 ns 
// Fig. 41: neutron nhits  (ibgr = 2), NEW_PA, NO_NABS vs NEW_NABS, T>700 ns 
// Fig. 42: photon  nhits  (ibgr = 3), NEW_PA, NO_NABS vs NEW_NABS, T>700 ns 
// Fig. 43: DIO     nhits  (ibgr = 0), NEW_PA, NO_NABS vs NEW_NABS, E>0.1 MeV
// Fig. 44: proton  nhits  (ibgr = 1), NEW_PA, NO_NABS vs NEW_NABS, E>0.1 MeV
// Fig. 45: neutron nhits  (ibgr = 2), NEW_PA, NO_NABS vs NEW_NABS, E>0.1 MeV
// Fig. 46: photon  nhits  (ibgr = 3), NEW_PA, NO_NABS vs NEW_NABS, E>0.1 MeV
// Fig. 47: DIO     nhits  (ibgr = 0), NEW_PA, NO_NABS vs NEW_NABS, E>1.0 MeV
// Fig. 48: proton  nhits  (ibgr = 1), NEW_PA, NO_NABS vs NEW_NABS, E>1.0 MeV
// Fig. 49: neutron nhits  (ibgr = 2), NEW_PA, NO_NABS vs NEW_NABS, E>1.0 MeV
// Fig. 50: photon  nhits  (ibgr = 3), NEW_PA, NO_NABS vs NEW_NABS, E>1.0 MeV
//
// Fig.102: emitted protons  301         : ZR dist of the origins of particles hitting the calorimeter
// Fig.103: emitted neutrons 301         : ZR dist of the origins of particles hitting the calorimeter
// Fig.112: emitted protons  dev         : ZR dist of the origins of particles hitting the calorimeter
// Fig.113: emitted neutrons dev         : ZR dist of the origins of particles hitting the calorimeter
// Fig.123: emitted neutrons dev/no_nabs : ZR dist of the origins of particles hitting the calorimeter
//
// Fig.201: protons - ratio of the number of hits - NEW_PABS/OLD_PABS for disk#1 vs radius
// Fig.202: protons - ratio of the number of hits - NEW_PABS/OLD_PABS for disk#1 vs energy
// Fig.203: protons - PDG code  NEW_PABS/OLD_PABS 
// Fig.204: protons - number of straw hits 9NEW_PABS vs OLD_PABS)
// 
////////////////////////////////////////////////////////////////////////////////
#include "plot/TMu2e3280.hh"

#include "TPaveLabel.h"
#include "TArrow.h"
#include "TObjArray.h"
#include "TGraphErrors.h"
#include "TObjString.h"
#include <iostream>

ClassImp(TMu2e3280)

//_____________________________________________________________________________
TMu2e3280::TMu2e3280(int PlotMode, int BlessingMode): TPlotNote(PlotMode,BlessingMode) {

  const char* hist_dir;
  char         res_dir[500];
  
  fFiguresDir   = Form("%s/mu2e_3280",gEnv->GetValue("mu2e.Figures","."));
  fWorkDir      = gSystem->Getenv("WORK_DIR");

  sprintf(res_dir,"%s/results",fWorkDir.Data());
  hist_dir      = gEnv->GetValue("mu2e.HistDir",res_dir);

  bgr_tcalm002_301[0] = Form("%s/v3_0_1_cal/dio_000_105_001_disk_670_330_700_Mau8_tcalm002.hist", hist_dir); 
  bgr_tcalm002_301[1] = Form("%s/v3_0_1_cal/epr_000_300_001_disk_670_330_700_Mau8_tcalm002.hist", hist_dir); 
  bgr_tcalm002_301[2] = Form("%s/v3_0_1_cal/enu_000_100_001_disk_670_330_700_Mau8_tcalm002.hist", hist_dir); 
  bgr_tcalm002_301[3] = Form("%s/v3_0_1_cal/eph_000_007_001_disk_670_330_700_Mau8_tcalm002.hist", hist_dir); 

  bgr_tcalm002_dev[0] = Form("%s/v3_0_1_plus/dio_000_105_001_disk_670_330_700_Mau8_tcalm002.hist", hist_dir); 
  bgr_tcalm002_dev[1] = Form("%s/v3_0_1_plus/epr_000_300_001_disk_670_330_700_Mau8_tcalm002.hist", hist_dir); 
  bgr_tcalm002_dev[2] = Form("%s/v3_0_1_plus/enu_000_100_001_disk_670_330_700_Mau8_tcalm002.hist", hist_dir); 
  bgr_tcalm002_dev[3] = Form("%s/v3_0_1_plus/eph_000_007_001_disk_670_330_700_Mau8_tcalm002.hist", hist_dir); 

  bgr_tcalm002_nab[0] = Form("%s/v3_0_1_plus/dio_000_105_001_no_nabs_disk_670_330_700_Mau8_tcalm002.hist", hist_dir); 
  bgr_tcalm002_nab[1] = Form("%s/v3_0_1_plus/epr_000_300_001_no_nabs_disk_670_330_700_Mau8_tcalm002.hist", hist_dir); 
  bgr_tcalm002_nab[2] = Form("%s/v3_0_1_plus/enu_000_100_001_no_nabs_disk_670_330_700_Mau8_tcalm002.hist", hist_dir); 
  bgr_tcalm002_nab[3] = Form("%s/v3_0_1_plus/eph_000_007_001_no_nabs_disk_670_330_700_Mau8_tcalm002.hist", hist_dir); 

  bgr_tcalm005_301[0] = Form("%s/v3_0_1_cal/dio_000_105_001_disk_670_330_700_Mau8_tcalm005.hist", hist_dir); 
  bgr_tcalm005_301[1] = Form("%s/v3_0_1_cal/epr_000_300_001_disk_670_330_700_Mau8_tcalm005.hist", hist_dir); 
  bgr_tcalm005_301[2] = Form("%s/v3_0_1_cal/enu_000_100_001_disk_670_330_700_Mau8_tcalm005.hist", hist_dir); 
  bgr_tcalm005_301[3] = Form("%s/v3_0_1_cal/eph_000_007_001_disk_670_330_700_Mau8_tcalm005.hist", hist_dir); 

  bgr_tcalm005_dev[0] = Form("%s/v3_0_1_plus/dio_000_105_001_disk_670_330_700_Mau8_tcalm005.hist", hist_dir); 
  bgr_tcalm005_dev[1] = Form("%s/v3_0_1_plus/epr_000_300_001_disk_670_330_700_Mau8_tcalm005.hist", hist_dir); 
  bgr_tcalm005_dev[2] = Form("%s/v3_0_1_plus/enu_000_100_001_disk_670_330_700_Mau8_tcalm005.hist", hist_dir); 
  bgr_tcalm005_dev[3] = Form("%s/v3_0_1_plus/eph_000_007_001_disk_670_330_700_Mau8_tcalm005.hist", hist_dir); 

  bgr_tcalm005_nab[0] = Form("%s/v3_0_1_plus/dio_000_105_001_no_nabs_disk_670_330_700_Mau8_tcalm005.hist", hist_dir); 
  bgr_tcalm005_nab[1] = Form("%s/v3_0_1_plus/epr_000_300_001_no_nabs_disk_670_330_700_Mau8_tcalm005.hist", hist_dir); 
  bgr_tcalm005_nab[2] = Form("%s/v3_0_1_plus/enu_000_100_001_no_nabs_disk_670_330_700_Mau8_tcalm005.hist", hist_dir); 
  bgr_tcalm005_nab[3] = Form("%s/v3_0_1_plus/eph_000_007_001_no_nabs_disk_670_330_700_Mau8_tcalm005.hist", hist_dir); 

  sprintf(res_dir,"%s",fWorkDir.Data());
  hist_dir      = gEnv->GetValue("mu2e.HistDir",res_dir);
//-----------------------------------------------------------------------------
// different settings (colors etc) for talk, note and paper modes
//-----------------------------------------------------------------------------
  // if (PlotMode == kNoteMode) {
  // }
  // else if (PlotMode == kTalkMode) {
  // }
  // else if (PlotMode == kPaperMode) {
  // }

  fQExp[0]    = 19331.;  // DIO
  fQExp[1]    =  3117.;  // protons
  fQExp[2]    = 37407.;  // neutrons
  fQExp[3]    = 62343.;  // photons

  fQGen301[0] =  200000.;
  fQGen301[1] =  100000.;
  fQGen301[2] = 1000000.;
  fQGen301[3] = 1000000.;

  fQGenDev[0] =  200000.;
  fQGenDev[1] =  100000.;
  fQGenDev[2] = 1000000.;
  fQGenDev[3] = 1000000.;

  fQGenNab[0] =  200000.;
  fQGenNab[1] =  100000.;
  fQGenNab[2] = 1000000.;
  fQGenNab[3] = 1000000.;

  fBgrName[0] = "DIO";
  fBgrName[1] = "Protons";
  fBgrName[2] = "Neutrons";
  fBgrName[3] = "Photons";

};

//_____________________________________________________________________________
TMu2e3280::~TMu2e3280() {
};

//-----------------------------------------------------------------------------
const char* TMu2e3280::GetFiguresDir() {
  return fFiguresDir.Data();
}

void TMu2e3280::get_filename(int Figure, char* Filename) {

  if      (Figure ==   1) sprintf(Filename,"fig_%03i_%s",Figure,"dio_energy");
  else if (Figure ==   2) sprintf(Filename,"fig_%03i_%s",Figure,"prt_energy");
  else if (Figure ==   3) sprintf(Filename,"fig_%03i_%s",Figure,"ntr_energy");
  else if (Figure ==   4) sprintf(Filename,"fig_%03i_%s",Figure,"pht_energy");
  else if (Figure ==   5) sprintf(Filename,"fig_%03i_%s",Figure,"dio_nhits");
  else if (Figure ==   6) sprintf(Filename,"fig_%03i_%s",Figure,"prt_nhits");
  else if (Figure ==   7) sprintf(Filename,"fig_%03i_%s",Figure,"ntr_nhits");
  else if (Figure ==   8) sprintf(Filename,"fig_%03i_%s",Figure,"pht_nhits");
  else if (Figure ==  28) sprintf(Filename,"fig_%03i_%s",Figure,"cal_1_hit_time");
  else if (Figure == 103) sprintf(Filename,"fig_%03i_%s",Figure,"r0z0_enu_301");
  else if (Figure == 113) sprintf(Filename,"fig_%03i_%s",Figure,"r0z0_enu_cd2");
  else if (Figure == 204) sprintf(Filename,"fig_%03i_%s",Figure,"n_straw_hits");
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
void TMu2e3280::remake_plots() {

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
// example: occupancy(bgr_tcalm002_dev, ..., kDio,"cal_1/energy","cal_1/energy")
//-----------------------------------------------------------------------------
TH1F* TMu2e3280::occupancy(TString HistSet[], const char* HistName, int IBgr, double* NGen, 
			   const char* HistRE, const char* HistR) {

  TH1F  *h1, *h_occ(0), *he;

  h1     = gh1(HistSet[IBgr],"TCalm002",HistRE);
  h_occ  = (TH1F*) h1->Clone(HistName);
  h_occ->Divide(h1,fHistRC[0],fQExp[IBgr]/NGen[IBgr],1.);
  
  he = gh1(HistSet[IBgr],"TCalm002",HistR);

  int nb = he->GetNbinsX();
  
  double qn, err;

  for (int i=1; i<=nb; i++) {
    qn  = he->GetBinContent(i);
    err = h_occ->GetBinContent(i)/sqrt(qn+1.e-12);
    
    h_occ->SetBinError(i,err);
  }

  return h_occ;
}



//-----------------------------------------------------------------------------
double TMu2e3280::q_background(TString* HistSet, int IBgr) {

  if      (HistSet == bgr_tcalm002_301) return fQGen301[IBgr];
  else if (HistSet == bgr_tcalm002_dev) return fQGenDev[IBgr];
  else if (HistSet == bgr_tcalm002_nab) return fQGenNab[IBgr];
  else {
    printf(">>> TMu2e3280::q_background ERROR : undefined background\n");
    return 0.;
  }

}

//-----------------------------------------------------------------------------
// plot energy per crystal vs R for various background sources
//-----------------------------------------------------------------------------
void TMu2e3280::plot_energy_vs_r(TPad* Pad, TString* HistSet1, TString* HistSet2, 
				 int IBgr, const char* HistRE, const char* HistR, double EMax) {
  TH1F  *he_301, *he_dev;
//-----------------------------------------------------------------------------
// first disk, first - v3_0_1, then - development
//-----------------------------------------------------------------------------
  Pad->cd(1);

  double qbgr1 = q_background(HistSet1,IBgr);
  double qbgr2 = q_background(HistSet2,IBgr);

  fH1[10] = gh1(HistSet1[IBgr],"TCalm002",Form("%s_0",HistRE));
  fH1[11] = (TH1F*) fH1[10]->Clone("occupancy");
  fH1[11]->Divide(fH1[10],fHistRC[0],fQExp[IBgr]/qbgr1,1.);
  
  he_301 = gh1(HistSet1[IBgr],"TCalm002",Form("%s_0",HistR));

  int nb = he_301->GetNbinsX();
  
  double qn, err;

  for (int i=1; i<=nb; i++) {
    qn  = he_301->GetBinContent(i);
    //    ncr = fHistRC[0]->GetBinContent(i);
    err = fH1[11]->GetBinContent(i)/sqrt(qn+1.e-12);
    
    fH1[11]->SetBinError(i,err);
  }
					// development: 

  fH1[12] = gh1(HistSet2[IBgr],"TCalm002",Form("%s_0",HistRE));
  fH1[13] = (TH1F*) fH1[12]->Clone("occupancy");
  fH1[13]->Divide(fH1[12],fHistRC[0],fQExp[IBgr]/qbgr2,1.);

  he_dev = gh1(HistSet2[IBgr],"TCalm002",Form("%s_0",HistR));

  for (int i=1; i<=nb; i++) {
    qn  = he_dev->GetBinContent(i);
    //    ncr = fHistRC[0]->GetBinContent(i);
    err = fH1[13]->GetBinContent(i)/sqrt(qn+1.e-12);
    
    fH1[13]->SetBinError(i,err);
  }

  fH1[11]->SetMarkerStyle(20);
  fH1[11]->SetMarkerSize (1);
  fH1[11]->GetXaxis()->SetRangeUser(300.,699.);
  fH1[11]->SetTitle(Form("%s: <E/crystal> vs radius, per #mubunch, Disk 0",fBgrName[IBgr].Data()));
  fH1[11]->GetXaxis()->SetTitle("Radius, cm");
  fH1[11]->SetStats(0);
  if (EMax > 0)  fH1[11]->SetMaximum(EMax);
  fH1[11]->Draw("pe");
  
  fH1[13]->SetMarkerStyle(24);
  fH1[13]->SetMarkerSize (1);
  fH1[13]->Draw("pe,same");

  TLegend* leg = new TLegend(0.6,0.7,0.9,0.8,"");

  if      (HistSet1 == bgr_tcalm002_301) leg->AddEntry(fH1[11],"old PA, new INA","ep");
  else if (HistSet1 == bgr_tcalm002_nab) leg->AddEntry(fH1[11],"new PA, no  INA","ep");

  leg->AddEntry(fH1[13],"new PA, new INA","ep");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
//-----------------------------------------------------------------------------
// right subpad: second disk
//-----------------------------------------------------------------------------
  Pad->cd(2);

  fH1[10] = gh1(HistSet1[IBgr],"TCalm002",Form("%s_1",HistRE));
  fH1[11] = (TH1F*) fH1[10]->Clone("occupancy");
  fH1[11]->Divide(fH1[10],fHistRC[1],fQExp[IBgr]/qbgr1,1.);
  
  he_301 = gh1(HistSet1[IBgr],"TCalm002",Form("%s_1",HistR));
  
  nb = he_301->GetNbinsX();

  for (int i=1; i<=nb; i++) {
    qn  = he_301->GetBinContent(i);
    //    ncr = fHistRC[1]->GetBinContent(i);
    err = fH1[11]->GetBinContent(i)/sqrt(qn+1.e-12);
    
    fH1[11]->SetBinError(i,err);
  }
					// development: 

  fH1[12] = gh1(HistSet2[IBgr],"TCalm002",Form("%s_1",HistRE));
  fH1[13] = (TH1F*) fH1[12]->Clone("occupancy");
  fH1[13]->Divide(fH1[12],fHistRC[1],fQExp[IBgr]/qbgr2,1.);
  
  he_dev = gh1(HistSet2[IBgr],"TCalm002",Form("%s_1",HistR));
  
  for (int i=1; i<=nb; i++) {
    qn  = he_dev->GetBinContent(i);
    //    ncr = fHistRC[0]->GetBinContent(i);
    err = fH1[13]->GetBinContent(i)/sqrt(qn+1.e-12);
    
    fH1[13]->SetBinError(i,err);
  }

  fH1[11]->SetMarkerStyle(20);
  fH1[11]->SetMarkerSize (1);
  fH1[11]->GetXaxis()->SetRangeUser(300.,699.);
  fH1[11]->SetTitle(Form("%s: <E/crystal> vs radius, per #mubunch, Disk 1",fBgrName[IBgr].Data()));
  fH1[11]->GetXaxis()->SetTitle("Radius, cm");
  fH1[11]->SetStats(0);
  if (EMax > 0)  fH1[11]->SetMaximum(EMax);
  fH1[11]->Draw("pe");
  
  fH1[13]->SetMarkerStyle(24);
  fH1[13]->SetMarkerSize (1);
  fH1[13]->Draw("pe,same");
  
  leg = new TLegend(0.6,0.7,0.9,0.8,"");

  if      (HistSet1 == bgr_tcalm002_301) leg->AddEntry(fH1[11],"old PA, new INA","ep");
  else if (HistSet1 == bgr_tcalm002_nab) leg->AddEntry(fH1[11],"new PA, no  INA","ep");

  leg->AddEntry(fH1[13],"new PA, new INA","ep");

  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

}


//-----------------------------------------------------------------------------
// plot energy per crystal vs R for various background sources
//-----------------------------------------------------------------------------
void TMu2e3280::plot_nhits_vs_r(TPad* Pad, TString HistSet1[], TString HistSet2[], 
				int IBgr, const char* HistRN, const char* HistR, 
				double QMax) {

  TH1F  *he_301, *he_dev;
//-----------------------------------------------------------------------------
// first disk, first - v3_0_1, then - development
//-----------------------------------------------------------------------------
  double qbgr1 = q_background(HistSet1,IBgr);
  double qbgr2 = q_background(HistSet2,IBgr);

  Pad->cd(1);

  fH1[10] = gh1(HistSet1[IBgr],"TCalm002",Form("%s_0",HistRN));
  fH1[11] = (TH1F*) fH1[10]->Clone("occupancy");
  fH1[11]->Divide(fH1[10],fHistRC[0],fQExp[IBgr]/qbgr1,1.);
  
  he_301 = gh1(HistSet1[IBgr],"TCalm002",Form("%s_0",HistR));

  int nb = he_301->GetNbinsX();
  
  double qn, err;

  for (int i=1; i<=nb; i++) {
    qn  = he_301->GetBinContent(i);
    //    ncr = fHistRC[0]->GetBinContent(i);
    err = fH1[11]->GetBinContent(i)/sqrt(qn+1.e-12);
    
    fH1[11]->SetBinError(i,err);
  }
					// development: 

  fH1[12] = gh1(HistSet2[IBgr],"TCalm002",Form("%s_0",HistRN));
  fH1[13] = (TH1F*) fH1[12]->Clone("occupancy");
  fH1[13]->Divide(fH1[12],fHistRC[0],fQExp[IBgr]/qbgr2,1.);

  he_dev = gh1(HistSet2[IBgr],"TCalm002",Form("%s_0",HistR));

  for (int i=1; i<=nb; i++) {
    qn  = he_dev->GetBinContent(i);
    //    ncr = fHistRC[0]->GetBinContent(i);
    err = fH1[13]->GetBinContent(i)/sqrt(qn+1.e-12);
    
    fH1[13]->SetBinError(i,err);
  }

  fH1[11]->SetMarkerStyle(20);
  fH1[11]->SetMarkerSize (1);
  fH1[11]->GetXaxis()->SetRangeUser(300.,699.);
  fH1[11]->SetTitle(Form("%s: <N(hits)/crystal> vs R per #mubunch, Disk 0",fBgrName[IBgr].Data()));
  fH1[11]->GetXaxis()->SetTitle("Radius, cm");
  fH1[11]->SetStats(0);
  if (QMax > 0)  fH1[11]->SetMaximum(QMax);
  fH1[11]->Draw("pe");
  
  fH1[13]->SetMarkerStyle(24);
  fH1[13]->SetMarkerSize (1);
  fH1[13]->Draw("pe,same");

  TLegend* leg = new TLegend(0.6,0.7,0.9,0.8,"");
  if      (HistSet1 == bgr_tcalm002_301) leg->AddEntry(fH1[11],"old PA, new INA","ep");
  else if (HistSet1 == bgr_tcalm002_nab) leg->AddEntry(fH1[11],"new PA, no  INA","ep");

  leg->AddEntry(fH1[13],"new PA, new INA","ep");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
//-----------------------------------------------------------------------------
// right subpad: second disk
//-----------------------------------------------------------------------------
  Pad->cd(2);

  fH1[10] = gh1(HistSet1[IBgr],"TCalm002",Form("%s_1",HistRN));
  fH1[11] = (TH1F*) fH1[10]->Clone("occupancy");
  fH1[11]->Divide(fH1[10],fHistRC[1],fQExp[IBgr]/qbgr1,1.);
  
  he_301 = gh1(HistSet1[IBgr],"TCalm002",Form("%s_1",HistR));
  
  nb = he_301->GetNbinsX();

  for (int i=1; i<=nb; i++) {
    qn  = he_301->GetBinContent(i);
    //    ncr = fHistRC[1]->GetBinContent(i);
    err = fH1[11]->GetBinContent(i)/sqrt(qn+1.e-12);
    
    fH1[11]->SetBinError(i,err);
  }
					// development: 

  fH1[12] = gh1(HistSet2[IBgr],"TCalm002",Form("%s_1",HistRN));
  fH1[13] = (TH1F*) fH1[12]->Clone("occupancy");
  fH1[13]->Divide(fH1[12],fHistRC[1],fQExp[IBgr]/qbgr2,1.);
  
  he_dev = gh1(HistSet2[IBgr],"TCalm002",Form("%s_1",HistR));
  
  for (int i=1; i<=nb; i++) {
    qn  = he_dev->GetBinContent(i);
    //    ncr = fHistRC[0]->GetBinContent(i);
    err = fH1[13]->GetBinContent(i)/sqrt(qn+1.e-12);
    
    fH1[13]->SetBinError(i,err);
  }

  fH1[11]->SetMarkerStyle(20);
  fH1[11]->SetMarkerSize (1);
  fH1[11]->GetXaxis()->SetRangeUser(300.,699.);
  fH1[11]->SetTitle(Form("%s: <N(hits)/crystal> vs R per #mubunch, Disk 1",fBgrName[IBgr].Data()));
  fH1[11]->GetXaxis()->SetTitle("Radius, cm");
  fH1[11]->SetStats(0);
  if (QMax > 0)  fH1[11]->SetMaximum(QMax);
  fH1[11]->Draw("pe");
  
  fH1[13]->SetMarkerStyle(24);
  fH1[13]->SetMarkerSize (1);
  fH1[13]->Draw("pe,same");
  
  leg = new TLegend(0.6,0.7,0.9,0.8,"");

  if      (HistSet1 == bgr_tcalm002_301) leg->AddEntry(fH1[11],"old PA, new INA","ep");
  else if (HistSet1 == bgr_tcalm002_nab) leg->AddEntry(fH1[11],"new PA, no  INA","ep");

  leg->AddEntry(fH1[13],"new PA, new INA","ep");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

}


//-----------------------------------------------------------------------------
void TMu2e3280::plot(Int_t Figure, const char* CanvasName) {
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
					// should be the same for everybody

  fHistRC[0] = gh1(bgr_tcalm002_301[0],"TCalm002", "rc_0");
  fHistRC[1] = gh1(bgr_tcalm002_301[0],"TCalm002", "rc_1");
//-----------------------------------------------------------------------------
// Fig.  1: DIO energy (ibgr = 1)
//-----------------------------------------------------------------------------
  if (Figure == 1) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_energy_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,0,"cal_1/rwe","cal_1/r",30.);
  }
//-----------------------------------------------------------------------------
// Fig.  2: proton energy  (ibgr=1)
//-----------------------------------------------------------------------------
  if (Figure == 2) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_energy_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,1,"cal_1/rwe","cal_1/r",20.);
  }
//-----------------------------------------------------------------------------
// Fig.  3: neutron energy  (ibgr=2)
//-----------------------------------------------------------------------------
  if (Figure == 3) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_energy_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,2,"cal_1/rwe","cal_1/r",2.);
  }
//-----------------------------------------------------------------------------
// Fig.  4: photon energy  (ibgr=3)
//-----------------------------------------------------------------------------
  if (Figure == 4) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_energy_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,3,"cal_1/rwe","cal_1/r",0.8);
  }
//-----------------------------------------------------------------------------
// Fig.  5: DIO number of hits  (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 5) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,0,"cal_1/r","cal_1/r",3.);
  }
//-----------------------------------------------------------------------------
// Fig.  6: Protons number of hits  (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 6) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,1,"cal_1/r","cal_1/r",3.);
  }
//-----------------------------------------------------------------------------
// Fig.  7: Neutrons number of hits  (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 7) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,2,"cal_1/r","cal_1/r",2.5);
  }
//-----------------------------------------------------------------------------
// Fig.  8: Photons number of hits  (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 8) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,3,"cal_1/r","cal_1/r",0.6);
  }
//-----------------------------------------------------------------------------
// Fig.  9: DIO number of hits  T>700 (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 9) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,kDio,"cal_1/r700","cal_1/r700",3);
  }
//-----------------------------------------------------------------------------
// Fig. 10: Protons number of hits  T>700 (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 10) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,kEpr,"cal_1/r700","cal_1/r700",3.);
  }
//-----------------------------------------------------------------------------
// Fig. 11: Neutrons number of hits  T>700 (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 11) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,kEnu,"cal_1/r700","cal_1/r700",2.5);
  }
//-----------------------------------------------------------------------------
// Fig. 12: Photons number of hits  T>700 (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 12) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,kEph,"cal_1/r700","cal_1/r700",0.6);
  }
//-----------------------------------------------------------------------------
// Fig. 13: DIO number of hits E > 0.1 (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 13) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,kDio,"cal_2/r","cal_2/r",3.);
  }
//-----------------------------------------------------------------------------
// Fig.  14: Protons number of hits E > 0.1 (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 14) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,kEpr,"cal_2/r","cal_2/r",3.);
  }
//-----------------------------------------------------------------------------
// Fig.  15: Neutrons number of hits E > 0.1 (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 15) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,kEnu,"cal_2/r","cal_2/r",2.5);
  }
//-----------------------------------------------------------------------------
// Fig.  16: Photons number of hits E > 0.1 (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 16) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,kEph,"cal_2/r","cal_2/r",0.6);
  }
//-----------------------------------------------------------------------------
// Fig. 17: DIO number of hits E > 1 (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 17) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,kDio,"cal_3/r","cal_3/r",3.);
  }
//-----------------------------------------------------------------------------
// Fig.  18: Protons number of hits E > 1 (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 18) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,kEpr,"cal_3/r","cal_3/r",3.0);
  }
//-----------------------------------------------------------------------------
// Fig.  19: Neutrons number of hits E > 1 (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 19) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,kEnu,"cal_3/r","cal_3/r",2.5);
  }
//-----------------------------------------------------------------------------
// Fig.  20: Photons number of hits E > 1 (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 20) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_301,bgr_tcalm002_dev,kEph,"cal_3/r","cal_3/r",0.6);
  }
//-----------------------------------------------------------------------------
// Fig.  21: Hit energy distributions E > 0
//-----------------------------------------------------------------------------
  if (Figure == 21) {
    TPaveLabel*  label;
    TLegend*     leg;

    fCanvas  = new_slide(name,title,2,2,1200,800);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fP1->cd(1);
    gPad->SetLogy(1);

    fH1[0] = gh1( bgr_tcalm002_301[kDio],"TCalm002","cal_1/energy_0");
    fH1[1] = gh1( bgr_tcalm002_301[kDio],"TCalm002","cal_1/energy_1");

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    fH1[0]->GetXaxis()->SetRangeUser(0,59.9);
    fH1[0]->GetXaxis()->SetTitle("Energy,  MeV");
    fH1[0]->SetTitle("");

    fH1[0]->Draw("ep");

    fH1[1]->SetMarkerStyle(24);
    fH1[1]->SetMarkerSize (1);
    fH1[1]->Draw("ep,same");

    DrawPaveLabelNDC(label,"DIO",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
			
    leg = new TLegend(0.4,0.7,0.7,0.8,"");
    leg->AddEntry(fH1[0],"First Disk" ,"ep");
    leg->AddEntry(fH1[1],"Second Disk","ep");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    fP1->cd(2);

    fH1[0] = gh1( bgr_tcalm002_301[kEpr],"TCalm002","cal_1/energy_0");
    fH1[1] = gh1( bgr_tcalm002_301[kEpr],"TCalm002","cal_1/energy_1");

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    fH1[0]->GetXaxis()->SetRangeUser(0,59.9);
    fH1[0]->GetXaxis()->SetTitle("Energy,  MeV");
    fH1[0]->SetTitle("");

    fH1[0]->Draw("ep");

    fH1[1]->SetMarkerStyle(24);
    fH1[1]->SetMarkerSize (1);
    fH1[1]->Draw("ep,same");
			
    DrawPaveLabelNDC(label,"Protons",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
			
    fP1->cd(3);
    gPad->SetLogy(1);

    fH1[0] = gh1( bgr_tcalm002_301[kEnu],"TCalm002","cal_1/energy_0");
    fH1[1] = gh1( bgr_tcalm002_301[kEnu],"TCalm002","cal_1/energy_1");

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    fH1[0]->GetXaxis()->SetRangeUser(0,59.9);
    fH1[0]->GetXaxis()->SetTitle("Energy,  MeV");
    fH1[0]->SetTitle("");

    fH1[0]->Draw("ep");

    fH1[1]->SetMarkerStyle(24);
    fH1[1]->SetMarkerSize (1);
    fH1[1]->Draw("ep,same");
			
    DrawPaveLabelNDC(label,"Neutrons",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
			
    fP1->cd(4);
    gPad->SetLogy(1);

    fH1[0] = gh1( bgr_tcalm002_301[kEph],"TCalm002","cal_1/energy_0");
    fH1[1] = gh1( bgr_tcalm002_301[kEph],"TCalm002","cal_1/energy_1");

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    fH1[0]->GetXaxis()->SetRangeUser(0,59.9);
    fH1[0]->GetXaxis()->SetTitle("Energy,  MeV");
    fH1[0]->SetTitle("");

    fH1[0]->Draw("ep");

    fH1[1]->SetMarkerStyle(24);
    fH1[1]->SetMarkerSize (1);
    fH1[1]->Draw("ep,same");

    DrawPaveLabelNDC(label,"Photons",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
  }
//-----------------------------------------------------------------------------
// Fig.  22: Hit energy distributions E > 0.1 MeV 
//-----------------------------------------------------------------------------
  if (Figure == 22) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,2,2,1200,800);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fP1->cd(1);
    gPad->SetLogy(1);

    fH1[0] = gh1( bgr_tcalm002_301[kDio],"TCalm002","cal_2/energy_0");
    fH1[1] = gh1( bgr_tcalm002_301[kDio],"TCalm002","cal_2/energy_1");

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    fH1[0]->GetXaxis()->SetRangeUser(0,59.9);
    fH1[0]->GetXaxis()->SetTitle("Energy,  MeV");
    fH1[0]->SetTitle("");

    fH1[0]->Draw("ep");

    fH1[1]->SetMarkerStyle(24);
    fH1[1]->SetMarkerSize (1);
    fH1[1]->Draw("ep,same");

    DrawPaveLabelNDC(label,"DIO",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
			
    leg = new TLegend(0.4,0.7,0.7,0.8,"");
    leg->AddEntry(fH1[0],"First Disk" ,"ep");
    leg->AddEntry(fH1[1],"Second Disk","ep");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    fP1->cd(2);
    fH1[0] = gh1( bgr_tcalm002_301[kEpr],"TCalm002","cal_2/energy_0");
    fH1[1] = gh1( bgr_tcalm002_301[kEpr],"TCalm002","cal_2/energy_1");

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    fH1[0]->GetXaxis()->SetRangeUser(0,59.9);
    fH1[0]->GetXaxis()->SetTitle("Energy,  MeV");
    fH1[0]->SetTitle("");

    fH1[0]->Draw("ep");

    fH1[1]->SetMarkerStyle(24);
    fH1[1]->SetMarkerSize (1);
    fH1[1]->Draw("ep,same");

    DrawPaveLabelNDC(label,"Protons",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
			
    fP1->cd(3);
    gPad->SetLogy(1);

    fH1[0] = gh1( bgr_tcalm002_301[kEnu],"TCalm002","cal_2/energy_0");
    fH1[1] = gh1( bgr_tcalm002_301[kEnu],"TCalm002","cal_2/energy_1");

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    fH1[0]->GetXaxis()->SetRangeUser(0,59.9);
    fH1[0]->GetXaxis()->SetTitle("Energy,  MeV");
    fH1[0]->SetTitle("");

    fH1[0]->Draw("ep");

    fH1[1]->SetMarkerStyle(24);
    fH1[1]->SetMarkerSize (1);
    fH1[1]->Draw("ep,same");

    DrawPaveLabelNDC(label,"Neutrons",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
			
    fP1->cd(4);
    gPad->SetLogy(1);

    fH1[0] = gh1( bgr_tcalm002_301[kEph],"TCalm002","cal_2/energy_0");
    fH1[1] = gh1( bgr_tcalm002_301[kEph],"TCalm002","cal_2/energy_1");

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    fH1[0]->GetXaxis()->SetRangeUser(0,59.9);
    fH1[0]->GetXaxis()->SetTitle("Energy,  MeV");
    fH1[0]->SetTitle("");

    fH1[0]->Draw("ep");

    fH1[1]->SetMarkerStyle(24);
    fH1[1]->SetMarkerSize (1);
    fH1[1]->Draw("ep,same");
			
    DrawPaveLabelNDC(label,"Photons",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);			
  }
//-----------------------------------------------------------------------------
// Fig.  23: total occupancies - E>0 vs E>1 MeV
//-----------------------------------------------------------------------------
  if (Figure == 23) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new TCanvas(name,title,1100,700);
    fP1      = (TPad*) fCanvas->GetPad(0);
    
    fP1->cd(1);

    fH1[0] = occupancy(bgr_tcalm002_dev,"h_dio_23",kDio, fQGenDev, "cal_1/r_0", "cal_1/r_0");
    fH1[1] = occupancy(bgr_tcalm002_dev,"h_epr_23",kEpr, fQGenDev, "cal_1/r_0", "cal_1/r_0");
    fH1[2] = occupancy(bgr_tcalm002_dev,"h_enu_23",kEnu, fQGenDev, "cal_1/r_0", "cal_1/r_0");
    fH1[3] = occupancy(bgr_tcalm002_dev,"h_eph_23",kEph, fQGenDev, "cal_1/r_0", "cal_1/r_0");


    fH1[0]->Add(fH1[1],1.);
    fH1[0]->Add(fH1[2],1.);
    fH1[0]->Add(fH1[3],1.);

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    fH1[0]->GetXaxis()->SetRangeUser(300.,699);
    fH1[0]->GetXaxis()->SetTitle("Radius,  cm");
    fH1[0]->GetYaxis()->SetTitle("N(hits)/crystal/#mubunch");
    fH1[0]->SetTitle("");

    fH1[0]->SetStats(0);

    fH1[0]->Draw();


    fH1[10] = occupancy(bgr_tcalm002_dev,"h_dio_23",kDio, fQGenDev, "cal_3/r_0", "cal_3/r_0");
    fH1[11] = occupancy(bgr_tcalm002_dev,"h_epr_23",kEpr, fQGenDev, "cal_3/r_0", "cal_3/r_0");
    fH1[12] = occupancy(bgr_tcalm002_dev,"h_enu_23",kEnu, fQGenDev, "cal_3/r_0", "cal_3/r_0");
    fH1[13] = occupancy(bgr_tcalm002_dev,"h_eph_23",kEph, fQGenDev, "cal_3/r_0", "cal_3/r_0");

    fH1[10]->Add(fH1[11],1.);
    fH1[10]->Add(fH1[12],1.);
    fH1[10]->Add(fH1[13],1.);

    fH1[10]->SetMarkerStyle(24);
    fH1[10]->SetMarkerSize (1);
    fH1[10]->Draw("pe,same");

    leg = new TLegend(0.4,0.7,0.7,0.8,"");
    leg->AddEntry(fH1[0],"all hits" ,"ep");
    leg->AddEntry(fH1[10],"E > 1 MeV","ep");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    if (fPlotMode == kTalkMode) {
      DrawPaveLabelNDC(label,"expected occupancy per #mubunch",0.35,0.85,0.8,0.95,52);
    }
    else {
      DrawPaveLabelNDC(label,"total expected occupancy per #mubunch, new PA, new INA",0.35,0.85,0.8,0.95,52);
    }
    label->SetTextSize(0.4);
  }			
//-----------------------------------------------------------------------------
// Fig. 24: total occupancies: new PA vs old PA
//-----------------------------------------------------------------------------
  if (Figure == 24) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,1,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fP1->cd(1);

    fH1[0] = occupancy(bgr_tcalm002_301,"h_dio_24",kDio, fQGen301, "cal_1/r_0", "cal_1/r_0");
    fH1[1] = occupancy(bgr_tcalm002_301,"h_epr_24",kEpr, fQGen301, "cal_1/r_0", "cal_1/r_0");
    fH1[2] = occupancy(bgr_tcalm002_301,"h_enu_24",kEnu, fQGen301, "cal_1/r_0", "cal_1/r_0");
    fH1[3] = occupancy(bgr_tcalm002_301,"h_eph_24",kEph, fQGen301, "cal_1/r_0", "cal_1/r_0");


    fH1[0]->Add(fH1[1],1.);
    fH1[0]->Add(fH1[2],1.);
    fH1[0]->Add(fH1[3],1.);

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    fH1[0]->GetXaxis()->SetRangeUser(300.,699);
    fH1[0]->GetXaxis()->SetTitle("Radius,  cm");
    fH1[0]->SetTitle("");
    fH1[0]->GetYaxis()->SetRangeUser(0.,6);

    fH1[0]->SetStats(0);

    fH1[0]->Draw();


    fH1[10] = occupancy(bgr_tcalm002_dev,"h_dio_24",kDio, fQGenDev, "cal_1/r_0", "cal_1/r_0");
    fH1[11] = occupancy(bgr_tcalm002_dev,"h_epr_24",kEpr, fQGenDev, "cal_1/r_0", "cal_1/r_0");
    fH1[12] = occupancy(bgr_tcalm002_dev,"h_enu_24",kEnu, fQGenDev, "cal_1/r_0", "cal_1/r_0");
    fH1[13] = occupancy(bgr_tcalm002_dev,"h_eph_24",kEph, fQGenDev, "cal_1/r_0", "cal_1/r_0");

    fH1[10]->Add(fH1[11],1.);
    fH1[10]->Add(fH1[12],1.);
    fH1[10]->Add(fH1[13],1.);

    fH1[10]->SetMarkerStyle(24);
    fH1[10]->SetMarkerSize (1);
    fH1[10]->Draw("pe,same");

    leg = new TLegend(0.4,0.7,0.7,0.8,"");
    leg->AddEntry(fH1[ 0],"old PA" ,"ep");
    leg->AddEntry(fH1[10],"new PA","ep");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"total expected occupancy per microbunch, E_{HIT} > 0",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
  }			
//-----------------------------------------------------------------------------
// Fig. 25: total occupancies, E > 0, new PA, NO_NABS vs New_Abs
//-----------------------------------------------------------------------------
  if (Figure == 25) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,1,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fP1->cd(1);

    fH1[0] = occupancy(bgr_tcalm002_nab,"h0_dio_25",kDio, fQGenNab, "cal_1/r_0", "cal_1/r_0");
    fH1[1] = occupancy(bgr_tcalm002_nab,"h1_epr_25",kEpr, fQGenNab, "cal_1/r_0", "cal_1/r_0");
    fH1[2] = occupancy(bgr_tcalm002_nab,"h2_enu_25",kEnu, fQGenNab, "cal_1/r_0", "cal_1/r_0");
    fH1[3] = occupancy(bgr_tcalm002_nab,"h3_eph_25",kEph, fQGenNab, "cal_1/r_0", "cal_1/r_0");


    fH1[0]->Add(fH1[1],1.);
    fH1[0]->Add(fH1[2],1.);
    fH1[0]->Add(fH1[3],1.);

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    fH1[0]->GetXaxis()->SetRangeUser(300.,699);
    fH1[0]->GetXaxis()->SetTitle("Radius,  cm");
    fH1[0]->GetYaxis()->SetRangeUser(0.,6);
    fH1[0]->SetTitle("");

    fH1[0]->SetStats(0);

    fH1[0]->Draw();


    fH1[10] = occupancy(bgr_tcalm002_dev,"h10_dio_25",kDio, fQGenDev, "cal_1/r_0", "cal_1/r_0");
    fH1[11] = occupancy(bgr_tcalm002_dev,"h11_epr_25",kEpr, fQGenDev, "cal_1/r_0", "cal_1/r_0");
    fH1[12] = occupancy(bgr_tcalm002_dev,"h12_enu_25",kEnu, fQGenDev, "cal_1/r_0", "cal_1/r_0");
    fH1[13] = occupancy(bgr_tcalm002_dev,"h13_eph_25",kEph, fQGenDev, "cal_1/r_0", "cal_1/r_0");

    fH1[10]->Add(fH1[11],1.);
    fH1[10]->Add(fH1[12],1.);
    fH1[10]->Add(fH1[13],1.);

    fH1[10]->SetMarkerStyle(24);
    fH1[10]->SetMarkerSize (1);
    fH1[10]->Draw("pe,same");

    leg = new TLegend(0.4,0.7,0.7,0.8,"");
    leg->AddEntry(fH1[ 0],"no neutron absorber" ,"ep");
    leg->AddEntry(fH1[10],"new neutron absorber","ep");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"total expected occupancy per microbunch, E_{HIT} > 0",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
  }			
//-----------------------------------------------------------------------------
// Fig. 26: total occupancies: old PA vs new PA for hit E > 1 MeV : 
//-----------------------------------------------------------------------------
  if (Figure == 26) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,1,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fP1->cd(1);

    fH1[0] = occupancy(bgr_tcalm002_301,"h_dio_26",kDio, fQGen301, "cal_3/r_0", "cal_3/r_0");
    fH1[1] = occupancy(bgr_tcalm002_301,"h_epr_26",kEpr, fQGen301, "cal_3/r_0", "cal_3/r_0");
    fH1[2] = occupancy(bgr_tcalm002_301,"h_enu_26",kEnu, fQGen301, "cal_3/r_0", "cal_3/r_0");
    fH1[3] = occupancy(bgr_tcalm002_301,"h_eph_26",kEph, fQGen301, "cal_3/r_0", "cal_3/r_0");


    fH1[0]->Add(fH1[1],1.);
    fH1[0]->Add(fH1[2],1.);
    fH1[0]->Add(fH1[3],1.);

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    fH1[0]->GetXaxis()->SetRangeUser(300.,699);
    fH1[0]->GetXaxis()->SetTitle("Radius,  cm");
    fH1[0]->SetTitle("");
    fH1[0]->GetYaxis()->SetRangeUser(0.,6);

    fH1[0]->SetStats(0);

    fH1[0]->Draw();


    fH1[10] = occupancy(bgr_tcalm002_dev,"h_dio_26",kDio, fQGenDev, "cal_3/r_0", "cal_3/r_0");
    fH1[11] = occupancy(bgr_tcalm002_dev,"h_epr_26",kEpr, fQGenDev, "cal_3/r_0", "cal_3/r_0");
    fH1[12] = occupancy(bgr_tcalm002_dev,"h_enu_26",kEnu, fQGenDev, "cal_3/r_0", "cal_3/r_0");
    fH1[13] = occupancy(bgr_tcalm002_dev,"h_eph_26",kEph, fQGenDev, "cal_3/r_0", "cal_3/r_0");

    fH1[10]->Add(fH1[11],1.);
    fH1[10]->Add(fH1[12],1.);
    fH1[10]->Add(fH1[13],1.);

    fH1[10]->SetMarkerStyle(24);
    fH1[10]->SetMarkerSize (1);
    fH1[10]->Draw("pe,same");

    leg = new TLegend(0.4,0.7,0.7,0.8,"");
    leg->AddEntry(fH1[ 0],"old PA" ,"ep");
    leg->AddEntry(fH1[10],"new PA","ep");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"total expected occupancy per microbunch, E_{HIT} > 1 MeV",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
  }			
//-----------------------------------------------------------------------------
// Fig. 27: total occupancies, E > 1 MeV, new PA, NO_NABS vs NEW_NABS
//-----------------------------------------------------------------------------
  if (Figure == 27) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,1,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fP1->cd(1);

    fH1[0] = occupancy(bgr_tcalm002_nab,"h_dio_26",kDio, fQGenNab, "cal_3/r_0", "cal_3/r_0");
    fH1[1] = occupancy(bgr_tcalm002_nab,"h_epr_26",kEpr, fQGenNab, "cal_3/r_0", "cal_3/r_0");
    fH1[2] = occupancy(bgr_tcalm002_nab,"h_enu_26",kEnu, fQGenNab, "cal_3/r_0", "cal_3/r_0");
    fH1[3] = occupancy(bgr_tcalm002_nab,"h_eph_26",kEph, fQGenNab, "cal_3/r_0", "cal_3/r_0");


    fH1[0]->Add(fH1[1],1.);
    fH1[0]->Add(fH1[2],1.);
    fH1[0]->Add(fH1[3],1.);

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    fH1[0]->GetXaxis()->SetRangeUser(300.,699);
    fH1[0]->GetXaxis()->SetTitle("Radius,  cm");
    fH1[0]->SetTitle("");
    fH1[0]->GetYaxis()->SetRangeUser(0.,6);

    fH1[0]->SetStats(0);

    fH1[0]->Draw();


    fH1[10] = occupancy(bgr_tcalm002_dev,"h_dio_26",kDio, fQGenDev, "cal_3/r_0", "cal_3/r_0");
    fH1[11] = occupancy(bgr_tcalm002_dev,"h_epr_26",kEpr, fQGenDev, "cal_3/r_0", "cal_3/r_0");
    fH1[12] = occupancy(bgr_tcalm002_dev,"h_enu_26",kEnu, fQGenDev, "cal_3/r_0", "cal_3/r_0");
    fH1[13] = occupancy(bgr_tcalm002_dev,"h_eph_26",kEph, fQGenDev, "cal_3/r_0", "cal_3/r_0");

    fH1[10]->Add(fH1[11],1.);
    fH1[10]->Add(fH1[12],1.);
    fH1[10]->Add(fH1[13],1.);

    fH1[10]->SetMarkerStyle(24);
    fH1[10]->SetMarkerSize (1);
    fH1[10]->Draw("pe,same");

    leg = new TLegend(0.4,0.7,0.7,0.8,"");
    leg->AddEntry(fH1[ 0],"no neutron absorber ","ep");
    leg->AddEntry(fH1[10],"new neutron absorber","ep");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"total expected occupancy per microbunch, E_{HIT} > 1 MeV",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
  }			
//-----------------------------------------------------------------------------
// Fig. 28: Hit timing distributions E > 0
//-----------------------------------------------------------------------------
  if (Figure == 28) {
    TPaveLabel*  label;
    TLegend*     leg;

    fCanvas  = new_slide(name,title,2,2,1200,800);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fP1->cd(1);

    fH1[0] = (TH1F*) gh1( bgr_tcalm002_301[kDio],"TCalm002","cal_1/time_0")->Clone("h1_0_fig_28_pad_1");
    fH1[1] = (TH1F*) gh1( bgr_tcalm002_301[kDio],"TCalm002","cal_1/time_1")->Clone("h1_1_fig_28_pad_1");

    fH1[0]->Rebin(5);
    fH1[1]->Rebin(5);

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    //    fH1[0]->GetXaxis()->SetRangeUser(0,59.9);
    fH1[0]->GetXaxis()->SetTitle("Time, ns");
    fH1[0]->SetTitle("");

    fH1[0]->SetMaximum(250);
    fH1[0]->Draw("ep");

    fH1[1]->SetMarkerStyle(24);
    fH1[1]->SetMarkerSize (1);
    fH1[1]->Draw("ep,same");

    DrawPaveLabelNDC(label,"DIO",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
			
    leg = new TLegend(0.4,0.7,0.7,0.8,"");
    leg->AddEntry(fH1[0],"First Disk" ,"ep");
    leg->AddEntry(fH1[1],"Second Disk","ep");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    fP1->cd(2);

    fH1[0] = (TH1F*) gh1( bgr_tcalm002_301[kEpr],"TCalm002","cal_1/time_0")->Clone("h1_0_fig_28_pad_2");
    fH1[1] = (TH1F*) gh1( bgr_tcalm002_301[kEpr],"TCalm002","cal_1/time_1")->Clone("h1_1_fig_28_pad_2");

    fH1[0]->Rebin(5);
    fH1[1]->Rebin(5);

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    //    fH1[0]->GetXaxis()->SetRangeUser(0,59.9);
    fH1[0]->GetXaxis()->SetTitle("Time, ns");
    fH1[0]->SetTitle("");

    fH1[0]->SetMaximum(100.);
    fH1[0]->Draw("ep");

    fH1[1]->SetMarkerStyle(24);
    fH1[1]->SetMarkerSize (1);
    fH1[1]->Draw("ep,same");
			
    DrawPaveLabelNDC(label,"Protons",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
			
    fP1->cd(3);

    fH1[0] = (TH1F*) gh1( bgr_tcalm002_301[kEnu],"TCalm002","cal_1/time_0")->Clone("h1_0_fig_28_pad_3");
    fH1[1] = (TH1F*) gh1( bgr_tcalm002_301[kEnu],"TCalm002","cal_1/time_1")->Clone("h1_1_fig_28_pad_3");

    fH1[0]->Rebin(5);
    fH1[1]->Rebin(5);

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    //    fH1[0]->GetXaxis()->SetRangeUser(0,59.9);
    fH1[0]->GetXaxis()->SetTitle("Time, ns");
    fH1[0]->SetTitle("");

    fH1[0]->Draw("ep");

    fH1[1]->SetMarkerStyle(24);
    fH1[1]->SetMarkerSize (1);
    fH1[1]->Draw("ep,same");
			
    DrawPaveLabelNDC(label,"Neutrons",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
			
    fP1->cd(4);

    fH1[0] = (TH1F*) gh1( bgr_tcalm002_301[kEph],"TCalm002","cal_1/time_0")->Clone("h1_0_fig_28_pad_4");
    fH1[1] = (TH1F*) gh1( bgr_tcalm002_301[kEph],"TCalm002","cal_1/time_1")->Clone("h1_1_fig_28_pad_4");

    fH1[0]->Rebin(5);
    fH1[1]->Rebin(5);

    fH1[0]->SetMarkerStyle(20);
    fH1[0]->SetMarkerSize (1);
    //    fH1[0]->GetXaxis()->SetRangeUser(0,59.9);
    fH1[0]->GetXaxis()->SetTitle("time, ns");
    fH1[0]->SetTitle("");

    fH1[0]->Draw("ep");

    fH1[1]->SetMarkerStyle(24);
    fH1[1]->SetMarkerSize (1);
    fH1[1]->Draw("ep,same");

    DrawPaveLabelNDC(label,"Photons",0.35,0.85,0.8,0.95,52);
    label->SetTextSize(0.4);
  }
//-----------------------------------------------------------------------------
// Fig. 31: DIO energy (ibgr = 0), NEW_PA, NO_NABS vs NEW_NABS
//-----------------------------------------------------------------------------
  if (Figure == 31) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_energy_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,0,"cal_1/rwe","cal_1/r",30.);
  }
//-----------------------------------------------------------------------------
// Fig. 32: proton energy  (ibgr=1), NEW_PA, NO_NABS vs NEW_NABS
//-----------------------------------------------------------------------------
  if (Figure == 32) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_energy_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,1,"cal_1/rwe","cal_1/r",20.);
  }
//-----------------------------------------------------------------------------
// Fig. 33: neutron energy  (ibgr=2), NEW_PA, NO_NABS vs NEW_NABS
//-----------------------------------------------------------------------------
  if (Figure == 33) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_energy_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,2,"cal_1/rwe","cal_1/r",2.);
  }
//-----------------------------------------------------------------------------
// Fig. 34: photon energy  (ibgr=3), NEW_PA, NO_NABS vs NEW_NABS
//-----------------------------------------------------------------------------
  if (Figure == 34) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_energy_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,3,"cal_1/rwe","cal_1/r",0.8);
  }
//-----------------------------------------------------------------------------
// Fig. 35: DIO number of hits  (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 35) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,0,"cal_1/r","cal_1/r",3.);
  }
//-----------------------------------------------------------------------------
// Fig. 36: Protons number of hits  (ibgr=1)
//-----------------------------------------------------------------------------
  if (Figure == 36) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,1,"cal_1/r","cal_1/r",3.);
  }
//-----------------------------------------------------------------------------
// Fig. 37: Neutrons number of hits  (ibgr=2)
//-----------------------------------------------------------------------------
  if (Figure == 37) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,2,"cal_1/r","cal_1/r",2.5);
  }
//-----------------------------------------------------------------------------
// Fig. 38: Photons number of hits  (ibgr=3)
//-----------------------------------------------------------------------------
  if (Figure == 38) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,3,"cal_1/r","cal_1/r",0.6);
  }
//-----------------------------------------------------------------------------
// Fig. 39: DIO number of hits  T>700 (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 39) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,kDio,"cal_1/r700","cal_1/r700",3);
  }
//-----------------------------------------------------------------------------
// Fig. 40: Protons number of hits  T>700 (ibgr=1)
//-----------------------------------------------------------------------------
  if (Figure == 40) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,kEpr,"cal_1/r700","cal_1/r700",3.);
  }
//-----------------------------------------------------------------------------
// Fig. 41: Neutrons number of hits  T>700 (ibgr=2)
//-----------------------------------------------------------------------------
  if (Figure == 41) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,kEnu,"cal_1/r700","cal_1/r700",2.5);
  }
//-----------------------------------------------------------------------------
// Fig. 42: Photons number of hits  T>700 (ibgr=3)
//-----------------------------------------------------------------------------
  if (Figure == 42) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,kEph,"cal_1/r700","cal_1/r700",0.6);
  }
//-----------------------------------------------------------------------------
// Fig. 43: DIO number of hits E > 0.1 MeV (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 43) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,kDio,"cal_2/r","cal_2/r",3.);
  }
//-----------------------------------------------------------------------------
// Fig. 44: Protons number of hits E > 0.1 MeV (ibgr=1)
//-----------------------------------------------------------------------------
  if (Figure == 44) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,kEpr,"cal_2/r","cal_2/r",3.);
  }
//-----------------------------------------------------------------------------
// Fig.  45: Neutrons number of hits E > 0.1 MeV (ibgr=2)
//-----------------------------------------------------------------------------
  if (Figure == 45) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,kEnu,"cal_2/r","cal_2/r",2.5);
  }
//-----------------------------------------------------------------------------
// Fig.  46: Photons number of hits E > 0.1 MeV (ibgr=3)
//-----------------------------------------------------------------------------
  if (Figure == 46) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,kEph,"cal_2/r","cal_2/r",0.6);
  }
//-----------------------------------------------------------------------------
// Fig. 47: DIO number of hits E > 1 MeV (ibgr=0)
//-----------------------------------------------------------------------------
  if (Figure == 47) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,kDio,"cal_3/r","cal_3/r",3.);
  }
//-----------------------------------------------------------------------------
// Fig.  48: Protons number of hits E > 1 MeV (ibgr=1)
//-----------------------------------------------------------------------------
  if (Figure == 48) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,kEpr,"cal_3/r","cal_3/r",3.0);
  }
//-----------------------------------------------------------------------------
// Fig.  49: Neutrons number of hits E > 1 MeV (ibgr=2)
//-----------------------------------------------------------------------------
  if (Figure == 49) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,kEnu,"cal_3/r","cal_3/r",2.5);
  }
//-----------------------------------------------------------------------------
// Fig.  50: Photons number of hits E > 1 MeV (ibgr=3)
//-----------------------------------------------------------------------------
  if (Figure == 50) {
    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
   
    plot_nhits_vs_r(fP1,bgr_tcalm002_nab,bgr_tcalm002_dev,kEph,"cal_3/r","cal_3/r",0.6);
  }
//-----------------------------------------------------------------------------
// Fig.102: emitted protons 301 : ZR dist of the origins of particles hitting the calorimeter
//-----------------------------------------------------------------------------
  if (Figure == 102) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    fH2[0] = gh2(bgr_tcalm005_301[kEpr],"TCalm005","stp_11/r0z0");
    fH2[0]->SetMarkerColor(2);
    fH2[1] = gh2(bgr_tcalm005_301[kEpr],"TCalm005","stp_12/r0z0");
    fH2[1]->SetMarkerColor(4);

    fH2[0]->GetXaxis()->SetRangeUser(2000.,14999.);
    fH2[0]->GetYaxis()->SetRangeUser(0.   , 1499.);

    fH2[0]->GetXaxis()->SetTitle("Z,  mm");
    fH2[0]->GetYaxis()->SetTitle("R,  mm");
    fH2[0]->SetTitle("");

    fH2[0]->SetStats(0);

    fH2[0]->Draw();

    fH2[1]->Draw("same");

    leg = new TLegend(0.6,0.82,0.9,0.89,"");
    leg->AddEntry(fH2[0],"red : e+/e-/photons","");
    leg->AddEntry(fH2[1],"blue: neutrons     ","");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"Protons: sources of particles hitting the calorimeter, old PABS",0.10,0.85,0.55,0.95,52);
    label->SetTextSize(0.3);
  }			
//-----------------------------------------------------------------------------
// Fig.103: emitted neutrons 301 : ZR dist of the origins of particles hitting the calorimeter
//-----------------------------------------------------------------------------
  if (Figure == 103) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    fH2[0] = gh2(bgr_tcalm005_301[kEnu],"TCalm005","stp_11/r0z0");
    fH2[0]->SetMarkerColor(2);
    fH2[1] = gh2(bgr_tcalm005_301[kEnu],"TCalm005","stp_12/r0z0");
    fH2[1]->SetMarkerColor(4);

    fH2[0]->GetXaxis()->SetRangeUser(2000.,14999.);
    fH2[0]->GetYaxis()->SetRangeUser(0.   , 1499.);

    fH2[0]->GetXaxis()->SetTitle("Z,  mm");
    fH2[0]->GetYaxis()->SetTitle("R,  mm");
    fH2[0]->SetTitle("");

    fH2[0]->SetStats(0);

    fH2[0]->Draw();

    fH2[1]->Draw("same");

    leg = new TLegend(0.6,0.82,0.9,0.89,"");
    leg->AddEntry(fH2[0],"red : e+/e-/photons","");
    leg->AddEntry(fH2[1],"blue: neutrons     ","");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"Ejected neutrons: sources of particles hitting the calorimeter",0.10,0.85,0.55,0.95,52);
    label->SetTextSize(0.3);
  }			
//-----------------------------------------------------------------------------
// Fig.112: emitted protons DEV: ZR dist of the origins of particles hitting the calorimeter
//-----------------------------------------------------------------------------
  if (Figure == 112) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    fH2[0] = gh2(bgr_tcalm005_dev[kEpr],"TCalm005","stp_11/r0z0");
    fH2[0]->SetMarkerColor(2);
    fH2[1] = gh2(bgr_tcalm005_dev[kEpr],"TCalm005","stp_12/r0z0");
    fH2[1]->SetMarkerColor(4);

    fH2[0]->GetXaxis()->SetRangeUser(2000.,14999.);
    fH2[0]->GetYaxis()->SetRangeUser(0.   , 1499.);

    fH2[0]->GetXaxis()->SetTitle("Z,  mm");
    fH2[0]->GetYaxis()->SetTitle("R,  mm");
    fH2[0]->SetTitle("");

    fH2[0]->SetStats(0);

    fH2[0]->Draw();

    fH2[1]->Draw("same");

    leg = new TLegend(0.6,0.82,0.9,0.89,"");
    leg->AddEntry(fH2[0],"red : e+/e-/photons","");
    leg->AddEntry(fH2[1],"blue: neutrons     ","");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"Protons: sources of particles hitting the calorimeter, new PABS",0.10,0.85,0.55,0.95,52);
    label->SetTextSize(0.3);
  }			
//-----------------------------------------------------------------------------
// Fig.113: emitted neutrons DEV: ZR dist of the origins of particles hitting the calorimeter
//-----------------------------------------------------------------------------
  if (Figure == 113) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    fH2[0] = gh2(bgr_tcalm005_dev[kEnu],"TCalm005","stp_11/r0z0");
    fH2[0]->SetMarkerColor(2);
    fH2[1] = gh2(bgr_tcalm005_dev[kEnu],"TCalm005","stp_12/r0z0");
    fH2[1]->SetMarkerColor(4);

    fH2[0]->GetXaxis()->SetRangeUser(2000.,14999.);
    fH2[0]->GetYaxis()->SetRangeUser(0.   , 1499.);

    fH2[0]->GetXaxis()->SetTitle("Z,  mm");
    fH2[0]->GetYaxis()->SetTitle("R,  mm");
    fH2[0]->SetTitle("");

    fH2[0]->SetStats(0);

    fH2[0]->Draw();

    fH2[1]->Draw("same");

    leg = new TLegend(0.6,0.83,0.9,0.89,"");
    leg->AddEntry(fH2[0],"red : e+/e-/photons","");
    leg->AddEntry(fH2[1],"blue: neutrons     ","");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"Ejected neutrons: sources of particles hitting the calorimeter",0.10,0.85,0.55,0.95,52);
    label->SetTextSize(0.3);
  }			
//-----------------------------------------------------------------------------
// Fig.123: emitted neutrons DEV/NO_NABS: ZR dist of the origins of particles hitting the calorimeter
//-----------------------------------------------------------------------------
  if (Figure == 123) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    fH2[0] = gh2(bgr_tcalm005_nab[kEnu],"TCalm005","stp_11/r0z0");
    fH2[0]->SetMarkerColor(2);
    fH2[1] = gh2(bgr_tcalm005_nab[kEnu],"TCalm005","stp_12/r0z0");
    fH2[1]->SetMarkerColor(4);

    fH2[0]->GetXaxis()->SetRangeUser(2000.,14999.);
    fH2[0]->GetYaxis()->SetRangeUser(0.   , 1499.);

    fH2[0]->GetXaxis()->SetTitle("Z,  mm");
    fH2[0]->GetYaxis()->SetTitle("R,  mm");
    fH2[0]->SetTitle("");

    fH2[0]->SetStats(0);

    fH2[0]->Draw();

    fH2[1]->Draw("same");

    leg = new TLegend(0.6,0.83,0.9,0.89,"");
    leg->AddEntry(fH2[0],"red : e+/e-/photons","");
    leg->AddEntry(fH2[1],"blue: neutrons     ","");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    DrawPaveLabelNDC(label,"Ejected neutrons: sources of particles hitting the calorimeter",0.10,0.85,0.55,0.95,52);
    label->SetTextSize(0.3);
  }			
//-----------------------------------------------------------------------------
// Fig.201: protons - ratio of the number of hits - NEW_PABS/OLD_PABS for disk#1 vs radius
//-----------------------------------------------------------------------------
  if (Figure == 201) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fH1[0] = (TH1F*) gh1(bgr_tcalm002_301[kEpr],"TCalm002","cal_1/r_0")->Clone();
    fH1[1] = (TH1F*) gh2(bgr_tcalm002_dev[kEpr],"TCalm002","cal_1/r_0")->Clone();

    fH1[2] = (TH1F*) fH1[0]->Clone("h1_2_fig_201");

    fH1[0]->Sumw2();
    fH1[1]->Sumw2();
					// Pad#1: the hit energy distributions
    fP1->cd(1);

    fH1[1]->SetTitle("");
    fH1[1]->GetXaxis()->SetRangeUser(300.,699.);
    fH1[1]->GetXaxis()->SetTitle("Hit Crystal Radius, mm");
    fH1[1]->Draw("h");

    fH1[0]->SetFillColor(kBlue-7);
    fH1[0]->SetFillStyle(3003);

    fH1[0]->Draw("h,same");	

    leg = new TLegend(0.6,0.83,0.9,0.89,"");
    leg->AddEntry(fH1[0],"old PABS","f");
    leg->AddEntry(fH1[1],"new PABS","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
				// the histogram ratio 
    fP1->cd(2);

    fH1[2]->Divide(fH1[1],fH1[0],1,1);

    fH1[2]->SetMarkerStyle(20);
    fH1[2]->SetMarkerSize(1);
    fH1[2]->GetXaxis()->SetRangeUser(300.,699.);
    fH1[2]->GetXaxis()->SetTitle("Hit Crystal Radius, mm");
    fH1[2]->Draw("e,p");

    DrawPaveLabelNDC(label,"Protons: N(NEW_PABS)/N(OLD_PABS) vs R",0.10,0.85,0.55,0.95,52);
    label->SetTextSize(0.3);
  }			
//-----------------------------------------------------------------------------
// Fig.202: protons - ratio of the number of hits - NEW_PABS/OLD_PABS for disk#1 vs energy
//-----------------------------------------------------------------------------
  if (Figure == 202) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fH1[0] = (TH1F*) gh1(bgr_tcalm002_301[kEpr],"TCalm002","cal_1/energy_0")->Clone();
    fH1[1] = (TH1F*) gh2(bgr_tcalm002_dev[kEpr],"TCalm002","cal_1/energy_0")->Clone();

    fH1[0]->Sumw2();
    fH1[1]->Sumw2();
					// Pad#1: the hit energy distributions
    fP1->cd(1);

    fH1[1]->SetTitle("");
    fH1[1]->GetXaxis()->SetRangeUser(0.,49.9);
    fH1[1]->GetXaxis()->SetTitle("Hit Energy, MeV");
    fH1[1]->Draw("h");

    fH1[0]->SetFillColor(kBlue-7);
    fH1[0]->SetFillStyle(3003);

    fH1[0]->Draw("h,same");

    leg = new TLegend(0.6,0.83,0.9,0.89,"");
    leg->AddEntry(fH1[0],"old PABS","f");
    leg->AddEntry(fH1[1],"new PABS","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
					// the histogram ratio 
    fP1->cd(2);

    fH1[2] = (TH1F*) fH1[0]->Clone("h1_2_fig_201");

    fH1[2]->Divide(fH1[1],fH1[0],1,1);

    fH1[2]->SetMarkerStyle(20);
    fH1[2]->SetMarkerSize(1);
    fH1[2]->GetXaxis()->SetRangeUser(0,49.9);
    fH1[2]->GetXaxis()->SetTitle("Hit Energy, MeV");
    fH1[2]->Draw("e,p");

    DrawPaveLabelNDC(label,"Protons: N(NEW_PABS)/N(OLD_PABS) vs hit energy",0.10,0.85,0.55,0.95,52);
    label->SetTextSize(0.3);
  }			
//-----------------------------------------------------------------------------
// Fig.203: protons - PDG code for particles hitting the calorimeter
//-----------------------------------------------------------------------------
  if (Figure == 203) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,2,1,1250,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fH1[0] = (TH1F*) gh1(bgr_tcalm005_301[kEpr],"TCalm005","stp_10/pdg_1")->Clone();
    fH1[1] = (TH1F*) gh2(bgr_tcalm005_dev[kEpr],"TCalm005","stp_10/pdg_1")->Clone();

    fH1[0]->Sumw2();
    fH1[1]->Sumw2();
					// Pad#1: the hit energy distributions
    fP1->cd(1);
    gPad->SetLogx(1);

    fH1[0]->SetTitle("PDG code of the particle hitting the first disk");
    fH1[0]->GetXaxis()->SetTitle("PDG code");
    fH1[0]->Draw("h");

    fH1[0]->Draw("h");	
    DrawPaveLabelNDC(label,"Protons: old PABS",0.10,0.85,0.55,0.95,52);
    label->SetTextSize(0.3);

    fP1->cd(2);
    gPad->SetLogx(1);

    fH1[1]->SetTitle("PDG code of the particle hitting the first disk");
    fH1[1]->GetXaxis()->SetTitle("PDG code");
    fH1[1]->Draw("h");

    fH1[1]->Draw("h");	

    DrawPaveLabelNDC(label,"Protons: new PABS",0.10,0.85,0.55,0.95,52);
    label->SetTextSize(0.3);
  }			
//-----------------------------------------------------------------------------
// Fig.204: protons - distributions for the number of straw hits per event
//-----------------------------------------------------------------------------
  if (Figure == 204) {
    TPaveLabel* label;
    TLegend*    leg;

    fCanvas  = new_slide(name,title,1,1,1200,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fH1[0] = (TH1F*) gh1(bgr_tcalm002_301[kEpr],"TCalm002","evt_0/nsh_0")->Clone();
    fH1[1] = (TH1F*) gh2(bgr_tcalm002_dev[kEpr],"TCalm002","evt_0/nsh_0")->Clone();

    // fH1[0]->Sumw2();
    // fH1[1]->Sumw2();
					// Pad#1: the hit energy distributions
    fP1->cd(1);
    gPad->SetLogy(1);

    fH1[1]->SetTitle("Number of straw hits per event");
    fH1[1]->GetXaxis()->SetRangeUser(0,119);
    fH1[1]->GetXaxis()->SetTitle("N(hits)");
    fH1[1]->Draw("h");

    fH1[0]->SetFillColor(kBlue-7);
    fH1[0]->SetFillStyle(3003);
    fH1[0]->Draw("sames,h");	

    leg = new TLegend(0.5,0.73,0.7,0.79,"");
    leg->AddEntry(fH1[1],"new PABS","f");
    leg->AddEntry(fH1[0],"old PABS","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    // DrawPaveLabelNDC(label,"Protons: new PABS",0.10,0.85,0.55,0.95,52);
    // label->SetTextSize(0.3);
  }			

}
