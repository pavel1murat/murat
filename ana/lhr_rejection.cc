//

#include "TH1.h"
#include "TEnv.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "Stntuple/val/stntuple_val_functions.hh"
#include "Stntuple/alg/TEmuLogLH.hh"
#include <math.h>

#include "murat/ana/lhr_rejection.hh"

namespace {

  TCanvas* c(NULL);
  TVirtualPad  *p1(0), *p2(0);

  TH1F  *h_dt_e, *h_dt_m, *h_dt_es, *h_dt_ms;
  TH2F  *h_ep_vs_s_e, *h_ep_vs_s_m, *h_ep_vs_s_e1, *h_ep_vs_s_m1, *h_ep_vs_s_es, *h_ep_vs_s_ms;
  
  TH1F  *h_llhr_cal_e, *h_llhr_cal_m;
  TH1F  *h_prob_e, *h_prob_m;
};

//-----------------------------------------------------------------------------
// smear E, but loop over S bins as well, E/P is Y
//-----------------------------------------------------------------------------
void smear_ep_vs_s_hist(TH2F* H1, TH2F* Hs, double SigEE) {

  TRandom3* rn = new TRandom3();

  int nbx = H1->GetNbinsX();
  int nby = H1->GetNbinsY();

  double r, s, ep, ep_new ;

  int q = 100;

  for (int ix=0; ix<nbx; ix++) {
    for (int iy=0; iy<nby; iy++) {
      int ne = H1->GetCellContent(ix+1,iy+1);
    
      s  = H1->GetXaxis()->GetBinCenter(ix+1);
      ep = H1->GetYaxis()->GetBinCenter(iy+1);
    
      if (ne > 0) {
	int nn = q*ne;
      
	for (int ii=0; ii<nn; ii++) {
	  r       = rn->Gaus(0.,SigEE);
//-----------------------------------------------------------------------------
// take into account that de/e = dee*sqrt(e);
//-----------------------------------------------------------------------------
	  ep_new  = ep+r*sqrt(ep);
	  Hs->Fill(s,ep_new);
	}
      }
    }
  }
				// preserve the normalization
  Hs->Scale(1./q);
}

//-----------------------------------------------------------------------------
void smear_dt_hist(TH1F* H1, TH1F* Hs, double SigT) {

  TRandom3* rn = new TRandom3();

  int nb = H1->GetNbinsX();

  double r, dt, dt_new ;

  int q = 100;

  for (int i=0; i<nb; i++) {
    int ne = H1->GetBinContent(i+1);
    
    dt = H1->GetBinCenter(i+1);
    
    if (ne > 0) {
      int nn = q*ne;
      
      for (int ii=0; ii<nn; ii++) {
	r = rn->Gaus(0.,SigT);
//-----------------------------------------------------------------------------
// take into account that de/e = dee*sqrt(e);
//-----------------------------------------------------------------------------
	dt_new  = dt+r;
	Hs->Fill(dt_new);
      }
    }
  }

  Hs->Scale(1./q);
}

//-----------------------------------------------------------------------------
void lhr_rejection(int HistSet, double SigEE, double SigT, int NEvents) {
//-----------------------------------------------------------------------------
// from murat/ana/TTrackAnaModule.cc: 
// -----------------------
// TRK_25 : Set C tracks with P > 100 , E/P > 0, chi2(match) < 100
// TRK_19 : Set C tracks with E/P > 0 ; // and chi2(match) < 100
//-----------------------------------------------------------------------------
  const char* folder = "trk_25";

  char   fn_e[200], fn_m[200];
  
  const char* hist_dir = gEnv->GetValue("mu2e.HistDir","");

  if (HistSet == 1412) {
    sprintf(fn_e,"%s/v4_2_1/e00s1412.track_ana.hist",hist_dir);
    sprintf(fn_m,"%s/v4_2_1/m00s1412.track_ana.hist",hist_dir);
  }
  else if (HistSet == 0041) {
    sprintf(fn_e,"%s/v4_2_4/e00s0041.track_ana.hist",hist_dir);
    sprintf(fn_m,"%s/v4_2_4/m00s0041.track_ana.hist",hist_dir);
  }
  else if (HistSet == 571) {
    sprintf(fn_e,"%s/v5_7_0/e00s5710.track_ana.hist",hist_dir);
    sprintf(fn_m,"%s/v5_7_0/m00s5710.track_ana.hist",hist_dir);
  }
//-----------------------------------------------------------------------------
// allow multiple histograms with the same name in the same directory
//-----------------------------------------------------------------------------
  TH1::AddDirectory(kFALSE);
//-----------------------------------------------------------------------------
// smear DT histograms
// h_dt_es, h_dt_ms are the smeared ones
//-----------------------------------------------------------------------------
  h_dt_e = (TH1F*) gh1(fn_e,"TrackAna",Form("%s/dt",folder));
  h_dt_m = (TH1F*) gh1(fn_m,"TrackAna",Form("%s/dt",folder));

  h_dt_es = (TH1F*) h_dt_e->Clone("h_dt_es");
  h_dt_ms = (TH1F*) h_dt_m->Clone("h_dt_ms");

  h_dt_es->Reset();
  h_dt_ms->Reset();

  smear_dt_hist(h_dt_e,h_dt_es,SigT);
  smear_dt_hist(h_dt_m,h_dt_ms,SigT);
//-----------------------------------------------------------------------------
// smear EP_vs_s histograms
//-----------------------------------------------------------------------------
  h_ep_vs_s_e = (TH2F*) gh2(fn_e,"TrackAna",Form("%s/ep_vs_path",folder));
  h_ep_vs_s_m = (TH2F*) gh2(fn_m,"TrackAna",Form("%s/ep_vs_path",folder));

  h_ep_vs_s_e1 = (TH2F*) ((TH1F*) h_ep_vs_s_e->Clone("h_ep_vs_s_e1")); // ->Rebin(2,2);
  h_ep_vs_s_m1 = (TH2F*) ((TH1F*) h_ep_vs_s_m->Clone("h_ep_vs_s_m1")); // ->Rebin(2,2);

  h_ep_vs_s_es = (TH2F*) h_ep_vs_s_e->Clone("h_ep_vs_s_es");
  h_ep_vs_s_ms = (TH2F*) h_ep_vs_s_m->Clone("h_ep_vs_s_ms");

  h_ep_vs_s_es->Reset();
  h_ep_vs_s_ms->Reset();

  smear_ep_vs_s_hist(h_ep_vs_s_e1,h_ep_vs_s_es,SigEE);
  smear_ep_vs_s_hist(h_ep_vs_s_m1,h_ep_vs_s_ms,SigEE);
//-----------------------------------------------------------------------------
// build the likelihood
//-----------------------------------------------------------------------------
  TEmuLogLH* llh = new TEmuLogLH();

  llh->InitEleDtHist(h_dt_es);
  llh->InitEleEpHist(h_ep_vs_s_es);

  llh->InitMuoDtHist(h_dt_ms);
  llh->InitMuoEpHist(h_ep_vs_s_ms);
//-----------------------------------------------------------------------------
// sample electron distributions and build LLHR distribution for electrons
//-----------------------------------------------------------------------------
  h_llhr_cal_e = new TH1F("h_llhr_cal_e","LLHR(cal) Electrons",500,-100,100);

  double llhr_ep, llhr_dt, llhr_cal;

  TEmuLogLH::PidData_t data;

  for (int i=0; i<NEvents; i++) {
    h_ep_vs_s_es->GetRandom2(data.fPath,data.fEp);
    data.fDt = h_dt_es->GetRandom();

    llhr_ep  = llh->LogLHREp(&data);
    llhr_dt  = llh->LogLHRDt(&data);
    llhr_cal = llhr_ep+llhr_dt;
					// debug
    if (llhr_cal > 80) {
      printf(" ------------- i = %i10 Electron LLHR_CAL = %12.5e\n",i,llhr_cal);
      printf("llhr_ep = %12.5e llhr_dt = %12.5e\n", llhr_ep,llhr_dt);
      printf("data.fPath, data.fEp, data.fDt = %12.5e %12.5e %12.5e\n",
	     data.fPath,data.fEp,data.fDt);
    }
    h_llhr_cal_e->Fill(llhr_cal);
  }
  h_llhr_cal_e->Scale(1./NEvents);
//-----------------------------------------------------------------------------
// sample muon distributions and build LLHR distribution for muons
//-----------------------------------------------------------------------------
  h_llhr_cal_m = new TH1F("h_llhr_cal_m","LLHR(cal) Muons    ",500,-100,100);
  
  for (int i=0; i<NEvents; i++) {
    h_ep_vs_s_ms->GetRandom2(data.fPath,data.fEp);
    data.fDt = h_dt_ms->GetRandom();

    llhr_ep  = llh->LogLHREp(&data);
    llhr_dt  = llh->LogLHRDt(&data);
    llhr_cal = llhr_ep+llhr_dt;

    h_llhr_cal_m->Fill(llhr_cal);
  }
  h_llhr_cal_m->Scale(1./NEvents);
//-----------------------------------------------------------------------------
// normalize likelihood distributions 
//-----------------------------------------------------------------------------
//  h_llhr_cal_e->DrawNormalized("",1);
  h_llhr_cal_e->Draw("");

  h_llhr_cal_m->SetFillColor(kBlue-7);
  h_llhr_cal_m->SetFillStyle(3002);

  //  h_llhr_cal_m->DrawNormalized("same",1);
  h_llhr_cal_m->Draw("same");

  h_prob_e = (TH1F*) h_llhr_cal_e->Clone("h_prob_e");
  h_prob_m = (TH1F*) h_llhr_cal_m->Clone("h_prob_m");

  h_prob_e->Reset();
  h_prob_m->Reset();

  int nb = h_prob_e->GetNbinsX();

  for (int i=0; i<nb; i++) {
    h_prob_e->SetBinContent(i+1,h_llhr_cal_e->Integral(i+1,nb));
    h_prob_m->SetBinContent(i+1,h_llhr_cal_m->Integral(i+1,nb));
  }
//-----------------------------------------------------------------------------
// plot histograms
//-----------------------------------------------------------------------------
  if (c == NULL) {
    c = new TCanvas("c","c",1400,800);
    c->Divide(2,1);

    c->cd(1);

    p1= gPad;
    p1->Divide(1,2);

    c->cd(2);
    p2 = gPad;

    p2->Divide(1,2);
  }
//-----------------------------------------------------------------------------
// 1. E/P distributions for electrons and muons after the cuts
//-----------------------------------------------------------------------------
  p1->cd(1);

  h_ep_vs_s_ms->SetStats(0);
  h_ep_vs_s_ms->SetTitle("");
  h_ep_vs_s_ms->GetXaxis()->SetTitle("E/P");
  h_ep_vs_s_ms->SetFillColor(kBlue-7);
  h_ep_vs_s_ms->SetFillStyle(3002);
  h_ep_vs_s_ms->DrawNormalized("",1);

  h_ep_vs_s_es->DrawNormalized("same",1);

  TLegend  /* *leg1,*/ *leg2, *leg3;

  // leg1 = new TLegend(0.6,0.75,0.8,0.85);
  // leg1->AddEntry(h_ep_ms,"muons","f");
  // leg1->AddEntry(h_ep_es,"electrons","f");
  // leg1->SetFillStyle(0);
  // leg1->SetFillColor(0);
  // leg1->SetBorderSize(0);
  // leg1->Draw();
//-----------------------------------------------------------------------------
// 3. distributions in Delta(T) for electrons and muons after the cuts
//-----------------------------------------------------------------------------
  p1->cd(2);
  h_dt_es->SetStats(0);
  h_dt_es->SetTitle("#Delta T = T_{trk}-T_{cal}");
  h_dt_es->GetXaxis()->SetTitle("#Delta T, ns");
  h_dt_es->DrawNormalized("",1);

  h_dt_ms->SetFillColor(kBlue-7);
  h_dt_ms->SetFillStyle(3002);
  h_dt_ms->DrawNormalized("same",1);

  leg2 = new TLegend(0.6,0.75,0.8,0.85);
  leg2->AddEntry(h_dt_es,"electrons","f");
  leg2->AddEntry(h_dt_ms,"muons","f");
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->Draw();
//-----------------------------------------------------------------------------
// 2. probability distributions for electrons and muons
//    why the electron efficiency doesn't saturate at 1 ?
//-----------------------------------------------------------------------------
  p2->cd(1);

  h_prob_e->SetTitle("logLHR = log(LH_{cal}(e)/LH_{cal}(#mu)) ");
  h_prob_e->GetXaxis()->SetTitle("log(LH_{cal}(e)/LH_{cal}(#mu)) ");
  h_prob_e->GetXaxis()->SetRangeUser(-20,19.9);
  h_prob_e->SetStats(0);
  h_prob_e->Draw();

  h_prob_m->SetStats(0);
  h_prob_m->Draw("same");

  leg3 = new TLegend(0.7,0.75,0.9,0.85);
  leg3->AddEntry(h_prob_e,"electrons","f");
  leg3->AddEntry(h_prob_m,"muons"    ,"f");
  leg3->SetFillStyle(0);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->Draw();
//-----------------------------------------------------------------------------
// 4. electron efficiency vs muon rejection
//-----------------------------------------------------------------------------
  p2->cd(2);

  gPad->SetLogy(1);
  gPad->SetGridy(1);

  double prob_e[1000], prob_m[1000];
//-----------------------------------------------------------------------------
// scale factor - N(trk_13)/N(trk_25)
//-----------------------------------------------------------------------------
  double qnm_25, qnm_13, qne_25, qne_13;

  qne_13 = gh1(fn_e,"TrackAna","trk_13/llhr_cal")->GetEntries();
  qne_25 = gh1(fn_e,"TrackAna",Form("%s/llhr_cal",folder))->GetEntries();

  qnm_13 = gh1(fn_m,"TrackAna","trk_13/llhr_cal")->GetEntries();
  qnm_25 = gh1(fn_m,"TrackAna",Form("%s/llhr_cal",folder))->GetEntries();

  for (int i=0; i<nb; i++) {
    prob_e[i] = h_prob_e->GetBinContent(i+1)*(qne_25/qne_13);
    prob_m[i] = 1./(h_prob_m->GetBinContent(i+1)+1.e-6)*(qnm_13/qnm_25);
  }

  TGraph* gr = new TGraph(nb,prob_e,prob_m);

  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1);
  gr->Draw("ALP");

}


