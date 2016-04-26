///////////////////////////////////////////////////////////////////////////////
// to run:
// ------
//        x = new lhr_rejection; x->run(571,0.01,0.2,100000);
//
///////////////////////////////////////////////////////////////////////////////

#include "Stntuple/val/stntuple_val_functions.hh"
#include <math.h>

#include "murat/ana/lhr_rejection.hh"


//-----------------------------------------------------------------------------
lhr_rejection::lhr_rejection() : TObject() {
  c_ep_vs_s    = NULL;
  c_dt         = NULL;
  c_prob       = NULL;
  c_eff_vs_rej = NULL;

}

//-----------------------------------------------------------------------------
lhr_rejection::~lhr_rejection() {
}

//-----------------------------------------------------------------------------
// smear E, but loop over S bins as well, E/P is Y
//-----------------------------------------------------------------------------
void lhr_rejection::smear_ep_vs_s_hist(TH2F* H1, TH2F* Hs, double SigEE) {

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
void lhr_rejection::smear_dt_hist(TH1F* H1, TH1F* Hs, double SigT) {

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
int lhr_rejection::build_cal_llhr_histogram(TH2F* HistEP, TH1F* HistDT, int NEvents, TH1F* HistLLHR) {

  double llhr_ep, llhr_dt, llhr_cal;

  TEmuLogLH::PidData_t data;

  if (NEvents <= 0) {
    printf(" ERROR in build_cal_llhr_histogram: NEvents = %10i. BAIL OUT\n",NEvents);
    return -1;
  }

  for (int i=0; i<NEvents; i++) {
    HistEP->GetRandom2(data.fPath,data.fEp);
    data.fDt = HistDT->GetRandom();

    llhr_ep  = llh->LogLHREp(&data);
    llhr_dt  = llh->LogLHRDt(&data);
    llhr_cal = llhr_ep+llhr_dt;

    HistLLHR->Fill(llhr_cal);
  }
  
  HistLLHR->Scale(1./NEvents);

  return 0;
}

//-----------------------------------------------------------------------------
void lhr_rejection::run(int HistSet, double SigEE, double SigT, int NEvents) {
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
  else if (HistSet == 572) {
    sprintf(fn_e,"%s/v5_7_0/e40s5720.track_ana.hist",hist_dir);
    sprintf(fn_m,"%s/v5_7_0/m40s5720.track_ana.hist",hist_dir);
  }
//-----------------------------------------------------------------------------
// allow multiple histograms with the same name in the same directory
//-----------------------------------------------------------------------------
  TH1::AddDirectory(kFALSE);
//-----------------------------------------------------------------------------
// smear DT histograms, h_dt_es, h_dt_ms are the smeared ones
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

  h_ep_vs_s_e1 = (TH2F*) ((TH1F*) h_ep_vs_s_e->Clone("h_ep_vs_s_e1"));
  h_ep_vs_s_m1 = (TH2F*) ((TH1F*) h_ep_vs_s_m->Clone("h_ep_vs_s_m1"));

  h_ep_vs_s_es = (TH2F*) h_ep_vs_s_e->Clone("h_ep_vs_s_es");
  h_ep_vs_s_ms = (TH2F*) h_ep_vs_s_m->Clone("h_ep_vs_s_ms");

  h_ep_vs_s_es->Reset();
  h_ep_vs_s_ms->Reset();

  smear_ep_vs_s_hist(h_ep_vs_s_e1,h_ep_vs_s_es,SigEE);
  smear_ep_vs_s_hist(h_ep_vs_s_m1,h_ep_vs_s_ms,SigEE);
//-----------------------------------------------------------------------------
// build the likelihood and replace DT and EP histograms with the smeared ones
//-----------------------------------------------------------------------------
  llh = new TEmuLogLH();

  llh->InitEleDtHist(h_dt_es);
  llh->InitEleEpHist(h_ep_vs_s_es);

  llh->InitMuoDtHist(h_dt_ms);
  llh->InitMuoEpHist(h_ep_vs_s_ms);
//-----------------------------------------------------------------------------
// sample electron distributions and build LLHR distribution for electrons
//-----------------------------------------------------------------------------
  h_llhr_cal_e = new TH1F("h_llhr_cal_e","LLHR(cal) Electrons",500,-100,100);

  build_cal_llhr_histogram(h_ep_vs_s_es,h_dt_es,NEvents,h_llhr_cal_e);
//-----------------------------------------------------------------------------
// sample muon distributions and build LLHR distribution for muons
//-----------------------------------------------------------------------------
  h_llhr_cal_m = new TH1F("h_llhr_cal_m","LLHR(cal) Muons    ",500,-100,100);
  
  build_cal_llhr_histogram(h_ep_vs_s_ms,h_dt_ms,NEvents,h_llhr_cal_m);
//-----------------------------------------------------------------------------
// running integrals
//-----------------------------------------------------------------------------
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
// 0. likelihood distributions 
//-----------------------------------------------------------------------------
  c_llhr_cal = new TCanvas("c_llhr_cal","LLHR(CAL) : electrons vs muons",1400,900);
  c_llhr_cal->cd();

  h_llhr_cal_e->SetTitle("LLHR(CAL) : electrons vs muons");
  h_llhr_cal_e->GetXaxis()->SetRangeUser(-50,49.99);
  h_llhr_cal_e->Draw("");

  h_llhr_cal_m->SetFillColor(kBlue-7);
  h_llhr_cal_m->SetFillStyle(3002);
  h_llhr_cal_m->Draw("same");

  TLegend* leg0 = new TLegend(0.2,0.75,0.4,0.85);
  leg0->AddEntry(h_llhr_cal_e,"electrons","f");
  leg0->AddEntry(h_llhr_cal_m,"muons"    ,"f");
  leg0->SetFillStyle(0);
  leg0->SetFillColor(0);
  leg0->SetBorderSize(0);
  leg0->Draw();
  

// 1. E/P distributions for electrons and muons after the cuts
//-----------------------------------------------------------------------------
  c_ep_vs_s = new TCanvas("c_ep_vs_s","EP vs SPath",1400,900);
  c_ep_vs_s->cd();
  
  h_ep_vs_s_ms->SetStats(0);
  h_ep_vs_s_ms->SetTitle("");
  h_ep_vs_s_ms->GetXaxis()->SetRangeUser(0,450);
  h_ep_vs_s_ms->GetXaxis()->SetTitle("E/P");
  h_ep_vs_s_ms->GetYaxis()->SetRangeUser(0,1.1);
  h_ep_vs_s_ms->SetMarkerColor(kBlue+3);
  h_ep_vs_s_ms->SetFillStyle(3002);
  h_ep_vs_s_ms->DrawNormalized("",1);

  h_ep_vs_s_es->SetMarkerColor(kRed+1);
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
// 2. distributions in Delta(T) for electrons and muons after the cuts
//-----------------------------------------------------------------------------
  c_dt = new TCanvas("c_dt","DT",1400,900);
  c_dt->cd();
  
  h_dt_es->SetStats(0);
  h_dt_es->SetTitle("#Delta T = T_{trk}-T_{cal}");
  h_dt_es->GetXaxis()->SetTitle("#Delta T, ns");
  h_dt_es->GetXaxis()->SetRangeUser(-10,10);
  //  h_dt_es->GetYaxis()->SetRangeUser(0,0.08);
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
// 3. probability distributions for electrons and muons
//    why the electron efficiency doesn't saturate at 1 ?
//-----------------------------------------------------------------------------
  c_prob = new TCanvas("c_prob","prob",1400,900);
  c_prob->cd();
  
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
  c_eff_vs_rej = new TCanvas("c_eff_vs_rej","Electron efficiency vs muon rejection",1400,900);
  c_eff_vs_rej->cd();

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

    if (prob_e[i] > 0.8) {
      printf(" i, prob_e[i], prob_m[i]: %5i %10.3f %10.3f\n",i, prob_e[i], prob_m[i]);
    }
  }

  TGraph* gr = new TGraph(nb,prob_e,prob_m);

  gr->SetTitle("Electron ID efficiency vs muon rejection");
    
  gr->GetXaxis()->SetTitle("electron ID efficiency");
  gr->GetYaxis()->SetTitle("muon rejection factor" );

  gr->GetXaxis()->SetRangeUser(0.8,1);
  gr->GetYaxis()->SetRangeUser(10,2.e4);

  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1);
  gr->Draw("ALP");

}


