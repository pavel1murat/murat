//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
/*
auto* f  = TFile::Open("hst.mu2e.kpp00s000r104n013.make_det_003.125470_000001.root")
auto* h1 = (TH2F*) f->Get("//Ana/PanelDtAna/tc_01/edep_vs_pnl");
auto* hp = h1->ProjectionY("py",1,1);
*/
#include <cmath>

#include <iostream>
#include <format>
#include "TFile.h"
#include "TF1.h"
#include "TH2.h"
#include "TFitResult.h"
#include "TCanvas.h"

#include "murat/util/fit_edep.hh"


//-----------------------------------------------------------------------------
double fit_edep::fit_function(double* X, double* P) {
    double f = 0;
    
    double x = X[0];
    
    if (x < P[1]) {
      f = 0;
    }
    else {
      f = P[0]*(1-exp(-(x-P[1])/P[2]))*exp(-(x-P[1])/P[3]);
    }
  
    return f;
  }
  
//-----------------------------------------------------------------------------
// different from the caller, channel number if within the panel , 0 <= channel < 96
// 'Panel = -1 : all panels
//-----------------------------------------------------------------------------
void fit_edep::fit_histogram(TH1* Hist, const char* FitOpt, double XMin, double XMax, fit_result_t* Fr) {

  f_fit = new TF1("f_fit",fit_function,0,0.01,4);
//-----------------------------------------------------------------------------
// perform fit
//-----------------------------------------------------------------------------
  f_fit->SetParameter(0,2500);
  f_fit->SetParameter(1,0.0002);
  f_fit->SetParameter(2,0.0003);
  f_fit->SetParameter(3,0.0013);

  TFitResultPtr tfr = Hist->Fit(f_fit,FitOpt,"",XMin,XMax);
  
  if (tfr->IsEmpty() or (! tfr->IsValid())) {
    Fr->msg = std::format("# FIT ERROR");
    return;
  }
//-----------------------------------------------------------------------------
// fit succeeded, extract parameters
//-----------------------------------------------------------------------------
  if (tfr->Ndf() > 0) {
    double sf2 = tfr->Chi2()/tfr->Ndf();
    double sf  = sqrt(sf2);
    
    for (int i=0; i<kNPar; i++) {
      Fr->p [i]  = tfr->Parameter(i);
      Fr->e [i]  = tfr->Error(i)*sf;
    }          
    
    Fr->chi2d = sf2;
    Fr->ndf   = tfr->Ndf();
  }
  else {
    Fr->msg = std::format(" TROUBLE : tfr->Ndf = %i",tfr->Ndf());
    
  }

  return;
}

//-----------------------------------------------------------------------------
int fit_edep::fit_panels(int IMin, int IMax, const char* FileName) {

  TFile* f  = TFile::Open(FileName);
  TH2F*  h1 = (TH2F*) f->Get("//Ana/PanelDtAna/tc_01/edep_vs_pnl");

  std::cout << " panel      p[0]      p[1]       p[2]       p[3]       chi2d       ndf\n";
  std::cout << " ---------------------------------------------------------------------\n";
  for (int i=IMin; i<IMax+1; i++) {
    std::string hname = std::format("py_{:03d}",i);
    fit_result_t* fr = &frr[i];
    h_edep[i] = h1->ProjectionY(hname.data(),i+1,i+1);
    std::string title = std::format("run {:6d}: panel {:3} eDep",fRunNumber,i);
    h_edep[i]->SetTitle(title.data());

    float qn = h_edep[i]->Integral();

    if (qn > 100) {
      fit_histogram(h_edep[i],"sql",0.0004,0.003,fr);
      if (fr->chi2d > 0) {
        std::cout << std::format(" {:03d} {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:5d}\n",
                                 i,
                                 fr->p[0],
                                 fr->p[1],
                                 fr->p[2],
                                 fr->p[3],
                                 // Fr->p[4],
                                 fr->chi2d,
                                 fr->ndf);
      }
      else {
        std::cout << std::format("panel {:03d}: {}\n",i,fr->msg);
      }
    }
    else {
      std::cout << std::format("panel {:03d}: ERROR not enough sttistics: integral:{}\n",i,qn);
    }
  }
  return 0;
}


//-----------------------------------------------------------------------------
void fit_edep::plot_summary_histograms() {
//-----------------------------------------------------------------------------
// and plot summary histogram in the end
//-----------------------------------------------------------------------------
  std::string title;

  title   = std::format("run {:6d}: <eDep> vs panel",fRunNumber);
  TH1F* h = new TH1F("h_sum",title.data(),216,0,216);

  title   = std::format("run {:6d} eDep_{{MP}} from the fit",fRunNumber);
  TH1F* h_all = new TH1F("h_all",title.data(),100,0,1.e-3);

  for (int i=0; i<216; i++) {
    float xmax(0);

    if (frr[i].chi2d > 0) {
      float p3 = frr[i].p[3];
      float p1 = frr[i].p[1];
      float p2 = frr[i].p[2];
      
      xmax = p1 + p2*std::log((p2+p3)/p2);
    }
    h->SetBinContent(i,xmax);
    h_all->Fill(xmax);
  }

  TCanvas* c = new TCanvas("c_all","c_all",1200,800);
  c->Divide(2,1);
  c->cd(1);
  h->Draw();
  c->cd(2);
  h_all->Draw();
}
