//
#ifndef __lhr_rejection_hh__
#define __lhr_rejection_hh__

// histograms for E/Mu + MIXP-x1
// histset: 1412 for offline v4_2_1
//          0041 for offline v4_2_4 (new dataset naming conventions)

#include "TH1.h"
#include "TEnv.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "Stntuple/alg/TEmuLogLH.hh"

class lhr_rejection : public TObject {
public:
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
  TCanvas* c;
  TVirtualPad  *p1, *p2;

  TH1F  *h_dt_e, *h_dt_m, *h_dt_es, *h_dt_ms;
  TH2F  *h_ep_vs_s_e, *h_ep_vs_s_m, *h_ep_vs_s_e1, *h_ep_vs_s_m1, *h_ep_vs_s_es, *h_ep_vs_s_ms;
  
  TH1F  *h_llhr_cal_e, *h_llhr_cal_m;
  TH1F  *h_prob_e, *h_prob_m;
  
  TEmuLogLH* llh;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
  lhr_rejection();
  ~lhr_rejection();

  void smear_ep_vs_s_hist(TH2F* H1, TH2F* Hs, double SigEE);
  void smear_dt_hist     (TH1F* H1, TH1F* Hs, double SigT );
  
  int  build_cal_llhr_histogram(TH2F* HistEP, TH1F* HistDT, int NEvents, TH1F* HistLLHR);
  
  void run(int HistSet, double SigEE, double SigT, int NEvents = 10000);
};


#endif
