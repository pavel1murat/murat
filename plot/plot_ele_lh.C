/////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////
#ifndef __CINT__

#include "TCanvas.h"
#include "TH1.h"
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"

#endif

TH1* get_mu2e_hist(const char* Filename, const char* Histname) {
  TFile  *f;
  TH1    *hist(0);

  f = gROOT->GetFile(Filename);
  if (f == 0) f = TFile::Open(Filename);

  if (f) {
    hist = (TH1*) f->Get(Histname);
  }

  return hist;
}

namespace {

  struct Hist_t {
    TH1F* fEleDt;
    TH1F* fEleEp;
    TH1F* fMuoDt;
    TH1F* fMuoEp;
  
    TH1F* fMuoLH[2];
    TH1F* fEleLH[2];
    TH1F* fDelLH[2];
  } Hist;
  
  TEmuLogLH*  gLh;
}

//-----------------------------------------------------------------------------
plot_ele_lh(int NEvents = 1000) {
  
  TH1F* h_ele_dt = (TH1F*) get_mu2e_hist("../v2_1_6/hist/e0000001.tcalm002.hist","//TCalm002/Hist/trk_1/dt");
  TH1F* h_ele_ep = (TH1F*) get_mu2e_hist("../v2_1_6/hist/e0000001.tcalm002.hist","//TCalm002/Hist/trk_1/ep");

  TH1F* h_muo_dt = (TH1F*) get_mu2e_hist("../v2_1_6/hist/m0000001.tcalm002.hist","//TCalm002/Hist/trk_1/dt");
  TH1F* h_muo_ep = (TH1F*) get_mu2e_hist("../v2_1_6/hist/m0000001.tcalm002.hist","//TCalm002/Hist/trk_1/ep");

  gLh = new TEmuLogLH();

  gLh->SetEleDtHist(h_ele_dt);
  gLh->SetEleEpHist(h_ele_ep);
  gLh->SetMuoDtHist(h_muo_dt);
  gLh->SetMuoEpHist(h_muo_ep);

  Hist.fEleLH[0] = new TH1F("h_ele_lh_0","ele LH[0] - electron sample",500,-50,0);
  Hist.fMuoLH[0] = new TH1F("h_muo_lh_0","muo LH[0] - electron sample",500,-50,0);

  Hist.fEleLH[1] = new TH1F("h_ele_lh_1","ele LH[1] - muon sample",500,-50,0);
  Hist.fMuoLH[1] = new TH1F("h_muo_lh_1","muo LH[1] - muon sample",500,-50,0);

  Hist.fDelLH[0] = new TH1F("h_del_lh_0","del LH[0] - electron sample",1000,-50,50);
  Hist.fDelLH[1] = new TH1F("h_del_lh_1","del LH[1] - muon sample",1000,-50,50);

  double   ele_lh, muo_lh;

  TEmuLogLH::Data_t dat;

  Hist.fEleLH[0]->Reset();
  Hist.fEleLH[1]->Reset();
  Hist.fMuoLH[0]->Reset();
  Hist.fMuoLH[1]->Reset();
  Hist.fDelLH[0]->Reset();
  Hist.fDelLH[1]->Reset();

  for (int i=0; i<NEvents; i++) {
//-----------------------------------------------------------------------------
// simulate an "electron event"
//-----------------------------------------------------------------------------
    dat.fDt = h_ele_dt->GetRandom();
    dat.fEp = h_ele_ep->GetRandom();

    // calculate muon and electron likelihoods 

    ele_lh = gLh->LogLH(&dat,11);
    muo_lh = gLh->LogLH(&dat,13);

    Hist.fEleLH[0]->Fill(ele_lh);
    Hist.fMuoLH[0]->Fill(muo_lh);

    Hist.fDelLH[0]->Fill(ele_lh-muo_lh);
//-----------------------------------------------------------------------------
// simulate a "muon event"
//-----------------------------------------------------------------------------
    dat.fDt = h_muo_dt->GetRandom();
    dat.fEp = h_muo_ep->GetRandom();

    ele_lh = gLh->LogLH(&dat,11);
    muo_lh = gLh->LogLH(&dat,13);

    Hist.fEleLH[1]->Fill(ele_lh);
    Hist.fMuoLH[1]->Fill(muo_lh);

    Hist.fDelLH[1]->Fill(ele_lh-muo_lh);

    //    printf("r1, r2, i1, i2, p1, p2: = %10.4f %10.4f %5i %5i %12.4g %12.4g \n",r1, r2, i1,i2, p1, p2);

    //   //    printf("r1, r2, i1, i2, p1, p2: = %10.4f %10.4f %5i %5i %12.4g %12.4g \n",r1, r2, i1,i2, p1, p2);
  }

}
