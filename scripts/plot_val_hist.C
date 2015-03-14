//-----------------------------------------------------------------------------
// hist name should be qualified with the full path within the file, for example
// "//neutronMixer/hNEvents"
//
// MU2E histogram files contain ntuples
//
// integral for R>36cm starts from bin=37
// 
//-----------------------------------------------------------------------------
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



//-----------------------------------------------------------------------------
// multiplicities of the background events
//-----------------------------------------------------------------------------
void plot_01(const char* Filename) {
  TH1* h;
  TCanvas* c = new TCanvas("c_plot_01","c_plot_01",0,0,800,600);
  c->Divide(2,2);

  c->cd(1);
  h = get_mu2e_hist(Filename,"//dioMixer/hNEvents");
  if (h) {
    h->SetTitle("number of mixed in DIO events");
    h->Draw();
  }

  c->cd(2);
  h = get_mu2e_hist(Filename,"//neutronMixer/hNEvents");
  if (h) {
    h->SetTitle("number of mixed in neutron events");
    h->Draw();
  }

  c->cd(3);
  h = get_mu2e_hist(Filename,"//protonMixer/hNEvents");
  if (h) {
    h->SetTitle("number of mixed in proton events");
    h->Draw();
  }

  c->cd(4);
  h = get_mu2e_hist(Filename,"//photonMixer/hNEvents");
  if (h) {
    h->SetTitle("number of mixed in photon events");
    h->Draw();
  }
}

//-----------------------------------------------------------------------------
// tracking efficiency vs cos(th)
//-----------------------------------------------------------------------------
void plot_02(const char* Filename) {
  TH1   *h1, *h2, *h3;

  TCanvas* c = new TCanvas("c_plot_02","Tracking eff vs eta, no cuts",0,0,800,800);
  c->Divide(1,3);

  c->cd(1);
  h1 = get_mu2e_hist(Filename,"//TCalm001/Hist/evt_0/ce_costh");
  if (h1) {
    h1->SetTitle("number of mixed in DIO events");
    h1->Draw();
  }
					// events with found tracks
  c->cd(2);
  h2 = get_mu2e_hist(Filename,"//TCalm001/Hist/evt_1/ce_costh");
  if (h2) {
    h2->Draw();
  }

  c->cd(3);
  h3 = (TH1*) h1->Clone("h3");
  h3->Reset();
  h3->Divide(h2,h1);
  if (h3) {
    h3->SetMarkerStyle(20);
    h3->SetMarkerSize(1);
    h3->Draw("pe");
  }
}

//-----------------------------------------------------------------------------
//  efficiency of the R > 35 cm cut
//
// 1. position the first disk
//-----------------------------------------------------------------------------
void plot_tcalm003_001(const char* Filename = "me000001.tcalm003.hist") {
  TH1   *h_qtot, *h1, *h2, *h3;

  TCanvas* c = new TCanvas("c_plot_tcalm003_001","Efficiency vs disk position, no cuts",0,0,1200,500);
  c->Divide(2,1);

  const int np = 40;
  float   eff [np], err[np], x[np], ex[np], eff2[np], err2[np];

  h_qtot  = get_mu2e_hist(Filename,"//TCalm003/Hist/evt_0/ce_costh");
  double qtot = h_qtot->GetEntries();

  double step = 5.; // in cm

  double zmax = np*step;

  h_eff = new TH1F("h_eff","eff",1,0,150);
  

  for (int i=0; i<np; i++) {
    h1     = get_mu2e_hist(Filename,Form("//TCalm003/Hist/evt_0/rt_%02i",i));
    x  [i] = 5*i;
    ex[i]  = 0;
    eff[i] = h1->Integral(37,100)/qtot;
    err[i] = sqrt(h1->Integral(37,100))/qtot;

//     h_eff->SetBinContent(i+1,eff[i]);
//     h_eff->SetBinError  (i+1,err[i]);
  }

  c->cd(1);
  h_eff->GetYaxis()->SetRangeUser(0.15,0.5);
  //  h_eff->Draw();

  TH2F* h2f = new TH2F("h2f","h2f",1,0.,zmax,1,0.,1.);

  h2f->Draw();
  TGraphErrors* gre = new TGraphErrors(np, x,eff,ex,err);
  gre->Draw("APe,same");

  TF1 *f1 = new TF1("f1","[0]+[1]*sin([2]*x+[3])",0,zmax);
  f1->SetParameters(0.3,0.05,0.,0.);
  gre->Fit("f1","","ALP",0,zmax-10);

  double p2 = f1->GetParameter(2);
  double p3 = f1->GetParameter(3);

  double offset = -(p3-TMath::ACos(0))/p2;

  printf(" disk#1: offset = %10.3f\n",offset);

//-----------------------------------------------------------------------------
// position the second disc
//-----------------------------------------------------------------------------
  c->cd(2);

  for (int i=0; i<np; i++) {
    h1     = get_mu2e_hist(Filename,Form("//TCalm003/Hist/evt_1/rt_%02i",i));
    x   [i] = 5*i;
    ex  [i]  = 0;
    eff2[i] = h1->Integral(37,100)/qtot;
    err2[i] = sqrt(h1->Integral(37,100))/qtot;
  }
  
  TGraphErrors* gre2 = new TGraphErrors(np, x,eff2,ex,err2);
  gre2->Draw("APe,same");

  f1->SetParameters(0.1,0.1,0.,0.);
  gre2->Fit("f1","","ALP",0,zmax-10.);

  p2 = f1->GetParameter(2);
  p3 = f1->GetParameter(3);

  offset = -(p3-TMath::ACos(0))/p2;

  printf(" disk#2: offset = %10.3f\n",offset);
}

//-----------------------------------------------------------------------------
//  efficiency of the R > 35 cm cut
//
// 1. position the first disk
//-----------------------------------------------------------------------------
void plot_tcalm003_002(const char* Filename = "me000001.tcalm003.hist") {
  TH1   *h_qtot, *h1, *h2, *h3;

  TCanvas* c = new TCanvas("c_plot_tcalm003_001","Efficiency vs disk position, no cuts",0,0,1200,500);
  c->Divide(2,1);

  const int np = 40;
  float   eff [np], err[np], x[np], ex[np], eff2[np], err2[np];

  h_qtot  = get_mu2e_hist(Filename,"//TCalm003/Hist/evt_0/ce_costh");
  double qtot = h_qtot->GetEntries();

  double step = 5.; // in cm

  double zmax = np*step;

  h_eff = new TH1F("h_eff","eff",1,0,150);
  

  for (int i=0; i<np; i++) {
    h1     = get_mu2e_hist(Filename,Form("//TCalm003/Hist/evt_0/rt_%02i",i));
    x  [i] = 5*i;
    ex[i]  = 0;
    eff[i] = h1->Integral(37,100)/qtot;
    err[i] = sqrt(h1->Integral(37,100))/qtot;

//     h_eff->SetBinContent(i+1,eff[i]);
//     h_eff->SetBinError  (i+1,err[i]);
  }

  c->cd(1);
  h_eff->GetYaxis()->SetRangeUser(0.15,0.5);
  //  h_eff->Draw();

  TH2F* h2f = new TH2F("h2f","h2f",1,0.,zmax,1,0.,1.);

  h2f->Draw();
  TGraphErrors* gre = new TGraphErrors(np, x,eff,ex,err);
  gre->Draw("APe,same");

  TF1 *f1 = new TF1("f1","[0]+[1]*sin([2]*x+[3])",0,zmax);
  f1->SetParameters(0.3,0.05,0.,0.);
  gre->Fit("f1","","ALP",0,zmax-10);

  double p2 = f1->GetParameter(2);
  double p3 = f1->GetParameter(3);

  double offset = -(p3-TMath::ACos(0))/p2;

  printf(" disk#1: offset = %10.3f\n",offset);

//-----------------------------------------------------------------------------
// position the second disc
//-----------------------------------------------------------------------------
  c->cd(2);

  for (int i=0; i<np; i++) {
    h1     = get_mu2e_hist(Filename,Form("//TCalm003/Hist/evt_1/rt_%02i",i));
    x   [i] = 5*i;
    ex  [i]  = 0;
    eff2[i] = h1->Integral(37,100)/qtot;
    err2[i] = sqrt(h1->Integral(37,100))/qtot;
  }
  
  TGraphErrors* gre2 = new TGraphErrors(np, x,eff2,ex,err2);
  gre2->Draw("APe,same");

  f1->SetParameters(0.1,0.1,0.,0.);
  gre2->Fit("f1","","ALP",0,zmax-10.);

  p2 = f1->GetParameter(2);
  p3 = f1->GetParameter(3);

  offset = -(p3-TMath::ACos(0))/p2;

  printf(" disk#2: offset = %10.3f\n",offset);
}

