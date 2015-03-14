//////////////////////////////////////////////////////////////////////////////
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
//  efficiency of the R > 35 cm cut
//
// 1. position the first disk
//-----------------------------------------------------------------------------
void plot_tcalm003_001(const char* Filename = "me000001.tcalm003.hist") {
  TH1   *h_qtot, *h1, *h2, *h3;

  TCanvas* c = new TCanvas("c_plot_tcalm003_001","Efficiency vs disk position, no cuts",0,0,1200,500);
  c->Divide(2,1);

  const int np = 60;
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
// for each plane plot efficiency for its misses 
//-----------------------------------------------------------------------------
void plot_tcalm003_002(const char* Filename = "me000001.tcalm003.hist", int IPL1) {
  TH1   *h_qtot, *h1, *h2, *h3;

  TCanvas* c = new TCanvas("c_plot_tcalm003_002","Efficiency vs disk position, no cuts",0,0,1200,500);
  c->Divide(2,1);

  const int np = 60;
  float   eff [np], err[np], x[np], ex[np], eff2[np], err2[np];
//-----------------------------------------------------------------------------
// normalize to the number of tracks - 1 track/event at max...
//-----------------------------------------------------------------------------
  h_qtot  = get_mu2e_hist(Filename,"//TCalm003/Hist/trk_0/tdip");

  double qtot = h_qtot->GetEntries();

  double step = 5.; // in cm

  double zmax = np*step;

 //  for (int i1=0; i1<np; i1++) {

  h0     = get_mu2e_hist(Filename,Form("//TCalm003/Hist/evt_0/rt_%02i",IPL1));

  double eff1 = h0->Integral(37,100)/qtot;

  TH1F* h_eff = (TH1F*) h0->Clone("h_eff");

  for (int i2=IPL1+1; i2<np; i2++) {
    h1     = get_mu2e_hist(Filename,Form("//TCalm003/Hist/evt_0/rt0_%02i_%02i",IPL1,i2));
    x [i2] = 5*i2;
    ex[i2]  = 0;
    eff[i2] = h1->Integral(37,100)/qtot;
    err[i2] = sqrt(h1->Integral(37,100))/qtot;

    eff2[i2] = h1->Integral(37,100)/qtot  + eff1;
    err2[i2] = sqrt(h1->Integral(37,100)+h0->Integral(37,100))/qtot;
    
  }

  c->cd(1);

  TH2F* h2f = new TH2F("h2f","h2f",1,0.,zmax,1,0.,1.);
  h2f->Draw();

  TGraphErrors* gre = new TGraphErrors(np, x,eff,ex,err);
  gre->SetTitle("Disk 2 acceptance vs Z_2, Z_1 = 50 cm");
  gre->Draw("APe,same");

  gre->Fit("gaus","","ALP,same",IPL1*step+20,IPL1*step+110);

  c->cd(2);
  TGraphErrors* gre2 = new TGraphErrors(np, x,eff2,ex,err2);
  gre2->SetTitle("Total calorimeter acceptance vs Z_2, Z_1 = 50 cm");
  gre2->Draw("APe,same");


  printf("eff(plane=%2i) = %8.3f\n",IPL1,eff1);
}
//-----------------------------------------------------------------------------
// for each plane plot efficiency for its misses 
//-----------------------------------------------------------------------------
void plot_tcalm003_003(const char* Filename = "me000001.tcalm003.hist") {
  TH1   *h_qtot, *h1, *h2, *h3;

  TCanvas* c = new TCanvas("c_plot_tcalm003_003","Efficiency vs disk position, no cuts",0,0,1200,500);
  c->Divide(2,1);

  const int np = 60;
  float   eff [np], err[np], x[np], ex[np], eff2[np], err2[np], d2[np];

  double  eff1, eff2_max, eff_2, d2_max;

//-----------------------------------------------------------------------------
// normalize to the number of tracks - 1 track/event at max...
//-----------------------------------------------------------------------------
  h_qtot  = get_mu2e_hist(Filename,"//TCalm003/Hist/trk_0/tdip");

  double qtot = h_qtot->GetEntries();

  double step = 5.; // in cm

  double zmax = np*step;

 //  for (int i1=0; i1<np; i1++) {

  for (int i1=0; i1<40; i1++) {

    h0 = get_mu2e_hist(Filename,Form("//TCalm003/Hist/evt_0/rt_%02i",i1));

    double eff1 = h0->Integral(37,100)/qtot;

    eff2_max = 0;

    for (int i2=i1+1; i2<np; i2++) {
      h1    = get_mu2e_hist(Filename,Form("//TCalm003/Hist/evt_0/rt0_%02i_%02i",i1,i2));
      eff_2 = h1->Integral(37,100)/qtot;

      if (eff_2 > eff2_max) {
	eff2_max  = eff_2;
	d2_max    = (i2-i1)*5;
      }
    }
    eff[i1] = eff1+eff2_max;
    err[i1] = 0.01;
    x  [i1] = step*i1;
    ex [i1] = 0;
    d2 [i1] = d2_max;
  }

  c->cd(1);

  TH2F* h2f = new TH2F("h2f","h2f",1,0.,zmax,1,0.,1.);
  h2f->Draw();

  TGraphErrors* gre = new TGraphErrors(np, x,eff,ex,err);
  gre->SetTitle("Efficiency vs Z_{1}");
  gre->GetYaxis()->SetRangeUser(0.8,1.1);
  gre->SetMarkerStyle(20);
  gre->SetMarkerSize(1);
  gre->Draw("APe,same");

  c->cd(2);

  h2f = new TH2F("h2f2","h2f2",1,0.,zmax,1,0.,100.);
  h2f->Draw();

  TGraphErrors* gre2 = new TGraphErrors(np, x,d2,ex,err);
  gre2->SetTitle("Distance between the calorimeter disks");
  gre2->GetYaxis()->SetRangeUser(0.,100.);
  gre2->SetMarkerStyle(20);
  gre2->SetMarkerSize(1);
  gre2->Draw("APe,same");

}
