//-----------------------------------------------------------------------------
// DiskNumber : 0 or 1
//-----------------------------------------------------------------------------
void plot_e_vs_r(const char* Filename, const char* Folder = "evt_0") {

  TFile* f = TFile::Open(Filename);
					// N(crystals vs radius)

  TH1F* h_r  = gh1(Filename,"CalAna","rc_0");

  TH1D* h_er = (TH1D*) gh1(Filename,"CalAna",Form("%s/ecr_vs_r_0",Folder));
  TH1D* h_er1 = (TH1D*) gh1(Filename,"CalAna",Form("%s/ecr_vs_r_1",Folder));

  TH1F* h    = (TH1F*) h_r->Clone("ecr_vs_r_0");
  TH1F* h1   = (TH1F*) h_r->Clone("ecr_vs_r_1");

  int  nev = gh1(Filename,"CalAna","evt_0/rv")->GetEntries();

  h->Reset();

  //  h->Divide(h_er,h_r);

  double y, ey, y1, ey1;
  int nb = h_er->GetNbinsX();

  for (int i=1; i<=nb; i++) {
    if (h_r->GetBinContent(i) == 0) {
      y  = 0;
      ey = 0;
      y1  = 0;
      ey1 = 0;
    }
    else {
      y  = h_er->GetBinContent(i)/(h_r->GetBinContent(i)+1.e-12)/(nev+1.e-12);
      ey = h_er->GetBinError  (i)/(h_r->GetBinContent(i)+1.e-12)/(nev+1.e-12);

      y1  = h_er1->GetBinContent(i)/(h_r->GetBinContent(i)+1.e-12)/(nev+1.e-12);
      ey1 = h_er1->GetBinError  (i)/(h_r->GetBinContent(i)+1.e-12)/(nev+1.e-12);
    }

    h->SetBinContent(i, y);
    h->SetBinError  (i,ey);

    h1->SetBinContent(i, y1);
    h1->SetBinError  (i,ey1);
  }

  h->SetTitle("<E(crystal)> vs radius");
  h->GetXaxis()->SetTitle("R, mm");
  h->GetYaxis()->SetTitle("<E>, MeV");
  h->GetXaxis()->SetRangeUser(300.,699.);
  h->SetMarkerStyle(20);
  h->Draw();

  h1->SetMarkerStyle(24);
  h1->Draw("same");
//-----------------------------------------------------------------------------
// plot legend
//-----------------------------------------------------------------------------
  TLegend* leg = new TLegend(0.5,0.7,0.7,0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h,"disk #0","pe");
  leg->AddEntry(h1,"disk #1","pe");
  leg->Draw();
}

//-----------------------------------------------------------------------------
// DiskNumber : 0 or 1, compare occupancies on two disks
//-----------------------------------------------------------------------------
void plot_nhits_vs_r(const char* Filename, const char* Folder = "evt_0") {

  TString title;

  TFile* f = TFile::Open(Filename);
					// N(crystals vs radius)

  TH1F* h_r = gh1(Filename,"CalAna","rc_0");

  TH1D* h_er0 = (TH1D*) gh1(Filename,"CalAna",Form("%s/nh_vs_r_0",Folder));
  TH1D* h_er1 = (TH1D*) gh1(Filename,"CalAna",Form("%s/nh_vs_r_1",Folder));

  TH1F* h0    = (TH1F*) h_r->Clone("nhits_per_crystal_vs_r_0");
  TH1F* h1    = (TH1F*) h_r->Clone("nhits_per_crystal_vs_r_1");

  h0->Reset();
  h1->Reset();

  if      (strcmp(Folder,"evt_0") == 0) title = "<N hits per crystal> vs radius, E > 0";
  else if (strcmp(Folder,"evt_1") == 0) title = "<N hits per crystal> vs radius, E > 0.1 MeV";
  else if (strcmp(Folder,"evt_2") == 0) title = "<N hits per crystal> vs radius, E > 0.5 MeV";
  else if (strcmp(Folder,"evt_3") == 0) title = "<N hits per crystal> vs radius, E > 1.0 MeV";

  int  nev = gh1(Filename,"CalAna","evt_0/rv")->GetEntries();

  float y0, y1, ey0, ey1;

 //  h->Divide(h_er,h_r);

  int nb = h_er0->GetNbinsX();

  for (int i=1; i<=nb; i++) {
    if (h_r->GetBinContent(i) == 0) {
      y0  = 0;
      ey0 = 0;
      y1  = 0;
      ey1 = 0;
    }
    else {
      y0  = h_er0->GetBinContent(i)/(h_r->GetBinContent(i)+1.e-12)/(nev+1.e-12);
      ey0 = h_er0->GetBinError  (i)/(h_r->GetBinContent(i)+1.e-12)/(nev+1.e-12);

      y1  = h_er1->GetBinContent(i)/(h_r->GetBinContent(i)+1.e-12)/(nev+1.e-12);
      ey1 = h_er1->GetBinError  (i)/(h_r->GetBinContent(i)+1.e-12)/(nev+1.e-12);
    }

    h0->SetBinContent(i, y0);
    h0->SetBinError  (i,ey0);

    h1->SetBinContent(i, y1);
    h1->SetBinError  (i,ey1);
  }

  h0->GetXaxis()->SetTitle("R, mm");
  h0->GetXaxis()->SetRangeUser(300.,699.);
  h0->SetTitle(title);
  h0->SetStats(0);
  h0->SetMarkerSize(1);
  h0->SetMarkerStyle(20);
  h0->SetMaximum(10);
  h0->Draw("pe");

  h1->SetMarkerSize(1);
  h1->SetMarkerStyle(24);
  h1->Draw("same,pe");
//-----------------------------------------------------------------------------
// plot legend
//-----------------------------------------------------------------------------
  TLegend* leg = new TLegend(0.5,0.7,0.7,0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h0,"disk #0","pe");
  leg->AddEntry(h1,"disk #1","pe");
  leg->Draw();
  
}


//-----------------------------------------------------------------------------
// total number of hits vs the readout threshold
//-----------------------------------------------------------------------------
void plot_nhits_vs_thr(const char* Filename) {
  
  TString title;

  TFile* f = TFile::Open(Filename);

  enum {kNHits = 6};

  double thr   [kNHits] = { 0., 0.1, 0.5, 1.0, 1.5, 2.0};

  double qnhits[kNHits];
					// N(crystals vs radius)

  qnhits[0] = gh1(Filename,"CalAna","evt_0/nhits_0")->GetMean()+gh1(Filename,"CalAna","evt_0/nhits_1")->GetMean();
  qnhits[1] = gh1(Filename,"CalAna","evt_1/nhits_0")->GetMean()+gh1(Filename,"CalAna","evt_1/nhits_1")->GetMean();
  qnhits[2] = gh1(Filename,"CalAna","evt_2/nhits_0")->GetMean()+gh1(Filename,"CalAna","evt_2/nhits_1")->GetMean();
  qnhits[3] = gh1(Filename,"CalAna","evt_3/nhits_0")->GetMean()+gh1(Filename,"CalAna","evt_3/nhits_1")->GetMean();
  qnhits[4] = gh1(Filename,"CalAna","evt_4/nhits_0")->GetMean()+gh1(Filename,"CalAna","evt_4/nhits_1")->GetMean();
  qnhits[5] = gh1(Filename,"CalAna","evt_5/nhits_0")->GetMean()+gh1(Filename,"CalAna","evt_5/nhits_1")->GetMean();

  TGraph* gr = new TGraph(kNHits,thr,qnhits);

  gr->GetXaxis()->SetTitle("Readout Threshold, MeV");
  gr->SetMarkerSize(1);
  gr->SetMarkerStyle(20);
  gr->SetTitle("Total Number of hits above the threshold"); 
  gr->Draw("ALP");
}

//-----------------------------------------------------------------------------
// total number of clusters vs the energy threshold
//-----------------------------------------------------------------------------
void plot_ncl_vs_thr(const char* Filename) {
  
  TString title;

  TFile* f = TFile::Open(Filename);

  enum {kNPoints = 4};

  double thr   [kNPoints] = { 0., 20, 50, 70};

  double qncl[kNPoints];
					// N(crystals vs radius)

  qncl[0] = gh1(Filename,"CalAna","evt_0/ncl")->GetMean();
  qncl[1] = gh1(Filename,"CalAna","evt_0/ncl20")->GetMean();
  qncl[2] = gh1(Filename,"CalAna","evt_0/ncl50")->GetMean();
  qncl[3] = gh1(Filename,"CalAna","evt_0/ncl70")->GetMean();

  TGraph* gr = new TGraph(kNPoints,thr,qncl);

  gr->GetXaxis()->SetTitle("Cluster Energy, MeV");
  gr->SetMarkerSize(1);
  gr->SetMarkerStyle(20);
  gr->SetTitle("Number of clusters above the threshold"); 
  gr->Draw("ALP");
}

//-----------------------------------------------------------------------------
// total number of clusters vs the energy threshold
//-----------------------------------------------------------------------------
void plot_etot(const char* Filename) {
  
  TString title;

  TFile* f = TFile::Open(Filename);

  TH1F* h0 = gh1(Filename,"CalAna","evt_0/etot_0");
  TH1F* h1 = gh1(Filename,"CalAna","evt_0/etot_1");


  h0->SetMaximum(h1->GetMaximum()*1.1);
  h0->GetXaxis()->SetTitle("Energy per disk, MeV");
  h0->Draw();

  h1->SetFillStyle(3002);
  h1->SetFillColor(kBlue-7);
  h1->Draw("sames");

//-----------------------------------------------------------------------------
// plot legend
//-----------------------------------------------------------------------------
  TLegend* leg = new TLegend(0.5,0.7,0.7,0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h0,"disk #0","f");
  leg->AddEntry(h1,"disk #1","f");
  leg->Draw();
}


//-----------------------------------------------------------------------------
void plot_cal_occupancies(const char* Filename) {
  
}
