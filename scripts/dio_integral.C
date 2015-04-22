//
//-----------------------------------------------------------------------------
struct Data_t {
  double  fCnvNGen;
  TH1F*   fCnvMom;
  TH1F*   fCnvMomNorm;
  double  fDioNGen;
  TH1F*   fDioMom;
  TH1F*   fDioMomNorm;
  TH1F*   fDioIntegral;
};

Data_t  d5331;

TCanvas     *c_sens (0), *c_dio_bgr (0);

double mu_capture_prob(0.609);
					// total number of 8 GeV protons on target in 3 years

double n_prot      (1.2e20 * 3.);

					// assumed, for illustration purposes, B(mu->e) conversion 
double br_mu2e_conv(1.e-16);

double NDio(0), NCnv(0);

double NStoppedMuons, bin_width;

					// N(DIO) in [100 < E < 105]
double qn_dio, qn_cnv;

double PMin;

//-----------------------------------------------------------------------------
void check_dio_norm() {

  int nbins = 100;

  TH1D* h = new TH1D("h_dio", "DIO end spectrum",1100,0,110);


  double me(0.511);
  double p, e, wt;

  for (int i=10; i<50; i++) {

    p = 100+(0.5+i)*0.1;

    e = sqrt(p*p+me*me);

    wt = TStntuple::DioWeightAl(e);

    h->Fill(p,wt);
    //    h->Fill(p);

    printf(" p = %12.4f wt = %12.4f \n",p,wt);
  }

  TFile* f = TFile::Open("/grid/fermiapp/mu2e/users/murat/hist/mu2e/Czarnecki.root");

  TH1D* h_dio = (TH1D*) f->Get("Spectrum");

  h_dio->Draw();

  h->SetFillStyle(3001);
  h->SetFillColor(kBlue-7);

  // h->Draw();

  gPad->SetLogy(1);

   h->Draw("sames");
}


//-----------------------------------------------------------------------------
void plot_sensitivity() {
  while (1) {
    c_sens = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("mu2e_sens");
    if (c_sens == NULL) break;
    delete c_sens;
  }

  c_sens = new TCanvas("mu2e_sens","Mu2e Sensitivity",1200,800);
  c_sens->cd();

  TH2F* h2 = new TH2F("h2","",1,101,106,1,0,0.6);

  h2->SetStats(0);
  h2->GetYaxis()->SetTitle(Form("N /%5.3f MeV/c",bin_width));
  h2->GetXaxis()->SetTitle(Form("Momentum, MeV/c"));
  h2->Draw();
//-----------------------------------------------------------------------------
// draw "signal" box"
//-----------------------------------------------------------------------------
  PMin= 103.60; // for BGR 

  TPave* pave = new TPave(0.1,0.1,0.9,0.9,1);
  pave->SetX1(PMin);
  pave->SetX2(105.0);
  pave->SetY1(0.0);
  pave->SetY2(0.4);
  pave->SetLineWidth(1);

  pave->SetFillStyle(3003);
  pave->SetFillColor(kGreen-10);
  
  pave->Draw();

  d5331.fDioMomNorm->SetStats(0);
  d5331.fDioMomNorm->GetXaxis()->SetRangeUser(101,105.9999);
  d5331.fDioMomNorm->GetYaxis()->SetRangeUser(0,0.34999);

  d5331.fDioMomNorm->SetLineColor(kBlue+3);
  d5331.fDioMomNorm->SetLineWidth(2);
  d5331.fDioMomNorm->SetMarkerStyle(20);
  d5331.fDioMomNorm->SetMarkerSize(1);
  d5331.fDioMomNorm->SetMarkerColor(kBlue+3);

  d5331.fDioMomNorm->Draw("same,pe");
//-----------------------------------------------------------------------------
// now the signal part
//-----------------------------------------------------------------------------
  double sf_conv = qn_cnv/d5331.fCnvNGen;

  d5331.fCnvMomNorm = (TH1F*) d5331.fCnvMom->Clone("d5331_fCnvMomNorm");

  int nbins = d5331.fCnvMomNorm->GetNbinsX();

  for (int i=0; i<nbins; i++) {
    double cont = d5331.fCnvMom->GetBinContent(i+1)*sf_conv;
    double err  = d5331.fCnvMom->GetBinError  (i+1)*sf_conv;
    d5331.fCnvMomNorm->SetBinContent(i+1,cont);
    d5331.fCnvMomNorm->SetBinError  (i+1,err );
  }

  d5331.fCnvMomNorm->SetStats(0);
  d5331.fCnvMomNorm->SetLineColor(2);
  d5331.fCnvMomNorm->SetLineWidth(2);
  d5331.fCnvMomNorm->Draw("h,sames");

  int ixmin = -1;

  for (int i=0; i<nbins; i++) {
    double p =  d5331.fCnvMomNorm->GetBinLowEdge(i);
    if (p <= PMin) {
      ixmin = i;
    }
    else {
      break;
    }
  }

  printf("  IXMIN = %5i\n",ixmin);

  double bgr = d5331.fDioMomNorm->Integral(ixmin,300);
  printf(" Expected DIO bgr = %12.4e\n",bgr);

  double sig = d5331.fCnvMomNorm->Integral(ixmin,300);
  printf(" Expected Signal = %12.4e\n",sig);

  TLatex* lat = new TLatex(101.2,0.37,Form("N(protons)      : %8.2e",n_prot));
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->Draw();

  TLatex* lat2 = new TLatex(101.2,0.34,Form("N(stopped #mu^{-}) : %8.2e",NStoppedMuons));
  lat2->SetTextFont(42);
  lat2->SetTextSize(0.04);
  lat2->Draw();

  TLatex* lat3 = new TLatex(101.2,0.31,Form("R(#mu #rightarrow e)        : %8.1e",br_mu2e_conv));
  lat3->SetTextFont(42);
  lat3->SetTextSize(0.04);
  lat3->Draw();

  TLatex* lat4 = new TLatex(101.2,0.28,Form("CE signal        : %8.2e",sig));
  lat4->SetTextFont(42);
  lat4->SetTextSize(0.04);
  lat4->Draw();

  TLatex* lat5 = new TLatex(101.2,0.25,Form("DIO                 : %8.2e",bgr));
  lat5->SetTextFont(42);
  lat5->SetTextSize(0.04);
  lat5->Draw();

  TLatex* lat6 = new TLatex(101.2,0.22,Form("SES                : %8.2e",1./sig*1e-16));
  lat6->SetTextFont(42);
  lat6->SetTextSize(0.04);
  lat6->Draw();
}


//-----------------------------------------------------------------------------
void plot_background_vs_threshold() {

  while (1) {
    c_dio_bgr = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("c_dio_bgr");
    if (c_dio_bgr == NULL) break;
    delete c_dio_bgr;
  }

  c_dio_bgr = new TCanvas("c_dio_bgr","DIO background",1200,800);
  c_dio_bgr->cd();

  TH2F* h2 = new TH2F("h2","Integral of the DIO background over [x,105] MeV/c",1,103,105,1,0,1);

  h2->SetStats(0);
  h2->GetYaxis()->SetTitle(Form("N /%5.3f MeV/c",bin_width));
  h2->GetXaxis()->SetTitle(Form("Momentum, MeV/c"));
  h2->Draw();


  d5331.fDioIntegral = (TH1F*) d5331.fDioMom->Clone("d5331_fDioIntegral");
  d5331.fDioIntegral->Reset();

  int nbins = d5331.fDioIntegral->GetNbinsX();

  
  double sum;
  for (int i=0; i<nbins; i++) {
    if (i < 300) {
      sum = d5331.fDioMomNorm->Integral(i+1,300);
    }
    else {
      sum = 0;
    }
    d5331.fDioIntegral->SetBinContent(i+1,sum);
  }
 
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  d5331.fDioIntegral->Draw("same");
  
  TArrow* arr = new TArrow(PMin,0.5,PMin,0.1,0.015);
  arr->SetLineWidth(3);
  arr->SetLineColor(2);
  arr->Draw();
}


//-----------------------------------------------------------------------------
// calculate total expected number of DIO events 
//-----------------------------------------------------------------------------
void dio_integral(double& NConv, double& NDio) {

					// muon (transport+stopping) efficiency, TDR number : 1.9e-3
  NStoppedMuons = n_prot*1.9e-3;
					// track reco eff for DIO's is accounted for lated 

  NDio  = NStoppedMuons*(1-mu_capture_prob);

  NConv = NStoppedMuons*mu_capture_prob*br_mu2e_conv;

  printf("NStoppedMuons   = %12.4e\n",NStoppedMuons);
  printf("NDIO            = %12.3e\n",NDio);
  printf("NC)nv           = %12.3e\n",NConv);


  TFile* f = TFile::Open("/grid/fermiapp/mu2e/users/murat/hist/mu2e/Czarnecki.root");

  TH1D* h_dio = (TH1D*) f->Get("Spectrum");

  double f_101_105 = h_dio->Integral(1011,1050)/h_dio->Integral();

  printf("f_101_105       = %12.3e\n",f_101_105);

  double n_101_105 = NDio*f_101_105;
  
  printf("n_101_105       = %12.3e\n",n_101_105);
}

//-----------------------------------------------------------------------------
void dio(const char* Folder = "trk_1") {

  double  cont, err;
					// N events generated in [100,105]
  d5331.fDioNGen = 980000.; 
  d5331.fCnvNGen = 980000.;

  //  d5331.fCnvMom = (TH1F*) gh1("hist/cnvs5331.track_ana.hist","TrackAna",Form("%s/p"   ,Folder))->Clone("d5331_fCnvMom");
  //  d5331.fDioMom = (TH1F*) gh1("hist/dios5331.track_ana.hist","TrackAna",Form("%s/pdio",Folder))->Clone("d5331_fDioMom");
  d5331.fCnvMom = (TH1F*) gh1("hist/cnvs5331.track_comp.hist","TrackComp",Form("%s/p"   ,Folder))->Clone("d5331_fCnvMom");
  d5331.fDioMom = (TH1F*) gh1("hist/dios5331.track_comp.hist","TrackComp",Form("%s/pdio",Folder))->Clone("d5331_fDioMom");

  //  d5331.fDioMom->Rebin(2);

  bin_width = d5331.fDioMom->GetBinWidth(1);

  printf("bin width: %12.4e\n",bin_width);

  dio_integral(qn_cnv,qn_dio);

  double s2 = 4*qn_dio/d5331.fDioNGen;

  printf("DIO scale factor: %12.4e \n",s2);

  //  dios5331.fDioMom->Draw();

  d5331.fDioMomNorm = (TH1F*) d5331.fDioMom->Clone("d5331_fDioMomNorm");

  int nbins = d5331.fDioMomNorm->GetNbinsX();

  for (int i=0; i<nbins; i++) {
    cont = d5331.fDioMom->GetBinContent(i+1)*s2;
    err  = d5331.fDioMom->GetBinError  (i+1)*s2;
    d5331.fDioMomNorm->SetBinContent(i+1,cont);
    d5331.fDioMomNorm->SetBinError  (i+1,err );
  }

  plot_sensitivity();

  plot_background_vs_threshold();
}


//-----------------------------------------------------------------------------
void test_001() {
  double n = dio_integral();
}

//-----------------------------------------------------------------------------
void test_002() {
  check_dio_norm();
}
