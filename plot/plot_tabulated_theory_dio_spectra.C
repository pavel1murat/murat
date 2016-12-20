///////////////////////////////////////////////////////////////////////////////
// read tabulated DIO spectra from difeerent theoretical papers
// and compare them
//-----------------------------------------------------------------------------
#include "Stntuple/val/stntuple_val_functions.hh"

TGraph   *gr_w, *gr_h;
TSpline3 *s_w , *s_h;
TH1F*     h1, *h2;

//-----------------------------------------------------------------------------
void plot_czarnecki_vs_watanabe_al27 (const char* Opt) {
  TNtuple* nt_h, *nt_w;

  create_ntuple("ConditionsService/data/czarnecki_Al.tbl",nt_h);
  create_ntuple("murat/data/dio_spectrum_al27_watanabe.txt",nt_w);

  float energy, weight;

  nt_w->SetBranchAddress("var0",&energy);
  nt_w->SetBranchAddress("var1",&weight);

  int nw = nt_w->GetEntries();

  float x[10000], y[10000];
  
  nt_h->SetBranchAddress("var0",&energy);
  nt_h->SetBranchAddress("var1",&weight);

  int nh = nt_h->GetEntries();

  for (int i=0; i<nh; i++) {
    nt_h->GetEntry(i);

    x[i] = energy;
    y[i] = weight;
  }

  gr_h = new TGraph(nh,x,y);

  gr_h->SetLineColor(kBlue+2);
  gr_h->SetLineWidth(1);

  gr_h->GetYaxis()->SetRangeUser(1.e-20,0.05);
  gr_h->Draw("alp");
 
  for (int i=0; i<nw; i++) {
    nt_w->GetEntry(i);

    x[i] = energy;
    y[i] = weight;
  }
  
  gr_w = new TGraph(nw,x,y);
  gr_w->SetLineColor(2);
  gr_w->SetLineWidth(1);

  gr_w->Draw("lp,same");
  
}

//-----------------------------------------------------------------------------
void plot_watanabe_vs_hertzog_ca40 (const char* Opt) {
  TNtuple* nt_h, *nt_w;

  create_ntuple("murat/data/dio_spectrum_ca40_watanabe.txt",nt_w);
  create_ntuple("murat/data/dio_spectrum_ca40_hertzog_alder.txt",nt_h);

  float energy, weight;

  nt_w->SetBranchAddress("var0",&energy);
  nt_w->SetBranchAddress("var1",&weight);

  int nw = nt_w->GetEntries();

  float x[100], y[100];
  
  for (int i=0; i<nw; i++) {
    nt_w->GetEntry(i);

    x[i] = energy;
    y[i] = weight*100;
  }
  
  gr_w = new TGraph(nw,x,y);
  //  float q = gr_w->Integral();

  // for (int i=0; i<nw; i++) {
  //   y[i] = y[i]/q/1.00258;
  // }

  // delete gr_w; gr_w = new TGraph(nw,x,y);

  s_w = new TSpline3("spline_w",gr_w);

  gr_w->SetLineColor(2);
  gr_w->SetLineWidth(2);

  if ((*Opt == 'a') || (*Opt == 'w')) {
    gr_w->Draw("alp");
    s_w->Draw("same");
  }
  
  nt_h->SetBranchAddress("var0",&energy);
  nt_h->SetBranchAddress("var1",&weight);

  int nh = nt_h->GetEntries();

  for (int i=0; i<nh; i++) {
    nt_h->GetEntry(i);

    x[i] = energy*0.511;
    y[i] = weight;
  }

  gr_h = new TGraph(nh,x,y);

  // q = gr_h->Integral();

  // for (int i=0; i<nh; i++) {
  //   y[i] = y[i]/q/1.07814;
  // }

  // delete gr_h; gr_h = new TGraph(nh,x,y);
  
  s_h = new TSpline3("spline_h",gr_h);

  
  gr_h->SetLineColor(2);
  gr_h->SetLineWidth(2);

  if      (*Opt == 'a') {
    gr_h->Draw("lp,same");
    s_h->Draw("same");
  }
  else if (*Opt == 'h') {
    gr_h->Draw("alp");
    s_h->Draw("same");
  }


  h1 = new TH1F("h1","h1",1000,0,100);
  h2 = new TH1F("h2","h2",1000,0,100);

  for (int i=1; i<=1000; i++) {
    float x1 = (i-0.5)*0.1;
    float yw = s_w->Eval(x1);
    float yh = s_h->Eval(x1);

    h1->SetBinContent(i,yw);
    h2->SetBinContent(i,yh);
  }
  
}

