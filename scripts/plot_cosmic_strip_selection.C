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
// plot cosmic strip selection cuts
//-----------------------------------------------------------------------------
void plot_01(const char* FnCosmics, const char* FnConversion) {
  TH1  *h1, *h2, *h3, *h4;
  TArrow* arr;
  TCanvas* c = new TCanvas("c_plot_01","c_plot_01",0,0,1100,500);
  c->Divide(2,1);

//-----------------------------------------------------------------------------
// D0 selection
//-----------------------------------------------------------------------------
  c->cd(1);
  h1 = get_mu2e_hist(FnConversion,"//TCalm002/Hist/trk_0/d0");
  h1->SetFillColor(603);
  h1->SetFillStyle(3001);
  h1->SetTitle("D0 selection");
  h1->GetXaxis()->SetTitle("D0, mm");
  h2 = get_mu2e_hist(FnCosmics   ,"//TCalm002/Hist/trk_0/d0");

  h1->DrawNormalized();
  h2->DrawNormalized("same");

  arr = new TArrow(-200,0.04,-200,0.025,0.015);
  arr->Draw();
  arr = new TArrow( 200,0.04, 200,0.025,0.015);
  arr->Draw();

//-----------------------------------------------------------------------------
// Z0 selection
//-----------------------------------------------------------------------------
  c->cd(2);

  h3 = get_mu2e_hist(FnConversion,"//TCalm002/Hist/trk_0/z0");

  h3->GetXaxis()->SetRangeUser(-2000,2000);
  h3->SetFillColor(603);
  h3->SetFillStyle(3001);
  h3->SetTitle("Z0 selection");
  h3->GetXaxis()->SetTitle("Z0, mm");

  h4 = get_mu2e_hist(FnCosmics   ,"//TCalm002/Hist/trk_0/z0");

  h3->DrawNormalized();
  h4->DrawNormalized("same");

  arr = new TArrow(-1000,0.1,-1000,0.05,0.015);
  arr->Draw();
  arr = new TArrow( 1000,0.1, 1000,0.05,0.015);
  arr->Draw();


}

