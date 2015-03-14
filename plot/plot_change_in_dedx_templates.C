///////////////////////////////////////////////////////////////////////////////
// assume Vadim's file format - 10 slices in path length, from 0 to 10 mm
// Particle: "e" or "m"
//
// old histogram files: 
// --------------------
// ParticleID/test/electrontemplates.root and ParticleID/test/muontemplates.hist
// 
// new histogram files: 
// --------------------
// v4_2_7/hist/electrontemplates_v427.hist and v4_2_7/hist/muontemplates_v427.hist
// 
// new templates have been produced with murat/plot/plot_ehit_vs_path_slices.C ,
// see make_ehit_vs_path_templates routine there, and using as input 2D histograms 
// from 
//
// v4_2_7/hist/trackRecoCheck_ele.hist and v4_2_7/hist/trackRecoCheck_muo.hist
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
void print_canvas_with_date(TCanvas* C, const char* Name) {
  TDatime d;
  C->Print(Form("%i-%02i-%02i-%s.eps",d.GetYear(),d.GetMonth(),d.GetDay(),Name));
}


void plot_change_in_dedx_templates(const char* Fn1, const char* Fn2, const char* Particle) {

  TCanvas* c = new TCanvas("c","c",800,1000);

  c->Divide(2,3);

  TFile* f1 = TFile::Open(Fn1);
  TFile* f2 = TFile::Open(Fn2);

  char hname[100];

  for (int i=1;i<7; i++) {
    c->cd(i);
  
    sprintf(hname,"htemp%s%i",Particle,i);

    TH1F* h1 = (TH1F*) f1->Get(hname);
    TH1F* h2 = (TH1F*) f2->Get(hname);

    h2->DrawNormalized("",3./5.);

    h1->SetFillColor(623);
    h1->SetFillStyle(3001);
    h1->SetLineColor(1);
    h1->DrawNormalized("same");

    printf("bin1, bin2: %10.5f %10.5f\n",
	   h1->GetXaxis()->GetBinWidth(1),
	   h2->GetXaxis()->GetBinWidth(1));
  }

  print_canvas_with_date(c,Form("%s_path_slices_01",Particle));
 

  TCanvas* c2 = new TCanvas("c2","c2",800,1000);

  c2->Divide(2,3);

  TFile* f1 = TFile::Open(Fn1);
  TFile* f2 = TFile::Open(Fn2);

  char hname[100];

  for (int i=7;i<11; i++) {
    c2->cd(i-6);
  
    sprintf(hname,"htemp%s%i",Particle,i);

    TH1F* h1 = (TH1F*) f1->Get(hname);
    TH1F* h2 = (TH1F*) f2->Get(hname);

    h2->DrawNormalized("",3./5.);
    h1->SetFillColor(623);
    h1->SetFillStyle(3001);
    h1->SetLineColor(1);
    h1->DrawNormalized("same");
  }

  print_canvas_with_date(c2,Form("%s_path_slices_02",Particle));
}
