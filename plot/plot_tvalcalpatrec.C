////////////////////////////////////////////////////////////////////////////////////
// List of plot:
//   1. 2-dim plot for showing TrkPatRec efficiency Vs calPatRec efficiency
//   2. overlapped momentum distribution with no quality cuts
//   3. overlapped momentum distribution applying set C cuts
//   4. ...
///////////////////////////////////////////////////////////////////////////////////


plot_box_eff() {

  TCanvas* c1 = new TCanvas("c1","c1",900,600);
  
  TString hist_dir="/mu2e/app/users/murat/2014-03-05/hist/cnvs0302.val_cpr.hist.500";
  TString hist_name;
  hist_name = hist_dir;
  hist_name += ;
  
  TH2F* h_ef  = gh2(hist_dir.Data(),"ValCalPatRec","evt_0/nt2_vs_nt1");
     
  
  //  h_ef->GetXaxis()->SetRangeUser(-15,14.99);
  h_ef->GetXaxis()->SetTitle("TrkPatRec");
  h_ef->GetYaxis()->SetTitle("CalPatRec");
  h_ef->SetTitle("All tracks");
  h_ef->Draw("box");
  //  c1->Draw();cnvs0302_reco_eff_all
}

plot_box_eff_good() {

  TCanvas* c1 = new TCanvas("c1","c1",900,600);
  
  TString hist_dir="/mu2e/app/users/murat/2014-03-05/hist/cnvs0302.val_cpr.hist.500  ";
  TString hist_name;
  hist_name = hist_dir;
  hist_name += ;
  
  TH2F* h_ef = gh2(hist_dir.Data(),"ValCalPatRec","evt_0/ngt2_vs_ngt1");
    
  
  //  h_ef->GetXaxis()->SetRangeUser(-15,14.99);
  h_ef->GetXaxis()->SetTitle("TrkPatRec");
  h_ef->GetYaxis()->SetTitle("CalPatRec");
  h_ef->SetTitle("Set C tracks");
  h_ef->Draw("box");
  //  c1->Draw();cnvs0302_reco_eff_setC
}

//-----------------------------------------------------------------------------
void plot_box_nexp_vs_nhits(const char* Fol = "trk1_0") {

  TCanvas* c1 = new TCanvas("c1","c1",900,900);
  
  TString hist_dir="/mu2e/app/users/murat/2014-03-05/hist/cnvs0302.val_cpr.hist";
  TString hist_name;
  hist_name = hist_dir;
  hist_name += ;
  
  TH2F* h1 = gh2(hist_dir.Data(),"ValCalPatRec",Form("%s/nep_vs_nhp",Fol));
    
  
  //  h_ef->GetXaxis()->SetRangeUser(-15,14.99);

  h1->GetXaxis()->SetRangeUser(0,59.9);
  h1->GetYaxis()->SetRangeUser(0,59.9);

  h1->SetTitle(Form("%s: N(expected) vs N(reconstructed)",Fol));
  h1->Draw("box");
}


plot_p_all() {

  TCanvas* c = new TCanvas("c","c",900,600);
  
  TString hist_dir="/mu2e/app/users/murat/2014-03-05/hist/cnvs0302.val_cpr.hist.500";
  TString hist_name;
  hist_name = hist_dir;
  hist_name += ;
  
  TH1F* h_pc  = gh1(hist_dir.Data(),"ValCalPatRec","trk2_1/p");
  TH1F* h_pt  = gh1(hist_dir.Data(),"ValCalPatRec","trk1_0/p");
  
  //  h_ef->GetXaxis()->SetRangeUser(-15,14.99);
  h_pt->GetXaxis()->SetTitle("P [MeV/c]");
  h_pt->GetXaxis()->SetRangeUser(90.,110.);
  double binsize=h_pt->GetBinWidth(1);
  h_pt->GetYaxis()->SetTitle(Form("Entries [#]/( %3.2f MeV/c)", binsize));
  h_pt->SetName("TrkPatRec");
  h_pt->SetTitle("Momentum distribution, all tracks");
  
  h_pc->SetName("CalPatRec");
  h_pc->SetFillStyle(3002);
  h_pc->SetFillColor(kRed);
  h_pc->SetLineColor(kRed);

  
  h_pt->Draw();
  h_pc->Draw("sames");//cnvs0302_p_all
  TLegend *leg = new TLegend(0.2, 0.72, 0.39, 0.85, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(h_pt,"TrkPatRec","L");
  leg->AddEntry(h_pc,"CalPatRec","L");
  leg->Draw();
}

plot_p_good() {

  TCanvas* c = new TCanvas("c","c",900,600);
  
  TString hist_dir="/mu2e/app/users/murat/2014-03-05/hist/cnvs0302.val_cpr.hist.500";
  TString hist_name;
  hist_name = hist_dir;
  hist_name += ;
  
  TH1F* h_pt  = gh1(hist_dir.Data(),"ValCalPatRec","trk1_4/p");
  TH1F* h_pc  = gh1(hist_dir.Data(),"ValCalPatRec","trk2_5/p");
  
  //  h_ef->GetXaxis()->SetRangeUser(-15,14.99);
  h_pt->GetXaxis()->SetTitle("P [MeV/c]");
  h_pt->GetXaxis()->SetRangeUser(90.,110.);
  double binsize=h_pt->GetBinWidth(1);
  h_pt->GetYaxis()->SetTitle(Form("Entries [#]/( %3.2f MeV/c)", binsize));
  h_pt->SetName("TrkPatRec");
  h_pt->SetTitle("Momentum distribution with cut-set C, all tracks");
  
  h_pc->SetName("CalPatRec");
  h_pc->SetFillStyle(3002);
  h_pc->SetFillColor(kRed);
  h_pc->SetLineColor(kRed);

  h_pt->Draw();
  h_pc->Draw("sames");
 TLegend *leg = new TLegend(0.2, 0.72, 0.39, 0.85, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(h_pt,"TrkPatRec","L");
  leg->AddEntry(h_pc,"CalPatRec","L");
  leg->Draw();
}






plot_chi_all() {

  TCanvas* c = new TCanvas("c","c",900,600);
  
  TString hist_dir="/mu2e/app/users/murat/2014-03-05/hist/cnvs0302.val_cpr.hist.500";
  TString hist_name;
  hist_name = hist_dir;
  hist_name += ;
  
  TH1F* h_pt  = gh1(hist_dir.Data(),"ValCalPatRec","trk1_0/chi2d");
  TH1F* h_pc  = gh1(hist_dir.Data(),"ValCalPatRec","trk2_1/chi2d");
  
  //  h_ef->GetXaxis()->SetRangeUser(-15,14.99);
  h_pt->GetXaxis()->SetTitle("#chi^{2}/ndof");
  //  h_pt->GetXaxis()->SetRangeUser(90.,110.);
  //  double binsize=h_pt->GetBinWidth(1);
  h_pt->GetYaxis()->SetTitle("Entries [#]");
  h_pt->SetName("TrkPatRec");
  h_pt->SetTitle("#chi^{2}/ndof distribution, all tracks");
  
  h_pc->SetName("CalPatRec");
  h_pc->SetFillStyle(3002);
  h_pc->SetFillColor(kRed);
  h_pc->SetLineColor(kRed);

  h_pt->Draw();
  h_pc->Draw("sames");

TLegend *leg = new TLegend(0.37, 0.72, 0.56, 0.85, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(h_pt,"TrkPatRec","L");
  leg->AddEntry(h_pc,"CalPatRec","L");
  leg->Draw();
}

plot_chi_good() {

  TCanvas* c = new TCanvas("c","c",900,600);
  
  TString hist_dir="/mu2e/app/users/murat/2014-03-05/hist/cnvs0302.val_cpr.hist.500";
  TString hist_name;
  hist_name = hist_dir;
  hist_name += ;
  
  TH1F* h_pt  = gh1(hist_dir.Data(),"ValCalPatRec","trk1_4/chi2d");
  TH1F* h_pc  = gh1(hist_dir.Data(),"ValCalPatRec","trk2_5/chi2d");
  
  //  h_ef->GetXaxis()->SetRangeUser(-15,14.99);
  h_pt->GetXaxis()->SetTitle("#chi^{2}/ndof");
  //  h_pt->GetXaxis()->SetRangeUser(90.,110.);
  double binsize=h_pt->GetBinWidth(1);
  h_pt->GetYaxis()->SetTitle("Entries [#]");
  h_pt->SetName("TrkPatRec");
  h_pt->SetTitle("#chi^{2}/ndof distribution with cut-set C, all tracks");
  
  h_pc->SetName("CalPatRec");
  h_pc->SetFillStyle(3002);
  h_pc->SetFillColor(kRed);
  h_pc->SetLineColor(kRed);

  h_pt->Draw();
  h_pc->Draw("sames");
  
TLegend *leg = new TLegend(0.37, 0.72, 0.56, 0.85, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(h_pt,"TrkPatRec","L");
  leg->AddEntry(h_pc,"CalPatRec","L");
  leg->Draw();
c->Draw();
}

plot_nactv_all() {

  TCanvas* c = new TCanvas("c","c",900,600);
  
  TString hist_dir="/mu2e/app/users/murat/2014-03-05/hist/cnvs0302.val_cpr.hist.500";
  TString hist_name;
  hist_name = hist_dir;
  hist_name += ;
  
  TH1F* h_pt  = gh1(hist_dir.Data(),"ValCalPatRec","trk1_0/nactv");
  TH1F* h_pc  = gh1(hist_dir.Data(),"ValCalPatRec","trk2_1/nactv");
  
  //  h_ef->GetXaxis()->SetRangeUser(-15,14.99);
  h_pt->GetXaxis()->SetTitle("nactvivated");
  h_pt->GetXaxis()->SetRangeUser(0.,100.);
  //  double binsize=h_pt->GetBinWidth(1);
  h_pt->GetYaxis()->SetTitle("Entries [#]");
  h_pt->SetName("TrkPatRec");
  h_pt->SetTitle("nactvivated distribution, all tracks");
  
  h_pc->SetName("CalPatRec");
  h_pc->SetFillStyle(3002);
  h_pc->SetFillColor(kRed);
  h_pc->SetLineColor(kRed);

  h_pt->Draw();
  h_pc->Draw("sames");
  
TLegend *leg = new TLegend(0.2, 0.72, 0.39, 0.85, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(h_pt,"TrkPatRec","L");
  leg->AddEntry(h_pc,"CalPatRec","L");
  leg->Draw();
c->Draw();
}

plot_nactv_good() {

  TCanvas* c = new TCanvas("c","c",900,600);
  
  TString hist_dir="/mu2e/app/users/murat/2014-03-05/hist/cnvs0302.val_cpr.hist.500";
  TString hist_name;
  hist_name = hist_dir;
  hist_name += ;
  
  TH1F* h_pt  = gh1(hist_dir.Data(),"ValCalPatRec","trk1_4/nactv");
  TH1F* h_pc  = gh1(hist_dir.Data(),"ValCalPatRec","trk2_5/nactv");
  
  //  h_ef->GetXaxis()->SetRangeUser(-15,14.99);
  h_pt->GetXaxis()->SetTitle("nactvivated");
  h_pt->GetXaxis()->SetRangeUser(0.,100.);
  double binsize=h_pt->GetBinWidth(1);
  h_pt->GetYaxis()->SetTitle("Entries [#]");
  h_pt->SetName("TrkPatRec");
  h_pt->SetTitle("nactvivated distribution with cut-set C, all tracks");
  
  h_pc->SetName("CalPatRec");
  h_pc->SetFillStyle(3002);
  h_pc->SetFillColor(kRed);
  h_pc->SetLineColor(kRed);

  h_pt->Draw();
  h_pc->Draw("sames");
  
TLegend *leg = new TLegend(0.2, 0.72, 0.39, 0.85, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);
  leg->AddEntry(h_pt,"TrkPatRec","L");
  leg->AddEntry(h_pc,"CalPatRec","L");
  leg->Draw();
c->Draw();
}

