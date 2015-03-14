// CalPatRec
//-----------------------------------------------------------------------------
void print_canvas_with_date(TCanvas* C, const char* Name) {
  TDatime d;
  C->Print(Form("%i-%02i-%02i-%s.eps",d.GetYear(),d.GetMonth(),d.GetDay(),Name));
}


void plot_calpatrec_timing(const char* Filename, int PrintFlag) {

  TCanvas* c1 = new TCanvas("c_kalrep","KalRep::update timing",0,0,1000,700);

  c1->cd();

  int n1(3), n2(11);

  double dat_kalrep_update_cpr[] = {
    52.38, 4.26, 4.06
  };

  double dat_kalrep_update_tpr[] = {
    11.87, 4.11, 3.94
  };

  double dat_kalbend_update_cpr[] = {
    39.23, 2.01, 2.00, 2.00, 1.89, 1.79, 0.67, 0.62, 0.52, 0.51, 0.50
  };

  double dat_kalbend_update_tpr[] = {
    8.99,  0.46, 0.46, 0.42, 0.35, 0.34, 0.18, 0.16, 0.13, 0.11, 0.11
  };

  char* kalrep_call[] = {
    "KalBend::update",
    "KalMaterial::update",
    "KalHit::update"
  };

  char* kalbend_call[] = {
    "BFieldIntegrator::deltaMomentum",
    "std::_Vector_base",
    "CLHEP::HepSymMatrix::similarity",
    "HelixTraj::derivDeflect",
    "HelixTraj::derivDisplace",
    "KalSite::setTraj",
    "CLHEP::operator+",
    "HelixTraj::direction",
    "CLHEP::operator*",
    "KalSite::reset",
    "KalParams::KalParams"
  };

  TH1F* h1_kpr = new TH1F("h1_kpr","KalRep::update timing (CalPatRec)",n1,0,n1);
  TH1F* h1_tpr = new TH1F("h1_tpr","KalRep::update timing (TrkPatRec)",n1,0,n1);

  for (int i=1; i<=n1; i++) {
    h1_kpr->GetXaxis()->SetBinLabel(i,kalrep_call[i-1]);
    h1_kpr->SetBinContent(i,dat_kalrep_update_cpr[i-1]);

    h1_tpr->GetXaxis()->SetBinLabel(i,kalrep_call[i-1]);
    h1_tpr->SetBinContent(i,dat_kalrep_update_tpr[i-1]);
  }

  TH1F* h2_kpr = new TH1F("h2_kpr","KalBend::update timing (CalPatRec)",n2,0,n2);
  TH1F* h2_tpr = new TH1F("h2_tpr","KalBend::update timing (TrkPatRec)",n2,0,n2);

  for (int i=1; i<=n2; i++) {
    h2_kpr->GetXaxis()->SetBinLabel(i,kalbend_call[i-1]);
    h2_kpr->SetBinContent(i,dat_kalbend_update_cpr[i-1]);

    h2_tpr->GetXaxis()->SetBinLabel(i,kalbend_call[i-1]);
    h2_tpr->SetBinContent(i,dat_kalbend_update_tpr[i-1]);
  }

  h1_kpr->SetStats(0);
  h1_kpr->SetMinimum(0.);
  h1_kpr->Draw();
  //  h1_kpr->Draw("same,text45");

  h1_tpr->SetFillColor(kRed-3);
  h1_tpr->SetFillStyle(3001);
  h1_tpr->Draw("same");
  //  h1_tpr->Draw("same,text45");

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  //  leg->AddEntry(h_mpr,"TrkPatRec+CalPatRec","f");
  leg->AddEntry(h1_kpr,"CalPatRec","f");
  leg->AddEntry(h1_tpr,"TrkPatRec","f");

  leg->Draw();
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  TCanvas* c2 = new TCanvas("c_kalbend","KalBend::update timing",0,0,1000,700);

  c2->cd();

  h2_kpr->SetStats(0);
  h2_kpr->SetMinimum(0.);
  h2_kpr->DrawNormalized();
  //  h1_kpr->Draw("same,text45");

  h2_tpr->SetFillColor(kRed-3);
  h2_tpr->SetFillStyle(3001);
  h2_tpr->DrawNormalized("same");
  //  h1_tpr->Draw("same,text45");

  leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(h1_kpr,"CalPatRec","f");
  leg->AddEntry(h1_tpr,"TrkPatRec","f");

  leg->Draw();

  if (PrintFlag) {
    print_canvas_with_date(c1,"kalrep_timing");
    print_canvas_with_date(c2,"kalbend_timing");
  }
}




// CalPatRec
//-----------------------------------------------------------------------------
void plot_mu2e_mergepatrec_timing() {

  TCanvas* c1 = new TCanvas("c_mpr","Mu2e MergePatRec timing",0,0,1000,700);

  c1->cd();

  int n1(26);

  double dat_mu2e_time[] = {
    0.002507 , 0.118637 , 0.134234 , 0.167593 , 0.302655 ,
    0.146095 , 0.298217 , 1.467009 , 0.047263 , 0.037459 ,
    1.639404 , 0.014748 , 2.318028 , 0.033918 , 0.435585 ,
    0.000674 , 0.033595 , 0.001077 , 0.167130 , 0.335949 ,
    0.003421 , 0.229353 , 0.000258 , 0.011618 , 0.000194 ,
    0.001814 
  };

  char* mu2e_call[] = {
    "generate",
    "g4run",
    "dioMixer",
    "protonMixer",
    "neutronMixer",
    "photonMixer",
    "ootMixer",
    "flashMixer",
    "protonTimeMap",
    "muonTimeMap",
    "makeSD",
    "makeSH",
    "CaloReadoutHitsMaker",
    "CaloCrystalHitsMaker",
    "MakeCaloCluster",
    "FSHPreStereo",
    "MakeStereoHits",
    "FlagStrawHits",
    "FlagBkgHits",
    "TrkPatRec",
    "MakeStrawHitPositions",
    "CalPatRec",
    "MergePatRec",
    "TrkExtrapol",
    "CaloMatching",
    "ParticleID"
  };

  TH1F* h1_mu2e = new TH1F("h1_mu2e","Mu2e timing (MergePatRec)",n1,0,n1);

  for (int i=1; i<=n1; i++) {
    h1_mu2e->GetXaxis()->SetBinLabel(i,mu2e_call[i-1]);
    h1_mu2e->SetBinContent(i,dat_mu2e_time[i-1]);

  }

  h1_mu2e->SetStats(0);
  h1_mu2e->SetMinimum(0.);
  h1_mu2e->Draw();
  //  h1_kpr->Draw("same,text45");

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  //  leg->AddEntry(h_mpr,"TrkPatRec+CalPatRec","f");
  leg->AddEntry(h1_mu2e,"Mu2e executable","f");
  //  leg->AddEntry(h1_tpr,"TrkPatRec","f");

  leg->Draw();

  print_canvas_with_date(c1,"mu2e_timing");
}
