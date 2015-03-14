// 2014-01-22: Andrei's files in /mu2e/data/tdr/beam/g4s3p3/mergedMuonStops
// why the first foil is 

void plot_muon_stops(const char* Filename = "/mu2e/data/tdr/beam/g4s3p3/mergedMuonStops/mustops.1025a_1025a_1316a.14278089.root") {

  TFile* f = TFile::Open(Filename);

  TDirectory* d = f->GetDirectory("stoppedMuonDumper");

  TTree* nt = (TTree*) d->Get("stops"); 

  nt->Draw("x");

  TH1F* h_z = new TH1F("h_z","Z",500,5400,6400);

  nt->Draw("z>>h_z");

  h_z->Draw();
  
}
