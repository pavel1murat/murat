// example of use of TGenPhaseSpace
//Author: Valerio Filippini

{
  if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");

  TLorentzVector target(0.0, 0.0, 0.0, 0.1056);
  TLorentzVector beam  (0.0, 0.0, 0.,   .0);
  TLorentzVector W = beam + target;

  //(Momentum, Energy units are Gev/C, GeV)
  Double_t masses[3] = { 0., 0., 0.000511} ;

  TGenPhaseSpace event;
  event.SetDecay(W, 3, masses);

  TH1F *h2 = new TH1F("h1","h1", 100,0,0.1);

  for (Int_t n=0;n<100000;n++) {
    Double_t weight = event.Generate();
    
    TLorentzVector *nu1  = event.GetDecay(0);
    
    TLorentzVector *nu2  = event.GetDecay(1);
    TLorentzVector *ele  = event.GetDecay(2);
    
    h1->Fill(ele.Energy(), weight);
  }
  
  h1->Draw();
}
