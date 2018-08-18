//------------------------------------------------------------------------------
// consider attenuation of the light emitted by a charged particle in GaAs sensor
// for given attenuation length (varied from 4mm to 4 cm)
//
// conclusion: light attenuation in the sensor can't be responsible for the
// width of the energy distribution
//-----------------------------------------------------------------------------

#include "TRandom3.h"
#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"


TH1F  *hist(0), *hx1, *hx2;

//-----------------------------------------------------------------------------
void attenuation_a(double AttLength = 0.4, int NEvents =10000) {
  // x coordinate of the alpha-particle varies from -4mm to 4mm

  TRandom3 rn3;

  if (hist) {
    delete hist;
    delete hx1;
    delete hx2;
  }

  TCanvas* c = new TCanvas("c","c",1200,800);
  c->Divide(2,2);
  
  hist = new TH1F("h",
		  Form("measured #alpha-particle energy, #lambda = %4.1f cm",AttLength),
		  500,0.0,5.0);

  hx1 = new TH1F("hx1","dx1",500,0.,5.);
  hx2 = new TH1F("hx2","dx2",500,0.,5.);
  
  for (int i=0; i<NEvents; i++) {
    
    double x = -0.4+0.8*rn3.Rndm();

    double dx1 = 0.4-x;
    double dx2 = 0.4+x+0.8;
    
    double att1 = TMath::Exp(-dx1/AttLength);
    double att2 = TMath::Exp(-dx2/AttLength);

    double e = 4.48*(att1+att2)/2.;

    hist->Fill(e);
    hx1->Fill(dx1);
    hx2->Fill(dx2);
  }

  c->Draw();
  
  c->cd(1);
  hist->SetStats(0);
  hist->GetXaxis()->SetTitle("E, MeV");
  hist->Draw();
  
  c->cd(2);
  hx1->Draw();

  c->cd(3);
  hx2->Draw();


  c->Modified();
  c->Update();
}
