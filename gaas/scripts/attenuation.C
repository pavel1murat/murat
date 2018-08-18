//

#include "TRandom3.h"
#include "TMath.h"
#include "TH1.h"

//------------------------------------------------------------------------------
// consider attenuation of the light emitted by a charged particle in GaAs sensor
// for given attenuation length (varied from 4mm to 4 cm)
//
// conclusion: light attenuation in the sensor can't be responsible for the
// width of the energy distribution
//-----------------------------------------------------------------------------

void attenuation(double AttLength = 0.4, int NEvents =10000) {
  // x coordinate of the alpha-particle varies from -4mm to 4mm

  TRandom3 rn3;


  TH1F* hist = new TH1F("h",
			Form("measured #alpha-particle energy, #lambda = %4.1f cm",AttLength),
			500,0.0,5.0);
  
  for (int i=0; i<NEvents; i++) {
    double x = -0.4+0.8*rn3.Rndm();

    double dx = 0.4-x;
    double att = TMath::Exp(-dx/AttLength);

    double e = 4.48*att;

    hist->Fill(e);
  }

  hist->SetStats(0);
  hist->GetXaxis()->SetTitle("E, MeV");
  hist->Draw();
}
