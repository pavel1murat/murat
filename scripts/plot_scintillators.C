//-------------------------------------------------------------------------------
// from file://home/murat/library/clfv/PhysRevD.66.096002.kitano_mu2e_rate_nuclei.pdf
// ------------------------------------------------------------------------------

void plot_scintillator_properties() {
  
  struct scint_data_t {
    const char* fName;
    double      fDensity;    // g/cm^3
    double      fLightYield; // photons/MeV
    double      fDecayTime;  // ns
  };


  scint_data_t data = {
    {"BaF2", 4.88, 2000., 0.6},
    {0     , 4.88, 2000., 0.6}
  }

  double x[1000], y[1000];


  TH2F* h2 = new TH2F("h2","h2",1000,0,200000,10000,0,1000);


  h2->Draw();

  TEllipse* el;

  int np(0);

  for (int i=0; data[i].fName != 0; i++) {

    double x = data[i].fLightYield;
    double y = data[i].fDecayTime;
    double r = data[i].fDensity*10;

    el = new TEllipse(x,y,r);
    
    el->Draw();
  }

}
 
