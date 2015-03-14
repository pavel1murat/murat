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


  TCanvas* c = new TCanvas("c","C",600,800);

  c->Draw();


  TH2F* h2 = new TH2F("h2","h2",1000,0,100000,1000,0,100);

  h2->GetXaxis()->SetTitle("Light Yield, Nph/MeV");
  h2->GetYaxis()->SetTitle("Decay Time, ns");

  h2->SetStats(0);
  
  h2->Draw();

  TMarker*  m;

  int np(0), done (0);

  double x, y, r, r1, r2, density, light_yield, decay_time;

  char c[1000], name[100];

  FILE* f  = fopen("murat/scripts/scintillator_properties.dat","r");

  while ( ((c[0]=getc(f)) != EOF) && !done) {
					// skip comment lines
    if (c[0] != '#') {
      ungetc(c[0],f);
					// read scintillator data
      fscanf(f,"%s" ,name        );
      fscanf(f,"%lf" ,&density    );
      fscanf(f,"%lf" ,&light_yield);
      fscanf(f,"%lf",&decay_time );

      printf("%20s  %lf %lf %lf\n",name,density,light_yield,decay_time);

      r2 = density/2;

      el = new TEllipse(light_yield,decay_time,r1,r2);

      m = new TMarker(light_yield,decay_time,20);
      m->SetMarkerSize(r2);
      m->SetMarkerStyle(3001);
      m->SetMarkerColor(42);

      m->Draw();

      t = new TText(light_yield,decay_time,name);
      t->SetTextFont(42);
      t->SetTextSize(0.035);
      t->Draw();
    }
					// skip line
    fgets(c,100,f);
  }

  fclose(f);
}
 
