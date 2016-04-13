//

//-----------------------------------------------------------------------------
void plot_dio_spectrum(const char* Fn, const char* Hist = "trk_1/pdio") {

  h_norm = gh1(Fn,"TrackAna","evt_0/ce_costh");

					// account for [-0.6, 0.8] region

  double q_events = h_norm->Integral()*20./14.;

  TH1D* h_dio  = (TH1D*) gh1(Fn,"TrackAna",Hist)->Clone("h_dio");
//-----------------------------------------------------------------------------
// now - normalization:
//-----------------------------------------------------------------------------
  double mu_capture_prob(0.609);
					// roughly speaking, fraction of captures within the 
					// active window [700,1700] ns  (0.53 = 2594./4889.)
  double live_fraction(0.53);
					// total number of 8 GeV protons on target in 3 years
  double n_prot      (1.2e20 * 3.);
					// muon (transport+stopping) efficiency : 1.6e-3

  double n_stopped_muons = n_prot*1.6e-3;

					// track reco eff for DIO's is accounted for later

  double n_dio = n_stopped_muons*(1-mu_capture_prob)*live_fraction;

  printf("n_stopped_muons = %12.4e\n",n_stopped_muons);
  printf("n_dio           = %12.4e\n",n_dio);

  double sf = n_dio/q_events;
//-----------------------------------------------------------------------------
// at this point h_dio is just the track reco efficiency
//-----------------------------------------------------------------------------
  h_dio->Scale(sf);
//-----------------------------------------------------------------------------
// next: calculate expected number of the DIO events 
//-----------------------------------------------------------------------------
  double xmin = 103.2;
  double xmax = 106.0;

  TH1D* h2 = (TH1D*) h_dio->Clone("h2");
  h2->Reset();

  int n = 10000;
  double step  = (xmax-xmin)/n;
  
  double qint = 0;
  for (int i=0; i<n; i++) {
    double x = xmin+(i+0.5)*step;
    qint += step*TStntuple::DioWeightAl(x); 
  }

  printf(" QINT[%8.3f,%8.f] = %12.5e\n",xmin,xmax,qint);

  h_dio->Draw();
//-----------------------------------------------------------------------------
// calculate integral from 
//-----------------------------------------------------------------------------
  int    min_bin(-1);

  int    nbx  = h_dio->GetNbinsX();
  double x1, xmin1(103.2), xmax1(106.0);

  for (int i=1; i<=nbx; i++) {
    x1 = h_dio->GetBinLowEdge(i);
    if (x1 >= xmin1) {
      min_bin = i;
      break;
    }
  }

  double dio_integral = h_dio->Integral(min_bin,nbx);
  printf("dio_integral = %12.4e\n", dio_integral);
}
