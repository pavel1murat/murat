//

#include "petr/SimPulse.hh"


//-----------------------------------------------------------------------------
SimPulse::SimPulse() {

  fNPhotonsPerMeV = 2000.;		// BaF2;
  fPDE            = 0.2;
					// assume 511KeV = 2000/2*0.2

  n_phot    = fNPhotonsPerMeV*0.511*fPDE;

  printf ("nphotons per MeV = %10.3f\n",fNPhotonsPerMeV);
  printf ("PDE              = %10.3f\n",fPDE);
  printf ("<nphotons>       = %10i\n"  ,n_phot);

  t_decay   = 0.8;

  sigma_res = 0.2;

  f_res     = new TF1("f",SimPulse::resolution,-10,10,3);

  f_res->SetParameter(1,sigma_res);
  f_res->SetParameter(2,1);

  fRn   = new TRandom3();

  for (int i=0; i<1000; i++) {
    fEvent[i].fWaveform = 0;
  }

  fDt = new TH1F("dt","dt",200,-1,1);
}


//-----------------------------------------------------------------------------
SimPulse::~SimPulse() {
}


//-----------------------------------------------------------------------------
// single photon time resolution (SPTR)
//-----------------------------------------------------------------------------
double  SimPulse::resolution(double* x, double* p) {
  double f, sig;

  double dt = x[0]-p[0];

  sig = p[1];

  f = p[2]*TMath::Exp(-(dt*dt)/(2*sig*sig));

  return f;
}

//-----------------------------------------------------------------------------
// unused
//-----------------------------------------------------------------------------
double  SimPulse::pulse(double* x, double* p) {
  double f, t, td, tf;

  t = x[0];
  
  if (t < 0) return 0;
  
  td = p[0];
  tf = p[1];
  
  f = TMath::Exp(-t/td)*(1-TMath::Exp(-t/tf));

  return f;

}

//-----------------------------------------------------------------------------
// for now, use histogram as a waveform storage
//-----------------------------------------------------------------------------
int  SimPulse::SimulateWaveform(Event_t* Event) {

  double     t, xmin, xmax;

  printf(">>> simulating event %3i\n",Event->fNumber);
    
  if (Event->fWaveform) delete Event->fWaveform;

  Event->fWaveform = new TH1F(Form("h_wf_%03i",Event->fNumber),
			      Form("waveform event %03i",Event->fNumber),
			      kNBins,-5,15);

  Event->fT0 = 0.2*fRn->Rndm(); // 0;
					// reset waveform
  for (int j=0; j<kNBins; j++) {
    fCharge[j] = 0;
  }
					// simulate collected charge
  for (int i=0; i<n_phot; i++) {
    t = fRn->Exp(t_decay);
    f_res->SetParameter(0,t);
    
    for (int j=0; j<kNBins; j++) {
      xmin = j*0.2-5;
      xmax = xmin+0.2;
      fCharge[j] += f_res->Integral(xmin,xmax);
    }
  }

  for (int i=0; i<kNBins; i++) {
    //      printf(" bin = %4i charge = %12.5e \n",i,fCharge[i]);
    Event->fWaveform->SetBinContent(i+1,fCharge[i]);
  }

  return 0;
}


//-----------------------------------------------------------------------------
// reconstruction
//-----------------------------------------------------------------------------
int SimPulse::ReconstructWaveform(Event_t* Event) {
  TH1F*  hist;

  double v, vmax, v1, v2, t1; //, t2;

  printf(">>>      reconstructing event %3i\n",Event->fNumber);

  hist = Event->fWaveform;
  double bin = hist->GetBinWidth(1);

  vmax = -1.e6;
  for (int i=0; i<kNBins; i++) {
    v = hist->GetBinContent(i+1);
    if ( v > vmax) {
      vmax = v;
    }
  }
//-----------------------------------------------------------------------------
// consider 2 points 
//-----------------------------------------------------------------------------
  for (int i=0; i<kNBins; i++) {
    
    v1 = hist->GetBinContent(i+1);

    if (v1 > 0.1*vmax) {

      t1 = hist->GetBinCenter(i);

      v2 =  hist->GetBinContent(i+2);
      //      t2 = t1+bin;

      // v = v1+(v2-v1)/bin*(t-t1)

      Event->fT01 = (vmax*0.1-v1)/(v2-v1)*bin+t1;
					// done
      break;
    }
  }

  double dt = Event->fT01-Event->fT0;
  fDt->Fill(dt);

  return 0;
}


//-----------------------------------------------------------------------------
void SimPulse::Run(int Nev) {

  Event_t* ev;

//-----------------------------------------------------------------------------
// digitization with with 200 ps bin
//-----------------------------------------------------------------------------
  for (int iev=0; iev<Nev; iev++) {
    ev = &fEvent[iev];

    ev->fNumber = iev;

    SimulateWaveform(ev);
//-----------------------------------------------------------------------------
// reconstruct waveform
//-----------------------------------------------------------------------------
    ReconstructWaveform(ev);
  } 

  // double qm = h_wf[0]->GetBinContent(32);

  // h_wf[0]->SetLineColor(1);
  // h_wf[0]->GetYaxis()->SetRangeUser(0,qm*1.5);
  // h_wf[0]->Draw();
  
  // h_wf[1]->SetLineColor(2);
  // h_wf[1]->Draw("same");

  // h_wf[2]->SetLineColor(4);
  // h_wf[2]->Draw("same");

}


  
