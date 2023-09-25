// vincenzo : 13.5ns, 40ns

#include "ana/TCaloTimeRes.hh"


namespace {
  double Tau1;
  double Tau2;
}

//-----------------------------------------------------------------------------
TCaloTimeRes::TCaloTimeRes() {

  fTau1     = 20.;			// decay time
  fTau2     = 10.;			// rise time 

  fNSamples = 20;
  
  fHist.fWaveform = new TH1F("wform","waveform",40,-10,190);
  fHist.fRes      = new TH1F("tres" ,"resolution",1000,-2.5,2.5);
  fHist.fPull     = new TH1F("pull" ,"Pull"      ,200,-10,10);

  //  fPulse          = new TF1("f_pulse",TCaloTimeRes::pulse,-10,200,2);
  fPulse          = new TF1("f_pulse",TCaloTimeRes::pulse2,-10,200,2);

  fClusterEnergy  = 50.; // MeV
  fNoise          = 0.3/sqrt(20)/sqrt(2); // MeV, assume averaging

  fNPePerMeV      = 30.; // 30 "pe" per MeV per APD channel
  fTimeBin        = 5;    // 5ns
  fMaxFitTime     = 40.;
}

//-----------------------------------------------------------------------------
TCaloTimeRes::~TCaloTimeRes() {
  delete fPulse;
}


//-----------------------------------------------------------------------------
// parameterization: P[0]*(exp((t-P[1])/Tau1) - exp((t-P[0])/Tau2))
// assume the pulse shape is constant
//-----------------------------------------------------------------------------
double TCaloTimeRes::pulse(double* X, double* P) {

  double x = X[0]-P[1];

  if (x < 0) return 0;

  double f = P[0]*(-TMath::Exp(-x/Tau1)+TMath::Exp(-x/Tau2))/(Tau1-Tau2);

  return f;
}

//-----------------------------------------------------------------------------
// pulse1: integral
// parameterization: P[0]*(exp((t-P[1])/Tau1) - exp((t-P[0])/Tau2))
// assume the pulse shape is constant
// assume histogram is being fit
//-----------------------------------------------------------------------------
double TCaloTimeRes::pulse2(double* X, double* P) {

  double f, r1, r2;

  double t1 = X[0]-P[1]-2.5;
  if (t1 < 0) t1 = 0;

  double t2 = X[0]+2.5-P[1];

  if (t2 < 0) return 0;

  r1 = Tau1*(exp(-t1/Tau1)-exp(-t2/Tau1));
  r2 = Tau2*(exp(-t1/Tau2)-exp(-t2/Tau2));

  f = P[0]*(r1-r2)/(Tau1-Tau2);

  return f;
}

//-----------------------------------------------------------------------------
// calculate integral of the pulse shape normalized to 1, over the bin
//-----------------------------------------------------------------------------
double TCaloTimeRes::func_integral(double t0, double tmin, double tmax) {

  double t1, t2, r1, r2, result;

  t1 = tmin-t0;
  if (t1 < 0) t1 = 0;

  t2 = tmax-t0;

  r1 = Tau1*(exp(-t1/Tau1)-exp(-t2/Tau1));
  r2 = Tau2*(exp(-t1/Tau2)-exp(-t2/Tau2));

  result = (r1-r2)/(Tau1-Tau2);

  return result;
}

//-----------------------------------------------------------------------------
void TCaloTimeRes::time_resolution(int NEvents) {
  
  // simulate the pulse

  double   q[1000], qe, qnpe, err;
  double   e, t0, noise, tmin, tmax, bin, sum_q, tfit, efit, dt;

  fNEvents    = NEvents;
  Tau1        = fTau1;
  Tau2        = fTau2;
  qnpe        = fNPePerMeV*fClusterEnergy;
  fSigmaNoise = fNPePerMeV*fNoise;

  fHist.fRes->Reset();
  fHist.fPull->Reset();

  for (int ievent=0; ievent<NEvents; ievent++) {

				// pulse "energy", in units of "photoelectrons"
				// pulse starts in the first bin 
    e  = fRn.Poisson(qnpe); 
    t0 = fRn.Rndm(0)*fTimeBin;

    // now need to integrate the charges in each bin and reconstruct the waveform

    sum_q = 0;

    fHist.fWaveform->Reset();

    for (int it=0; it<20; it++) {
      tmin = fTimeBin*it;
      tmax = tmin+fTimeBin;
//-----------------------------------------------------------------------------
// simulate the bin charges. The integral is normalized to ~1 (the tail 
// is not accounted for)
//-----------------------------------------------------------------------------
      qe = func_integral(t0,tmin,tmax)*e;
//-----------------------------------------------------------------------------
// for the fast scintillator, do not need to account for statistical fluctuations 
// of the number of photons - those define the amplitude only
//-----------------------------------------------------------------------------
      //      q[it] = fRn.Poisson(qe);
      q[it] = qe;

				// add the electronic noise

      noise  = fRn.Gaus(0,fSigmaNoise);
      q[it]  += noise;
      sum_q  += q[it];

      bin = it+3;
      fHist.fWaveform->SetBinContent(bin,q[it]);
      //      err = sqrt(fSigmaNoise*fSigmaNoise+qe);
      err = fSigmaNoise;
      fHist.fWaveform->SetBinError  (bin,err);
    }
//-----------------------------------------------------------------------------
// the 'sampled waveform' is simulated, plot it first
// for fitting purpose - assume T0 at center of the first bin 
//-----------------------------------------------------------------------------
    fPulse->SetParameter(0,sum_q);
    fPulse->SetParameter(1,2.5);

    fHist.fWaveform->Fit(fPulse,"Q","",0,fMaxFitTime);

    tfit = fPulse->GetParameter(1);
    efit = fPulse->GetParError(1);

    //    printf(" t0, tfit, efit = %10.3f %10.3f %10.3f \n",t0, tfit, efit);
    
    dt = tfit-t0;
    fHist.fRes->Fill(dt);
    fHist.fPull->Fill(dt/efit);
  }

  fHist.fRes->Fit("gaus");
}
