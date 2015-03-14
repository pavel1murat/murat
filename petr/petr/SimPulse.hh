#ifndef pet_SimPulse_hh
#define pet_SimPulse_hh

#include "TMath.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TF1.h"

class SimPulse {
public:

  enum { kNBins = 100};

  struct Event_t {
    int    fNumber;
    TH1F*  fWaveform;
    double fT0;				// simulated T0
    double fT01;			// reconstructed T0: CFD_01;
  };

  Event_t    fEvent[1000];

  TF1*       f_res;
  TRandom3*  fRn;

  double     fNPhotonsPerMeV;
  double     fPDE;
  int        n_phot;

  double     fCharge[kNBins];

  double     t_decay;
  double     sigma_res;


  TH1F*      fDt;

//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  SimPulse();
  ~SimPulse();

  static double  resolution(double* x, double* p);
  static double  pulse     (double* x, double* p);

  int  SimulateWaveform   (Event_t* Event);
  int  ReconstructWaveform(Event_t* Event);

  void Run(int Nev);

};

#endif
