#ifndef __murat_ana_TrackPar_t__
#define __murat_ana_TrackPar_t__

#include "TrackParBase_t.hh"

struct TrackPar_t : public TrackParBase_t {
  int     fNHPl;
  int     fNEPl;
  int     fNDPl;
  int     fIDWord[20];

  double  fMVAOut[20];			// outputs of different MVA classifiers
  double  fProb;
  
  double  fP   ; 			// momentum, corrected to set DPF at 0
  float   fDpF ;			// tracker-only resolution
  float   fDp0 ;
  float   fDp2 ;
  float   fDpFSt;
  double  fDioWt;			// DIO LO weight
  double  fDioWtRC;			// DIO LO weight with rad corrections
  double  fLumWt;			// luminosity weight
  double  fTotWt;			// total weight
  double  fTotWtRC;			// total weight with rad corrections

  double  fDtZ0;			// delta(T) at z=0
  double  fDtBack;			// delta(T) at z=Z(TT_Back)
  
  double  fEcl;
  double  fEp;
  double  fDx;
  double  fDy;
  double  fDz;
  double  fDt;
  double  fDu;			// rotated residuals
  double  fDv;
  double  fChi2Tcm;
  double  fChi2XY;
  double  fChi2T;
  double  fPath;
  
  double  fSinTC;			// angle between the cluster and the track
  double  fDrTC;
  double  fSInt;

  double  fLogLHDedm;			// downstram electron (DE) vs downstream muon (DM)

  double  fTMean;			// mean time over the track hits
};
#endif
