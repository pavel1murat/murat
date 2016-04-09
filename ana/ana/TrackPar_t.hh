#ifndef __murat_ana_TrackPar_t__
#define __murat_ana_TrackPar_t__

#include "TrackParBase_t.hh"

struct TrackPar_t : public TrackParBase_t {
  int     fNHPl;
  int     fNEPl;
  int     fNDPl;
  int     fIDWord[20];
  
  double  fP   ; 			// momentum, corrected to set DPF at 0
  float   fDpF ;			// tracker-only resolution
  float   fDp0 ;
  float   fDp2 ;
  float   fDpFSt;
  double  fDioWt;

  double  fDtZ0;
  
  double  fEcl;
  double  fEp;
  double  fDx;
  double  fDy;
  double  fDz;
  double  fDt;
  double  fDu;			// rotated residuals
  double  fDv;
  double  fChi2Match;
  double  fChi2XY;
  double  fChi2T;
  double  fPath;
  
  double  fSinTC;			// angle between the cluster and the track
  double  fDrTC;
  double  fSInt;
};
#endif
