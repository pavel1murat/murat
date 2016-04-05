  struct TrackPar_t {
    int     fNHPl;
    int     fNEPl;
    int     fNDPl;
    int     fIDWord;
    int     fIDWord_A;
    int     fIDWord_01;

    float   fDpF ;    // tracker-only resolution
    float   fDp0 ;
    float   fDp2 ;
    float   fDpFSt;
    double  fDioWt;

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
