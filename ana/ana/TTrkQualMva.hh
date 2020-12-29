#ifndef __murat_ana_TTrkQualMva_hh__
#define __murat_ana_TTrkQualMva_hh__

class TBranch;

namespace murat {
  struct TmvaTrainData_t {
    float    fP;
    float    fPMC;
    float    fTanDip;
    float    fNActive;
    float    fNaFract;
    float    fNDoublets;                // total number of doublets
    float    fNDa;                      // number of doublets with all hits active
    float    fChi2Dof;
    float    fFitCons;
    float    fMomErr;
    float    fT0Err;
    float    fD0;
    float    fRMax;
    float    fNdaOverNa;
    float    fNzaOverNa;
    float    fNmaOverNm;
    float    fNdaOverNd;
    float    fZ1;			// Z-coordinate of the first hit
    float    fWeight;
  };

  struct TmvaTrainBranches_t {
    TBranch*  fP;
    TBranch*  fPMC;
    TBranch*  fTanDip;
    TBranch*  fNActive;
    TBranch*  fNaFract;
    TBranch*  fNDoublets;
    TBranch*  fNDa;
    TBranch*  fChi2Dof;
    TBranch*  fFitCons;
    TBranch*  fMomErr;
    TBranch*  fT0Err;
    TBranch*  fD0;
    TBranch*  fRMax;
    TBranch*  fNdaOverNa;
    TBranch*  fNzaOverNa;
    TBranch*  fNmaOverNm;
    TBranch*  fNdaOverNd;
    TBranch*  fZ1;			       // Z-coordinate of the first hit
    TBranch*  fWeight;			       // for background only
  };
}
#endif
