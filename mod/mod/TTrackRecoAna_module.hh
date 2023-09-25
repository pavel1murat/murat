//
#ifndef __murat_inc_TTrackRecoAna_module_hh__
#define __murat_inc_TTrackRecoAna_module_hh__

#include "TObject.h"
#include "TObjArray.h"
#include "TString.h"

// #include "Offline/Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#ifndef __CINT__

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"

#include "Offline/BTrkData/inc/TrkStrawHit.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/alg/TEmuLogLH.hh"
#include "Stntuple/obj/TStnTrack.hh"
#include "Stntuple/obj/AbsEvent.hh"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#else

namespace art {
  class Event;
}

#endif

class TStnTrackBlock;
class TEmuLogLH;
class KalRep;

namespace mu2e {
  class StrawHit;
  class CaloCluster;
  class TrkToCaloExtrapol;
  class StepPointMC;
  class GenParticle;
  class SimParticle;
  class CalTimePeak;
  class TrackClusterMatch;
}

namespace mu2e {

class TTrackRecoAna : public art::EDAnalyzer {

  struct Doublet_t {
    const mu2e::TrkStrawHit*  fHit  [10];
    int                       fIndex[10];
    int                       fNHits;
  };

  struct DoubletRecoPar_t {
    const KalRep*             fTrack;
    int                       fIDWord;
    int                       fNHits;
    double                    fChi2T;
    int                       fNDoublets;
    int                       fNTriplets;
    int                       fNQuadruplets;
    const mu2e::TrkStrawHit*  fStrawHit  [200];
    Doublet_t                 fDoublet   [200];
    Doublet_t*                fDoubletPtr[200];
    double                    fChi2;
    double                    fSlope;
  };

  struct Hist_t {
    TH1F*  fNDoublets;
    TH1F*  fNHits;
    
    TH1F*  fChi2T;
    TH1F*  fChi2D;
    TH1F*  fSlope;
  };
 
  enum { fNHistSets = 2 };

protected:

  const art::Event*   fEvent;

  int                 fDiagLevel;
  std::string         fTrkPatRecLabel;

  DoubletRecoPar_t    fEle;

  Hist_t              fHist[fNHistSets];

public:

  TTrackRecoAna(fhicl::ParameterSet const& pset);
  ~TTrackRecoAna();
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  void analyze(const art::Event& e) override;

  void beginRun ( const art::Run& r) override;

  int  MakeDoublets(const KalRep* Track, DoubletRecoPar_t* Tp);

  void doubletmaker(const KalRep* EleTrk);

  void FillHistograms(Hist_t* Hist);
  void BookHistograms();

  virtual void beginJob();

};

}
#endif
