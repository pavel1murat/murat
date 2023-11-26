//
// Read the tracks added to the event by the track finding and fitting code.
//
// $Id: TCalm001_module.cc,v 1.6 2014/08/11 23:09:52 murat Exp $
// $Author: murat $
// $Date: 2014/08/11 23:09:52 $
//
// Contact person Rob Kutschke
//

// Mu2e includes.
// #include "CLHEP/Geometry/HepPoint.h"
#include "BTrk/BbrGeom/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"

#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"
// storable objects (data products)
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// C++ includes.
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"

#include "Stntuple/obj/TStnTrack.hh"
#include "Stntuple/mod/THistModule.hh"

namespace mu2e {

  class TCalm001 : public THistModule {
  public:

    struct ClusterHist_t {
      TH1F*    fVaneID;
      TH1F*    fEnergy;
      TH1F*    fT0;
      TH1F*    fRow;
      TH1F*    fCol;
      TH1F*    fX;
      TH1F*    fY;
      TH1F*    fZ;
      TH1F*    fFrE12;
    };

    struct TrackHist_t {
      TH1F*    fPTot;
      TH1F*    fPt;
      TH1F*    fCosTh;
      TH1F*    fChi2;
      TH1F*    fNDof;
      TH1F*    fChi2Dof;
      TH1F*    fNActive;
      TH1F*    fT0;
      TH1F*    fQ;
      TH1F*    fFitCons;		// fit consistency (0 to 1)
      TH1F*    fD0;
      TH1F*    fZ0;
      TH1F*    fTanDip;
    };

    struct EventHist_t {
      TH1F*    fEleCosTh;
      TH1F*    fNClusters;
      TH1F*    fNTracks;
      TH1F*    fNStrawHits[2];
      TH1F*    fNGoodSH;
      TH1F*    fDtClT;
      TH1F*    fEMax;			// energy of the first reco cluster
      TH1F*    fDtClS;
      TH1F*    fNHyp;
      TH1F*    fBestHyp[2];
    };

    enum { kNEventHistSets   = 100 };
    enum { kNTrackHistSets   = 400 };
    enum { kNClusterHistSets = 100 };

    struct Hist_t {
      EventHist_t*   fEvent  [kNEventHistSets];
      TrackHist_t*   fTrack  [kNTrackHistSets];
      ClusterHist_t* fCluster[kNClusterHistSets];
    };

    Hist_t                 fHist;
//-----------------------------------------------------------------------------
// parameters
//-----------------------------------------------------------------------------
    std::string            fHistFileName;
    double                 fMinTActive;

    CaloClusterCollection* fListOfClusters;
    KalRepPtrCollection*   fListOfTracks[4]; // hypotheses: de-, ue+, dmu-, umu+
    GenParticleCollection* fListOfGenParticles;
    StrawHitCollection*    fListOfStrawHits;

    const KalRep*          fTrack;
    CaloCluster*           fCluster;

    GenParticle*           fEle;
    int                    fNClusters;
    int                    fNTracks[4];
    int                    fNGoodTracks; // for the downstream e- hypothesis
    int                    fNStrawHits;

    int                    fNHyp;
    int                    fBestHyp[2]; // [0]: by chi2, [1]: by fit consistency
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
    explicit TCalm001(fhicl::ParameterSet const& pset);
    ~TCalm001();

    void BookClusterHistograms(ClusterHist_t* Hist, const char* Folder);
    void BookEventHistograms  (EventHist_t*   Hist, const char* Folder);
    void BookTrackHistograms  (TrackHist_t*   Hist, const char* Folder);
    void BookHistograms();

    void FillClusterHistograms(ClusterHist_t* Hist, CaloCluster* Cls);
    void FillEventHistograms  (EventHist_t*   Hist);
    void FillTrackHistograms  (TrackHist_t*   Hist, KalRep*      Trk);
    void FillHistograms();

    void getData(const art::Event& event);
    void Init   ();

    void DefineBestHypothesis();
//-----------------------------------------------------------------------------
// overloaded framework module methods
//-----------------------------------------------------------------------------
    virtual bool filter  (art::Event& event);
    virtual void beginJob();
    virtual void endJob  ();

  };


//-----------------------------------------------------------------------------
  TCalm001::TCalm001(fhicl::ParameterSet const& pset): 
    THistModule(pset,"TCalm001"),
    fHistFileName(pset.get<std::string>   ("histFileName",   "tcalm001.hist")),
    fMinTActive  (pset.get<double>        ("minTActive"  ,   710.           ))
  {
					// reset all histogram pointers

    for (int i=0; i<kNEventHistSets; i++) {
      fHist.fEvent[i] = 0;
    }
  }


//-----------------------------------------------------------------------------
  TCalm001::~TCalm001() {
//     fFolder->Delete();
//     delete fFolder;
  }

//-----------------------------------------------------------------------------
  void TCalm001::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
//     char name [200];
//     char title[200];
//-----------------------------------------------------------------------------
//  
//-----------------------------------------------------------------------------
    HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
    HBook1F(Hist->fNClusters ,"ncl"      ,Form("%s: Number of Reconstructed Clusters",Folder),200,0,200,Folder);
    HBook1F(Hist->fNTracks   ,"ntrk"     ,Form("%s: Number of Reconstructed Tracks"  ,Folder),100,0,100,Folder);
    HBook1F(Hist->fNStrawHits[0],"nsh_0" ,Form("%s: Number of Straw Hits [0]"        ,Folder),250,0,250,Folder);
    HBook1F(Hist->fNStrawHits[1],"nsh_1" ,Form("%s: Number of Straw Hits [1]"        ,Folder),250,0,5000,Folder);
    HBook1F(Hist->fNGoodSH   ,"nsh50"    ,Form("%s: N(SH) +/-50"                     ,Folder),300,0,1500,Folder);
    HBook1F(Hist->fDtClT     ,"dt_clt"   ,Form("%s: DT(cluster-track)"               ,Folder),100,-100,100,Folder);
    HBook1F(Hist->fDtClS     ,"dt_cls"   ,Form("%s: DT(cluster-straw hit)"           ,Folder),200,-200,200,Folder);
    HBook1F(Hist->fEMax      ,"emax"     ,Form("%s: Max cluster energy"              ,Folder),150,0,150,Folder);
    HBook1F(Hist->fNHyp      ,"nhyp"     ,Form("%s: N(fit hypotheses)"               ,Folder),5,0,5,Folder);
    HBook1F(Hist->fBestHyp[0],"bfh_0"    ,Form("%s: Best Fit Hyp (e-,e+,mu-,mu+) chi2",Folder),5,0,5,Folder);
    HBook1F(Hist->fBestHyp[1],"bfh_1"    ,Form("%s: Best Fit Hyp (e-,e+,mu-,mu+) fcon",Folder),5,0,5,Folder);
  }

//-----------------------------------------------------------------------------
  void TCalm001::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
//     char name [200];
//     char title[200];
//-----------------------------------------------------------------------------
//  
//-----------------------------------------------------------------------------
    HBook1F(Hist->fPTot      ,"ptot" ,Form("%s: Track P(total)"   ,Folder),600,0,120,Folder);
    HBook1F(Hist->fPt        ,"pt"   ,Form("%s: Track Pt"         ,Folder),600,0,120,Folder);
    HBook1F(Hist->fCosTh     ,"costh",Form("%s: Track cos(theta)" ,Folder),100,-1,1,Folder);
    HBook1F(Hist->fChi2      ,"chi2" ,Form("%s: Track chi2 total" ,Folder),200, 0,200,Folder);
    HBook1F(Hist->fNDof      ,"ndof" ,Form("%s: Number of DOF"    ,Folder),200, 0,200,Folder);
    HBook1F(Hist->fChi2Dof   ,"chi2d",Form("%s: track chi2/N(dof)",Folder),100, 0, 10,Folder);
    HBook1F(Hist->fNActive   ,"nactv",Form("%s: N(active)"        ,Folder),200, 0,200,Folder);
    HBook1F(Hist->fT0        ,"t0"   ,Form("%s: track T0"         ,Folder),200, 0,2000,Folder);
    HBook1F(Hist->fQ         ,"q"    ,Form("%s: track Q"          ,Folder),  4,-2,   2,Folder);
    HBook1F(Hist->fFitCons   ,"fcon" ,Form("%s: track fit cons"   ,Folder),200, 0,   1,Folder);
    HBook1F(Hist->fD0        ,"d0"   ,Form("%s: track D0      "   ,Folder),200,-500, 500,Folder);
    HBook1F(Hist->fZ0        ,"z0"   ,Form("%s: track Z0      "   ,Folder),200,-20000,20000,Folder);
    HBook1F(Hist->fTanDip    ,"tdip" ,Form("%s: track tan(dip)"   ,Folder),100,-2.5 ,2.5,Folder);
  }

//-----------------------------------------------------------------------------
  void TCalm001::BookClusterHistograms(ClusterHist_t* Hist, const char* Folder) {
//     char name [200];
//     char title[200];
//-----------------------------------------------------------------------------
//  
//-----------------------------------------------------------------------------
    HBook1F(Hist->fVaneID ,"vane_id",Form("%s: Vane ID"       ,Folder), 10, 0,  10,Folder);
    HBook1F(Hist->fEnergy ,"energy" ,Form("%s: Cluster Energy",Folder),150, 0, 150,Folder);
    HBook1F(Hist->fT0     ,"t0"     ,Form("%s: cluster T0"    ,Folder),200, 0,2000,Folder);
    HBook1F(Hist->fRow    ,"row"    ,Form("%s: cluster Row"   ,Folder),200, 0, 200,Folder);
    HBook1F(Hist->fCol    ,"col"    ,Form("%s: cluster column",Folder),200, 0, 200,Folder);
    HBook1F(Hist->fX      ,"x"      ,Form("%s: cluster X"     ,Folder),200, 0,1000,Folder);
    HBook1F(Hist->fY      ,"y"      ,Form("%s: cluster Y"     ,Folder),200, 0,1000,Folder);
    HBook1F(Hist->fZ      ,"z"      ,Form("%s: cluster Z"     ,Folder),200, 0,1000,Folder);
    HBook1F(Hist->fFrE12  ,"fre12"  ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);
  }

//-----------------------------------------------------------------------------
  void TCalm001::BookHistograms() {
    TFolder  *hist_folder, *fol;
    char     folder_name[200];

    hist_folder = (TFolder*) fFolder->FindObject("Hist");

//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[0] = 1;		// all events
    book_event_histset[1] = 1;	        // events with a reconstructed track
    book_event_histset[2] = 1;	        // events without reconstructed tracks
    book_event_histset[3] = 1;	        // events with a reconstructed cluster
    book_event_histset[4] = 1;	        // events without reconstructed clusters
    book_event_histset[5] = 1;	        // events w/o reconstructed tracks, |costh|<0.4
    book_event_histset[6] = 1;	        // events with tracks passing "Set C" cuts

    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
	sprintf(folder_name,"evt_%i",i);
	fol = (TFolder*) hist_folder->FindObject(folder_name);
	if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
	fHist.fEvent[i] = new EventHist_t;
	BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
      }
    }
//-----------------------------------------------------------------------------
// book track histograms
//-----------------------------------------------------------------------------
    int book_track_histset[kNTrackHistSets];
    for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

    book_track_histset[  0] = 1;		// all tracks e- fit hyp
    book_track_histset[  1] = 1;		// tracks, e- fit hyp, events with clusters
    book_track_histset[  2] = 1;		// tracks, e- fit hyp, events without clusters
    book_track_histset[  3] = 1;		// tracks, e- fit hyp, tracks passing "C" cuts
    
    book_track_histset[100] = 1;		// all tracks e+
    book_track_histset[200] = 1;		// all tracks mu-
    book_track_histset[300] = 1;		// all tracks mu+

    for (int i=0; i<kNTrackHistSets; i++) {
      if (book_track_histset[i] != 0) {
	sprintf(folder_name,"trk_%i",i);
	fol = (TFolder*) hist_folder->FindObject(folder_name);
	if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
	fHist.fTrack[i] = new TrackHist_t;
	BookTrackHistograms(fHist.fTrack[i],Form("Hist/%s",folder_name));
      }
    }
//-----------------------------------------------------------------------------
// book cluster histograms
//-----------------------------------------------------------------------------
    int book_cluster_histset[kNClusterHistSets];
    for (int i=0; i<kNClusterHistSets; i++) book_cluster_histset[i] = 0;

    book_cluster_histset[0] = 1;		// all clusters
    book_cluster_histset[1] = 1;		// clusters in events with the reconstructed e-

    for (int i=0; i<kNClusterHistSets; i++) {
      if (book_cluster_histset[i] != 0) {
	sprintf(folder_name,"cls_%i",i);
	fol = (TFolder*) hist_folder->FindObject(folder_name);
	if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
	fHist.fCluster[i] = new ClusterHist_t;
	BookClusterHistograms(fHist.fCluster[i],Form("Hist/%s",folder_name));
      }
    }
  }

//-----------------------------------------------------------------------------
// begin job - book histograms etc
//-----------------------------------------------------------------------------
  void TCalm001::beginJob() {

    //    TFolder* fol;
    //    TFolder* hist_folder;
    
    //    char folder_name[200];

    BookHistograms();
  }


//-----------------------------------------------------------------------------
  void TCalm001::endJob() {
    TDirectory* old_dir = gDirectory;
    
    TFile* f = TFile::Open(fHistFileName.data(),"recreate");

    SaveFolder(fFolder,f);
    f->Write();
    f->Close();
    
    old_dir->cd();
  }

//-----------------------------------------------------------------------------
  void TCalm001::FillClusterHistograms(ClusterHist_t* Hist, CaloCluster* Cluster) {

    int   row, col;
    float  x, y, z;

    row = -999; // Cluster->cogRow();
    col = -999; // Cluster->cogColumn();
    x   = Cluster->cog3Vector().x();
    y   = Cluster->cog3Vector().y();
    z   = Cluster->cog3Vector().z();

//     if ((row < 0) || (row > 9999)) row = -9999;
//     if ((col < 0) || (col > 9999)) col = -9999;

    Hist->fVaneID->Fill(Cluster->diskID());
    Hist->fEnergy->Fill(Cluster->energyDep());
    Hist->fT0->Fill(Cluster->time());
    Hist->fRow->Fill(row);
    Hist->fCol->Fill(col);
    Hist->fX->Fill(x);
    Hist->fY->Fill(y);
    Hist->fZ->Fill(z);
  }

//-----------------------------------------------------------------------------
  void TCalm001::FillTrackHistograms(TrackHist_t* Hist, KalRep* Trk) {

    //    HelixParams h (Trk-helix(0.));

    Hist->fPTot->Fill(Trk->momentum().mag());
    Hist->fPt->Fill(Trk->momentum().perp());
    Hist->fCosTh->Fill(Trk->momentum().cosTheta());
    Hist->fChi2->Fill(Trk->chisq());
    Hist->fNDof->Fill(Trk->nDof());
    Hist->fChi2Dof->Fill(Trk->chisq()/Trk->nDof());
    Hist->fNActive->Fill(Trk->nActive());
    Hist->fT0->Fill(Trk->t0().t0());
    Hist->fQ->Fill(Trk->charge());
    Hist->fFitCons->Fill(Trk->chisqConsistency().consistency());

    Hist->fD0->Fill(Trk->helix(0).d0());
    Hist->fZ0->Fill(Trk->helix(0).z0());
    Hist->fTanDip->Fill(Trk->helix(0).tanDip());
  }

//-----------------------------------------------------------------------------
  void TCalm001::FillEventHistograms(EventHist_t* Hist) {

    double   cos_th;

    cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();

    Hist->fEleCosTh->Fill(cos_th);

    Hist->fNClusters->Fill(fNClusters);
    Hist->fNTracks->Fill(fNTracks[0]);
    Hist->fNStrawHits[0]->Fill(fNStrawHits);
    Hist->fNStrawHits[1]->Fill(fNStrawHits);

    double emax   = -1;
    double t0_cls = -1;
    double dt     = 9999.;
    if (fCluster) {
      emax   = fCluster->energyDep();
      t0_cls = fCluster->time();
    }

    double t0_trk = -1;
    if (fTrack) {
      t0_trk = fTrack->t0().t0();
    }

    if (fTrack && fCluster) {
      dt = t0_cls-t0_trk;
    }

    Hist->fDtClT->Fill(dt);
    Hist->fEMax->Fill(emax);

    StrawHit* sh;
    int n_good_hits = 0;
    for (int i=0; i<fNStrawHits; i++ ) {
      sh = &fListOfStrawHits->at(i);
      dt = t0_cls-sh->time() + 15;
      Hist->fDtClS->Fill(dt);

      if (fabs(dt+15.)< 50) n_good_hits += 1;
    }

    Hist->fNGoodSH->Fill(n_good_hits);

    Hist->fNHyp->Fill(fNHyp);
    Hist->fBestHyp[0]->Fill(fBestHyp[0]);
    Hist->fBestHyp[1]->Fill(fBestHyp[1]);
  }

//-----------------------------------------------------------------------------
  void TCalm001::DefineBestHypothesis() {
    double   chi2, chi2min(1.e6);
    double   fcon, fconmax(-1.);

    KalRep*  trk;

    fNHyp    = 0;
    fBestHyp[0] = -1;
    fBestHyp[1] = -1;

    for (int ih=0; ih<4; ih++) {

      if (fNTracks[ih] > 0) {
	fNHyp += 1;
	trk = (KalRep*) &fListOfTracks[ih]->at(0);

	chi2 = trk->chisq()/trk->nDof();
	if (chi2 < chi2min) {
	  chi2min     = chi2;
	  fBestHyp[0] = ih;
	}

	fcon = trk->chisqConsistency().consistency();
	if (fcon > fconmax) {
	  fconmax = fcon;
	  fBestHyp[1] = ih;
	}

      }
    }
  }

//-----------------------------------------------------------------------------
  void TCalm001::FillHistograms() {
    double cos_th;

    cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();
//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
    FillEventHistograms(fHist.fEvent[0]);

    if (fNTracks[0] > 0) FillEventHistograms(fHist.fEvent[1]);
    else                 FillEventHistograms(fHist.fEvent[2]);

    if (fNClusters  > 0) FillEventHistograms(fHist.fEvent[3]);
    else                 FillEventHistograms(fHist.fEvent[4]);

    if ((fNTracks[0] == 0) && (fabs(cos_th) < 0.4)) {
      FillEventHistograms(fHist.fEvent[5]); 
    }
					// fill event histograms for events with good tracks
    if (fNGoodTracks > 0) {
      FillEventHistograms(fHist.fEvent[6]); 
    }
//-----------------------------------------------------------------------------
// track histograms
//-----------------------------------------------------------------------------
    KalRep*      trk;

    int loc;
    for (int i=0; i<4; i++) {
      loc = 100*i;
      for (int it=0; it<fNTracks[i]; it++) {
	trk = (KalRep*) &fListOfTracks[i]->at(it);
	FillTrackHistograms(fHist.fTrack[loc],trk);
      }
    }
//-----------------------------------------------------------------------------
// track histograms ID=1: events with    clusters
// track histograms ID=2: events without clusters
// these are for downstream electron ("e-d") fit only
//-----------------------------------------------------------------------------
    for (int it=0; it<fNTracks[0]; it++) {
      trk = (KalRep*) &fListOfTracks[0]->at(it);
      if (fNClusters > 0) FillTrackHistograms(fHist.fTrack[1],trk);
      else                FillTrackHistograms(fHist.fTrack[2],trk);
      //      if ( ### for tracks passing "C" cuts
    }

    //    HelixParams  helix;
//-----------------------------------------------------------------------------
// cluster histograms, assume one list of tracks
//-----------------------------------------------------------------------------
    CaloCluster* cl;
    for (int i=0; i<fNClusters; ++i ) {
      cl = /*(CaloCluster*)*/ &fListOfClusters->at(i);
      FillClusterHistograms(fHist.fCluster[0],cl);
      if (fNTracks[0] > 0) {
	FillClusterHistograms(fHist.fCluster[1],cl);
      }
    }
  }


//-----------------------------------------------------------------------------
  void TCalm001::getData(const art::Event& event) {
//-----------------------------------------------------------------------------
// generated particles
//-----------------------------------------------------------------------------
    art::Handle<GenParticleCollection> genpHandle;
    event.getByLabel("generate",genpHandle);
    fListOfGenParticles = (GenParticleCollection*) &(*genpHandle);
    fEle                = (GenParticle*) &fListOfGenParticles->at(0);
//-----------------------------------------------------------------------------
// reconstructed tracks
//-----------------------------------------------------------------------------
    art::Handle<KalRepPtrCollection> krepsHandle;
    event.getByLabel("trkPatRec1","DownstreameMinus", krepsHandle);
    fListOfTracks[0] = (KalRepPtrCollection*) krepsHandle.product();
    fNTracks     [0] = fListOfTracks[0]->size();

    fTrack   = 0;
    if (fNTracks[0]) fTrack = fListOfTracks[0]->at(0).get();
//-----------------------------------------------------------------------------
// alternative hypotheses
//-----------------------------------------------------------------------------
    art::Handle<KalRepPtrCollection> krepsHandle2;
    event.getByLabel("trkPatRec2","UpstreamePlus", krepsHandle2);
    fListOfTracks[1] = (KalRepPtrCollection*) krepsHandle2.product();
    fNTracks     [1] = fListOfTracks[1]->size();

    art::Handle<KalRepPtrCollection> krepsHandle3;
    event.getByLabel("trkPatRec3","DownstreammuMinus", krepsHandle3);
    fListOfTracks[2] = (KalRepPtrCollection*) krepsHandle3.product();
    fNTracks     [2] = fListOfTracks[2]->size();

    art::Handle<KalRepPtrCollection> krepsHandle4;
    event.getByLabel("trkPatRec4","UpstreammuPlus", krepsHandle4);
    fListOfTracks[3] = (KalRepPtrCollection*) krepsHandle4.product();
    fNTracks     [3] = fListOfTracks[3]->size();
//-----------------------------------------------------------------------------
// reconstructed calorimeter clusters
//-----------------------------------------------------------------------------
    art::Handle<CaloClusterCollection> calo_cluster_handle;
    event.getByLabel("makeCaloCluster","AlgoCLOSESTSeededByENERGY",calo_cluster_handle);
    fListOfClusters = (CaloClusterCollection*) &(*calo_cluster_handle);
    fNClusters      = fListOfClusters->size();

    fCluster        = 0;
    if (fNClusters > 0) fCluster = &fListOfClusters->at(0);
//-----------------------------------------------------------------------------
// straw hits 
//-----------------------------------------------------------------------------
    art::Handle<StrawHitCollection> shHandle;
    event.getByLabel("makeSH",shHandle);
    fListOfStrawHits = (StrawHitCollection*) &(*shHandle);
    fNStrawHits      = fListOfStrawHits->size();
  }

//-----------------------------------------------------------------------------
  void TCalm001::Init() {
    double   fitmom_err, tan_dip;
    const KalRep*  trk;

    fNGoodTracks = 0;
    for (int i=0; i<fNTracks[0]; ++i ) {
      trk = fListOfTracks[0]->at(i).get();

      // get the fit at the first hit
      //      CLHEP::Hep3Vector entpos = det->toDetector(vdg->getGlobal(VirtualDetectorId::TT_FrontPA));
      //      double zent           = entpos.z();
      double firsthitfltlen = trk->firstHit()->kalHit()->hit()->fltLen() - 10;
      double lasthitfltlen  = trk->lastHit ()->kalHit()->hit()->fltLen() - 10;
      double entlen         = std::min(firsthitfltlen,lasthitfltlen);

      //      TrkHelixUtils::findZFltlen(trk->traj(),zent,entlen,0.1); 

      CLHEP::Hep3Vector fitmom = trk->momentum(entlen);
      //      double     ptot   = fitmom.mag();
      CLHEP::Hep3Vector momdir = fitmom.unit();

      BbrVectorErr momerr = trk->momentumErr(entlen);

      CLHEP::HepVector momvec(3);

      for(int icor=0;icor<3;icor++) momvec[icor] = momdir[icor];

      fitmom_err = sqrt(momerr.covMatrix().similarity(momvec));
	
      tan_dip = trk->helix(0).tanDip();

      if ((trk->chisqConsistency().consistency() > 0.001      ) &&
	  (trk->t0().t0()                        > fMinTActive) &&
	  (trk->nActive()                        > 25         ) &&
	  (trk->t0().t0Err()                     < 1.0        ) &&
          (fitmom_err                            < 0.18       ) && 
	  (fabs(tan_dip)                         > 0.577      ) && 
	  (fabs(tan_dip)                         < 1.0        )    ) {
	fNGoodTracks += 1;
      }
    }

    DefineBestHypothesis();
  }
//-----------------------------------------------------------------------------
  bool TCalm001::filter(art::Event& event) {
    const char* name = "MyExo::analyze";

    printf(" >>>>>>> [%s] EVENT : %10i\n",name,event.event());

    getData(event);

    Init();
//-----------------------------------------------------------------------------
// calculate number of good tracks
//-----------------------------------------------------------------------------

    FillHistograms();

    return true;
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::TCalm001)
