//-----------------------------------------------------------------------------
// 1. search for downstream tracks
// debug bits:
// -----------
//  0: DEM (pt > 80.) && (t->fNActive > 35) && (t->NClusters() == 0)
// 11: DEM tracks passing all ID cuts, but not electron
// 12: for all events with DEM tracks print track IDWORD, E/P and DZ
//-----------------------------------------------------------------------------
#include "TEnv.h"

#include "murat/mod/TCalm004_module.hh"

#include "Stntuple/val/stntuple_val_functions.hh"
#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"

#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "BTrk/KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "TrackCaloMatching/inc/TrkToCaloExtrapol.hh"

namespace mu2e {

//-----------------------------------------------------------------------------
  TCalm004::TCalm004(fhicl::ParameterSet const& pset): 
    THistModule(pset,"TCalm004"),
    fHistFileName (pset.get<std::string> ("histFileName" , "tcalm004.hist")),
    fStrawHitMaker(pset.get<std::string> ("strawHitMaker", "undefined"    )),
    fClusterMaker (pset.get<std::string> ("clusterMaker" , "undefined"    )),

    fTrkPatRecDem (pset.get<std::string> ("trkPatRecDem" , "trkPatRecDem" )),
    fTrkPatRecUem (pset.get<std::string> ("trkPatRecUem" , "trkPatRecUem" )),

    fTrkExtrapolDem (pset.get<std::string> ("trkExtrapolDem" , "tkExtrapolDem")),
    fTrkCalMatchDem (pset.get<std::string> ("caloMatchingDem", "undefined"    )),

    fTrkExtrapolUem (pset.get<std::string> ("trkExtrapolUem" , "tkExtrapolUem")),
    fTrkCalMatchUem (pset.get<std::string> ("caloMatchingUem", "undefined"    )),
    fPidDem         (pset.get<std::string> ("pidDem"         , "undefined"    )),

    fMinTActive   (pset.get<double>      ("minTActive"   ,   710.         )),
    fMinECrystal  (pset.get<double>      ("minECrystal"  ,    0.1         ))
  {
					// reset all histogram pointers

    for (int i=0; i<kNEventHistSets; i++) {
      fHist.fEvent[i] = 0;
    }

    fTrackBlock[0]  = new TStnTrackBlock  ();
    fTrackBlock[1]  = new TStnTrackBlock  ();
    fClusterBlock   = new TStnClusterBlock();

    fTrackIDSetC    = new TStnTrackID();
    fTrackIDSetB    = new TStnTrackID();

    fTrackIDSetB->SetMinNActive  (20   );
    fTrackIDSetB->SetMaxT0Err    (1.0  );
    fTrackIDSetB->SetMaxFitMomErr(0.2  );
    fTrackIDSetB->SetMinFitCons  (1.e-4);

    fLogLH          = new TEmuLogLH();
  }



//-----------------------------------------------------------------------------
  TCalm004::~TCalm004() {
//     fFolder->Delete();
//     delete fFolder;

    delete fTrackBlock[0];
    delete fTrackBlock[1];
    delete fClusterBlock;
    delete fTrackIDSetB;
    delete fTrackIDSetC;
    delete fLogLH;
  }

//-----------------------------------------------------------------------------
  void TCalm004::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
//     char name [200];
//     char title[200];
//-----------------------------------------------------------------------------
//  
//-----------------------------------------------------------------------------
    HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
    HBook1F(Hist->fRv        ,"rv"       ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 1000,Folder);
    HBook1F(Hist->fZv        ,"zv"       ,Form("%s: Z(Vertex)"                       ,Folder), 300, 0,15000,Folder);
    HBook1F(Hist->fNClusters ,"ncl"      ,Form("%s: Number of Reconstructed Clusters",Folder),200,0,200,Folder);
    HBook1F(Hist->fNTracksDem,"ntrk_dem" ,Form("%s: N(Reconstructed Tracks) DEM"     ,Folder),100,0,100,Folder);
    HBook1F(Hist->fNTracksUem,"ntrk_uem" ,Form("%s: N(Reconstructed Tracks) UEM"     ,Folder),100,0,100,Folder);
    HBook1F(Hist->fNStrawHits[0],"nsh_0" ,Form("%s: Number of Straw Hits [0]"        ,Folder),250,0,250,Folder);
    HBook1F(Hist->fNStrawHits[1],"nsh_1" ,Form("%s: Number of Straw Hits [1]"        ,Folder),250,0,5000,Folder);
    HBook1F(Hist->fNGoodSH   ,"nsh50"    ,Form("%s: N(SH) +/-50"                     ,Folder),300,0,1500,Folder);
    HBook1F(Hist->fDtClT     ,"dt_clt"   ,Form("%s: DT(cluster-track)"               ,Folder),100,-100,100,Folder);
    HBook1F(Hist->fDtClS     ,"dt_cls"   ,Form("%s: DT(cluster-straw hit)"           ,Folder),200,-200,200,Folder);
    HBook1F(Hist->fEMax      ,"emax"     ,Form("%s: Max cluster energy"              ,Folder),150,0,150,Folder);
    HBook1F(Hist->fNHyp      ,"nhyp"     ,Form("%s: N(fit hypotheses)"               ,Folder),5,0,5,Folder);
    HBook1F(Hist->fBestHyp[0],"bfh0"     ,Form("%s: Best Fit Hyp[0](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
    HBook1F(Hist->fBestHyp[1],"bfh1"     ,Form("%s: Best Fit Hyp[1](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  }

//-----------------------------------------------------------------------------
  void TCalm004::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
//     char name [200];
//     char title[200];
//-----------------------------------------------------------------------------
//  
//-----------------------------------------------------------------------------
    HBook1F(Hist->fP          ,"p"        ,Form("%s: Track P(total)"    ,Folder),1000, 0,200,Folder);
    HBook1F(Hist->fPt         ,"pt"       ,Form("%s: Track Pt"          ,Folder), 600, 0,120,Folder);
    HBook1F(Hist->fCosTh      ,"costh"    ,Form("%s: Track cos(theta)"  ,Folder), 100,-1,1,Folder);
    HBook1F(Hist->fChi2       ,"chi2"     ,Form("%s: Track chi2 total"  ,Folder), 200, 0,200,Folder);
    HBook1F(Hist->fNDof       ,"ndof"     ,Form("%s: Number of DOF"     ,Folder), 200, 0,200,Folder);
    HBook1F(Hist->fChi2Dof    ,"chi2d"    ,Form("%s: track chi2/N(dof)" ,Folder), 100, 0, 10,Folder);
    HBook1F(Hist->fNActive    ,"nactv"    ,Form("%s: N(active)"         ,Folder), 200, 0,200,Folder);
    HBook1F(Hist->fT0         ,"t0"       ,Form("%s: track T0"          ,Folder), 200, 0,2000,Folder);
    HBook1F(Hist->fQ          ,"q"        ,Form("%s: track Q"           ,Folder),   4,-2,   2,Folder);
    HBook1F(Hist->fFitCons    ,"fcon"     ,Form("%s: track fit cons"    ,Folder), 200, 0,   1,Folder);
    HBook1F(Hist->fD0         ,"d0"       ,Form("%s: track D0      "    ,Folder), 200,-500, 500,Folder);
    HBook1F(Hist->fZ0         ,"z0"       ,Form("%s: track Z0      "    ,Folder), 200,-20000,20000,Folder);
    HBook1F(Hist->fTanDip     ,"tdip"     ,Form("%s: track tan(dip)"    ,Folder), 100,-2.5 ,2.5,Folder);
    HBook1F(Hist->fDt         ,"dt"       ,Form("%s: track delta(T)"    ,Folder), 200,-20  ,20 ,Folder);
    HBook1F(Hist->fDy         ,"dy"       ,Form("%s: track delta(Y)"    ,Folder), 100,-500 ,500,Folder);
    HBook1F(Hist->fDz         ,"dz"       ,Form("%s: track delta(Z)"    ,Folder), 100,-500 ,500,Folder);
    HBook1F(Hist->fNClusters  ,"ncl"      ,Form("%s: track N(clusters)" ,Folder),  10, 0   , 10,Folder);
    HBook1F(Hist->fVaneID     ,"vid"      ,Form("%s: track vane ID"     ,Folder),  10,-5   ,  5,Folder);
    HBook1F(Hist->fXCal       ,"xcal"     ,Form("%s: track XCal"        ,Folder), 200,-1000,1000,Folder);
    HBook1F(Hist->fYCal       ,"ycal"     ,Form("%s: track YCal"        ,Folder), 200,-1000,1000,Folder);
    HBook1F(Hist->fZCal       ,"zcal"     ,Form("%s: track ZCal"        ,Folder), 200, 1500,3500,Folder);
    HBook1F(Hist->fYTrk       ,"ytrk"     ,Form("%s: track YTrk"        ,Folder), 100,-250,250,Folder);
    HBook1F(Hist->fZTrk       ,"ztrk"     ,Form("%s: track ZTrk"        ,Folder), 200,-1000,1000,Folder);
    HBook1F(Hist->fEp         ,"ep"       ,Form("%s: track E/P"         ,Folder), 300, 0   ,1.5,Folder);
    HBook2F(Hist->fNHVsStation,"nh_vs_st" ,Form("%s: N(hits) Vs Station",Folder),  40, 0,40,10,-0.5,9.5,Folder);
    HBook2F(Hist->fNHVsNSt    ,"nh_vs_nst",Form("%s: N(hits) Vs NSt"    ,Folder),  10,-0.5,9.5,40,-0.5,39.5,Folder);

    HBook1F(Hist->fRSlope     ,"rslope"   ,Form("%s: Res Slope"         ,Folder), 200,-20 , 20,Folder);
    HBook1F(Hist->fXSlope     ,"xslope"   ,Form("%s: Res/Sig Slope"     ,Folder), 200,-20 , 20,Folder);

    HBook1F(Hist->fEleLogLHCal,"ele_llh_c",Form("%s: ELE Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
    HBook1F(Hist->fMuoLogLHCal,"muo_llh_c",Form("%s: MUO Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
    HBook1F(Hist->fLogLHRCal  ,"llhr_cal" ,Form("%s: LogLH(e/m) Cal"    ,Folder), 200,-100,100,Folder);
    //    HBook1F(Hist->fLogLHRSlope,"lhr_slope",Form("%s: LogLHR(e-m) Slope" ,Folder), 200,-20 , 20,Folder);
    HBook1F(Hist->fLogLHRDeDx ,"llhr_dedx",Form("%s: LogLH(e/m) De/Dx"  ,Folder), 200, -20, 20,Folder);
    HBook1F(Hist->fLogLHRXs   ,"llhr_xs"  ,Form("%s: LogLH(e/m) Xs"     ,Folder), 200, -20, 20,Folder);
    HBook1F(Hist->fLogLHRTrk  ,"llhr_trk" ,Form("%s: LogLH(e/m) Trk"    ,Folder), 200, -20, 20,Folder);

    HBook1F(Hist->fPdgCode    ,"pdg"      ,Form("%s: track PDG code"    ,Folder), 100,-50,50,Folder);
    HBook1F(Hist->fFrGH       ,"fgh"      ,Form("%s: Fraction Goog Hits",Folder), 100, 0,1,Folder);
  }

//-----------------------------------------------------------------------------
  void TCalm004::BookClusterHistograms(ClusterHist_t* Hist, const char* Folder) {
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
    HBook1F(Hist->fX      ,"x"      ,Form("%s: cluster X"     ,Folder),200, -5000,5000,Folder);
    HBook1F(Hist->fY      ,"y"      ,Form("%s: cluster Y"     ,Folder),200,-1000,1000,Folder);
    HBook1F(Hist->fZ      ,"z"      ,Form("%s: cluster Z"     ,Folder),200, 0,20000,Folder);
    HBook1F(Hist->fR      ,"r"      ,Form("%s: cluster Radius",Folder),100, 0,  100,Folder);
    HBook1F(Hist->fYMean  ,"ymean"  ,Form("%s: cluster YMean" ,Folder),400,-200,200,Folder);
    HBook1F(Hist->fZMean  ,"zmean"  ,Form("%s: cluster ZMean" ,Folder),400,-200,200,Folder);
    HBook1F(Hist->fSigY   ,"sigy"   ,Form("%s: cluster SigY"  ,Folder),100, 0,100,Folder);
    HBook1F(Hist->fSigZ   ,"sigz"   ,Form("%s: cluster SigZ"  ,Folder),100, 0,100,Folder);
    HBook1F(Hist->fSigR   ,"sigr"   ,Form("%s: cluster SigR"  ,Folder),100, 0,100,Folder);
    HBook1F(Hist->fNCr0   ,"ncr0"   ,Form("%s: cluster NCR[0]",Folder),100, 0,100,Folder);
    HBook1F(Hist->fNCr1   ,"ncr1"   ,Form("%s: cluster NCR[1]",Folder),100, 0,100,Folder);
    HBook1F(Hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
    HBook1F(Hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);
    HBook1F(Hist->fSigE1  ,"sige1"   ,Form("%s: SigmaE/Etot"  ,Folder),200, 0, 10,Folder);
    HBook1F(Hist->fSigE2  ,"sige2"   ,Form("%s: SigmaE/Emean" ,Folder),200, 0, 10,Folder);
  }

//-----------------------------------------------------------------------------
  void TCalm004::BookHistograms() {
    TFolder  *hist_folder, *fol;
    char     folder_name[200];

    hist_folder = (TFolder*) fFolder->FindObject("Hist");

//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[0] = 1;		// all events
    book_event_histset[1] = 1;	        // events with an upstream track
    book_event_histset[2] = 1;	        // events without an upstream track

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

    book_track_histset[  0] = 1;		// all tracks e-
    book_track_histset[  1] = 1;		// all tracks e- passing Set C cuts 
    book_track_histset[  2] = 1;		// all tracks e- passing Set C cuts, events with clusters 
    book_track_histset[  3] = 1;		// all tracks e- passing Set C cuts, events w/o  clusters
    book_track_histset[  4] = 1;		// all tracks e- passing Set C cuts, events with clusters, no closest
    book_track_histset[  5] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns

    book_track_histset[ 20] = 1;		// all tracks e-
    book_track_histset[ 21] = 1;		// all tracks e- passing Set C cuts 
    book_track_histset[ 22] = 1;		// all tracks e- passing Set C cuts, events with clusters 
    book_track_histset[ 23] = 1;		// all tracks e- passing Set C cuts, events w/o  clusters
    book_track_histset[ 24] = 1;		// all tracks e- passing Set C cuts, events with clusters, no closest
    book_track_histset[ 25] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 26] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 27] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 28] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 29] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 30] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 31] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 32] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 33] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 34] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 35] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 36] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 37] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 38] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns

    book_track_histset[ 41] = 1;		// all tracks e- passing Set B cuts, |D0|<12cm 
    book_track_histset[ 42] = 1;		// all tracks e- passing Set B cuts, |D0|<12cm, cluster
    book_track_histset[ 43] = 1;		// all tracks e- passing Set B cuts, |D0|<12cm, no cluster
    book_track_histset[ 44] = 1;		// all tracks e- passing Set B cuts, |D0|<12cm, cluster, ID:ele
    book_track_histset[ 45] = 1;		// all tracks e- passing Set B cuts, |D0|<12cm, cluster, ID:muo
    book_track_histset[ 46] = 1;		// 
    book_track_histset[ 47] = 1;		// 
    book_track_histset[ 48] = 1;		// 
    book_track_histset[ 49] = 1;		// 
    book_track_histset[ 50] = 1;		// 
    book_track_histset[ 51] = 1;		// 
    book_track_histset[ 52] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 53] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 54] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 55] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 56] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 57] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
    book_track_histset[ 58] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns


    book_track_histset[100] = 1;		// all upstream tracks e-
    book_track_histset[101] = 1;		// all upstream tracks e- passing Set C cuts 
    book_track_histset[102] = 1;		// all upstream tracks e- passing Set C cuts, events with clusters 
    book_track_histset[103] = 1;		// all upstream tracks e- passing Set C cuts, events w/o  clusters
    book_track_histset[104] = 1;		// all upstream tracks e- passing Set C cuts, events with clusters, no closest
    book_track_histset[105] = 1;		// all upstream tracks e- passing Set C cuts, |dt| < 2.5ns

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
    book_cluster_histset[2] = 1;		// clusters in events with the track passing SetC cuts
    book_cluster_histset[3] = 1;		// clusters in events w/track passing SetC cuts and |dt|<2.5ns 

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
  void TCalm004::beginJob() {

    BookHistograms();

    fTrackIDSetC->SetMinT0(fMinTActive);
    fTrackIDSetB->SetMinT0(fMinTActive);
//-----------------------------------------------------------------------------
// define collection names to be used for initialization
//-----------------------------------------------------------------------------
    fClusterBlock->AddCollName("mu2e::CaloClusterCollection"          ,fClusterMaker.data(),"AlgoCLOSESTSeededByENERGY");
    fClusterBlock->AddCollName("mu2e::TrackClusterMatchCollection"    ,fTrkCalMatchDem.data(),"");

    fTrackBlock[0]->AddCollName("mu2e::KalRepCollection"              ,fTrkPatRecDem.data(),"DownstreameMinus");
    fTrackBlock[0]->AddCollName("mu2e::CaloClusterCollection"         ,fClusterMaker.data(),"AlgoCLOSESTSeededByENERGY");
    fTrackBlock[0]->AddCollName("mu2e::TrkToCaloExtrapolCollection"   ,fTrkExtrapolDem.data(),"");
    fTrackBlock[0]->AddCollName("mu2e::TrackClusterMatchCollection"   ,fTrkCalMatchDem.data(),"");
    fTrackBlock[0]->AddCollName("mu2e::StrawHitCollection"            ,fStrawHitMaker.data(),"");
    fTrackBlock[0]->AddCollName("mu2e::PtrStepPointMCVectorCollection",fStrawHitMaker.data(),"StrawHitMCPtr");
    fTrackBlock[0]->AddCollName("mu2e::PIDProductCollection"          ,fPidDem.data()       ,"");

    // no PID for upstream particles, should not crash

    fTrackBlock[1]->AddCollName("mu2e::KalRepCollection"              ,fTrkPatRecUem.data(),"UpstreameMinus");
    fTrackBlock[1]->AddCollName("mu2e::CaloClusterCollection"         ,fClusterMaker.data(),"AlgoCLOSESTSeededByENERGY");
    fTrackBlock[1]->AddCollName("mu2e::TrkToCaloExtrapolCollection"   ,fTrkExtrapolUem.data(),"");
    fTrackBlock[1]->AddCollName("mu2e::TrackClusterMatchCollection"   ,fTrkCalMatchUem.data(),"");
    fTrackBlock[1]->AddCollName("mu2e::StrawHitCollection"            ,fStrawHitMaker.data(),"");
    fTrackBlock[1]->AddCollName("mu2e::PtrStepPointMCVectorCollection",fStrawHitMaker.data(),"StrawHitMCPtr");
//-----------------------------------------------------------------------------
// initialize the likelihood histograms
// so far, the likelihood histograms are based on 
//-----------------------------------------------------------------------------
    const char*  mu2e_hist_dir;
    char         fn[200];

    mu2e_hist_dir = gEnv->GetValue("mu2e.HistDir","mu2e.HistDir.undefined");

    sprintf(fn,"%s/v2_1_0/e0000001.tcalm002.hist",mu2e_hist_dir);
    fLogLH->InitEleDtHist(fn);
    fLogLH->InitEleEpHist(fn);

    sprintf(fn,"%s/v2_1_0/m0000001.tcalm002.hist",mu2e_hist_dir);
    fLogLH->InitMuoDtHist(fn);
    fLogLH->InitMuoEpHist(fn);

    sprintf(fn,"%s/v3_0_0/e0000001_tcalm002_vane.hist",mu2e_hist_dir);
    TH1* he_xs = gh1(fn,"TCalm002","trk_1/xslope");
    //    TH1* he_de = gh1(fn,"TCalm002","trk_1/lhr_dedx");
    fLogLH->SetEleXsHist(he_xs);
    //    fLogLH->SetEleDeHist(he_de);

    sprintf(fn,"%s/v3_0_0/m0000001_tcalm002_vane.hist",mu2e_hist_dir);
    TH1* hm_xs = gh1(fn,"TCalm002","trk_1/xslope");
    //    TH1* hm_de = gh1(fn,"TCalm002","trk_1/lhr_dedx");
    fLogLH->SetMuoXsHist(hm_xs);
    //    fLogLH->SetMuoDeHist(hm_de);
  }


//-----------------------------------------------------------------------------
  void TCalm004::endJob() {
    TDirectory* old_dir = gDirectory;
    
    TFile* f = TFile::Open(fHistFileName.data(),"recreate");

    SaveFolder(fFolder,f);
    f->Write();
    f->Close();
    
    old_dir->cd();
  }

//-----------------------------------------------------------------------------
  void TCalm004::FillClusterHistograms(ClusterHist_t* Hist, TStnCluster* Cluster) {

    int   row, col;
    float  x, y, z, r;

    const CaloCluster* cl = Cluster->fCaloCluster;

    row = -999; // cl->cogRow();
    col = -999; // cl->cogColumn();
    x   = cl->cog3Vector().x();
    y   = cl->cog3Vector().y();
    z   = cl->cog3Vector().z();
    r   = sqrt(x*x+y*y);

    if ((row < 0) || (row > 9999)) row = -9999;
    if ((col < 0) || (col > 9999)) col = -9999;

    Hist->fVaneID->Fill(cl->sectionId());
    Hist->fEnergy->Fill(cl->energyDep());
    Hist->fT0->Fill(cl->time());
    Hist->fRow->Fill(row);
    Hist->fCol->Fill(col);
    Hist->fX->Fill(x);
    Hist->fY->Fill(y);
    Hist->fZ->Fill(z);
    Hist->fR->Fill(r);

    Hist->fYMean->Fill(Cluster->fYMean);
    Hist->fZMean->Fill(Cluster->fZMean);
    Hist->fSigY->Fill(Cluster->fSigY);
    Hist->fSigZ->Fill(Cluster->fSigZ);
    Hist->fSigR->Fill(Cluster->fSigR);
    Hist->fNCr0->Fill(Cluster->fNCrystals);
    Hist->fNCr1->Fill(Cluster->fNCr1);
    Hist->fFrE1->Fill(Cluster->fFrE1);
    Hist->fFrE2->Fill(Cluster->fFrE2);
    Hist->fSigE1->Fill(Cluster->fSigE1);
    Hist->fSigE2->Fill(Cluster->fSigE2);
  }

//-----------------------------------------------------------------------------
  void TCalm004::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track) {

    //    HelixParams h (Trk-helix(0.));

    KalRep* trk = Track->fKalRep[0];

    Hist->fP->Fill (Track->fP);
    Hist->fPt->Fill(Track->fPt);
    Hist->fCosTh->Fill(Track->Momentum()->CosTheta());
    Hist->fChi2->Fill (Track->fChi2);
    Hist->fNDof->Fill(Track->fNActive-5.);
    Hist->fChi2Dof->Fill(Track->fChi2/(Track->fNActive-5.));
    Hist->fNActive->Fill(Track->fNActive);
    Hist->fT0->Fill(Track->fT0);
    Hist->fQ->Fill(trk->charge());
    Hist->fFitCons->Fill(Track->fFitCons);

    Hist->fD0->Fill(trk->helix(0).d0());
    Hist->fZ0->Fill(trk->helix(0).z0());
    Hist->fTanDip->Fill(trk->helix(0).tanDip());

    int nh, nst_with_nh[10];

    for (int i=0; i<10; i++) nst_with_nh[i] = 0;

    for (int i=0; i<40; i++) {
      Hist->fNHVsStation->Fill(i,Track->fNHPerStation[i]);
      nh = Track->fNHPerStation[i];
      if (nh < 10) {
	nst_with_nh[nh] += 1;
      }
    }

    for (int i=0; i<10; i++) {
      Hist->fNHVsNSt->Fill(i,nst_with_nh[i]);
    }
//-----------------------------------------------------------------------------
// track-cluster matching part: 
// - for residuals, determine intersection with the most energetic cluster
// - for track -only parameters use intersection with lowest trajectory length
//-----------------------------------------------------------------------------
    TStnTrack::InterData_t*    vt = Track->fVMinS;  // track-only
    //    TStnTrack::InterData_t*    vr = Track->fVMaxEp; // residuals

    double  yt, zt;
    int     vid;

    if (vt) {
      vid = vt->fID;
      yt  = vt->fYTrk;
      zt  = vt->fZTrk;
    }
    else {
      vid = -1;
      yt  = -1.e12;
      zt  = -1.e12;
    }

    Hist->fVaneID->Fill(vid);
    Hist->fYTrk->Fill  (yt);
    Hist->fZTrk->Fill  (zt);

    Hist->fDt->Fill(Track->Dt());
    Hist->fEp->Fill(Track->Ep());
    Hist->fDy->Fill(Track->Dy());
    Hist->fDz->Fill(Track->Dz());

    int ncl = Track->NClusters();
    Hist->fNClusters->Fill(ncl);
    
    Hist->fRSlope->Fill(Track->RSlope());
    Hist->fXSlope->Fill(Track->XSlope());

    Hist->fEleLogLHCal->Fill(Track->EleLogLHCal());
    Hist->fMuoLogLHCal->Fill(Track->MuoLogLHCal());

    Hist->fLogLHRCal ->Fill(Track->LogLHRCal ());
    Hist->fLogLHRDeDx->Fill(Track->LogLHRDeDx());
    Hist->fLogLHRXs  ->Fill(Track->LogLHRXs  ());
    Hist->fLogLHRTrk ->Fill(Track->LogLHRTrk ());

    Hist->fPdgCode->Fill(Track->fPdgCode);
    Hist->fFrGH->Fill(Track->fNGoodMcHits/(Track->fNActive+1.e-5));

  }

//-----------------------------------------------------------------------------
  void TCalm004::FillEventHistograms(EventHist_t* Hist) {

    double   cos_th, xv, yv, rv, zv;

    cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();

    xv = fEle->position().x()+3904.;
    yv = fEle->position().y();
    rv = sqrt(xv*xv+yv*yv);
    zv = fEle->position().z();

    Hist->fEleCosTh->Fill(cos_th);
    Hist->fRv->Fill(rv);
    Hist->fZv->Fill(zv);

    int ntrk_dem = fTrackBlock[0]->NTracks();
    int ntrk_uem = fTrackBlock[1]->NTracks();

    Hist->fNClusters->Fill(fNClusters);

    Hist->fNTracksDem->Fill(ntrk_dem);
    Hist->fNTracksUem->Fill(ntrk_uem);

    Hist->fNStrawHits[0]->Fill(fNStrawHits);
    Hist->fNStrawHits[1]->Fill(fNStrawHits);

    double emax   = -1;
    double t0_cls = -1;
    double dt     = 9999.;
    if (fCluster) {
      emax   = fCluster->Energy();
      t0_cls = fCluster->Time();
    }

    double t0_trk = -1;
    if (fTrack) {
      t0_trk = fTrack->GetKalRep()->t0().t0();
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
// this analysis code is intended to run on the obsolete dataset, so it can 
// contain various kludges
//-----------------------------------------------------------------------------
  void TCalm004::FillHistograms() {
//     double cos_th;
//     cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();
//-----------------------------------------------------------------------------
// EVT[0]: all events
// EVT[1]: events with upstream negative tracks
// EVT[2]: events without upstream negative tracks
//-----------------------------------------------------------------------------
    fCosmicFlag  = 0;

    int ntrk_uem = fTrackBlock[1]->NTracks();

    FillEventHistograms(fHist.fEvent[0]);

    if (ntrk_uem   > 0) {
      FillEventHistograms(fHist.fEvent[1]);
      fCosmicFlag = 1;
    }
    else                FillEventHistograms(fHist.fEvent[2]);
//-----------------------------------------------------------------------------
// track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
    TStnTrack*   trk;

    for (int i=0; i<fNTracks[0]; ++i ) {

      trk                = fTrackBlock[0]->Track(i);
      fElectronIDWord[i] = 0;

      if (fabs(trk->D0())  > 120.)                               fElectronIDWord[i] |= kTrackD0Bit;
      if ((trk->Ep() < 0) || (trk->Ep() > 1.2))                  fElectronIDWord[i] |= kClusterBit;
      if ((fabs(trk->Dy()) > 150.) || (fabs(trk->Dz()) > 150.))  fElectronIDWord[i] |= kClusterDxBit;
      if (trk->LogLHRCal() <  0  )                               fElectronIDWord[i] |= kLHRCalBit;
//-----------------------------------------------------------------------------
// normalization of the track-only likelihood is such that requiring llhr_trk > 0
// zero may not be the right cut, cut at -2 gives efficiency of ~92% for electrons
//-----------------------------------------------------------------------------
      if (trk->LogLHRTrk() <  -2.)                               fElectronIDWord[i] |= kLHRTrkBit;

      if ((fElectronIDWord[i] & kTrackD0Bit) == 1) {
//-----------------------------------------------------------------------------
// there is a cosmics-like track in the "event"
//----------------------------------------------------------------------------- 
	fCosmicFlag = 1;
      }

      FillTrackHistograms(fHist.fTrack[0],trk);

      if (trk->fIDWord == 0) {
					// track passes selection "C" 

	FillTrackHistograms(fHist.fTrack[1],trk);
	
	if (fNClusters > 0) {
	  FillTrackHistograms(fHist.fTrack[2],trk);
	}
	else {
	  FillTrackHistograms(fHist.fTrack[3],trk);
	}
//-----------------------------------------------------------------------------
// events with a good track, reconstructed clusters but without a match
//-----------------------------------------------------------------------------
	if ((fNClusters > 0) && (trk->NClusters() == 0)) {
	  FillTrackHistograms(fHist.fTrack[4],trk);
	}
//-----------------------------------------------------------------------------
// TRK 5 : events with a good track, reconstructed clusters but without a match
//-----------------------------------------------------------------------------
	if ((trk->fVMaxEp) && (fabs(trk->fVMaxEp->fDt) < 2.5)) {
	  FillTrackHistograms(fHist.fTrack[5],trk);
	}
      }
//-----------------------------------------------------------------------------
// if an event has an upstream negative track - reject it
//-----------------------------------------------------------------------------
      if (ntrk_uem == 0) {
	FillTrackHistograms(fHist.fTrack[20],trk);

	if (trk->fIDWord == 0) {
//-----------------------------------------------------------------------------
// TRK 21: track passes selection "C" - this is denominator "C"
//-----------------------------------------------------------------------------
	  FillTrackHistograms(fHist.fTrack[21],trk);

	  if ((fElectronIDWord[i] & kTrackD0Bit) == 0) {
//-----------------------------------------------------------------------------
// TRK 22: track passes selection "C" and has |D0| < 12cm
//-----------------------------------------------------------------------------
	    FillTrackHistograms(fHist.fTrack[22],trk);
//-----------------------------------------------------------------------------
// TRK 37 : events with a class "C" track, electron by TRK_ONLY_ID, 
// TRK 38 : add 100<P<110
//-----------------------------------------------------------------------------
	    if ((fElectronIDWord[i] & kLHRTrkBit) == 0) {
	      FillTrackHistograms(fHist.fTrack[37],trk);
	      if ((trk->P() >= 100.) && (trk->P() < 110)) { 
		FillTrackHistograms(fHist.fTrack[38],trk);
	      }
	    }
//-----------------------------------------------------------------------------
// TRK 23: track passes selection "C", has |D0| < 12cm, and 0 < E/P < 1.2
//-----------------------------------------------------------------------------
	    if ((fElectronIDWord[i] & kClusterBit) == 0) {
	      FillTrackHistograms(fHist.fTrack[23],trk);

	      if ((fElectronIDWord[i] & kClusterDxBit) == 0) {
//-----------------------------------------------------------------------------
// TRK 24: track passes selection "C", has |D0| < 12cm, 0 < E/P < 1.2, and |DY,Z| < 150mm
//-----------------------------------------------------------------------------
		FillTrackHistograms(fHist.fTrack[24],trk);
		if ((fElectronIDWord[i] & kLHRCalBit) == 0) {
//-----------------------------------------------------------------------------
// TRK 25 : events with a class "C" track, good cluster, electron candidates
//-----------------------------------------------------------------------------
		  FillTrackHistograms(fHist.fTrack[25],trk);
		  if ((trk->P() >= 100.) && (trk->P() < 110)) { 
//-----------------------------------------------------------------------------
// TRK 34 : events with a class "C" track, good cluster, electron candidates, 100<P<110
//-----------------------------------------------------------------------------
		    FillTrackHistograms(fHist.fTrack[34],trk);
		  }
		}
		else {
//-----------------------------------------------------------------------------
// TRK 26 : identified muon candidates
//-----------------------------------------------------------------------------
		  FillTrackHistograms(fHist.fTrack[26],trk);
		}
	      }
	      else {
//-----------------------------------------------------------------------------
// TRK 27 : events with a class "C" track and the cluster far from the track:
//          rely on the tracker-only based ID
//-----------------------------------------------------------------------------
		FillTrackHistograms(fHist.fTrack[27],trk);
//-----------------------------------------------------------------------------
// TRK 28 : events with class "C" tracks, with faraway cluster, identified as electrons TRK_ONLY
// TRK 35 : add 100 < P < 110 MeV
//-----------------------------------------------------------------------------
		if ((fElectronIDWord[i] & kLHRTrkBit) == 0) {
		  FillTrackHistograms(fHist.fTrack[28],trk);
		  if ((trk->P() >= 100.) && (trk->P() < 110)) { 
		    FillTrackHistograms(fHist.fTrack[35],trk);
		  }
		}
	      }
	    }
	    else {
//-----------------------------------------------------------------------------
// TRK 29 : events with class "C" track, |D0| < 120mm, and cluster EP < 0 or EP > 1.2:
//-----------------------------------------------------------------------------
	      FillTrackHistograms(fHist.fTrack[29],trk);
	      if (trk->Ep() < 0) {
//-----------------------------------------------------------------------------
// TRK 30 : events with class "C" track, |D0| < 120mm, and no cluster
//          for them one can only rely on the straw tracker based ID
//-----------------------------------------------------------------------------
		FillTrackHistograms(fHist.fTrack[30],trk);
//-----------------------------------------------------------------------------
// TRK 31 : events with class "C" track, no cluster, identified as electrons
// TRK 36 : add 100 < P < 110 MeV
//-----------------------------------------------------------------------------
		if ((fElectronIDWord[i] & kLHRTrkBit) == 0) {
		  FillTrackHistograms(fHist.fTrack[31],trk);
		  if ((trk->P() >= 100.) && (trk->P() < 110)) { 
		    FillTrackHistograms(fHist.fTrack[36],trk);
		  }
		}
	      }
	      else {
//-----------------------------------------------------------------------------
// TRK 32 : events with class "C" track, |D0| < 120mm, and E/P > 1.2 : 
//          reject as cosmics
//-----------------------------------------------------------------------------
		FillTrackHistograms(fHist.fTrack[32],trk);
	      }
	    }
	  }
	  else {
//-----------------------------------------------------------------------------
// TRK 33 : events with class "C" track and |D0| > 120mm: identified cosmics
//-----------------------------------------------------------------------------
	    FillTrackHistograms(fHist.fTrack[33],trk);
	  }
	}
//-----------------------------------------------------------------------------
// Set B cuts
// ----------
// TRK 41: events with the downstream track with |D0|<12cm
// TRK 42: events with the downstream track, |D0|<12cm and matched cluster
// TRK 43: events with the downstream track, |D0|<12cm and without matched cluster
// TRK 44: events with the downstream track, |D0|<12cm, matched cluster, LHR:electron
// TRK 45: events with the downstream track, |D0|<12cm, matched cluster, LHR:muons
//-----------------------------------------------------------------------------
	if (fTrackIDWordSetB[i] == 0) {
//-----------------------------------------------------------------------------
// TRK 41: track passes selection "B" - this is denominator "B"
//-----------------------------------------------------------------------------
	  FillTrackHistograms(fHist.fTrack[41],trk);

	  if ((fElectronIDWord[i] & kTrackD0Bit) == 0) {
	    FillTrackHistograms(fHist.fTrack[42],trk);
//-----------------------------------------------------------------------------
// TRK 57 : events with a class "B" track, electron by TRK_ONLY_ID, 
// TRK 58 : add 100<P<110
//-----------------------------------------------------------------------------
	    if ((fElectronIDWord[i] & kLHRTrkBit) == 0) {
	      FillTrackHistograms(fHist.fTrack[57],trk);
	      if ((trk->P() >= 100.) && (trk->P() < 110)) { 
		FillTrackHistograms(fHist.fTrack[58],trk);
	      }
	    }
	    if ((fElectronIDWord[i] & kClusterBit) == 0) {
	      FillTrackHistograms(fHist.fTrack[43],trk);

	      if ((fElectronIDWord[i] & kClusterDxBit) == 0) {
		FillTrackHistograms(fHist.fTrack[44],trk);
//-----------------------------------------------------------------------------
// TRK 45 : events with a class "B" track, good cluster, electron candidates
//-----------------------------------------------------------------------------
		if ((fElectronIDWord[i] & kLHRCalBit) == 0) {
		  FillTrackHistograms(fHist.fTrack[45],trk);
		  if ((trk->P() >= 100.) && (trk->P() < 110)) { 
//-----------------------------------------------------------------------------
// TRK 54 : events with a class "B" track, good cluster, electron candidates, 100<P<110
//-----------------------------------------------------------------------------
		    FillTrackHistograms(fHist.fTrack[54],trk);
		  }
		}
		else {
//-----------------------------------------------------------------------------
// TRK 46 : identified muon candidates
//-----------------------------------------------------------------------------
		  FillTrackHistograms(fHist.fTrack[46],trk);
		}
	      }
	      else {
//-----------------------------------------------------------------------------
// TRK 47 : events with a class "B" track and the cluster far from the track:
//          rely on the tracker-only based ID
//-----------------------------------------------------------------------------
		FillTrackHistograms(fHist.fTrack[47],trk);
//-----------------------------------------------------------------------------
// TRK 48 : events with class "B" tracks, with away cluster, identified as electrons
// TRK 55 : add 100 < P < 110 MeV
//-----------------------------------------------------------------------------
		if ((fElectronIDWord[i] & kLHRTrkBit) == 0) {
		  FillTrackHistograms(fHist.fTrack[48],trk);
		  if ((trk->P() >= 100.) && (trk->P() < 110)) { 
		    FillTrackHistograms(fHist.fTrack[55],trk);
		  }
		}
	      }
	    }
	    else {
//-----------------------------------------------------------------------------
// TRK 49 : events with class "B" track, |D0| < 120mm, and  E/P < 0 or E/P > 1.2
//-----------------------------------------------------------------------------
	      FillTrackHistograms(fHist.fTrack[49],trk);
//-----------------------------------------------------------------------------
// TRK 50 : events with class "B" track, |D0| < 120mm, and no cluster:
//          for them one can only rely on the straw tracker based ID
//-----------------------------------------------------------------------------
	      if (trk->Ep() < 0) {
		FillTrackHistograms(fHist.fTrack[50],trk);
//-----------------------------------------------------------------------------
// TRK 51 : events with class "B" track, no cluster, identified as electrons
// TRK 56 : add 100 < P < 110 MeV
//-----------------------------------------------------------------------------
		if ((fElectronIDWord[i] & kLHRTrkBit) == 0) {
		  FillTrackHistograms(fHist.fTrack[51],trk);
		  if ((trk->P() >= 100.) && (trk->P() < 110)) { 
		    FillTrackHistograms(fHist.fTrack[56],trk);
		  }
		}
	      }
	      else {
//-----------------------------------------------------------------------------
// TRK 52 : events with class "B" track, |D0| < 120mm, and E/P > 1.2 : cosmics
//-----------------------------------------------------------------------------
		FillTrackHistograms(fHist.fTrack[52],trk);
	      }
	    }
	  }
	  else {
//-----------------------------------------------------------------------------
//  TRK 53 : tracks with large impact parameter - cosmics
//-----------------------------------------------------------------------------
	    FillTrackHistograms(fHist.fTrack[53],trk);
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// track histograms, for fTrackBlock[1] - upstream e- hypothesis
//-----------------------------------------------------------------------------
    for (int i=0; i<fNTracks[1]; ++i ) {
      trk = fTrackBlock[1]->Track(i);
      FillTrackHistograms(fHist.fTrack[100],trk);

      if (trk->fIDWord == 0) {
					// track passes selection "C" 

	FillTrackHistograms(fHist.fTrack[101],trk);
	
	if (fNClusters > 0) {
	  FillTrackHistograms(fHist.fTrack[102],trk);
	}
	else {
	  FillTrackHistograms(fHist.fTrack[103],trk);
	}
//-----------------------------------------------------------------------------
// events with a good track, reconstructed clusters but without a match
//-----------------------------------------------------------------------------
	if ((fNClusters > 0) && (trk->NClusters() == 0)) {
	  FillTrackHistograms(fHist.fTrack[104],trk);
	}
//-----------------------------------------------------------------------------
// TRK 5 : events with a good track, reconstructed clusters but without a match
//-----------------------------------------------------------------------------
	if ((trk->fVMaxEp) && (fabs(trk->fVMaxEp->fDt) < 2.5)) {
	  FillTrackHistograms(fHist.fTrack[105],trk);
	}
      }
    }

//-----------------------------------------------------------------------------
// cluster histograms
//-----------------------------------------------------------------------------
    TStnCluster* cl;
    for (int i=0; i<fNClusters; ++i ) {
      cl = fClusterBlock->Cluster(i);

      FillClusterHistograms(fHist.fCluster[0],cl);

      if (fNTracks        > 0) FillClusterHistograms(fHist.fCluster[1],cl);
      if (fNGoodTracks    > 0) FillClusterHistograms(fHist.fCluster[2],cl);
    }
  }

//-----------------------------------------------------------------------------
  void TCalm004::getData(const art::Event& Evt) {
//-----------------------------------------------------------------------------
// generated particles
//-----------------------------------------------------------------------------
    art::Handle<GenParticleCollection> genpHandle;
    Evt.getByLabel("generate",genpHandle);
    fListOfGenParticles = (GenParticleCollection*) &(*genpHandle);
    fEle                = (GenParticle*) &fListOfGenParticles->at(0);

    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    Evt.getByLabel(fStrawHitMaker,"StrawHitMCPtr",mcptrHandle);
    fListOfMcStrawHits = (PtrStepPointMCVectorCollection*) &(*mcptrHandle);
//-----------------------------------------------------------------------------
// straw hits 
//-----------------------------------------------------------------------------
    art::Handle<StrawHitCollection> shHandle;
    Evt.getByLabel("makeSH",shHandle);
    fListOfStrawHits = (StrawHitCollection*) &(*shHandle);
    fNStrawHits      = fListOfStrawHits->size();
  }



//-----------------------------------------------------------------------------
  void TCalm004::Init(art::Event& Evt) {
//    TStnCluster*    cluster;
    int             id_word; // , ntrk;
    TStnTrack*      track;

    TEmuLogLH::PidData_t  cdat;
    double                xs;

//-----------------------------------------------------------------------------
// initialize tracks and determine track quality
//-----------------------------------------------------------------------------
    StntupleInitMu2eTrackBlock(fTrackBlock[0],&Evt,0);

    fNGoodTracks    = 0;
    fNMatchedTracks = 0;
//-----------------------------------------------------------------------------
// determine number of DEM tracks matching to clusters in time
//-----------------------------------------------------------------------------
    fNTracks[0] = fTrackBlock[0]->NTracks();
    for (int i=0; i<fNTracks[0]; i++) {
      track               = fTrackBlock[0]->Track(i);
      fTrackIDWordSetB[i] = fTrackIDSetB->IDWord(track);
      id_word             = fTrackIDSetC->IDWord(track);
      track->fIDWord      = id_word;

      if (id_word == 0) {
	fNGoodTracks += 1;
	if ((track->fVMaxEp) && (fabs(track->fVMaxEp->fDt) < 2.5)) {
	  fNMatchedTracks += 1;
	}
      }

      cdat.fDt   = track->Dt();
      cdat.fEp   = track->Ep();
      cdat.fPath = -1.;
      if (track->fVMinS) cdat.fPath = track->fVMinS->fPath;
      
      track->fEleLogLHCal = fLogLH->LogLHCal(&cdat,11);
      track->fMuoLogLHCal = fLogLH->LogLHCal(&cdat,13);

      xs = track->XSlope  ();
      
      track->fLogLHRXs = fLogLH->LogLHRXs(xs);
    }
    if (fNTracks[0] > 0) fTrack = fTrackBlock[0]->Track(0);
//-----------------------------------------------------------------------------
// number of upstream negative tracks
//-----------------------------------------------------------------------------
    fNTracks[1] = fTrackBlock[1]->NTracks();
//-----------------------------------------------------------------------------
// initialize clusters
//-----------------------------------------------------------------------------
    StntupleInitMu2eClusterBlock(fClusterBlock,&Evt,0);
    fNClusters = fClusterBlock->NClusters();
    if (fNClusters > 0) fCluster = fClusterBlock->Cluster(0);
    else                fCluster = 0;
//-----------------------------------------------------------------------------
// then - clusters, tracks are supposed to be already initialized
//-----------------------------------------------------------------------------
//     const mu2e::CaloCluster       *cl;
//     unsigned int nm(0);

//     if (fTrkCalMap) nm = fTrkCalMap->size();

//     for (int i=0; i<fNClusters; i++) {
//       cluster               = fClusterBlock->NewCluster();

//       for(size_t i=0; i<nm; i++) {
// 	//	KalRepPtr const& trkPtr = fTrkCalMap->at(i).first->trk();
// 	//	const KalRep *  const &trk = *trkPtr;

// 	cl = &(*(fTrkCalMap->at(i).second));

// 	if (cl == cluster->fCaloCluster) { 
// 	  cluster->fClosestTrack = fTrackBlock->Track(i);
// 	  break;
// 	}
//       }
//     }
  }


//-----------------------------------------------------------------------------
  void TCalm004::Debug(const art::Event& Evt) {
    const char* oname = "TCalm004::Debug";
    double      pt;
    TStnTrack*  t;
    int         nt, nt1;

//-----------------------------------------------------------------------------
// bit 000: 
//-----------------------------------------------------------------------------
    if (TModule::fDebugBit[0] != 0) {
      nt = fTrackBlock[0]->NTracks();
      for (int i=0; i<nt; i++) {
	t = fTrackBlock[0]->Track(i);
	pt = t->Momentum()->Pt();
	if ((pt > 80.) && (t->fNActive > 35) && (t->NClusters() == 0)) {
	  printf(" >>>>>>> [%s] EVENT : %10i ERROR:000 TRACK Pt = %10.3f doesnt have a cluster\n",
		 oname,Evt.event(),pt);
	  //	  t->Print("");
	}
      }
    }
    
    if (TModule::fDebugBit[11] != 0) {
      nt  = fTrackBlock[0]->NTracks();
      nt1 = fTrackBlock[1]->NTracks();
      if (nt1 == 0) {
	for (int i=0; i<nt; i++) {
	  t = fTrackBlock[0]->Track(i);
	  if ((t->fIDWord == 0) && (fElectronIDWord[i] == 0) && (t->fPdgCode != 11)) {
	    printf(" >>>>>>> [%s] RUN:EVENT: %6i:%10i BIT:011 track ele ID=0, but not electron\n",
		   oname,Evt.run(),Evt.event());
	    //	  t->Print("");
	    fPassed = 1;
	  }
	}
      }
    }

    if (TModule::fDebugBit[12] != 0) {
      nt  = fTrackBlock[0]->NTracks();
      for (int i=0; i<nt; i++) {
	t = fTrackBlock[0]->Track(i);
	printf(" >>>>>>> [%s] RUN:EVENT: %6i:%10i BIT:012 track ele ID=0x%08x E/P = %10.3f DZ = %10.3f\n",
	       oname,Evt.run(),Evt.event(), t->fIDWord,t->Ep(), t->Dz());
      }
    }
  }

//-----------------------------------------------------------------------------
  bool TCalm004::filter(art::Event& Evt) {
    const char* oname = "TCalm004::analyze";

    bool rc(true);

    fPassed = 0;

    printf(" >>>>>>> [%s] Run:EVENT : %10i:%10i\n",oname,Evt.run(),Evt.event());

    getData(Evt);

    Init(Evt);

    FillHistograms();

    Debug(Evt);

    if (TModule::fDebugBit[11] != 0) {
      rc = fPassed;
    }

    TModule::filter(Evt);

    return rc;
  }

}  // end namespace mu2e

using mu2e::TCalm004;

DEFINE_ART_MODULE(TCalm004);
