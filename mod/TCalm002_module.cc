//-----------------------------------------------------------------------------
// General-purpose analysis module, fills multiple histograms
//
// Debug Printout:
// ---------------
// bit   0: tracks with IDWORD==0 and no cluster
// bit   1: tracks with IDWORD==0 and no cluster with E/P > 0.4
// bit 003: print downstream - tracks
// bit 004: events with a track with the fit consistency < 1e-5
// bit 005: for each event print NTracks, NClusters, E(max)
//
// filter:
// -------
// bit 51: n(matched tracks > 0
// bit 52: N(calo crystal hits > 0) (for DIO)
//-----------------------------------------------------------------------------
#include "TEnv.h"
#include "TSystem.h"

#include "murat/mod/TCalm002_module.hh"
#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/obj/TStnCrystal.hh"

// #include "murat/gui/TEvdCrystal.hh"

#include "Stntuple/val/stntuple_val_functions.hh"
#include "Stntuple/mod/TAnaDump.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "BTrk/KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalRepPtrCollection.hh"
#include "TrackCaloMatching/inc/TrkToCaloExtrapol.hh"

#include "Stntuple/obj/TDiskCalorimeter.hh"
#include "Stntuple/obj/TDisk.hh"
#include "Stntuple/base/TObjHandle.hh"

namespace mu2e {

//-----------------------------------------------------------------------------
  TCalm002::TCalm002(fhicl::ParameterSet const& pset): 
    THistModule(pset,"TCalm002"),
    fHistFileName   (pset.get<std::string> ("histFileName"   , "tcalm002.hist")),
    fG4ModuleLabel  (pset.get<std::string> ("g4ModuleLabel"  , "g4run"        )),
    fStrawHitMaker  (pset.get<std::string> ("strawHitMaker"  , "makeSH"       )),
    fTrkPatRecDem   (pset.get<std::string> ("trkPatRecDem"   , "trkPatRecDem" )),
    fTrkPatRecUem   (pset.get<std::string> ("trkPatRecUem"   , "trkPatRecUem" )),
    fCrystalHitMaker(pset.get<std::string> ("crystalHitMaker", "CaloCrystalHitsMaker")),
    fCaloClusterMaker(pset.get<std::string> ("caloClusterMaker", "makeCaloCluster")),
    fTrkExtrapol    (pset.get<std::string> ("trkExtrapol"    , "trkExtrapol"  )),
    fTrkCalMatch    (pset.get<std::string> ("trkCalMatch"    , "caloMatching" )),
    fPidDem         (pset.get<std::string> ("pidDem"         , "undefined"    )),
    fMinTActive     (pset.get<double>      ("minTActive"     ,   710.         )),
    fMinECrystal    (pset.get<double>      ("minECrystal"    ,    0.1         )),
    fMinCrystalFr   (pset.get<double>      ("minCrystalFr"   ,    1.0         ))
  {
					// reset all histogram pointers

    for (int i=0; i<kNEventHistSets; i++) {
      fHist.fEvent[i] = 0;
    }

    fCalDataBlock = new TCalDataBlock   ();
    fClusterBlock = new TStnClusterBlock();
    fTrackBlock   = new TStnTrackBlock  ();

    fTrackID      = new TStnTrackID("SetC");
    fLogLH        = new TEmuLogLH();

    fCalorimeterType = 0; // no default - make it undefined
  }



//-----------------------------------------------------------------------------
  TCalm002::~TCalm002() {
//     fFolder->Delete();
//     delete fFolder;

    delete fTrackBlock;
    delete fClusterBlock;
    delete fTrackID;
    delete fLogLH;
  }

//-----------------------------------------------------------------------------
  void TCalm002::BookCaloHistograms(CaloHist_t* Hist, const char* Folder) {
//     char name [200];
//     char title[200];
//-----------------------------------------------------------------------------
//  
//-----------------------------------------------------------------------------
    HBook1F(Hist->fVaneID ,"vane_id",Form("%s: Vane ID"       ,Folder), 10, 0,  10,Folder);

    for (int i=0; i<4; i++) {
      HBook1F(Hist->fEnergy  [i],Form("energy_%i",i),Form("%s: Hit Energy[%i]",Folder,i),200, 0, 100,Folder);
      HBook1F(Hist->fTime    [i],Form("time_%i"  ,i),Form("%s: Hit time  [%i]",Folder,i),200, 0,2000,Folder);
      HBook1F(Hist->fNHits   [i],Form("nhits_%i" ,i),Form("%s: NHits     [%i]",Folder,i), 50, 0,  50,Folder);
      HBook1F(Hist->fRadius  [i],Form("r_%i"     ,i),Form("%s: Radius    [%i]",Folder,i),100, 0,1000,Folder);
      HBook1F(Hist->fRadiusWE[i],Form("rwe_%i"   ,i),Form("%s: RadiusWE  [%i]",Folder,i),100, 0,1000,Folder);

      HBook1F(Hist->fE700    [i],Form("e700_%i",i),Form("%s: Hit Energy[%i] (T > 700ns)",Folder,i),200, 0, 100,Folder);
      HBook1F(Hist->fT700    [i],Form("t700_%i",i),Form("%s: Hit time  [%i] (T > 700ns)",Folder,i),200, 0,2000,Folder);
      HBook1F(Hist->fN700    [i],Form("n700_%i",i),Form("%s: NHits     [%i] (T > 700ns)",Folder,i), 50, 0,  50,Folder);

      HBook1F(Hist->fR700  [i],Form("r700_%i"  ,i),Form("%s: Radius (T>700) [%i]",Folder,i),100, 0,1000,Folder);
      HBook1F(Hist->fRWE700[i],Form("rwe700_%i",i),Form("%s: Radius*E(T>700)[%i]",Folder,i),100, 0,1000,Folder);
    }
  }

//-----------------------------------------------------------------------------
  void TCalm002::BookClusterHistograms(ClusterHist_t* Hist, const char* Folder) {
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
    HBook1F(Hist->fZ      ,"z"      ,Form("%s: cluster Z"     ,Folder),200, 11500,13500,Folder);
    HBook1F(Hist->fR      ,"r"      ,Form("%s: cluster Radius",Folder),100, 0,  1000,Folder);
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
  void TCalm002::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
//     char name [200];
//     char title[200];
//-----------------------------------------------------------------------------
//  
//-----------------------------------------------------------------------------
    HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
    HBook1F(Hist->fRv         ,"rv"       ,Form("%s: R(Vertex)"                      ,Folder), 100, 0, 1000,Folder);
    HBook1F(Hist->fZv         ,"zv"       ,Form("%s: Z(Vertex)"                      ,Folder), 300, 0,15000,Folder);
    HBook1F(Hist->fNClusters ,"ncl"      ,Form("%s: Number of Reconstructed Clusters",Folder),200,0,200,Folder);
    HBook1F(Hist->fNTracks   ,"ntrk"     ,Form("%s: Number of Reconstructed Tracks"  ,Folder),100,0,100,Folder);
    HBook1F(Hist->fNStrawHits[0],"nsh_0" ,Form("%s: Number of Straw Hits [0]"        ,Folder),250,0,250,Folder);
    HBook1F(Hist->fNStrawHits[1],"nsh_1" ,Form("%s: Number of Straw Hits [1]"        ,Folder),250,0,5000,Folder);
    HBook1F(Hist->fNGoodSH   ,"nsh50"    ,Form("%s: N(SH) +/-50"                     ,Folder),300,0,1500,Folder);
    HBook1F(Hist->fDtClT     ,"dt_clt"   ,Form("%s: DT(cluster-track)"               ,Folder),100,-100,100,Folder);
    HBook1F(Hist->fDtClS     ,"dt_cls"   ,Form("%s: DT(cluster-straw hit)"           ,Folder),200,-200,200,Folder);
    HBook1F(Hist->fSHTime    ,"shtime"   ,Form("%s: Straw Hit Time"                  ,Folder),400,0,2000,Folder);
    HBook1F(Hist->fEMax      ,"emax"     ,Form("%s: Max cluster energy"              ,Folder),150,0,150,Folder);
    HBook1F(Hist->fNHyp      ,"nhyp"     ,Form("%s: N(fit hypotheses)"               ,Folder),5,0,5,Folder);
    HBook1F(Hist->fBestHyp[0],"bfh0"     ,Form("%s: Best Fit Hyp[0](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
    HBook1F(Hist->fBestHyp[1],"bfh1"     ,Form("%s: Best Fit Hyp[1](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
    char  name[200];
    for (int i=0; i<2; i++) {
      sprintf(name,"ncch_%i",i);
      HBook1F(Hist->fNCaloCrystalHits[i],name,Form("%s: N(calo crystal hits) [%i]",Folder,i),500,0,1000,Folder);
      sprintf(name,"ncch_vs_vane_%i",i);
      HBook2F(Hist->fNCaloHitsVsVane[i],name,Form("%s: N(calo crystal hits) vs vane[%i]",Folder,i),4,0,4,200,0,200,Folder);
      sprintf(name,"ncch_vs_row_%i",i);
      HBook2F(Hist->fNCaloHitsVsRow[i],name,Form("%s: N(calo crystal hits) vs row [%i]",Folder,i),20,0,20,200,0,200,Folder);
      sprintf(name,"ncch_vs_col_%i",i);
      HBook2F(Hist->fNCaloHitsVsCol[i],name,Form("%s: N(calo crystal hits) vs col [%i]",Folder,i),50,0,50,200,0,200,Folder);
    }

    for (int i=0; i<4; i++) {
      HBook1F(Hist->fETot        [i],Form("etot_%i"    ,i),Form("%s: Etot[%i]",Folder,i), 300, 0,150,Folder);
      HBook2F(Hist->fECrVsR      [i],Form("ecr_vs_r_%i",i),Form("%s: E Cr Vs R [%i]"    ,Folder,i), 100, 0,1000,500,0,100,Folder);
      HBook2F(Hist->fNCrVsR      [i],Form("ncr_vs_r_%i",i),Form("%s: N Cr Vs R [%i]"    ,Folder,i), 100, 0,1000,100,0,100,Folder);

      HBook2F(Hist->fNCrystalHitsVsR[i],Form("ncrh_vs_r_%i",i),Form("%s: N Crystal Hits[%i] vs R",Folder,i), 100, 0, 1000,100,0,100,Folder);
      HBook2F(Hist->fNHitCrystalsVsR[i],Form("nhcr_vs_r_%i",i),Form("%s: N Hit Crystals[%i] vs R",Folder,i), 100, 0, 1000,100,0,100,Folder);
    }

    HBook1F(Hist->fNHitCrystalsTot,"nhcr_tot",Form("%s: NHit Crystals Tot",Folder), 100, 0,100,Folder);
  }


//-----------------------------------------------------------------------------
  void TCalm002::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
//     char name [200];
//     char title[200];
//-----------------------------------------------------------------------------
//  
//-----------------------------------------------------------------------------
    HBook1F(Hist->fP[0]       ,"p"        ,Form("%s: Track P(total)[0]" ,Folder),1000,  90  ,110. ,Folder);
    HBook1F(Hist->fP[1]       ,"p1"       ,Form("%s: Track P(total)[1]" ,Folder),1000, 104.5,105.5,Folder);
    HBook1F(Hist->fP[2]       ,"p2"       ,Form("%s: Track P(total)[1]" ,Folder),1000,   0  ,200. ,Folder);
    HBook1F(Hist->fPFront     ,"pf"       ,Form("%s: Track P(front)   " ,Folder),1000,  90  ,110. ,Folder);
    HBook1F(Hist->fDpFront    ,"dpf"      ,Form("%s: Track P-P(front) " ,Folder), 200,  -5. ,  5. ,Folder);
    HBook1F(Hist->fPStOut     ,"pstout"   ,Form("%s: Track P(ST_Out)  " ,Folder),1000,  90  ,110. ,Folder);
    HBook1F(Hist->fDpFSt      ,"dpfst"    ,Form("%s: Track Pf-Psto"     ,Folder), 200,  -5  ,  5. ,Folder);
    HBook1F(Hist->fPt         ,"pt"       ,Form("%s: Track Pt"          ,Folder), 600, 75,95,Folder);
    HBook1F(Hist->fCosTh      ,"costh"    ,Form("%s: Track cos(theta)"  ,Folder), 100,-1,1,Folder);
    HBook1F(Hist->fChi2       ,"chi2"     ,Form("%s: Track chi2 total"  ,Folder), 200, 0,200,Folder);
    HBook1F(Hist->fNDof       ,"ndof"     ,Form("%s: Number of DOF"     ,Folder), 200, 0,200,Folder);
    HBook1F(Hist->fChi2Dof    ,"chi2d"    ,Form("%s: track chi2/N(dof)" ,Folder), 500, 0, 10,Folder);
    HBook1F(Hist->fChi2DofC   ,"chi2dc"   ,Form("%s: track chi2/N calc" ,Folder), 500, 0, 10,Folder);
    HBook1F(Hist->fNActive    ,"nactv"    ,Form("%s: N(active)"         ,Folder), 200, 0,200,Folder);
    HBook1F(Hist->fT0         ,"t0"       ,Form("%s: track T0"          ,Folder), 200, 0,2000,Folder);
    HBook1F(Hist->fQ          ,"q"        ,Form("%s: track Q"           ,Folder),   4,-2,   2,Folder);
    HBook1F(Hist->fFitCons[0] ,"fcon"     ,Form("%s: track fit cons [0]",Folder), 200, 0,   1,Folder);
    HBook1F(Hist->fFitCons[1] ,"fcon1"    ,Form("%s: track fit cons [1]",Folder), 1000, 0,   0.1,Folder);
    HBook1F(Hist->fD0         ,"d0"       ,Form("%s: track D0      "    ,Folder), 200,-20, 20,Folder);
    HBook1F(Hist->fZ0         ,"z0"       ,Form("%s: track Z0      "    ,Folder), 200,-20000,20000,Folder);
    HBook1F(Hist->fTanDip     ,"tdip"     ,Form("%s: track tan(dip)"    ,Folder), 100,-2.5 ,2.5,Folder);
    HBook1F(Hist->fResid      ,"resid"    ,Form("%s: hit residuals"     ,Folder), 500,-0.5 ,0.5,Folder);
    HBook1F(Hist->fDt         ,"dt"       ,Form("%s: track delta(T)"    ,Folder), 200,-20  ,20 ,Folder);
    HBook1F(Hist->fDx         ,"dx"       ,Form("%s: track delta(X)"    ,Folder), 100,-500 ,500,Folder);
    HBook1F(Hist->fDy         ,"dy"       ,Form("%s: track delta(Y)"    ,Folder), 100,-500 ,500,Folder);
    HBook1F(Hist->fDz         ,"dz"       ,Form("%s: track delta(Z)"    ,Folder), 100,-500 ,500,Folder);
    HBook1F(Hist->fNClusters  ,"ncl"      ,Form("%s: track N(clusters)" ,Folder),  10, 0   , 10,Folder);
    HBook1F(Hist->fVaneID     ,"vid"      ,Form("%s: track vane ID"     ,Folder),  10,-5   ,  5,Folder);
    HBook1F(Hist->fXCal       ,"xcal"     ,Form("%s: track XCal"        ,Folder), 200,-1000,1000,Folder);
    HBook1F(Hist->fYCal       ,"ycal"     ,Form("%s: track YCal"        ,Folder), 200,-1000,1000,Folder);
    HBook1F(Hist->fZCal       ,"zcal"     ,Form("%s: track ZCal"        ,Folder), 200, 1500,3500,Folder);
    HBook1F(Hist->fXTrk       ,"xtrk"     ,Form("%s: track XTrk"        ,Folder), 200,-1000,1000,Folder);
    HBook1F(Hist->fYTrk       ,"ytrk"     ,Form("%s: track YTrk"        ,Folder), 200,-1000,1000,Folder);
    HBook1F(Hist->fRTrk       ,"rtrk"     ,Form("%s: track RTrk"        ,Folder), 200,-1000,1000,Folder);
    HBook1F(Hist->fZTrk       ,"ztrk"     ,Form("%s: track ZTrk"        ,Folder), 200,-1000,1000,Folder);
    HBook1F(Hist->fEp         ,"ep"       ,Form("%s: track E/P"         ,Folder), 300, 0   ,1.5,Folder);
    HBook2F(Hist->fNHVsStation,"nh_vs_st" ,Form("%s: N(hits) Vs Station",Folder),  40, 0,40,10,-0.5,9.5,Folder);
    HBook2F(Hist->fNHVsNSt    ,"nh_vs_nst",Form("%s: N(hits) Vs NSt"    ,Folder),  10,-0.5,9.5,40,-0.5,39.5,Folder);

    HBook1F(Hist->fRSlope     ,"rslope"   ,Form("%s: Res Slope"         ,Folder), 200,-20 , 20,Folder);
    HBook1F(Hist->fXSlope     ,"xslope"   ,Form("%s: Res/Sig Slope"     ,Folder), 200,-20 , 20,Folder);

    HBook1F(Hist->fEleLogLHCal,"ele_llh_c",Form("%s: ELE Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
    HBook1F(Hist->fMuoLogLHCal,"muo_llh_c",Form("%s: MUO Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
    HBook1F(Hist->fLogLHRCal  ,"llhr_cal" ,Form("%s: LogLH(e/m) Cal"    ,Folder), 200,-100,100,Folder);
    HBook1F(Hist->fLogLHRDeDx ,"llhr_dedx",Form("%s: LogLH(e/m) De/Dx"  ,Folder), 200,-20 , 20,Folder);
    HBook1F(Hist->fLogLHRXs   ,"llhr_xs"  ,Form("%s: LogLH(e/m) XSlope" ,Folder), 200,-20 , 20,Folder);
    HBook1F(Hist->fLogLHRTrk  ,"llhr_trk" ,Form("%s: LogLH(e/m) Trk"    ,Folder), 200,-20 , 20,Folder);

    HBook1F(Hist->fPdgCode    ,"pdg"      ,Form("%s: track PDG code"    ,Folder), 100,-50,50,Folder);
    HBook1F(Hist->fFrGH       ,"fgh"      ,Form("%s: Fraction Goog Hits",Folder), 100, 0,1,Folder);
  }

//-----------------------------------------------------------------------------
  void TCalm002::BookHistograms() {
    TFolder  *hist_folder, *fol;
    char     folder_name[200];

    hist_folder = (TFolder*) fFolder->FindObject("Hist");
//-----------------------------------------------------------------------------
// book crystal histograms
//-----------------------------------------------------------------------------
    HBook1F(fHist.fCrystalR[0],"rc_0"     ,Form("disk [0] crystal radius"),100,0,1000,"Hist");
    HBook1F(fHist.fCrystalR[1],"rc_1"     ,Form("disk [1] crystal radius"),100,0,1000,"Hist");

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
    book_event_histset[6] = 1;	        // events w/o tracks passing "Set C" cuts
    book_event_histset[7] = 1;	        // events with E(cluster) > 70 MeV
    book_event_histset[8] = 1;	        // events with the highest energy cluster on the 1st disk
    book_event_histset[9] = 1;	        // events with teh highest energy cluster on the 2nd diskk

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
    book_track_histset[  6] = 1;		// all tracks e-  E/P > 0.4
    book_track_histset[  7] = 1;		// all tracks e- passing Set C cuts, E/P > 0.4
    book_track_histset[  8] = 1;		// Set C tracks e- , |xslope| < 3.
    book_track_histset[  9] = 1;		// all tracks in the event where there aren't EM clusters with energy greater than 70 MeV
    book_track_histset[ 10] = 1;		// good tracks in the event where there aren't EM clusters with energy greater than 70 MeV
    book_track_histset[ 11] = 1;		// all tracks with P > 103.5
    book_track_histset[ 12] = 1;		// tracks with fcons < 1.e-4
    book_track_histset[ 13] = 1;		// "Set C" tracks with 100 <= P < 110 
    book_track_histset[ 14] = 1;		// tracks with fcons < 1.e-2
    book_track_histset[ 15] = 1;		// tracks intersecting the 1st disk
    book_track_histset[ 16] = 1;		// tracks intersecting the 2nd disk
    book_track_histset[ 17] = 1;		// tracks with no calorimeter intersections

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
    book_cluster_histset[4] = 1;		// clusters > 10 MeV
    book_cluster_histset[5] = 1;		// clusters > 70 MeV

    for (int i=0; i<kNClusterHistSets; i++) {
      if (book_cluster_histset[i] != 0) {
	sprintf(folder_name,"cls_%i",i);
	fol = (TFolder*) hist_folder->FindObject(folder_name);
	if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
	fHist.fCluster[i] = new ClusterHist_t;
	BookClusterHistograms(fHist.fCluster[i],Form("Hist/%s",folder_name));
      }
    }
//-----------------------------------------------------------------------------
// book calorimeter histograms
//-----------------------------------------------------------------------------
    int book_calo_histset[kNCaloHistSets];
    for (int i=0; i<kNCaloHistSets; i++) book_calo_histset[i] = 0;

    book_calo_histset[0] = 1;		// all crystals
    book_calo_histset[1] = 1;		// all crystals, e > 0
    book_calo_histset[2] = 1;		// all crystals, e > 0.1
    book_calo_histset[3] = 1;		// all crystals, e > 1.0

    for (int i=0; i<kNCaloHistSets; i++) {
      if (book_calo_histset[i] != 0) {
	sprintf(folder_name,"cal_%i",i);
	fol = (TFolder*) hist_folder->FindObject(folder_name);
	if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
	fHist.fCalo[i] = new CaloHist_t;
	BookCaloHistograms(fHist.fCalo[i],Form("Hist/%s",folder_name));
      }
    }
  }

//-----------------------------------------------------------------------------
// begin job - book histograms etc
//-----------------------------------------------------------------------------
  void TCalm002::beginJob() {

    //    TFolder* fol;
    //    TFolder* hist_folder;
    
    //    char folder_name[200];

    BookHistograms();

    fTrackID->SetMinT0(fMinTActive);
//-----------------------------------------------------------------------------
// define collection names to be used for initialization
//-----------------------------------------------------------------------------
    fClusterBlock->AddCollName("mu2e::CaloClusterCollection"       ,fCaloClusterMaker.data(),"");

    fTrackBlock->AddCollName("mu2e::KalRepCollection"              ,fTrkPatRecDem.data()    ,"DownstreameMinus");
    fTrackBlock->AddCollName("mu2e::KalRepCollection"              ,fTrkPatRecUem.data()    ,"UpstreameMinus");
    fTrackBlock->AddCollName("mu2e::CaloClusterCollection"         ,fCaloClusterMaker.data(),"AlgoCLOSESTSeededByENERGY");
    fTrackBlock->AddCollName("mu2e::TrkToCaloExtrapolCollection"   ,fTrkExtrapol.data()     ,"");
    //    fTrackBlock->AddCollName("mu2e::TrackClusterLink"              ,fTrkCalMatch.data()     ,"");
    fTrackBlock->AddCollName("mu2e::TrackClusterMatchCollection"   ,fTrkCalMatch.data()     ,"");
    fTrackBlock->AddCollName("mu2e::StrawHitCollection"            ,fStrawHitMaker.data()   ,"");
    fTrackBlock->AddCollName("mu2e::PtrStepPointMCVectorCollection",fStrawHitMaker.data()   ,"StrawHitMCPtr");
    fTrackBlock->AddCollName("mu2e::PIDProductCollection"          ,fPidDem.data()          ,"");
    fTrackBlock->AddCollName("mu2e::StepPointMCCollection"         ,fG4ModuleLabel.data()   ,"");

    fCalDataBlock->AddCollName("mu2e::CaloCrystalHitCollection"    ,fCrystalHitMaker.data(),"");

    TAnaDump::Instance()->AddObject("TCalm002::TrackBlock"  ,fTrackBlock  );
    TAnaDump::Instance()->AddObject("TCalm002::ClusterBlock",fClusterBlock);
//-----------------------------------------------------------------------------
// initialize the likelihood histograms
//-----------------------------------------------------------------------------
    const char*  mu2e_hist_dir;
    char         fn[200];

    mu2e_hist_dir = gEnv->GetValue("mu2e.HistDir",gSystem->Getenv("MU2E_HIST_DIR"));

    sprintf(fn,"%s/v2_1_0/e0000001.tcalm002.hist",mu2e_hist_dir);
    fLogLH->InitEleDtHist(fn);
    fLogLH->InitEleEpHist(fn);

    sprintf(fn,"%s/v2_1_0/m0000001.tcalm002.hist",mu2e_hist_dir);
    fLogLH->InitMuoDtHist(fn);
    fLogLH->InitMuoEpHist(fn);

    sprintf(fn,"%s/v3_0_0/e0000001_tcalm002_vane.hist",mu2e_hist_dir);
    TH1* he_xs = gh1(fn,"TCalm002","trk_1/xslope");
    fLogLH->SetEleXsHist(he_xs);

    sprintf(fn,"%s/v3_0_0/m0000001_tcalm002_vane.hist",mu2e_hist_dir);
    TH1* hm_xs = gh1(fn,"TCalm002","trk_1/xslope");
    fLogLH->SetMuoXsHist(hm_xs);
  }


//-----------------------------------------------------------------------------
  void TCalm002::endJob() {
    TDirectory* old_dir = gDirectory;
    
    TFile* f = TFile::Open(fHistFileName.data(),"recreate");

    SaveFolder(fFolder,f);
    f->Write();
    f->Close();
    
    old_dir->cd();
  }

//-----------------------------------------------------------------------------
// begin run - perform run-dependent and geometry-dependent initializations
//-----------------------------------------------------------------------------
  bool TCalm002::beginRun(art::Run& Run) {
    static int first_call (1);

    TDiskCalorimeter::GeomData_t data;

    const mu2e::DiskCalorimeter* cal;

    if (first_call == 1) {
					// geometry choice - disks or vanes

      art::ServiceHandle<mu2e::GeometryService> geom;

      if (geom->hasElement<mu2e::VaneCalorimeter>() ) {
	fCalorimeterType = 1;
      }
      else if (geom->hasElement<mu2e::DiskCalorimeter>() ) {
//-----------------------------------------------------------------------------
// create disk calorimeter and initialize it
//-----------------------------------------------------------------------------
	mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;
	const mu2e::Disk*                       disk;

	fCalorimeterType = 2;
	cal = dc.operator->();

	data.fNDisks = cal->nDisk();
	for (int i=0; i<data.fNDisks; i++) {
	  disk = &cal->disk(i);
	  data.fRMin[i] = disk->innerRadius();
	  data.fRMax[i] = disk->outerRadius();
	  data.fZ0  [i] = disk->origin().z();
	}
	data.fHexSize     = cal->caloGeomInfo().crystalHalfTrans();
	data.fMinFraction = fMinCrystalFr;

	fDiskCalorimeter  = new TDiskCalorimeter(&data);
      }

    }
    return true;
  }
//-----------------------------------------------------------------------------
  void TCalm002::FillCaloHistograms(CaloHist_t* Hist, TStnCrystal* Cr) {

    int                    nhits;
    float                  t, e, r, e700, n700;
    TObjHandle*            phit;
    mu2e::CaloCrystalHit*  hit;
					// determine crystal coordinates
    TDisk* disk = Cr->Disk();

    int idisk = disk->SectionID();
					// time needs to be defiend
    //    t  = Cr->Time();
    e     = Cr->Energy();
    r     = Cr->Radius();
    nhits = Cr->NHits();

    Hist->fVaneID->Fill(idisk);

    Hist->fEnergy  [idisk]->Fill(e);
    Hist->fNHits   [idisk]->Fill(nhits);
    //    Hist->fTime    [idisk]->Fill(t);
    Hist->fRadius  [idisk]->Fill(r);
    Hist->fRadiusWE[idisk]->Fill(r,e);
    
    e700 = 0;
    n700 = 0;
    for (int i=0; i<nhits; i++) {
      phit = (TObjHandle*) Cr->ListOfHits()->UncheckedAt(i);
      hit  = (mu2e::CaloCrystalHit*) phit->Object();
      t   = hit->time();
      Hist->fTime[idisk]->Fill(t);
      if (t > 700.) {
	n700 += 1;
	e700 += hit->energyDep();
	Hist->fT700[idisk]->Fill(t);
      }
    }

    Hist->fE700   [idisk]->Fill(e700);
    Hist->fN700   [idisk]->Fill(n700);

    if (n700 > 0) {
      Hist->fR700  [idisk]->Fill(r);
      Hist->fRWE700[idisk]->Fill(r,e700);
    }
  }

//-----------------------------------------------------------------------------
  void TCalm002::FillClusterHistograms(ClusterHist_t* Hist, TStnCluster* Cluster) {

    int   row, col;
    float  x, y, z, r;

    const CaloCluster* cl = Cluster->fCaloCluster;

    row = -999; // cl->cogRow();
    col = -999; // cl->cogColumn();
//-----------------------------------------------------------------------------
// starting from v5_3 the cluster coordinates are reconstructed in the local 
// coordinate system of the disk
//-----------------------------------------------------------------------------
    x   = cl->cog3Vector().x(); // +3904.;
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
  void TCalm002::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track) {

    //    HelixParams h (Trk-helix(0.));
    double dp, dpfst;

    KalRep* trk = Track->fKalRep[0];

    Hist->fP[0]->Fill (Track->fP);
    Hist->fP[1]->Fill (Track->fP);
    Hist->fP[2]->Fill (Track->fP);
    Hist->fPt->Fill(Track->fPt);
    Hist->fPFront->Fill(Track->fPFront);
    Hist->fPStOut->Fill(Track->fPStOut);

    dp    = Track->fP-Track->fPFront;
    dpfst = Track->fPFront-Track->fPStOut;

    Hist->fDpFront->Fill(dp);
    Hist->fDpFSt->Fill(dpfst);

    Hist->fCosTh->Fill(Track->Momentum()->CosTheta());
    Hist->fChi2->Fill (Track->fChi2);
    Hist->fNDof->Fill(Track->fNActive-5.);
    Hist->fChi2Dof->Fill(Track->fChi2/(Track->fNActive-5.));
    Hist->fNActive->Fill(Track->fNActive);
    Hist->fT0->Fill(Track->fT0);
    Hist->fQ->Fill(trk->charge());
    Hist->fFitCons[0]->Fill(Track->fFitCons);
    Hist->fFitCons[1]->Fill(Track->fFitCons);

    Hist->fD0->Fill(trk->helix(0).d0());
    Hist->fZ0->Fill(trk->helix(0).z0());
    Hist->fTanDip->Fill(trk->helix(0).tanDip());

    const mu2e::TrkStrawHit* hit;
    const TrkHotList*        hot_list = trk->hotList();

    double chi2c(0), dr;
    for(TrkHotList::hot_iterator it=hot_list->begin(); it<hot_list->end(); it++) {
      hit = (const mu2e::TrkStrawHit*) &(*it);
      dr  = hit->resid();
      Hist->fResid->Fill(dr);
					// normalize by 100 microns
      chi2c += dr*dr/0.01;
    }

    chi2c = chi2c/(Track->fNActive-5.);
    Hist->fChi2DofC->Fill(chi2c);

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
    TStnTrack::InterData_t*    vr = Track->fVMaxEp; // residuals
    double  ep, r;

    if (vt) {
      Hist->fVaneID->Fill(vt->fID  );
      Hist->fXTrk->Fill  (vt->fXTrk);
      Hist->fYTrk->Fill  (vt->fYTrk);

      r = sqrt(vt->fXTrk*vt->fXTrk+vt->fYTrk*vt->fYTrk);
      Hist->fRTrk->Fill  (r);

      Hist->fZTrk->Fill  (vt->fZTrk);
    }
    else {
					// fill histograms with numbers, easy to recognize as dummy
      Hist->fVaneID->Fill(-1.);
      Hist->fXTrk->Fill  (999.);
      Hist->fYTrk->Fill  (999.);
      Hist->fRTrk->Fill  (999.);
      Hist->fZTrk->Fill  (-1. );
    }

    if (vr) {
      Hist->fDt->Fill(vr->fDt);
      ep = vr->fCluster->energyDep()/Track->fP;
      Hist->fEp->Fill(ep);
      Hist->fDx->Fill(vr->fDx);
      Hist->fDy->Fill(vr->fDy);
      Hist->fDz->Fill(vr->fDz);
    }

    int ncl = Track->NClusters();
    Hist->fNClusters->Fill(ncl);

    Hist->fRSlope->Fill(Track->RSlope());
    Hist->fXSlope->Fill(Track->XSlope());

    Hist->fEleLogLHCal->Fill(Track->EleLogLHCal());
    Hist->fMuoLogLHCal->Fill(Track->MuoLogLHCal());
    Hist->fLogLHRCal->Fill(Track->LogLHRCal());

    double llhr_dedx, llhr_xs, llhr_trk;

    llhr_dedx = Track->LogLHRDeDx();
    llhr_xs   = Track->LogLHRXs();
    llhr_trk  = Track->LogLHRTrk();

    Hist->fLogLHRDeDx->Fill(llhr_dedx);
    Hist->fLogLHRXs->Fill(llhr_xs);
    Hist->fLogLHRTrk->Fill(llhr_trk);

    Hist->fPdgCode->Fill(Track->fPdgCode);
    Hist->fFrGH->Fill(Track->fNGoodMcHits/(Track->fNActive+1.e-5));

  }

//-----------------------------------------------------------------------------
  void TCalm002::FillEventHistograms(EventHist_t* Hist) {

    double   cos_th, xv, yv, rv, zv;

    cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();

    xv = fEle->position().x()+3904.;
    yv = fEle->position().y();
    
    rv = sqrt(xv*xv+yv*yv);
    zv = fEle->position().z();

    Hist->fEleCosTh->Fill(cos_th);
    Hist->fRv->Fill(rv);
    Hist->fZv->Fill(zv);

    Hist->fNClusters->Fill(fNClusters);
    Hist->fNTracks->Fill(fNTracks[0]);
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
      t0_trk = fTrack->fT0;
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
      Hist->fSHTime->Fill(sh->time());

      if (fabs(dt+15.)< 50) n_good_hits += 1;
    }

    Hist->fNGoodSH->Fill(n_good_hits);

    Hist->fNHyp->Fill(fNHyp);
    Hist->fBestHyp[0]->Fill(fBestHyp[0]);
    Hist->fBestHyp[1]->Fill(fBestHyp[1]);
//-----------------------------------------------------------------------------
// crystals - count crystals with E > 1MeV
//-----------------------------------------------------------------------------
    CaloCrystalHit* cch;

    int n_cch_1mev = 0;


    if (fCalorimeterType == 1) {
					// vane calorimeter

      int  nhits_vane[2][4], nhits_row [2][20], nhits_col[2][50];
      int  crystal_id, vane_id, local_id, vane_row, vane_col;

      for (int i=0; i<4; i++) {
	nhits_vane[0][i] = 0;
	nhits_vane[1][i] = 0;
      }
      
      for (int i=0; i<20; i++) {
	nhits_row[0][i] = 0;
	nhits_row[1][i] = 0;
      }

      for (int i=0; i<50; i++) {
	nhits_col[0][i] = 0;
	nhits_col[1][i] = 0;
      }
      
      for (int ic=0; ic<fNCaloCrystalHits; ic++) {
	cch        = &fListOfCaloCrystalHits->at(ic);
	crystal_id = cch->id();

	if (cch->energyDep() > 1.) {
	  n_cch_1mev += 1;
	}
					// for each crystal determine its row and column
					// the following is for vanes
	vane_id  = crystal_id/484.;
	local_id = crystal_id-vane_id*484;
	vane_row = local_id/44;
	vane_col = local_id-vane_row*44;
      
	nhits_vane[0][vane_id ] += 1;
	nhits_row [0][vane_row] += 1;
	nhits_col [0][vane_col] += 1;
      
	if (cch->energyDep() > 1.) {
	  nhits_row [1][vane_row] += 1;
	  nhits_col [1][vane_col] += 1;
	  nhits_vane[1][vane_id ] += 1;
	}
      }

      Hist->fNCaloCrystalHits[0]->Fill(fNCaloCrystalHits);
      Hist->fNCaloCrystalHits[1]->Fill(n_cch_1mev);

      for (int iv=0; iv<4; iv++) {
	Hist->fNCaloHitsVsVane[0]->Fill(iv,nhits_vane[0][iv]);
	Hist->fNCaloHitsVsVane[1]->Fill(iv,nhits_vane[1][iv]);
      }

      for (int ir=0; ir<20; ir++) {
	Hist->fNCaloHitsVsRow[0]->Fill(ir,nhits_row[0][ir]);
	Hist->fNCaloHitsVsRow[1]->Fill(ir,nhits_row[1][ir]);
      }

      for (int ic=0; ic<50; ic++) {
	Hist->fNCaloHitsVsCol[0]->Fill(ic,nhits_col[0][ic]);
	Hist->fNCaloHitsVsCol[1]->Fill(ic,nhits_col[1][ic]);
      }
    }
    else if (fCalorimeterType == 2) {
					// disk calorimeter

      int nc, ndisks, nhits, n_hit_crystals, n_hit_crystals_tot;
      double   etot;

      TDisk*       disk;
      TStnCrystal* cr;
      float        e, r;

      ndisks         = fDiskCalorimeter->NDisks();

      int   bin, nhits_r[100], n_hit_crystals_r[100];

      n_hit_crystals_tot = 0;

      for (int i=0; i<ndisks; i++) {
	disk = fDiskCalorimeter->Disk(i);

	nc             = disk->NCrystals();
	etot           = 0;
	n_hit_crystals = 0;
	nhits          = 0;

	for (int ib=0; ib<100; ib++) {
	  nhits_r[ib]         = 0;
	  n_hit_crystals_r[ib] = 0;
	}
					// loop over crystals
	for (int ic=0; ic<nc; ic++) {
	  cr   = disk->Crystal(ic);
	  r    = cr->Radius();
	  e    = cr->Energy();
	  etot = etot+e;

	  bin  = (int) (r/10.);

	  if (cr->Energy() > 0) {
	    nhits          += cr->NHits();
	    n_hit_crystals += 1;

	    nhits_r         [bin] += cr->NHits();
	    n_hit_crystals_r[bin] += 1;
	  }

	  Hist->fECrVsR[i]->Fill(r,e);
	  Hist->fNCrVsR[i]->Fill(r,1);
	}

	n_hit_crystals_tot += n_hit_crystals;

	Hist->fETot[i]->Fill(etot);
//-----------------------------------------------------------------------------
// 100 is fixed by the number of bins in the radial distributions
//-----------------------------------------------------------------------------
	for (int ib=0; ib<100; ib++) {
	  r = (ib+0.5)*10.;
	  Hist->fNCrystalHitsVsR[i]->Fill(r,nhits_r[ib]);
	  Hist->fNHitCrystalsVsR[i]->Fill(r,n_hit_crystals_r[ib]);
	}
      }

      Hist->fNHitCrystalsTot->Fill(n_hit_crystals_tot);
    }
  }
//-----------------------------------------------------------------------------
// so far I'm assuming that there is only one track
//-----------------------------------------------------------------------------
  void TCalm002::DefineBestHypothesis() {
    double   chi2, chi2min(1.e6), fcon, fconmax(-1.);

    KalRep*  trk;

    fNHyp       = 0;
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
	  fconmax     = fcon;
	  fBestHyp[1] = ih;
	}
      }
    }
  }

//-----------------------------------------------------------------------------
  void TCalm002::FillHistograms() {
    double     cos_th, cl_e(-1.);
    int        disk_id(-1);
    TStnCluster  *cl0;

    cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();

    if (fNClusters > 0) {
      cl0     = fClusterBlock->Cluster(0);
      cl_e    = cl0->Energy();
      disk_id = cl0->DiskID();
    }
//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
    FillEventHistograms(fHist.fEvent[0]);

    if (fNTracks[0]> 0) FillEventHistograms(fHist.fEvent[1]);
    else                FillEventHistograms(fHist.fEvent[2]);

    if (fNClusters > 0) FillEventHistograms(fHist.fEvent[3]);
    else                FillEventHistograms(fHist.fEvent[4]);

    if ((fNTracks[0] == 0) && (fabs(cos_th) < 0.4)) {
      FillEventHistograms(fHist.fEvent[5]); 
    }

    if (fNGoodTracks > 0) {
      FillEventHistograms(fHist.fEvent[6]); 
    }

    if (cl_e > 70.) {
      FillEventHistograms(fHist.fEvent[7]); 
    }

    if      (disk_id == 0) FillEventHistograms(fHist.fEvent[8]);
    else if (disk_id == 1) FillEventHistograms(fHist.fEvent[9]);
//-----------------------------------------------------------------------------
// track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
    TStnTrack*   trk;

    //    HelixParams  helix;

    for (int i=0; i<fNTracks[0]; ++i ) {
      trk = fTrackBlock->Track(i);
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
//-----------------------------------------------------------------------------
// TRK 8: good track, |xslope| < 3
//-----------------------------------------------------------------------------
	if (fabs(trk->XSlope()) < 3.) {
	  FillTrackHistograms(fHist.fTrack[8],trk);
	}
      }
//-----------------------------------------------------------------------------
// TRK 6 : events with a track and a cluster E/P > 0.4
// TRK 7 : events with a "Set C" track and a cluster E/P > 0.4
//-----------------------------------------------------------------------------
      if ((trk->Ep() > 0.4) && ( trk->Ep() < 1.2)) {
	FillTrackHistograms(fHist.fTrack[6],trk);
	if (trk->fIDWord == 0) {
	  FillTrackHistograms(fHist.fTrack[7],trk);
	}
      }
//----------------------------------------------------------------------------
//TRK  9: events with track and with no EM cluster with E < 70 MeV
//TRK 10: events with good track and with no EM cluster with E < 70 MeV
//TRK 10: events with good track and with no EM cluster with E < 70 MeV
//----------------------------------------------------------------------------
      if (cl_e > 70) {
	FillTrackHistograms(fHist.fTrack[9],trk);
	if (trk->fIDWord == 0) {
	  FillTrackHistograms(fHist.fTrack[10],trk);
	}
      }
//-----------------------------------------------------------------------------
// TRK 11: tracks with P > 103.5 MeV
//-----------------------------------------------------------------------------
      if (trk->P() > 103.5) FillTrackHistograms(fHist.fTrack[11],trk);
//-----------------------------------------------------------------------------
// TRK 12: tracks with fcon < 1e-4
// TRK 13: "Set C" tracks with 100 <= P < 110 
// TRK 14: tracks with fcon < 1e-2
//-----------------------------------------------------------------------------
      if (trk->fFitCons < 1.e-4) FillTrackHistograms(fHist.fTrack[12],trk);

      if ((trk->fIDWord == 0) && (trk->P() >= 100.) && (trk->P() < 110.)) {
	FillTrackHistograms(fHist.fTrack[13],trk);
      }

      if (trk->fFitCons < 1.e-2) FillTrackHistograms(fHist.fTrack[14],trk);

      TStnTrack::InterData_t*    vt = trk->fVMinS;  // track-only

//-----------------------------------------------------------------------------
// TRK 15: tracks which have intersection with the 1st disk
// TRK 16: tracks which have intersection with the 2nd disk
// TRK 17: tracks which do not have intersections with the calorimeter
//-----------------------------------------------------------------------------
      if (vt) {
	if      (vt->fID == 0) FillTrackHistograms(fHist.fTrack[15],trk);
	else if (vt->fID == 1) FillTrackHistograms(fHist.fTrack[16],trk);
      }
      else {
	FillTrackHistograms(fHist.fTrack[17],trk);
      }
    }
//-----------------------------------------------------------------------------
// cluster histograms
//-----------------------------------------------------------------------------
    TStnCluster* cl;
    for (int i=0; i<fNClusters; ++i ) {
      cl = fClusterBlock->Cluster(i);

      FillClusterHistograms(fHist.fCluster[0],cl);

      if (fNTracks[0]     >  0 ) FillClusterHistograms(fHist.fCluster[1],cl);
      if (fNGoodTracks    >  0 ) FillClusterHistograms(fHist.fCluster[2],cl);
      if (fNMatchedTracks >  0 ) FillClusterHistograms(fHist.fCluster[3],cl);
      if (cl->Energy()    > 10.) FillClusterHistograms(fHist.fCluster[4],cl);
      if (cl->Energy()    > 70.) FillClusterHistograms(fHist.fCluster[5],cl);
    }
//-----------------------------------------------------------------------------
// calorimeter histograms
//-----------------------------------------------------------------------------
    TDisk*         disk;
    TStnCrystal*   cr;

    if (fCalorimeterType == 2) {
      int nd = fDiskCalorimeter->NDisks();

      for (int i=0; i<nd; i++) {
	disk = fDiskCalorimeter->Disk(i);
	for (int ic=0; ic<disk->NCrystals(); ic++) {
	  cr = disk->Crystal(ic);
	  FillCaloHistograms(fHist.fCalo[0],cr);
	  if (cr->Energy() > 0) {
	    FillCaloHistograms(fHist.fCalo[1],cr);
	  }
	  if (cr->Energy() > 0.1) {
	    FillCaloHistograms(fHist.fCalo[2],cr);
	  }
	  if (cr->Energy() > 1.0) {
	    FillCaloHistograms(fHist.fCalo[3],cr);
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// radial distributions for crystals
//-----------------------------------------------------------------------------
    static int first_entry(1);

    if (first_entry == 1) {
      if (fCalorimeterType == 2) {
	int nd = fDiskCalorimeter->NDisks();
	
	for (int i=0; i<nd; i++) {
	  disk = fDiskCalorimeter->Disk(i);
	  for (int ic=0; ic<disk->NCrystals(); ic++) {
	    cr = disk->Crystal(ic);

	    fHist.fCrystalR[i]->Fill(cr->Radius());
	  }
	}
      }
    }

    first_entry = 0;
  }



//-----------------------------------------------------------------------------
  void TCalm002::getData(const art::Event& Evt) {
//-----------------------------------------------------------------------------
// generated particles
//-----------------------------------------------------------------------------
    art::Handle<GenParticleCollection> genpHandle;
    Evt.getByLabel("generate",genpHandle);
    fListOfGenParticles = (GenParticleCollection*) &(*genpHandle);
    fEle                = (GenParticle*) &fListOfGenParticles->at(0);

    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    Evt.getByLabel(fStrawHitMaker,"StrawHitMCPtr",mcptrHandle);

    fListOfMcStrawHits = 0;
    if (mcptrHandle.isValid()) {
      fListOfMcStrawHits = (PtrStepPointMCVectorCollection*) mcptrHandle.product();
    }
//-----------------------------------------------------------------------------
// reconstructed tracks
//-----------------------------------------------------------------------------
//     art::Handle<KalRepPtrCollection> krepsHandle;
//     Evt.getByLabel("TrkPatRec1","DownstreameMinus", krepsHandle);
//     fListOfTracks[0] = (KalRepPtrCollection*) &(*krepsHandle);
//     fNTracks     [0] = fListOfTracks[0]->size();
//-----------------------------------------------------------------------------
// alternative hypotheses
//-----------------------------------------------------------------------------
//     art::Handle<KalRepPtrCollection> krepsHandle2;
//     Evt.getByLabel("trkPatRec2","UpstreamePlus", krepsHandle2);
//     fListOfTracks[1] = (KalRepPtrCollection*) &(*krepsHandle2);
//     fNTracks[1]      = fListOfTracks[1]->size();

//     art::Handle<KalRepPtrCollection> krepsHandle3;
//     Evt.getByLabel("trkPatRec3","DownstreammuMinus", krepsHandle3);
//     fListOfTracks[2] = (KalRepPtrCollection*) &(*krepsHandle3);
//     fNTracks     [2] = fListOfTracks[2]->size();

//     art::Handle<KalRepPtrCollection> krepsHandle4;
//     Evt.getByLabel("trkPatRec4","UpstreammuPlus", krepsHandle4);
//     fListOfTracks[3] = (KalRepPtrCollection*) &(*krepsHandle4);
//     fNTracks     [3] = fListOfTracks[3]->size();
//-----------------------------------------------------------------------------
// reconstructed calorimeter clusters
//-----------------------------------------------------------------------------
//     art::Handle<CaloClusterCollection> calo_cluster_handle;
//     Evt.getByLabel("makeCaloCluster","AlgoCLOSESTSeededByENERGY",calo_cluster_handle);
//     fListOfClusters  = (CaloClusterCollection*) &(*calo_cluster_handle);
//     fNClusters       = fListOfClusters->size();
//-----------------------------------------------------------------------------
// straw hits 
//-----------------------------------------------------------------------------
    art::Handle<StrawHitCollection> shHandle;
    Evt.getByLabel("makeSH",shHandle);
    fListOfStrawHits = (StrawHitCollection*) &(*shHandle);
    fNStrawHits      = fListOfStrawHits->size();
//-----------------------------------------------------------------------------
// calorimeter clusters
//-----------------------------------------------------------------------------
    art::Handle<CaloCrystalHitCollection> ccHandle;
    Evt.getByLabel(fCrystalHitMaker.data(),ccHandle);
    fListOfCaloCrystalHits = (CaloCrystalHitCollection*) ccHandle.product();
    fNCaloCrystalHits      = fListOfCaloCrystalHits->size();
//-----------------------------------------------------------------------------
// track extrapolation results
//-----------------------------------------------------------------------------
//     art::Handle<TrkToCaloExtrapolCollection>  texHandle;
//     Evt.getByLabel("trkExtrapol",texHandle);
//     fListOfExtrapolatedTracks = (TrkToCaloExtrapolCollection*) &(*texHandle);

//     art::Handle<TrackClusterLink>  tcmapHandle;
//     Evt.getByLabel("caloMatching",tcmapHandle);
//     fTrkCalMap = tcmapHandle.operator->();
  }



//-----------------------------------------------------------------------------
  void TCalm002::Init(art::Event& Evt) {
//    TStnCluster*    cluster;

    int                id_word, ntrk;
    double             xs;
    TStnTrack*         track;
    TEmuLogLH::PidData_t  dat;
//-----------------------------------------------------------------------------
// initialize tracks and determine track quality
//-----------------------------------------------------------------------------
    StntupleInitMu2eTrackBlock(fTrackBlock,&Evt,0);

    fNGoodTracks    = 0;
    fNMatchedTracks = 0;

    fNTracks[0] = fTrackBlock->NTracks();
    if (fNTracks[0] == 0) fTrack = 0;
    else                  fTrack = fTrackBlock->Track(0);

    ntrk = fNTracks[0];

    for (int i=0; i<ntrk; i++) {
      track          = fTrackBlock->Track(i);
      id_word        = fTrackID->IDWord(track);
      track->fIDWord = id_word;
      if (id_word == 0) {
	fNGoodTracks += 1;
	if ((track->fVMaxEp != NULL) && (fabs(track->fVMaxEp->fDt) < 2.5)) {
	  fNMatchedTracks += 1;
	}
      }

      dat.fDt   = track->Dt();
      dat.fEp   = track->Ep();
      dat.fPath = -1.;
      if (track->fVMinS) dat.fPath = track->fVMinS->fPath;
      
      xs = track->XSlope();

      track->fEleLogLHCal = fLogLH->LogLHCal(&dat,11);
      track->fMuoLogLHCal = fLogLH->LogLHCal(&dat,13);
      track->fLogLHRXs    = fLogLH->LogLHRXs(xs);
    }
//-----------------------------------------------------------------------------
// initialize calorimeter and cluster data blocks
//-----------------------------------------------------------------------------
    StntupleInitMu2eCalDataBlock(fCalDataBlock,&Evt,0);
    StntupleInitMu2eClusterBlock(fClusterBlock,&Evt,0);


    fNClusters = fClusterBlock->NClusters();
    if (fNClusters == 0) fCluster = 0;
    else                 fCluster = fClusterBlock->Cluster(0);
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
//-----------------------------------------------------------------------------
// init calorimeter
//-----------------------------------------------------------------------------
    if      (fCalorimeterType == 1) {
    }
    else if (fCalorimeterType == 2) {
      fDiskCalorimeter->InitEvent(fCalDataBlock);
    }
  }

//-----------------------------------------------------------------------------
  void TCalm002::Debug(const art::Event& Evt) {
    const char* oname = "TCalm002::Debug";
    double      pt;
    TStnTrack*  t;
    int         nt;

 //-----------------------------------------------------------------------------
// bit 000: 
//-----------------------------------------------------------------------------
   if (TModule::fDebugBit[0] != 0) {
      nt = fTrackBlock->NTracks();
      for (int i=0; i<nt; i++) {
	t = fTrackBlock->Track(i);
	pt = t->Momentum()->Pt();
	if ((t->fIDWord == 0) && (t->NClusters() == 0)) {
	  printf(" >>>>>>> [%s] EVENT : %10i BIT:000 TRACK IDWORD==0 Pt = %10.3f doesnt have a cluster\n",
		 oname,Evt.event(),pt);
	  //	  t->Print("");
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 001: events with a  track, clusters, but no matching cluster with E/P > 0.4
//-----------------------------------------------------------------------------
    if (TModule::fDebugBit[1] != 0) {
      nt = fTrackBlock->NTracks();
      for (int i=0; i<nt; i++) {
	t = fTrackBlock->Track(i);
	if ((t->fIDWord == 0) && (fNClusters > 0) && (t->Ep() < 0.4)) {
	  printf(" >>>>>>> [%s] EVENT : %10i BIT:001 TRACK with ID_WORD = 0, E/P = %10.3f\n",
		 oname,Evt.event(),t->Ep());
	  //	  t->Print("");
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 002: events with a track, cluster, and DT < 3ns
//-----------------------------------------------------------------------------
    if (TModule::fDebugBit[2] != 0) {
      nt = fTrackBlock->NTracks();
      for (int i=0; i<nt; i++) {
	t = fTrackBlock->Track(i);
	if ((t->fIDWord == 0) && (t->fVMinS != NULL)) {

	  double dt = t->fVMinS->fDt;

	  if (fabs(dt) > 2.5) {
	    printf(" >>>>>>> [%s] EVENT : %10i BIT:002 TRACK with ID_WORD = 0 DT = %10.3f\n",
		   oname,Evt.event(),dt);
	  }
	  //	  t->Print("");
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 003: print downstream - tracks
//-----------------------------------------------------------------------------
    if (TModule::fDebugBit[3] != 0) {
      TModule::fDump->printKalRepCollection(fTrkPatRecDem.data(),"","",1);
    }
//-----------------------------------------------------------------------------
// bit 004: events with a track with the fit consistency < 1e-5
//-----------------------------------------------------------------------------
    if (TModule::fDebugBit[4] != 0) {
      nt = fTrackBlock->NTracks();
      for (int i=0; i<nt; i++) {
	t = fTrackBlock->Track(i);
	if (t->fFitCons < 1.e-5) {
	  printf(" >>>>>>> [%s] EVENT : %10i BIT:004 TRACK with FIT_CONS = %12.3e\n",
		 oname,Evt.event(),t->fFitCons);
	}
	if (t->fFitCons >= 1.0) {
	  printf(" >>>>>>> [%s] EVENT : %10i BIT:004 TRACK with FIT_CONS = %12.3e\n",
		 oname,Evt.event(),t->fFitCons);
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 005: for each event print NTracks, NClusters, E(max)
//-----------------------------------------------------------------------------
    if (TModule::fDebugBit[5] != 0) {
      int   nt, ncl;
      float pmax, emax;
      TStnCluster* cl;

      nt = fTrackBlock->NTracks();
      if (nt > 0) pmax = fTrackBlock->Track(0)->P();
      else        pmax = -1.;

      ncl = fClusterBlock->NClusters();
      if (ncl > 0) {
	cl   = fClusterBlock->Cluster(0);
	emax = cl->Energy();
      }
      else         emax = -1.;

      printf(" >>>>>>> [%s] EVENT : %10i BIT:005 NTRACKS = %2i Pmax = %10.3f ",
	     oname,Evt.event(),nt,pmax);
      printf(" NCLUSTERS = %2i  EMAX = %10.3f\n",ncl,emax);
    }

  }


//-----------------------------------------------------------------------------
  bool TCalm002::filter(art::Event& Evt) {
    const char* oname = "TCalm002::analyze";

    bool rc(true);
    
    if(Evt.event()% 50 == 0 ){
      printf(" >>>>>>> [%s] RUN:EVENT : %10i:%10i\n",
	     oname,Evt.run(),Evt.event());
    }
    fDump->SetEvent(Evt);

    getData(Evt);

    Init(Evt);

    Debug(Evt);

    FillHistograms();

    if (TModule::fDebugBit[51] != 0) {
      rc = (fNMatchedTracks > 0);
    }

    if (TModule::fDebugBit[52] != 0) {
      rc = (fNCaloCrystalHits > 0);
      if (rc != 0) {
	printf(" >>>>>>> [%s] BIT:052 RUN:EVENT : %10i:%10i N(CaloCrystalHits) : %i\n",
	       oname,Evt.run(),Evt.event(),fNCaloCrystalHits);
      
      }
    }

    TModule::filter(Evt);

    return rc;
  }

}  // end namespace mu2e

using mu2e::TCalm002;

DEFINE_ART_MODULE(TCalm002);
