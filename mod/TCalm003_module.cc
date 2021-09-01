//-----------------------------------------------------------------------------
// 2013-04-24 P.Murat: optimize Z positions of the calorimeter disks
// bit usage: 
// can run on MC generated with vane-based geometry
//-----------------------------------------------------------------------------

// #include "Stntuple/obj/AbsEvent.hh"
#include "murat/mod/TCalm003_module.hh"
#include "Stntuple/alg/TStnTrackID.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

#include "BTrk/KalmanTrack/KalRep.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

namespace mu2e {

//-----------------------------------------------------------------------------
  TCalm003::TCalm003(fhicl::ParameterSet const& pset): 
    THistModule(pset,"TCalm003"),
    fHistFileName (pset.get<std::string> ("histFileName" , "tcalm003.hist")),
    fStrawHitMaker(pset.get<std::string> ("strawHitMaker", "makeSH"       )),
    fTrkExtrapol  (pset.get<std::string> ("trkExtrapol"  , "trkExtrapol"  )),
    fTrkCalMatch  (pset.get<std::string> ("trkCalMatch"  , "caloMatching" )),
    fMinTActive   (pset.get<double>      ("minTActive"   ,   710.         ))
  {
					// reset all histogram pointers

    for (int i=0; i<kNEventHistSets; i++) {
      fHist.fEvent[i] = 0;
    }

    for (int i=0; i<kNTrackHistSets; i++) {
      fHist.fTrack[i] = 0;
    }

    fTrackBlock   = new TStnTrackBlock  ();
    fClusterBlock = new TStnClusterBlock();
    fTrackID      = new TStnTrackID     ();
  }



//-----------------------------------------------------------------------------
  TCalm003::~TCalm003() {
    //     fFolder->Delete();
    //     delete fFolder;

    delete fTrackBlock;
    //    delete fClusterBlock;
    delete fTrackID;
  }

//-----------------------------------------------------------------------------
  void TCalm003::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
//-----------------------------------------------------------------------------
//  
//-----------------------------------------------------------------------------
    HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
    HBook1F(Hist->fNTracks   ,"ntrk"     ,Form("%s: Number of Reconstructed Tracks"  ,Folder),100,0,100,Folder);
    
    for (int i=0; i<kNPlanes; i++) {
      HBook1F(Hist->fRTrack[i],
	      Form("rt_%02i",i),
	      Form("%s: R(track) in plane number %i, Z= %8.3f",Folder,i,fZPlane[i]),
	      100,0,1000,Folder);
    }

    char name[200], title[200];

    for (int i1=0; i1<kNPlanes; i1++) {
      for (int i2=0; i2<kNPlanes; i2++) {
	sprintf(name,"rt0_%02i_%02i",i1,i2);
	sprintf(title,"%s: R(track) in plane number %i, Z= %8.3f, MISS_36 in plane %i",Folder,i2,fZPlane[i2],i1);
	HBook1F(Hist->fRTrack0[i1][i2],name,title,100,0,1000,Folder);
	
	sprintf(name,"rt1_%02i_%02i",i1,i2);
	sprintf(title,"%s: R(track) in plane number %i, Z= %8.3f, HIT_36 in plane %i",Folder,i2,fZPlane[i2],i1);
	HBook1F(Hist->fRTrack1[i1][i2],name,title,100,0,1000,Folder);

	sprintf(name,"rt33_0_%02i_%02i",i1,i2);
	sprintf(title,"%s: R(track) in plane number %i, Z= %8.3f, MISS_33 in plane %i",Folder,i2,fZPlane[i2],i1);
	HBook1F(Hist->fRTrack330[i1][i2],name,title,100,0,1000,Folder);
	
	sprintf(name,"rt33_1_%02i_%02i",i1,i2);
	sprintf(title,"%s: R(track) in plane number %i, Z= %8.3f, HIT_33 in plane %i",Folder,i2,fZPlane[i2],i1);
	HBook1F(Hist->fRTrack331[i1][i2],name,title,100,0,1000,Folder);
      }
    }
  }

//-----------------------------------------------------------------------------
  void TCalm003::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
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
    HBook1F(Hist->fYTrk       ,"ytrk"     ,Form("%s: track YTrk"        ,Folder), 100,-250,250,Folder);
    HBook1F(Hist->fZTrk       ,"ztrk"     ,Form("%s: track ZTrk"        ,Folder), 200,-1000,1000,Folder);
    HBook1F(Hist->fEp         ,"ep"       ,Form("%s: track E/P"         ,Folder), 300, 0   ,1.5,Folder);
    HBook2F(Hist->fNHVsStation,"nh_vs_st" ,Form("%s: N(hits) Vs Station",Folder),  40, 0,40,10,-0.5,9.5,Folder);
    HBook2F(Hist->fNHVsNSt    ,"nh_vs_nst",Form("%s: N(hits) Vs NSt"    ,Folder),  10,-0.5,9.5,40,-0.5,39.5,Folder);
    HBook1F(Hist->fPdgCode    ,"pdg"      ,Form("%s: track PDG code"    ,Folder), 100,-50,50,Folder);
    HBook1F(Hist->fFrGH       ,"fgh"      ,Form("%s: Fraction Goog Hits",Folder), 100, 0,1,Folder);
  }

//-----------------------------------------------------------------------------
  void TCalm003::BookHistograms() {
    TFolder  *hist_folder, *fol;
    char     folder_name[200];

    hist_folder = (TFolder*) fFolder->FindObject("Hist");

//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[0] = 1;		// all events
    book_event_histset[1] = 1;	        // events with reconstructed track and R0 < 360 mm
    book_event_histset[2] = 1;	        // events with reconstructed track and R0 < 330 mm

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

    for (int i=0; i<kNTrackHistSets; i++) {
      if (book_track_histset[i] != 0) {
	sprintf(folder_name,"trk_%i",i);
	fol = (TFolder*) hist_folder->FindObject(folder_name);
	if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
	fHist.fTrack[i] = new TrackHist_t;
	BookTrackHistograms(fHist.fTrack[i],Form("Hist/%s",folder_name));
      }
    }
  }

//-----------------------------------------------------------------------------
// begin job - book histograms etc
//-----------------------------------------------------------------------------
  void TCalm003::beginJob() {

    fTrackID->SetMinT0(fMinTActive);
//-----------------------------------------------------------------------------
// initialize Z coordinates of the planes with respect to the tracker center 
//-----------------------------------------------------------------------------
    for (int i=0; i<kNPlanes; i++) {
      fZPlane[i] = 1500.+50.*i;
    }
					// optimized position of the first disk
    fZDisk1 = 1500.+230.;
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
    BookHistograms();
//-----------------------------------------------------------------------------
// define collection names
//-----------------------------------------------------------------------------
    fClusterBlock->AddCollName("mu2e::CaloClusterCollection",
  			       "makeCaloCluster","AlgoCLOSESTSeededByENERGY");

    fTrackBlock->AddCollName("mu2e::CaloClusterCollection", "makeCaloCluster","");
    fTrackBlock->AddCollName("mu2e::KalRepCollection","TrkPatRec1","DownstreameMinus");
    fTrackBlock->AddCollName("mu2e::TrkToCaloExtrapolCollection",fTrkExtrapol.data(),"");
    fTrackBlock->AddCollName("mu2e::TrackClusterMatchCollection",fTrkCalMatch.data(),"");
    fTrackBlock->AddCollName("mu2e::StrawHitCollection","makeSH","");
    fTrackBlock->AddCollName("mu2e::PtrStepPointMCVectorCollection",
			     "makeSH","StrawHitMCPtr");
  }


  //-----------------------------------------------------------------------------
  void TCalm003::endJob() {
    TDirectory* old_dir = gDirectory;
    
    TFile* f = TFile::Open(fHistFileName.data(),"recreate");

    SaveFolder(fFolder,f);
    f->Write();
    f->Close();
    
    old_dir->cd();
  }

//-----------------------------------------------------------------------------
  void TCalm003::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track) {

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

    Hist->fPdgCode->Fill(Track->fPdgCode);
    Hist->fFrGH->Fill(Track->fNGoodMcHits/(Track->fNActive+1.e-5));
  }

//-----------------------------------------------------------------------------
  void TCalm003::FillEventHistograms(EventHist_t* Hist) {
    
    double   cos_th, x1, y1, r1, x2, y2, r2;
    
    cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();
    
    Hist->fEleCosTh->Fill(cos_th);
    
    int ntrk = fTrackBlock->NTracks();

    Hist->fNTracks->Fill(ntrk);

    for (int ipl=0; ipl<kNPlanes; ipl++) {
      for (int it=0; it<ntrk; it++) {
	x1 = fPos[it][ipl].x();
	y1 = fPos[it][ipl].y();
	r1 = sqrt(x1*x1+y1*y1);
	Hist->fRTrack[ipl]->Fill(r1);
      }
    }

    for (int i1=0; i1<kNPlanes; i1++) {
      for (int it=0; it<ntrk; it++) {
	x1 = fPos[it][i1].x();
	y1 = fPos[it][i1].y();
	r1 = sqrt(x1*x1+y1*y1);
	
	for (int i2=0; i2<kNPlanes; i2++) {
	  x2 = fPos[it][i2].x();
	  y2 = fPos[it][i2].y();
	  r2 = sqrt(x2*x2+y2*y2);

	  if (r1 < 360.) Hist->fRTrack0[i1][i2]->Fill(r2);
	  else           Hist->fRTrack1[i1][i2]->Fill(r2);

	  if (r1 < 330.) Hist->fRTrack330[i1][i2]->Fill(r2);
	  else           Hist->fRTrack331[i1][i2]->Fill(r2);
	}
      }
    }
  }

//-----------------------------------------------------------------------------
  void TCalm003::FillHistograms() {
//-----------------------------------------------------------------------------
// event histograms, everything is in mm!!!
//-----------------------------------------------------------------------------
    FillEventHistograms(fHist.fEvent[0]);

    if ((fNTracks   > 0) && (fRDisk1[0] < 360.))  { 
      FillEventHistograms(fHist.fEvent[1]);
    }

    if ((fNTracks   > 0) && (fRDisk1[0] < 330.))  { 
      FillEventHistograms(fHist.fEvent[2]);
    }
//-----------------------------------------------------------------------------
// track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
    TStnTrack*   trk;
    int ntrk = fTrackBlock->NTracks();

    for (int i=0; i<ntrk; ++i ) {
      trk = fTrackBlock->Track(i);
      FillTrackHistograms(fHist.fTrack[0],trk);
    }
  }

//-----------------------------------------------------------------------------
  void TCalm003::getData(art::Event* Evt) {
//-----------------------------------------------------------------------------
// generated particles
//-----------------------------------------------------------------------------
    art::Handle<GenParticleCollection> genpHandle;
    Evt->getByLabel("generate",genpHandle);
    fListOfGenParticles = (GenParticleCollection*) &(*genpHandle);
    fEle                = (GenParticle*) &fListOfGenParticles->at(0);

    // art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    // Evt->getByLabel(fStrawHitMaker,"StrawHitMCPtr",mcptrHandle);
    // fListOfMcStrawHits = (PtrStepPointMCVectorCollection*) &(*mcptrHandle);
//-----------------------------------------------------------------------------
// reconstructed tracks
//-----------------------------------------------------------------------------
    art::Handle<KalRepPtrCollection> krepsHandle;
    Evt->getByLabel("TrkPatRec1","DownstreameMinus", krepsHandle);
    fListOfTracks[0] = (KalRepPtrCollection*) krepsHandle.product();
    fNTracks     [0] = fListOfTracks[0]->size();
    fTrack           = 0;
    if (fNTracks[0]) fTrack = (KalRep*) &fListOfTracks[0]->at(0);
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
    art::Handle<CaloClusterCollection> calo_cluster_handle;
    Evt->getByLabel("makeCaloCluster","AlgoCLOSESTSeededByENERGY",calo_cluster_handle);
    fListOfClusters  = (CaloClusterCollection*) &(*calo_cluster_handle);
    fNClusters       = fListOfClusters->size();
//     fCluster         = 0;
//     if (fNClusters > 0) fCluster = &fListOfClusters->at(0);
//-----------------------------------------------------------------------------
// straw hits 
//-----------------------------------------------------------------------------
//     art::Handle<StrawHitCollection> shHandle;
//     Evt->getByLabel("makeSH",shHandle);
//     fListOfStrawHits = (StrawHitCollection*) &(*shHandle);
//     fNStrawHits      = fListOfStrawHits->size();
//-----------------------------------------------------------------------------
// track extrapolation results
//-----------------------------------------------------------------------------
//     art::Handle<TrkToCaloExtrapolCollection>  texHandle;
//     Evt->getByLabel("trkExtrapol",texHandle);
//     fListOfExtrapolatedTracks = (TrkToCaloExtrapolCollection*) &(*texHandle);

    StntupleInitMu2eClusterBlock(fClusterBlock,Evt,0);
    StntupleInitMu2eTrackBlock  (fTrackBlock  ,Evt,0);
  }



//-----------------------------------------------------------------------------
  void TCalm003::Init(art::Event* Evt) {
    const char oname[] = "TCalm003::Init";

    //    TStnCluster*    cluster;
    int             id_word, ntrk;
    TStnTrack*      track;
    KalRep*         krep;
    double          slast, x, y, s1, s2, z1, z2, splane;
    TrkDifTraj*     traj;

    double          trj_len[50];
    HepPoint        pos    [50];
    TrkErrCode      trk_rc;
//-----------------------------------------------------------------------------
// initialize tracks and determine track quality
//-----------------------------------------------------------------------------
    fNGoodTracks    = 0;
    fNMatchedTracks = 0;

    ntrk = fTrackBlock->NTracks();

    for (int it=0; it<ntrk; it++) {
      track          = fTrackBlock->Track(it);
      id_word        = fTrackID->IDWord(track);
      track->fIDWord = id_word;
      krep           = track->fKalRep[0];
//-----------------------------------------------------------------------------
// optimize acceptance - use all reconstructed tracks
// extrapolate track to a given distance - 
//-----------------------------------------------------------------------------
      slast      = krep->lastHit ()->kalHit()->hit()->fltLen();
      traj       = &krep->traj();

      for (int j=0; j<50; j++) {
	trj_len[j] = slast+j*100.;
	trk_rc = krep->extendThrough(trj_len[j]);
	if (trk_rc.success() != 1) {
//-----------------------------------------------------------------------------
// apparently, an extrapolation error, handle it
//-----------------------------------------------------------------------------
	  printf(">>> ERROR in [%s]: extrapolation it = %i iplane = %i\n",
		 oname, it, j);
	  pos    [j] = HepPoint(0.,0.,0);
	}
	else {
	  pos    [j] = traj->position(trj_len[j]);
	}
      }
//-----------------------------------------------------------------------------
// now, for each plane find 2 points before and after it
//-----------------------------------------------------------------------------
      int jfirst = -1;
      for (int ipl=0; ipl<kNPlanes; ipl++) {
	for (int j=0; j<50; j++) {
	  if (fZPlane[ipl] < pos[j].z()) {
	    jfirst = j;
	    break;
	  }
	}
//-----------------------------------------------------------------------------
// hopefully, jfirst found
//-----------------------------------------------------------------------------
	z1 = pos[jfirst-1].z();
	z2 = pos[jfirst  ].z();

	s1 = trj_len[jfirst-1];
	s2 = trj_len[jfirst  ];

	splane = s1+(fZPlane[ipl]-z1)/(z2-z1)*(s2-s1);

	fPos[it][ipl] = traj->position(splane);
      }
//-----------------------------------------------------------------------------
// track position on the first disk
//-----------------------------------------------------------------------------
      for (int j=0; j<50; j++) {
	if (pos[j].z() > fZDisk1) {
	  jfirst = j;
	  break;
	}
      }
//-----------------------------------------------------------------------------
// hopefully, jfirst found
//-----------------------------------------------------------------------------
      z1 = pos[jfirst-1].z();
      z2 = pos[jfirst  ].z();

      s1 = trj_len[jfirst-1];
      s2 = trj_len[jfirst  ];

      splane        = s1+(fZDisk1-z1)/(z2-z1)*(s2-s1);

      fPosDisk1[it] = traj->position(splane);
      x             = fPosDisk1[it].x();
      y             = fPosDisk1[it].y();
      fRDisk1[it]   = sqrt(x*x+y*y);
    }
  }



//-----------------------------------------------------------------------------
  void TCalm003::Debug(art::Event* Evt) {
    const char* oname = "TCalm003::Debug";
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
	if ((pt > 80.) && (t->fNActive > 35) && (t->NClusters() == 0)) {
	  printf(" >>>>>>> [%s] EVENT : %10i ERROR:000 TRACK Pt = %10.3f doesnt have a cluster\n",
		 oname,Evt->event(),pt);
	  //	  t->Print("");
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 001: events with a  track, clusters, but no match
//-----------------------------------------------------------------------------
    if (TModule::fDebugBit[1] != 0) {
      nt = fTrackBlock->NTracks();
      for (int i=0; i<nt; i++) {
	t = fTrackBlock->Track(i);
	if ((t->fIDWord == 0) && (fNClusters > 0) && (t->NClusters() == 0)) {
	  printf(" >>>>>>> [%s] EVENT : %10i BIT:001 TRACK with ID_WORD = 0 doesnt have a cluster\n",
		 oname,Evt->event());
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
		   oname,Evt->event(),dt);
	  }
	  //	  t->Print("");
	}
      }
    }

  }

//-----------------------------------------------------------------------------
  bool TCalm003::filter(art::Event& Evt) {
    const char* oname = "TCalm003::analyze";

    bool rc(true);

    printf(" >>>>>>> [%s] EVENT : %10i\n",oname,Evt.event());

    getData(&Evt);
    Init   (&Evt);
    Debug  (&Evt);

    FillHistograms();

    if (TModule::fDebugBit[51] != 0) {
      rc = (fNMatchedTracks > 0);
    }

    TModule::filter(Evt);
    
    return rc;
  }
  
}  // end namespace mu2e

using mu2e::TCalm003;

DEFINE_ART_MODULE(TCalm003);
