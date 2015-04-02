///////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"

#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
#include "ana/TValCalPatRecModule.hh"

ClassImp(TValCalPatRecModule)
//-----------------------------------------------------------------------------
TValCalPatRecModule::TValCalPatRecModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fMinT0   = 0;                       // do not cut on time by default
  fTrackID = new TStnTrackID();
}

//-----------------------------------------------------------------------------
TValCalPatRecModule::~TValCalPatRecModule() {
}


//-----------------------------------------------------------------------------
void TValCalPatRecModule::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP          ,"p"        ,Form("%s: Track P(total)[1]" ,Folder),1000,   0  ,200. ,Folder);
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

  HBook1F(Hist->fPdgCode    ,"pdg"      ,Form("%s: track PDG code"    ,Folder), 100,-50,50,Folder);
  HBook1F(Hist->fFrGH       ,"fgh"      ,Form("%s: Fraction Goog Hits",Folder), 100, 0,1,Folder);

  HBook2F(Hist->fNEPlVsNHPl ,"nep_vs_nhp",Form("%s: Track NEXP vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
  HBook2F(Hist->fNDPlVsNHPl ,"ndp_vs_nhp",Form("%s: Track NDIF vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
  HBook2F(Hist->fChi2dVsNDPl,"chi2d_vs_ndp",Form("%s: Track Chi2/Dof vs NDP",Folder), 30, 0,30,100,0.,10,Folder);
  HBook2F(Hist->fDpFVsNDPl  ,"dpf_vs_ndp"  ,Form("%s: Track DpF vs NDP",Folder)     , 30, 0,30,100,-5,5,Folder);
}

//-----------------------------------------------------------------------------
void TValCalPatRecModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
  HBook1F(Hist->fRv         ,"rv"       ,Form("%s: R(Vertex)"                      ,Folder), 100, 0, 1000,Folder);
  HBook1F(Hist->fZv         ,"zv"       ,Form("%s: Z(Vertex)"                      ,Folder), 300, 0,15000,Folder);
  HBook1F(Hist->fNClusters ,"ncl"       ,Form("%s: Number of Reconstructed Clusters",Folder),200,0,200,Folder);
  HBook2F(Hist->fNT2VsNT1  ,"nt2_vs_nt1",Form("%s: NT2 vs NT1"  ,Folder),3,0,3,3,0,3,Folder);
  HBook2F(Hist->fNGT2VsNGT1  ,"ngt2_vs_ngt1",Form("%s: NGT2 vs NGT1"  ,Folder),3,0,3,3,0,3,Folder);
}

//_____________________________________________________________________________
void TValCalPatRecModule::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[0] = 1;		// all events

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
// book Track1 histograms
//-----------------------------------------------------------------------------
  int book_track1_histset[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) book_track1_histset[i] = 0;

  book_track1_histset[  0] = 1;		// all tracks
  book_track1_histset[  1] = 1;		// all tracks NTracks2      = 0
  book_track1_histset[  2] = 1;		// all tracks NTracks2      > 0
  book_track1_histset[  3] = 1;		// all tracks NClusters(70) = 0
  book_track1_histset[  4] = 1;		// all "Set C" tracks
  book_track1_histset[  5] = 1;		// all "Set C" tracks NTracks2      = 0
  book_track1_histset[  6] = 1;		// all "Set C" tracks NTracks2      = 1
  book_track1_histset[  7] = 1;		// all "Set C" tracks NClusters(70) = 0

  for (int i=0; i<kNTrackHistSets; i++) {
    if (book_track1_histset[i] != 0) {
      sprintf(folder_name,"trk1_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrack1[i] = new TrackHist_t;
      BookTrackHistograms(fHist.fTrack1[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book Track2 histograms
//-----------------------------------------------------------------------------
  int book_track2_histset[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) book_track2_histset[i] = 0;

  book_track2_histset[  0] = 1;		// all tracks
  book_track2_histset[  1] = 1;		// all tracks NTracks1      = 0
  book_track2_histset[  2] = 1;		// all tracks NTracks1      = 1
  book_track2_histset[  3] = 1;		// all tracks NClusters(70) = 0
  book_track2_histset[  4] = 1;		// all "Set C" tracks
  book_track2_histset[  5] = 1;		// all "Set C" tracks NTracks1      = 0
  book_track2_histset[  6] = 1;		// all "Set C" tracks NTracks1      = 1
  book_track2_histset[  7] = 1;		// all "Set C" tracks NClusters(70) = 0

  for (int i=0; i<kNTrackHistSets; i++) {
    if (book_track2_histset[i] != 0) {
      sprintf(folder_name,"trk2_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrack2[i] = new TrackHist_t;
      BookTrackHistograms(fHist.fTrack2[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
// need MC truth branch
//-----------------------------------------------------------------------------
void TValCalPatRecModule::FillEventHistograms(EventHist_t* Hist) {
  double   cos_th, xv, yv, rv, zv;

  cos_th = -1; // fEle->momentum().pz()/fEle->momentum().vect().mag();

  xv = -1.; // fEle->position().x()+3904.;
  yv = -1.; // fEle->position().y();
    
  rv = sqrt(xv*xv+yv*yv);
  zv = -1.; // fEle->position().z();

  Hist->fEleCosTh->Fill(cos_th);
  Hist->fRv->Fill(rv);
  Hist->fZv->Fill(zv);

  Hist->fNClusters->Fill(fNClusters);
  Hist->fNT2VsNT1->Fill(fNTracks[0],fNTracks[1]);
  Hist->fNGT2VsNGT1->Fill(fNGoodTracks[0],fNGoodTracks[1]);
}


//-----------------------------------------------------------------------------
void TValCalPatRecModule::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track, TrackPar_t* Tp) {

  double dp, dpfst, chi2c;

  Hist->fP->Fill (Track->fP);
  Hist->fPt->Fill(Track->fPt);
  Hist->fPFront->Fill(Track->fPFront);
  Hist->fPStOut->Fill(Track->fPStOut);

  dp    = Track->fP-Track->fPFront;
  dpfst = Track->fPFront-Track->fPStOut;

  Hist->fDpFront->Fill(dp);
  Hist->fDpFSt->Fill(dpfst);

  Hist->fCosTh->Fill(Track->Momentum()->CosTheta());
  Hist->fChi2->Fill (Track->fChi2);
  Hist->fNDof->Fill(Track->NActive()-5.);
  Hist->fChi2Dof->Fill(Track->fChi2/(Track->NActive()-5.));
  Hist->fNActive->Fill(Track->NActive());
  Hist->fT0->Fill(Track->fT0);
  printf("TValCalPatRecModule::FillTrackHistograms: track charge is not defined yet\n");
  Hist->fQ->Fill(-10);
  Hist->fFitCons[0]->Fill(Track->fFitCons);
  Hist->fFitCons[1]->Fill(Track->fFitCons);

  Hist->fD0->Fill(Track->fD0);
  Hist->fZ0->Fill(Track->fZ0);
  Hist->fTanDip->Fill(Track->fTanDip);

  chi2c = Track->fChi2C/(Track->NActive()-5.);
  Hist->fChi2DofC->Fill(chi2c);

  int nh, nst_with_nh[10];

  for (int i=0; i<10; i++) nst_with_nh[i] = 0;

  for (int i=0; i<40; i++) {
    Hist->fNHVsStation->Fill(i,Track->fNHPerStation[i]);
    nh = Track->fNHPerStation[i];
    if ((nh < 0) || (nh >= 10)) {
      //      GetHeaderBlock()->Print(Form("ERROR: nh = %i\n",nh));
    }
    else {
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
    ep = vr->fEnergy/Track->fP;
    Hist->fEp->Fill(ep);
    Hist->fDx->Fill(vr->fDx);
    Hist->fDy->Fill(vr->fDy);
    Hist->fDz->Fill(vr->fDz);
  }

  int ncl = Track->NClusters();
  Hist->fNClusters->Fill(ncl);

  Hist->fPdgCode->Fill(Track->fPdgCode);
  Hist->fFrGH->Fill(Track->fNGoodMcHits/(Track->NActive()+1.e-5));

  Hist->fNEPlVsNHPl ->Fill(Tp->fNEPl,Tp->fNHPl);
  Hist->fNDPlVsNHPl ->Fill(Tp->fNDPl,Tp->fNHPl);
  Hist->fChi2dVsNDPl->Fill(Tp->fNDPl,Track->Chi2Dof());
  Hist->fDpFVsNDPl  ->Fill(Tp->fNDPl,Tp->fDpF);
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TValCalPatRecModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("TrackBlock"    ,"TStnTrackBlock"   ,&fTrackBlock[0]);
  RegisterDataBlock("TrackBlock2"   ,"TStnTrackBlock"   ,&fTrackBlock[1]);
  RegisterDataBlock("ClusterBlock"  ,"TStnClusterBlock" ,&fClusterBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// initialize track identification
//-----------------------------------------------------------------------------
  fTrackID->SetMinT0(fMinT0);

  return 0;
}


//_____________________________________________________________________________
int TValCalPatRecModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
void TValCalPatRecModule::FillHistograms() {

  //  double       cos_th (-2.),  cl_e(-1.);
  //  int          disk_id(-1);
  //  TStnCluster  *cl0;
  TStnTrack    *trk;
  TrackPar_t*  tp;
//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);

//-----------------------------------------------------------------------------
// track1 histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<fNTracks[0]; ++i ) {
    trk = fTrackBlock[0]->Track(i);
    tp  = fTrackPar[0]+i;
    FillTrackHistograms(fHist.fTrack1[0],trk,tp);

    if (fNTracks[1] == 0) FillTrackHistograms(fHist.fTrack1[1],trk,tp);
    if (fNTracks[1] >  0) FillTrackHistograms(fHist.fTrack1[2],trk,tp);
    if (fNTracks[1] >  0) FillTrackHistograms(fHist.fTrack1[3],trk,tp);

    if (trk->fIDWord == 0) {
					// track passes selection "C" 

      FillTrackHistograms(fHist.fTrack1[4],trk,tp);

      if (fNTracks[1] == 0) FillTrackHistograms(fHist.fTrack1[5],trk,tp);
      if (fNTracks[1] >  0) FillTrackHistograms(fHist.fTrack1[6],trk,tp);
      if (fNTracks[1] >  0) FillTrackHistograms(fHist.fTrack1[7],trk,tp);
    }
  }
//-----------------------------------------------------------------------------
// track2 histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<fNTracks[1]; ++i ) {
    trk = fTrackBlock[1]->Track(i);
    tp  = fTrackPar[1]+i;
    FillTrackHistograms(fHist.fTrack2[0],trk,tp);

    if (fNTracks[0] == 0) FillTrackHistograms(fHist.fTrack2[1],trk,tp);
    if (fNTracks[0] >  0) FillTrackHistograms(fHist.fTrack2[2],trk,tp);
    if (fNTracks[0] >  0) FillTrackHistograms(fHist.fTrack2[3],trk,tp);

    if (trk->fIDWord == 0) {
					// track passes selection "C" 

      FillTrackHistograms(fHist.fTrack2[4],trk,tp);

      if (fNTracks[0] == 0) FillTrackHistograms(fHist.fTrack2[5],trk,tp);
      if (fNTracks[0] >  0) FillTrackHistograms(fHist.fTrack2[6],trk,tp);
      if (fNTracks[0] >  0) FillTrackHistograms(fHist.fTrack2[7],trk,tp);
    }
  }
}



//_____________________________________________________________________________
int TValCalPatRecModule::Event(int ientry) {

  //  double                xs;
  //  TEmuLogLH::CalData_t  dat;
  TStnTrack*            track;
  int                   id_word;
  TrackPar_t*  tp;

  //  TDiskCalorimeter::GeomData_t disk_geom;

  fTrackBlock[0]->GetEntry(ientry);
  fTrackBlock[1]->GetEntry(ientry);
  fClusterBlock ->GetEntry(ientry);


  fNTracks[0] = fTrackBlock[0]->NTracks();
  fNTracks[1] = fTrackBlock[1]->NTracks();

  fNClusters  = fClusterBlock->NClusters();

  for (int iblock=0; iblock<2; iblock++) {
    
    fNTracks       [iblock] = fTrackBlock[iblock]->NTracks();
    fNGoodTracks   [iblock] = 0;
    fNMatchedTracks[iblock] = 0;


    int ntrk = fNTracks[iblock];
    for (int i=0; i<ntrk; i++) {
      					// assume less 20 tracks

      tp             = fTrackPar[iblock]+i;
      track          = fTrackBlock[iblock]->Track(i);

      id_word        = fTrackID->IDWord(track);
      track->fIDWord = id_word;
      if (id_word == 0) {
	fNGoodTracks[iblock] += 1;
	if ((track->fVMaxEp != NULL) && (fabs(track->fVMaxEp->fDt) < 2.5)) {
	  fNMatchedTracks[iblock] += 1;
	}
      }
//-----------------------------------------------------------------------------
// process hit masks
//-----------------------------------------------------------------------------
      int i1, i2, n1(0) ,n2(0), ndiff(0);
      int nbits = track->fHitMask.GetNBits();
      for (int i=0; i<nbits; i++) {
	i1 = track->HitMask()->GetBit(i);
	i2 = track->ExpectedHitMask()->GetBit(i);
	n1 += i1;
	n2 += i2;
	if (i1 != i2) ndiff += 1;
      }
//-----------------------------------------------------------------------------
// define additional parameters
//-----------------------------------------------------------------------------
      tp->fNHPl = n1;
      tp->fNEPl = n2;
      tp->fNDPl = ndiff;

      tp->fDpF   = track->fP     -track->fPFront;
      tp->fDp0   = track->fP0    -track->fPFront;
      tp->fDp2   = track->fP2    -track->fPFront;
      tp->fDpFSt = track->fPFront-track->fPStOut;
    }
  }
    
  fNClusters = fClusterBlock->NClusters();
  
  FillHistograms();

  Debug();

  return 0;		       
}


//_____________________________________________________________________________
int TValCalPatRecModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TValCalPatRecModule::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}


//-----------------------------------------------------------------------------
void TValCalPatRecModule::Debug() {
  TrackPar_t*  tp;
//-----------------------------------------------------------------------------
// bit 5: events with N(tracks) > 1
//-----------------------------------------------------------------------------
  if (GetDebugBit(5) == 1) {
    if (fNTracks[0] > 1) {
      GetHeaderBlock()->Print(Form("NTracks = %5i",fNTracks[0]));
    }
  }
//-----------------------------------------------------------------------------
// bit 6: events with large mismatch between the Nexp and Nobs intersections
//        for  CalPatRec tracks
//-----------------------------------------------------------------------------
  if (GetDebugBit(6) == 1) {

    int iblock = 1;
    int ntrk    = fNTracks[iblock];

    for (int i=0; i<ntrk; i++) {
      tp             = fTrackPar[iblock]+i;

      if (tp->fNHPl - tp->fNEPl > 10) {
	GetHeaderBlock()->Print(Form("it = %2i N(hit planes) = %3i N(expected): %3i",
				     i, tp->fNHPl,tp->fNEPl));
      }
    }
  }
}
