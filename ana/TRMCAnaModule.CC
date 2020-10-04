//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
// 
//  3  : events with set C tracks and 70mm < |dx|  < 90 mm
//  4  : events with DpF > 1 MeV : obviously, misreconstructed ones
//  5  : events with N(tracks) > 1
//  6  : events trk_41 with 0.8< E/P < 1.1 - tracks missed by CalPatRec
//  7  : events (muo) with LogLHRCal >   20
//  8  : events (ele) with LogLHRCal < - 20
//  9  : events (muo) with 0.42 < E/P < 0.46
// 10  : events (muo) with Set C track with ECL > 80 MeV
// 28  : Set C DEM tracks with E/P > 1.1
// 29  : TRK_19 (Set C DEM tracks with a cluster) and LLHR(cal) < 0
// 31  : EVT_6 events with ce_costh > 0.8 
// 32  : TRK_1 events with chi2tcm > 100. 
// 33  : DU < -80mm - study edge effects
// 34  : EVT_7: events with E_CL > 60 and no tracks (makes sense only for single CE events)
// 35  : TRK_1: events with P > 106 MeV/c - misreconstruction
// 36  : TRK_23 events with P < 80: odd misidentified muons - turned out to be DIO electrons
// 37  : TRK_26 LLHR_CAL > 5
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
#include "ana/TRMCAnaModule.hh"

ClassImp(TRMCAnaModule)

//-----------------------------------------------------------------------------
TRMCAnaModule::TRMCAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fPtMin  = 1.;
  fTrackNumber.Set(100);

  fFillDioHist     = 1;

  fMinT0 = 0; // do not cut on time by default

  fTrackID = new TStnTrackID();
  fLogLH   = new TEmuLogLH();
//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
//-----------------------------------------------------------------------------
  fPdgCode       = 11;
  fGeneratorCode = 2;			// conversionGun, 28:StoppedParticleReactionGun
}

//-----------------------------------------------------------------------------
TRMCAnaModule::~TRMCAnaModule() {
}


//-----------------------------------------------------------------------------
void TRMCAnaModule::BookGenpHistograms(GenpHist_t* HistR, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP      ,"p"       ,Form("%s: Momentum"     ,Folder),1000,     0, 200,Folder);
  HBook1F(Hist->fPdgCode[0],"pdg_code_0",Form("%s: PDG Code[0]"     ,Folder),200, -100, 100,Folder);
  HBook1F(Hist->fPdgCode[1],"pdg_code_1",Form("%s: PDG Code[1]"     ,Folder),500, -2500, 2500,Folder);
  HBook1F(Hist->fGenID  ,"gen_id"  ,Form("%s: Generator ID" ,Folder), 100,     0, 100,Folder);
  HBook1F(Hist->fZ0     ,"z0"      ,Form("%s: Z0"           ,Folder), 500,  5400, 6400,Folder);
  HBook1F(Hist->fT0     ,"t0"      ,Form("%s: T0"           ,Folder), 200,     0, 2000,Folder);
  HBook1F(Hist->fR0     ,"r"       ,Form("%s: R0"           ,Folder), 100,     0,  100,Folder);
  HBook1F(Hist->fCosTh  ,"cos_th"  ,Form("%s: Cos(Theta)"   ,Folder), 200,   -1.,   1.,Folder);
}

//-----------------------------------------------------------------------------
void TRMCAnaModule::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP[0]       ,"p"        ,Form("%s: Track P(Z1)"       ,Folder), 800,  80  ,120. ,Folder);
  HBook1F(Hist->fP[1]       ,"p_1"      ,Form("%s: Track P(total)[1]" ,Folder), 100, 100  ,105. ,Folder);
  HBook1F(Hist->fP[2]       ,"p_2"      ,Form("%s: Track P(total)[1]" ,Folder),2000,   0  ,200. ,Folder);
  HBook1F(Hist->fP0         ,"p0"       ,Form("%s: Track P(Z0)"       ,Folder),1000,   0  ,200. ,Folder);
  HBook1F(Hist->fP2         ,"p2"       ,Form("%s: Track P(z=-1540)"  ,Folder),1000,   0  ,200. ,Folder);
//-----------------------------------------------------------------------------
  HBook1F(Hist->fFitMomErr  ,"momerr"   ,Form("%s: Track FitMomError" ,Folder), 200,   0  ,  1. ,Folder);
  HBook1F(Hist->fPFront     ,"pf"       ,Form("%s: Track P(front)   " ,Folder), 400,  90  ,110. ,Folder);
  HBook1F(Hist->fDpFront    ,"dpf"      ,Form("%s: Track P-P(front) " ,Folder),1000,  -5. ,  5. ,Folder);
  HBook1F(Hist->fXDpF       ,"xdpf"     ,Form("%s: DpF/momErr"        ,Folder),1000, -50. , 50. ,Folder);
  HBook1F(Hist->fDpFront0   ,"dp0f"     ,Form("%s: Track P0-P(front)" ,Folder),1000,  -5. ,  5. ,Folder);
  HBook1F(Hist->fDpFront2   ,"dp2f"     ,Form("%s: Track P2-P(front)" ,Folder),1000,  -5. ,  5. ,Folder);
  HBook1F(Hist->fPStOut     ,"pstout"   ,Form("%s: Track P(ST_Out)  " ,Folder), 400,  90. ,110. ,Folder);
  HBook1F(Hist->fDpFSt      ,"dpfst"    ,Form("%s: Track Pf-Psto"     ,Folder),1000,  -5  ,  5. ,Folder);
  HBook2F(Hist->fDpFVsZ1    ,"dpf_vs_z1",Form("%s: Track DPF Vs Z1"   ,Folder), 200, -2000.,0,200,-5.,5,Folder);

  HBook1F(Hist->fPt         ,"pt"       ,Form("%s: Track Pt"          ,Folder), 600, 75,95,Folder);
  HBook1F(Hist->fCosTh      ,"costh"    ,Form("%s: Track cos(theta)"  ,Folder), 100,-1,1  ,Folder);
  HBook1F(Hist->fChi2       ,"chi2"     ,Form("%s: Track chi2 total"  ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fChi2Dof    ,"chi2d"    ,Form("%s: track chi2/N(dof)" ,Folder), 500, 0, 10,Folder);

  HBook1F(Hist->fNActive    ,"nactv"    ,Form("%s: N(active)"         ,Folder), 200,  0  , 200 ,Folder);
  HBook1F(Hist->fNaFract    ,"nafr"     ,Form("%s: N(active fraction)",Folder), 110,  0.5,1.05 ,Folder);
  HBook1F(Hist->fDNa        ,"dna"      ,Form("%s: Nhits-Nactive"     ,Folder), 100, -0.5 ,99.5,Folder);
  HBook1F(Hist->fNWrong     ,"nwrng"    ,Form("%s: N(wrong drift sgn)",Folder), 100, 0,100,Folder);
  HBook1F(Hist->fNDoublets  ,"nd"       ,Form("%s: N(doublets)"       ,Folder),  50, 0, 50,Folder);
  HBook1F(Hist->fNadOverNd  ,"nad_nd"   ,Form("%s: Nad/N(doublets)"   ,Folder), 110, 0,1.1,Folder);
  HBook1F(Hist->fNSSD       ,"nssd"     ,Form("%s: N(SS doublets)"    ,Folder),  50, 0, 50,Folder);
  HBook1F(Hist->fNOSD       ,"nosd"     ,Form("%s: N(OS doublets)"    ,Folder),  50, 0, 50,Folder);
  HBook1F(Hist->fNdOverNa   ,"nd_na"    ,Form("%s: NDoublets/Nactive" ,Folder), 100, 0,0.5,Folder);
  HBook1F(Hist->fNssdOverNa ,"nssd_na"  ,Form("%s: NSSD/Nactive"      ,Folder), 100, 0,0.5,Folder);
  HBook1F(Hist->fNosdOverNa ,"nosd_na"  ,Form("%s: NOSD/Nactive"      ,Folder), 100, 0,0.5,Folder);
  HBook1F(Hist->fNZeroAmb   ,"nza"      ,Form("%s: N (Iamb = 0) hits" ,Folder), 100, 0,100,Folder);
  HBook1F(Hist->fNzaOverNa  ,"nza_na"   ,Form("%s: NZeroAmb/Nactive"  ,Folder), 100, 0,  1,Folder);
  HBook1F(Hist->fNMatActive ,"nma"      ,Form("%s: N (Mat Active"     ,Folder), 100, 0,100,Folder);
  HBook1F(Hist->fNmaOverNa  ,"nma_na"   ,Form("%s: NMatActive/Nactive",Folder), 200, 0,   2,Folder);
  HBook1F(Hist->fNBend      ,"nbend"    ,Form("%s: Nbend"             ,Folder), 100, 0,1000,Folder);

  HBook1F(Hist->fT0         ,"t0"       ,Form("%s: track T0"          ,Folder), 200, 0,2000,Folder);
  HBook1F(Hist->fT0Err      ,"t0err"    ,Form("%s: track T0Err"       ,Folder), 100, 0,  10,Folder);
  HBook1F(Hist->fQ          ,"q"        ,Form("%s: track Q"           ,Folder),   4,-2,   2,Folder);
  HBook1F(Hist->fFitCons[0] ,"fcon"     ,Form("%s: track fit cons [0]",Folder), 200, 0,   1,Folder);
  HBook1F(Hist->fFitCons[1] ,"fcon1"    ,Form("%s: track fit cons [1]",Folder), 1000, 0,   0.1,Folder);
  HBook1F(Hist->fD0         ,"d0"       ,Form("%s: track D0      "    ,Folder), 200,-200, 200,Folder);
  HBook1F(Hist->fZ0         ,"z0"       ,Form("%s: track Z0      "    ,Folder), 200,-2000,2000,Folder);
  HBook1F(Hist->fTanDip     ,"tdip"     ,Form("%s: track tan(dip)"    ,Folder), 200, 0.0 ,2.0,Folder);
  HBook1F(Hist->fRMax       ,"rmax"     ,Form("%s: track R(max)  "    ,Folder), 200, 0., 1000,Folder);
  HBook1F(Hist->fDtZ0       ,"dtz0"     ,Form("%s: T0_trk-T0_MC(Z=0)" ,Folder), 200, -10.0 ,10.0,Folder);
  HBook1F(Hist->fXtZ0       ,"xtz0"     ,Form("%s: DT(Z0)/sigT"       ,Folder), 200, -10.0 ,10.0,Folder);

  HBook1F(Hist->fResid      ,"resid"    ,Form("%s: hit residuals"     ,Folder), 500,-0.5 ,0.5,Folder);
  HBook1F(Hist->fAlgMask    ,"alg"      ,Form("%s: algorithm mask"    ,Folder),  10,  0, 10,Folder);

  HBook1F(Hist->fChi2Tcm  ,"chi2tcm"  ,Form("%s: chi2(t-c match)"   ,Folder), 250,  0  ,250 ,Folder);
  HBook1F(Hist->fChi2XY     ,"chi2xy"   ,Form("%s: chi2(t-c match) XY",Folder), 300,-50  ,250 ,Folder);
  HBook1F(Hist->fChi2T      ,"chi2t"    ,Form("%s: chi2(t-c match) T" ,Folder), 250,  0  ,250 ,Folder);

  HBook1F(Hist->fDt         ,"dt"       ,Form("%s: T(trk)-T(cl)"      ,Folder), 400,-20  ,20 ,Folder);
  HBook1F(Hist->fDx         ,"dx"       ,Form("%s: X(trk)-X(cl)"      ,Folder), 200,-500 ,500,Folder);
  HBook1F(Hist->fDy         ,"dy"       ,Form("%s: Y(trk)-Y(cl)"      ,Folder), 200,-500 ,500,Folder);
  HBook1F(Hist->fDz         ,"dz"       ,Form("%s: Z(trk)-Z(cl)"      ,Folder), 200,-250 ,250,Folder);
  HBook1F(Hist->fDu         ,"du"       ,Form("%s: track-cluster DU"  ,Folder), 250,-250 ,250,Folder);
  HBook1F(Hist->fDv         ,"dv"       ,Form("%s: track-cluster DV"  ,Folder), 200,-100 ,100,Folder);
  HBook1F(Hist->fPath       ,"path"     ,Form("%s: track sdisk"       ,Folder),  50,   0 ,500,Folder);

  HBook1F(Hist->fECl        ,"ecl"      ,Form("%s: cluster E"         ,Folder), 300, 0   ,150,Folder);
  HBook1F(Hist->fEClEKin    ,"ecl_ekin" ,Form("%s: cluster E/Ekin(mu)",Folder), 500, 0   ,5,Folder);
  HBook1F(Hist->fEp         ,"ep"       ,Form("%s: track E/P"         ,Folder), 300, 0   ,1.5,Folder);
  HBook1F(Hist->fDrDzCal    ,"drdzcal"  ,Form("%s: track dr/dz cal"   ,Folder), 200, -5  ,5  ,Folder);
  HBook1F(Hist->fDtClZ0     ,"dtclz0"   ,Form("%s: T(cl_z0)-T(Z0)"    ,Folder), 250, -5 , 5,Folder);
  HBook2F(Hist->fDtClZ0VsECl,"dtclz0_vs_ecl",Form("%s: DtClZ0 vs ECl" ,Folder), 100, 0 , 200, 250, -5 , 5,Folder);
  HBook2F(Hist->fDtClZ0VsP  ,"dtclz0_vs_p"  ,Form("%s: DtClZ0 vs p"   ,Folder), 100, 0 , 200, 250, -5 , 5,Folder);

  HBook2F(Hist->fFConsVsNActive,"fc_vs_na" ,Form("%s: FitCons vs NActive",Folder),  150, 0, 150, 200,0,1,Folder);
  HBook1F(Hist->fDaveTrkQual,"dtqual"   ,Form("%s:DaveTrkQual"        ,Folder), 200, -0.5, 1.5,Folder);
  HBook1F(Hist->fMVAOut     ,"mvaout"   ,Form("%s:MVAOut[0]"          ,Folder), 200, -0.5, 1.5,Folder);
  HBook1F(Hist->fDeltaMVA   ,"dmva"     ,Form("%s:MVAOut[0]-TrkQual"  ,Folder), 200, -1.0, 1.0,Folder);

}

//-----------------------------------------------------------------------------
void TRMCAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
  HBook1F(Hist->fEleMom    ,"ce_mom"   ,Form("%s: Conversion Electron Momentum"    ,Folder),1000,  0,200,Folder);
  HBook1D(Hist->fDioMom    ,"dio_mom"  ,Form("%s: DIO momentum"                    ,Folder),1000, 50,150,Folder);
  HBook1F(Hist->fRv         ,"rv"      ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 1000,Folder);
  HBook1F(Hist->fZv         ,"zv"      ,Form("%s: Z(Vertex)"                       ,Folder), 300, 0,15000,Folder);
  //  HBook1F(Hist->fNClusters ,"ncl"      ,Form("%s: Number of Reconstructed Clusters",Folder),200,0,200,Folder);
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
  HBook1F(Hist->fNGenp     ,"ngenp"    ,Form("%s: N(Gen Particles)"                ,Folder),500,0,500,Folder);

  // HBook1F(Hist->fNHitCrystalsTot,"nhcr_tot",Form("%s: NHit Crystals Tot",Folder), 100, 0,100,Folder);
  // HBook1F(Hist->fECal,"ecal",Form("%s: E(cal), sum over both disks",Folder), 500, 0,250,Folder);
  // HBook1F(Hist->fECalOverEKin,"ec_over_ek",Form("%s: E(cal)/E(kin)",Folder), 200, 0,2,Folder);
}

//-----------------------------------------------------------------------------
void TRMCAnaModule::BookSimpHistograms(SimpHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fPdgCode   ,"pdg"         ,Form("%s: PDG code"                     ,Folder),200,-100,100,Folder);
  HBook1F(Hist->fNStrawHits,"nsth"        ,Form("%s: n straw hits"                 ,Folder),200,   0,200,Folder);
  HBook1F(Hist->fMomTargetEnd    ,"ptarg" ,Form("%s: CE mom after Stopping Target" ,Folder),400,  90,110,Folder);
  HBook1F(Hist->fMomTrackerFront ,"pfront",Form("%s: CE mom at the Tracker Front"  ,Folder),400,  90,110,Folder);
}


//_____________________________________________________________________________
void TRMCAnaModule::BookHistograms() {

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

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;	        // events with a reconstructed track

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
// book simp histograms
//-----------------------------------------------------------------------------
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[ 0] = 1;		// all events

  for (int i=0; i<kNSimpHistSets; i++) {
    if (book_simp_histset[i] != 0) {
      sprintf(folder_name,"sim_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fSimp[i] = new SimpHist_t;
      BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track histograms
//-----------------------------------------------------------------------------
  int book_track_histset[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

  book_track_histset[  0] = 1;		// all tracks e-
  book_track_histset[  1] = 1;		// all tracks e- passing Set C cuts 

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
// book Genp histograms
//-----------------------------------------------------------------------------
  int book_genp_histset[kNGenpHistSets];
  for (int i=0; i<kNGenpHistSets; i++) book_genp_histset[i] = 0;

  book_genp_histset[0] = 1;		// all particles
//   book_genp_histset[1] = 1;		// all crystals, e > 0

  for (int i=0; i<kNGenpHistSets; i++) {
    if (book_genp_histset[i] != 0) {
      sprintf(folder_name,"gen_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fGenp[i] = new GenpHist_t;
      BookGenpHistograms(fHist.fGenp[i],Form("Hist/%s",folder_name));
    }
  }
}


//-----------------------------------------------------------------------------
// need MC truth branch
//-----------------------------------------------------------------------------
void TRMCAnaModule::FillEventHistograms(EventHist_t* Hist) {
  double            cos_th(-2), dio_wt(-1.), xv(-1.e6), yv(-1.e6), rv(-1.e6), zv(-1.e6), p(-1.);
  //  double            e, m, r;
  TLorentzVector    mom;

  if (fParticle) {
    fParticle->Momentum(mom);
    p      = mom.P();
    cos_th = mom.Pz()/p;
    xv     = fParticle->Vx()+3904.;
    yv     = fParticle->Vy();
    rv     = sqrt(xv*xv+yv*yv);
    zv     = fParticle->Vz();
    dio_wt = TStntuple::DioWeightAl(p);
  }

  Hist->fEleMom->Fill(p);
  Hist->fDioMom->Fill(p,dio_wt);
  Hist->fEleCosTh->Fill(cos_th);
  Hist->fRv->Fill(rv);
  Hist->fZv->Fill(zv);

  // Hist->fNClusters->Fill(fNClusters);
  Hist->fNTracks->Fill  (fNTracks[0]);
  // Hist->fNStrawHits[0]->Fill(fNStrawHits);
  // Hist->fNStrawHits[1]->Fill(fNStrawHits);

  double emax   = -1;
  double dt     = 9999.;

  Hist->fDtClT->Fill(dt);
  Hist->fEMax->Fill(emax);

  Hist->fNHyp->Fill(fNHyp);
  Hist->fBestHyp[0]->Fill(fBestHyp[0]);
  Hist->fBestHyp[1]->Fill(fBestHyp[1]);
  Hist->fNGenp->Fill(fNGenp);
}

//-----------------------------------------------------------------------------
void TRMCAnaModule::FillGenpHistograms(GenpHist_t* HistR, TGenParticle* Genp) {
  int    gen_id;
  float  p, cos_th, z0, t0, r0, x0, y0;

  TLorentzVector mom, v;

  Genp->Momentum(mom);
  //  Genp->ProductionVertex(v);

  p      = mom.P();
  cos_th = mom.CosTheta();

  x0     = Genp->Vx()+3904.;
  y0     = Genp->Vy();

  z0     = Genp->Vz();
  t0     = Genp->T();
  r0     = sqrt(x0*x0+y0*y0);
  gen_id = Genp->GetStatusCode();

  Hist->fPdgCode[0]->Fill(Genp->GetPdgCode());
  Hist->fPdgCode[1]->Fill(Genp->GetPdgCode());
  Hist->fGenID->Fill(gen_id);
  Hist->fZ0->Fill(z0);
  Hist->fT0->Fill(t0);
  Hist->fR0->Fill(r0);
  Hist->fP->Fill(p);
  Hist->fCosTh->Fill(cos_th);
}

//-----------------------------------------------------------------------------
void TRMCAnaModule::FillSimpHistograms(SimpHist_t* Hist, TSimParticle* Simp) {

  Hist->fPdgCode->Fill(Simp->fPdgCode);
  Hist->fMomTargetEnd->Fill(Simp->fMomTargetEnd);
  Hist->fMomTrackerFront->Fill(Simp->fMomTrackerFront);
  Hist->fNStrawHits->Fill(Simp->fNStrawHits);
}

//-----------------------------------------------------------------------------
// for DIO : ultimately, one would need to renormalize the distribution
//-----------------------------------------------------------------------------
void TRMCAnaModule::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track, TrackPar_t* Tp, double Weight) {

  TLorentzVector  mom;

  TrackHist_t* Hist = (TrackHist_t*) HistR;

					// Tp->fP - corrected momentum, fP0 and fP2 - not corrected
  Hist->fP[0]->Fill (Tp->fP,Weight);
  Hist->fP[1]->Fill (Tp->fP,Weight);
  Hist->fP[2]->Fill (Tp->fP,Weight);
					// fP0: momentum in the first point,  fP2 - in the last
  Hist->fP0->  Fill (Track->fP0,Weight);
  Hist->fP2->  Fill (Track->fP2,Weight);

  //  Hist->fPDio->Fill(Tp->fP,Tp->fDioWt);

  // Hist->fPlw->Fill   (Tp->fP, Tp->fLumWt);
  // Hist->fPDiolw->Fill(Tp->fP, Tp->fTotWt);

  Hist->fFitMomErr->Fill(Track->fFitMomErr,Weight);

  Hist->fPt    ->Fill(Track->fPt    , Weight);
  Hist->fPFront->Fill(Track->fPFront, Weight);
  Hist->fPStOut->Fill(Track->fPStOut, Weight);
//-----------------------------------------------------------------------------
// dp: Tracker-only resolution
//-----------------------------------------------------------------------------
  Hist->fDpFront ->Fill(Tp->fDpF   , Weight);
  Hist->fXDpF    ->Fill(Tp->fXDpF  , Weight);
  Hist->fDpFront0->Fill(Tp->fDp0   , Weight);
  Hist->fDpFront2->Fill(Tp->fDp2   , Weight);
  Hist->fDpFSt   ->Fill(Tp->fDpFSt , Weight);
  Hist->fDpFVsZ1 ->Fill(Track->fZ1 ,Tp->fDpF, Weight);

  Hist->fCosTh->Fill(Track->Momentum()->CosTheta(), Weight);
  Hist->fChi2->Fill (Track->fChi2, Weight);
  Hist->fChi2Dof->Fill(Track->fChi2/(Track->NActive()-5.), Weight);

  float na  = Track->NActive();
  float dna = Track->NHits()-na;

  Hist->fNActive->Fill(na, Weight);
  Hist->fNaFract->Fill(na/(Track->NHits()+0.), Weight);
  Hist->fDNa->Fill(dna, Weight);
  Hist->fNWrong->Fill(Track->NWrong(), Weight);

  float nd = Track->NDoublets();

  float nad = Track->NDoubletsAct();
  Hist->fNDoublets->Fill(nd, Weight);
  Hist->fNadOverNd->Fill(nad/nd, Weight);
  Hist->fNOSD->Fill(Track->NOSDoublets(), Weight);
  Hist->fNSSD->Fill(Track->NSSDoublets(), Weight);
  Hist->fNdOverNa->Fill(nd/na, Weight);
  Hist->fNosdOverNa->Fill(Track->NOSDoublets()/na , Weight);
  Hist->fNssdOverNa->Fill(Track->NSSDoublets()/na , Weight);
  Hist->fNZeroAmb  ->Fill(Track->NHitsAmbZero()   , Weight);
  Hist->fNzaOverNa ->Fill(Track->NHitsAmbZero()/na, Weight);

  int nma = Track->NMatActive();

  Hist->fNMatActive->Fill(nma   , Weight);
  Hist->fNmaOverNa ->Fill(nma/na, Weight);

  Hist->fNBend->Fill(Track->NBend(), Weight);

  Hist->fT0->Fill(Track->fT0, Weight);
  Hist->fT0Err->Fill(Track->fT0Err, Weight);
  Hist->fQ->Fill(Track->fCharge, Weight);
//-----------------------------------------------------------------------------
// the two histograms just have different limits
//-----------------------------------------------------------------------------
  Hist->fFitCons[0]->Fill(Track->fFitCons, Weight);
  Hist->fFitCons[1]->Fill(Track->fFitCons, Weight);

  Hist->fD0->Fill(Track->fD0, Weight);
  Hist->fZ0->Fill(Track->fZ0, Weight);
  Hist->fTanDip->Fill(Track->fTanDip, Weight);
  Hist->fDtZ0->Fill(Tp->fDtZ0, Weight);
  Hist->fXtZ0->Fill(Tp->fXtZ0, Weight);
  Hist->fRMax->Fill(Track->RMax(), Weight);
  
  Hist->fAlgMask->Fill(Track->AlgMask(), Weight);

  Hist->fChi2Tcm->Fill(Tp->fChi2Tcm, Weight);
  Hist->fChi2XY ->Fill(Tp->fChi2XY , Weight);
  Hist->fChi2T  ->Fill(Tp->fChi2T  , Weight);

  Hist->fDt->Fill(Tp->fDt, Weight);
  Hist->fDx->Fill(Tp->fDx, Weight);
  Hist->fDy->Fill(Tp->fDy, Weight);
  Hist->fDz->Fill(Tp->fDz, Weight);
  Hist->fDu->Fill(Tp->fDu, Weight);
  Hist->fDv->Fill(Tp->fDv, Weight);

  Hist->fPath->Fill(Tp->fPath, Weight);
  Hist->fECl ->Fill(Tp->fEcl , Weight);
//-----------------------------------------------------------------------------
// assume muon hypothesis
//-----------------------------------------------------------------------------
  double    ekin(-1.);
  if (fSimp) {
    double p, m;
    p    = Tp->fP;
    m    = 105.658; // muon mass, in MeV
    ekin = sqrt(p*p+m*m)-m;
  }

  Hist->fEClEKin->Fill(Tp->fEcl/ekin, Weight);
  Hist->fEp     ->Fill(Tp->fEp      , Weight);
  Hist->fDrDzCal->Fill(Tp->fDrDzCal , Weight);
  Hist->fDtClZ0 ->Fill(Tp->fDtClZ0  , Weight);

  Hist->fDtClZ0VsECl->Fill(Tp->fEcl,Tp->fDtClZ0, Weight);
  Hist->fDtClZ0VsP  ->Fill(Tp->fP  ,Tp->fDtClZ0, Weight);

  Hist->fFConsVsNActive->Fill(Track->NActive(),Track->fFitCons, Weight);
//-----------------------------------------------------------------------------
// MVA variables
//-----------------------------------------------------------------------------
  Hist->fDaveTrkQual->Fill(Track->DaveTrkQual(), Weight);
  Hist->fMVAOut     ->Fill(Tp->fMVAOut[0]      , Weight);

  float dmva=Tp->fMVAOut[0]-Track->DaveTrkQual();
  Hist->fDeltaMVA->Fill(dmva, Weight);
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TRMCAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("TrackBlockDar"       ,"TStnTrackBlock"      ,&fTrackBlock      );
  RegisterDataBlock("GenpBlock"           ,"TGenpBlock"          ,&fGenpBlock       );
  RegisterDataBlock("SimpBlock"           ,"TSimpBlock"          ,&fSimpBlock       );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// initialize likelihood histograms
//-----------------------------------------------------------------------------
  fTrackID->SetMinT0(fMinT0);
//-----------------------------------------------------------------------------
// PID initialization: read the likelihood templates
//-----------------------------------------------------------------------------
  fLogLH->Init("v5_7_0");

  return 0;
}


//_____________________________________________________________________________
int TRMCAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
void TRMCAnaModule::FillHistograms() {

  //  double       cos_th (-2.),  cl_e(-1.);
  // int          disk_id(-1); // , alg_mask, nsh, nactive;
  //  float        pfront, ce_pitch, reco_pitch, fcons, t0, sigt, sigp, p; 
  //  TStnCluster  *cl0;

  //  cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();

  // if (fNClusters > 0) {
  //   cl0     = fClusterBlock->Cluster(0);
  //   cl_e    = cl0->Energy();
  //   disk_id = cl0->DiskID();
  // }
  //-----------------------------------------------------------------------------
  // event histograms
  //-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);

  if (fNTracks[0]> 0) FillEventHistograms(fHist.fEvent[1]);
//-----------------------------------------------------------------------------
// Simp histograms
//-----------------------------------------------------------------------------
  if (fSimp) {
    FillSimpHistograms(fHist.fSimp[0],fSimp);
  }
//-----------------------------------------------------------------------------
// track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
  TStnTrack*   trk;
  //  TrackPar_t*  tp;

  for (int i=0; i<fNTracks[0]; ++i ) {
    trk = fTrackBlock->Track(i);
    //    tp  = fTrackPar+i;

    FillTrackHistograms(fHist.fTrack[0],trk);

    if (trk->fIDWord == 0) {
					// track passes selection "C" 

      FillTrackHistograms(fHist.fTrack[1],trk);
    }
  }
//-----------------------------------------------------------------------------
// fill GENP histograms
// GEN_0: all particles
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  for (int i=0; i<fNGenp; i++) {
    genp = fGenpBlock->Particle(i);
    FillGenpHistograms(fHist.fGenp[0],genp);
  }
  //  first_entry = 0;
}



//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TRMCAnaModule::Event(int ientry) {

  double                xs, p;
  TEmuLogLH::PidData_t  dat;
  TStnTrack*            track;
  int                   id_word;
  TLorentzVector        mom;

  //  TDiskCalorimeter::GeomData_t disk_geom;

  fTrackBlock  ->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fNGenp    = fGenpBlock->NParticles();

  TGenParticle* genp;
  int           pdg_code, generator_code;

  fParticle = NULL;
  for (int i=fNGenp-1; i>=0; i--) {
    genp           = fGenpBlock->Particle(i);
    pdg_code       = genp->GetPdgCode();
    generator_code = genp->GetStatusCode();
    if ((abs(pdg_code) == fPdgCode) && (generator_code == fGeneratorCode)) {
      fParticle = genp;
      break;
    }
  }
					// may want to revisit the definition of fSimp
  fSimp     = fSimpBlock->Particle(0);

  if (fParticle) {
    fParticle->Momentum(mom);
    p         = mom.P();
  }
  else p = 0.;

  fEleE     = sqrt(p*p+0.511*0.511);


  fNTracks[0] = fTrackBlock->NTracks();

  fNHyp       = -1;
  fBestHyp[0] = -1;
  fBestHyp[1] = -1;

  fNGoodTracks    = 0;
  fNMatchedTracks = 0;

  fNTracks[0] = fTrackBlock->NTracks();
  if (fNTracks[0] == 0) fTrack = 0;
  else                  fTrack = fTrackBlock->Track(0);

  int ntrk = fNTracks[0];

  TrackPar_t*   tp;

  for (int itrk=0; itrk<ntrk; itrk++) {
					// assume less 20 tracks
    tp             = fTrackPar+itrk;

    track          = fTrackBlock->Track(itrk);
    id_word        = fTrackID->IDWord(track);
    track->fIDWord = id_word;
    if (id_word == 0) {
      fNGoodTracks += 1;
      if ((track->fVMaxEp != NULL) && (fabs(track->fVMaxEp->fDt) < 2.5)) {
	fNMatchedTracks += 1;
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

    if (fFillDioHist == 0) tp->fDioWt = 1.;
    else                   tp->fDioWt = TStntuple::DioWeightAl(fEleE);
//-----------------------------------------------------------------------------
// track residuals
//-----------------------------------------------------------------------------
    TStnTrack::InterData_t*  vr = track->fVMaxEp; 
    double    nx, ny;

    tp->fEcl       = -1.e6;
    tp->fEp        = -1.e6;

    tp->fDu        = -1.e6;
    tp->fDv        = -1.e6;
    tp->fDx        = -1.e6;
    tp->fDy        = -1.e6;
    tp->fDz        = -1.e6;
    tp->fDt        = -1.e6;

    tp->fChi2Match = -1.e6;
    tp->fPath      = -1.e6;

    if (vr) {
      tp->fEcl = vr->fEnergy;
      tp->fEp  = tp->fEcl/track->fP;

      tp->fDx  = vr->fDx;
      tp->fDy  = vr->fDy;
      tp->fDz  = vr->fDz;
//-----------------------------------------------------------------------------
// v4_2_4: correct by additional 0.22 ns - track propagation by 6 cm
//-----------------------------------------------------------------------------
      tp->fDt  = vr->fDt - 0.22; // - 1.;

      nx  = vr->fNxTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);
      ny  = vr->fNyTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);

      tp->fDu        = vr->fDx*nx+vr->fDy*ny;
      tp->fDv        = vr->fDx*ny-vr->fDy*nx;
      tp->fChi2Match = vr->fChi2Match;
      tp->fPath      = vr->fPath;
    }

    //    if ((tp->fEp > 0) && (track->fEp > 0) && (fabs(tp->fEp-track->fEp) > 1.e-6)) {
    //      GetHeaderBlock()->Print(Form(" TRMCAnaModule ERROR: tp->fEp = %10.5f  track->fEp = %10.5f",tp->fEp,track->fEp));
    //    }
//-----------------------------------------------------------------------------
// PID likelihoods
//-----------------------------------------------------------------------------
    dat.fDt   = tp->fDt;
    dat.fEp   = tp->fEp;
    dat.fPath = tp->fPath;
      
    xs = track->XSlope();

    track->fEleLogLHCal = fLogLH->LogLHCal(&dat,11);
    track->fMuoLogLHCal = fLogLH->LogLHCal(&dat,13);

    double llhr_cal = track->fEleLogLHCal-track->fMuoLogLHCal;

    if (GetDebugBit(7)) {
      if ((id_word == 0) && (llhr_cal > 20)) {
	GetHeaderBlock()->Print(Form("bit:007: dt = %10.3f ep = %10.3f",track->Dt(),tp->fEp));
      }
    }

    if (GetDebugBit(8)) {
      if ((id_word == 0) && (llhr_cal < -20)) {
	GetHeaderBlock()->Print(Form("bit:008: p = %10.3f dt = %10.3f ep = %10.3f",
				     track->P(),track->Dt(),tp->fEp));
      }
    }

    track->fLogLHRXs    = fLogLH->LogLHRXs(xs);
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TRMCAnaModule::Debug() {

  TStnTrack* trk;
  //  TrackPar_t* tp;
  int ntrk = fTrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = fTrackBlock->Track(itrk);
    //    tp  = &fTrackPar[itrk];
//-----------------------------------------------------------------------------
// bit 3: Set C tracks with large DX : 70mm < |DX| < 90mm
//-----------------------------------------------------------------------------
    if (GetDebugBit(3) == 1) {
      if (trk->fIDWord == 0) {
	TStnTrack::InterData_t*    vr = trk->fVMaxEp; // residuals
	if ((vr && (fabs(vr->fDx) > 70) && (fabs(vr->fDx) < 90))) {
	  GetHeaderBlock()->Print(Form("large DX: %f",vr->fDx));
	}
      }
    }
  }
}

//_____________________________________________________________________________
int TRMCAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TRMCAnaModule::Test001() {
}
