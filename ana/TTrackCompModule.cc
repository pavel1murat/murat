//////////////////////////////////////////////////////////////////////////////
// Validation of CalPatREc: compare CalPatRec tracks to TrkPatRec ones
// assume running on stntuple with 3 tracking branches
//
// use of tmp:
//
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//  3  : events with CalPatRec tracks with DPF > 5
//  4  : events with TrkPatRec track and a 60+ MeV cluster, but with no CalPatRec track
//  5  : events with CalPatRec tracks with P > 106
//  6  : events with CalPatRec tracks with 1.5 < DPF < 5
//  7  : events with N(set "C" CalPatRec tracks)  > 0
//  8  : events with CalPatRec tracks with P > 105
//
// call: "track_comp(28,4)
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
#include "ana/TTrackCompModule.hh"

ClassImp(TTrackCompModule)
//-----------------------------------------------------------------------------
TTrackCompModule::TTrackCompModule(const char* name, const char* title):
  TStnModule    (name,title),
  fPdgCode      (11),                       // electron
  fGeneratorCode( 2)                        // 2:ConversionGun 28:StoppedParticleReactionGun
{
  fFillDioHist = 1;
//-----------------------------------------------------------------------------
// TrackID[0] : "SetC"
// i = 1..6 : cut on DaveTrkQual > 0.1*i instead
//-----------------------------------------------------------------------------
  fNID  = 7;
  for (int i=0; i<fNID; i++) {
    fTrackID[i] = new TStnTrackID();
    if (i > 0) {
      fTrackID[i]->SetMaxFitMomErr (100);
      fTrackID[i]->SetMaxT0Err     (100);
      //      fTrackID[i]->SetMinNActive   ( 10);
      fTrackID[i]->SetMinFitCons   (-1.);
      fTrackID[i]->SetMinTrkQual   (0.1*i);
    }
  }

  fBestTrackID = fTrackID[4];  // Dave's default: DaveTrkQual > 0.4
  fBestID      = 4;
  fLogLH       = new TEmuLogLH();
//-----------------------------------------------------------------------------
// debugging information
//-----------------------------------------------------------------------------
  fDebugCut[5].fXMin   = 106.;
  fDebugCut[5].fXMax   = 200.;

  fDebugCut[6].fXMin   = 1.5;
  fDebugCut[6].fXMax   = 10.0;
}

//-----------------------------------------------------------------------------
TTrackCompModule::~TTrackCompModule() {
  delete fLogLH;
  for (int i=0; i<fNID; i++) delete fTrackID[i];
}

//-----------------------------------------------------------------------------
void TTrackCompModule::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP[0]       ,"p"        ,Form("%s: Track P(Z1)"       ,Folder), 800,  80  ,120. ,Folder);
  HBook1F(Hist->fP[1]       ,"p_1"      ,Form("%s: Track P(total)[1]" ,Folder), 100, 100  ,105. ,Folder);
  HBook1F(Hist->fP[2]       ,"p_2"      ,Form("%s: Track P(total)[1]" ,Folder),2000,   0  ,200. ,Folder);
  HBook1F(Hist->fP0         ,"p0"       ,Form("%s: Track P(Z0)"       ,Folder),1000,   0  ,200. ,Folder);
  HBook1F(Hist->fP2         ,"p2"       ,Form("%s: Track P(z=-1540)"  ,Folder),1000,   0  ,200. ,Folder);

  HBook1D(Hist->fPDio       ,"pdio"     ,Form("%s: Track P(DIO WT)"   ,Folder), 800,  80  ,120. ,Folder);
  Hist->fPDio->Sumw2(kTRUE);
//-----------------------------------------------------------------------------
// luminosity-weighted distributions for signal and background
//-----------------------------------------------------------------------------
  HBook1D(Hist->fPlw        ,"plw"     ,Form("%s: Track P(Lumi-WT)"   ,Folder), 800,  80  ,120. ,Folder);
  Hist->fPlw->Sumw2(kTRUE);
  HBook1D(Hist->fPDiolw     ,"pdiolw"  ,Form("%s: Trk P WT(Lumi+DIO)" ,Folder), 800,  80  ,120. ,Folder);
  Hist->fPDiolw->Sumw2(kTRUE);
//-----------------------------------------------------------------------------
  HBook1F(Hist->fFitMomErr  ,"momerr"   ,Form("%s: Track FitMomError" ,Folder), 200,   0  ,  1. ,Folder);
  HBook1F(Hist->fPFront     ,"pf"       ,Form("%s: Track P(front)   " ,Folder), 400,  90  ,110. ,Folder);
  HBook1F(Hist->fDpFront    ,"dpf"      ,Form("%s: Track P-P(front) " ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(Hist->fXDpF       ,"xdpf"     ,Form("%s: DpF/momErr"        ,Folder),1000, -50. , 50. ,Folder);
  HBook1F(Hist->fDpFDio     ,"dpfdio"   ,Form("%s: Track DpF(DIO Wt)" ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(Hist->fDpFront0   ,"dp0f"     ,Form("%s: Track P0-P(front)" ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(Hist->fDpFront2   ,"dp2f"     ,Form("%s: Track P2-P(front)" ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(Hist->fPStOut     ,"pstout"   ,Form("%s: Track P(ST_Out)  " ,Folder), 400,  90. ,110. ,Folder);
  HBook1F(Hist->fDpFSt      ,"dpfst"    ,Form("%s: Track Pf-Psto"     ,Folder), 200,  -5  ,  5. ,Folder);
  HBook2F(Hist->fDpFVsZ1    ,"dpf_vs_z1",Form("%s: Track DPF Vs Z1"   ,Folder), 200, -2000.,0,200,-5.,5,Folder);

  HBook1F(Hist->fPt         ,"pt"       ,Form("%s: Track Pt"          ,Folder), 600, 75,95,Folder);
  HBook1F(Hist->fCosTh      ,"costh"    ,Form("%s: Track cos(theta)"  ,Folder), 100,-1,1,Folder);
  HBook1F(Hist->fChi2       ,"chi2"     ,Form("%s: Track chi2 total"  ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fNDof       ,"ndof"     ,Form("%s: Number of DOF"     ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fChi2Dof    ,"chi2d"    ,Form("%s: track chi2/N(dof)" ,Folder), 500, 0, 10,Folder);
  HBook1F(Hist->fChi2DofC   ,"chi2dc"   ,Form("%s: track chi2/N calc" ,Folder), 500, 0, 10,Folder);
  HBook1F(Hist->fNActive    ,"nactv"    ,Form("%s: N(active)"         ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fNWrong     ,"nwrong"   ,Form("%s: Nhits w wrong drft",Folder), 200, 0,200,Folder);
  HBook1F(Hist->fNDoublets  ,"nd"       ,Form("%s: N(Doublets)"       ,Folder),  50, 0, 50,Folder);
  HBook1F(Hist->fNOSDoublets,"nosd"     ,Form("%s: N(OS Doublets)"    ,Folder),  50, 0, 50,Folder);
  HBook1F(Hist->fNSSDoublets,"nssd"     ,Form("%s: N(SS Doublets)"    ,Folder),  50, 0, 50,Folder);
  HBook1F(Hist->fNAmb0      ,"namb0"    ,Form("%s: N(Amb=0)"          ,Folder),  50, 0, 50,Folder);
  HBook1F(Hist->fT0         ,"t0"       ,Form("%s: track T0"          ,Folder), 200, 0,2000,Folder);
  HBook1F(Hist->fT0Err      ,"t0err"    ,Form("%s: track T0Err"       ,Folder), 100, 0,  10,Folder);
  HBook1F(Hist->fQ          ,"q"        ,Form("%s: track Q"           ,Folder),   4,-2,   2,Folder);
  HBook1F(Hist->fFitCons[0] ,"fcon"     ,Form("%s: track fit cons [0]",Folder), 200, 0,   1,Folder);
  HBook1F(Hist->fFitCons[1] ,"fcon1"    ,Form("%s: track fit cons [1]",Folder), 1000, 0,   0.1,Folder);
  HBook1F(Hist->fD0         ,"d0"       ,Form("%s: track D0      "    ,Folder), 200,-200, 200,Folder);
  HBook1F(Hist->fZ0         ,"z0"       ,Form("%s: track Z0      "    ,Folder), 200,-2000,2000,Folder);
  HBook1F(Hist->fTanDip     ,"tdip"     ,Form("%s: track tan(dip)"    ,Folder), 200, 0.0 ,2.0,Folder);
  HBook1F(Hist->fDtZ0       ,"dtz0"     ,Form("%s: DT(Z0), MC"        ,Folder), 200, -10.0 ,10.0,Folder);

  HBook1F(Hist->fResid      ,"resid"    ,Form("%s: hit residuals"     ,Folder), 500,-0.5 ,0.5,Folder);
  HBook1F(Hist->fAlgMask    ,"alg"      ,Form("%s: algorithm mask"    ,Folder),  10,  0, 10,Folder);

  HBook1F(Hist->fChi2Match  ,"chi2tcm"  ,Form("%s: chi2(t-c match)"   ,Folder), 250,  0  ,250 ,Folder);
  HBook1F(Hist->fChi2XY     ,"chi2xy"   ,Form("%s: chi2(t-c match) XY",Folder), 300,-50  ,250 ,Folder);
  HBook1F(Hist->fChi2T      ,"chi2t"    ,Form("%s: chi2(t-c match) T" ,Folder), 250,  0  ,250 ,Folder);

  HBook1F(Hist->fDt         ,"dt"       ,Form("%s: track delta(T)"    ,Folder), 400,-20  ,20 ,Folder);
  HBook1F(Hist->fDx         ,"dx"       ,Form("%s: track delta(X)"    ,Folder), 200,-500 ,500,Folder);
  HBook1F(Hist->fDy         ,"dy"       ,Form("%s: track delta(Y)"    ,Folder), 200,-500 ,500,Folder);
  HBook1F(Hist->fDz         ,"dz"       ,Form("%s: track delta(Z)"    ,Folder), 200,-250 ,250,Folder);
  HBook1F(Hist->fDu         ,"du"       ,Form("%s: track-cluster DU)" ,Folder), 250,-250 ,250,Folder);
  HBook1F(Hist->fDv         ,"dv"       ,Form("%s: track-cluster DV)" ,Folder), 200,-100 ,100,Folder);
  HBook1F(Hist->fPath       ,"path"     ,Form("%s: track sdisk"       ,Folder),  50,   0 ,500,Folder);


  HBook2F(Hist->fFConsVsNActive,"fc_vs_na" ,Form("%s: FitCons vs NActive",Folder),  150, 0, 150, 200,0,1,Folder);
  HBook1F(Hist->fDaveTrkQual,"dtqual"   ,Form("%s:DaveTrkQual"        ,Folder), 200, -0.5, 1.5,Folder);

}

//-----------------------------------------------------------------------------
void TTrackCompModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1D(Hist->fLumWt     ,"lumwt"    ,Form("%s: Luminosity Weight"               ,Folder),200, 0,10,Folder);
  Hist->fLumWt->Sumw2(kTRUE);

  HBook1F(Hist->fRv        ,"rv"       ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 1000,Folder);
  HBook1F(Hist->fZv        ,"zv"       ,Form("%s: Z(Vertex)"                       ,Folder), 300, 0,15000,Folder);

  HBook1F(Hist->fPdgCode         ,"pdg"   ,Form("%s: PDG code"                     ,Folder),200,-100,100,Folder);
  HBook1F(Hist->fMomTargetEnd    ,"ptarg" ,Form("%s: CE mom after Stopping Target" ,Folder),400,  90,110,Folder);
  HBook1F(Hist->fMomTrackerFront ,"pfront",Form("%s: CE mom at the Tracker Front"  ,Folder),400,  90,110,Folder);
  HBook1F(Hist->fNshCE           ,"nsh_ce",Form("%s: CE Number of Straw Hits"      ,Folder),150,0,150,Folder);

  HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
  HBook1F(Hist->fEleMom    ,"ce_mom"   ,Form("%s: Conversion Electron Momentum"    ,Folder),1000,  0,200,Folder);
  HBook1F(Hist->fNTracks[0],"ntrk_0"   ,Form("%s: N(Reconstructed Tracks)[0]"      ,Folder),100,0,100,Folder);
  HBook1F(Hist->fNTracks[1],"ntrk_1"   ,Form("%s: N(Reconstructed Tracks)[1]"      ,Folder),100,0,100,Folder);
  HBook1F(Hist->fNshTot [0],"nshtot_0" ,Form("%s: Total Number of Straw Hits [0]"  ,Folder),250,0,250,Folder);
  HBook1F(Hist->fNshTot [1],"nshtot_1" ,Form("%s: Total Number of Straw Hits [1]"  ,Folder),250,0,5000,Folder);
  HBook1F(Hist->fNGoodSH   ,"nsh50"    ,Form("%s: N(SH) +/-50"                     ,Folder),300,0,1500,Folder);
  HBook1F(Hist->fDtClT     ,"dt_clt"   ,Form("%s: DT(cluster-track)"               ,Folder),100,-100,100,Folder);
  HBook1F(Hist->fDtClS     ,"dt_cls"   ,Form("%s: DT(cluster-straw hit)"           ,Folder),200,-200,200,Folder);
  HBook1F(Hist->fSHTime    ,"shtime"   ,Form("%s: Straw Hit Time"                  ,Folder),400,0,2000,Folder);
  HBook1F(Hist->fNHyp      ,"nhyp"     ,Form("%s: N(fit hypotheses)"               ,Folder),5,0,5,Folder);
  HBook1F(Hist->fBestHyp[0],"bfh0"     ,Form("%s: Best Fit Hyp[0](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(Hist->fBestHyp[1],"bfh1"     ,Form("%s: Best Fit Hyp[1](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(Hist->fNGenp     ,"ngenp"    ,Form("%s: N(Gen Particles)"                ,Folder),500,0,500,Folder);
  HBook1F(Hist->fNClusters ,"ncl"      ,Form("%s: N(Clusters)"                     ,Folder),100,0,100,Folder);
  HBook1F(Hist->fEClMax    ,"eclmax"   ,Form("%s: Max cluster energy"              ,Folder),150,0,150,Folder);
  HBook1F(Hist->fTClMax    ,"tclmax"   ,Form("%s: highest cluster time"            ,Folder),200,0,2000,Folder);
  HBook1F(Hist->fDp        ,"dp"       ,Form("%s: P(TPR)-P(CPR)"                   ,Folder),500,-2.5,2.5,Folder);
  HBook1F(Hist->fInstLumi  ,"dp"       ,Form("%s: Inst Luminosity"                 ,Folder),500,-2.5e6,2.5e6,Folder);
}

//_____________________________________________________________________________
void TTrackCompModule::BookHistograms() {

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
  book_event_histset[ 1] = 1;		// events with EclMax > 60 and TClMax > 550

					// TrkPatRec eff: histsets 10:19, CalpatRec efficiency:20-29

  for (int i=10; i<20; i++) book_event_histset[i] = 1;
  for (int i=20; i<30; i++) book_event_histset[i] = 1;

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

  book_track_histset[  0] = 1;          // good TrkPatRec tracks 
  book_track_histset[  1] = 1;          // good CalPatRec tracks in events where there is no good TrkPatRec tracks TrkQual>0.3
  book_track_histset[  2] = 1;          // good CalPatRec tracks in events where there is no good TrkPatRec tracks TrkQual>0.2
  book_track_histset[  3] = 1;          // best track


  book_track_histset[100] = 1;		// TrkPatRec all  tracks 
  book_track_histset[101] = 1;		// TrkPatRec BestTrackID
  book_track_histset[102] = 1;		// TrkPatRec BestTrackID no fitCons&momErr&t0Err tracks 
  book_track_histset[103] = 1;		// TrkPatRec BestTrackID, dpf>1 tracks 
  book_track_histset[104] = 1;          // TrkPatRec all  tracks events Ecl > 60
  book_track_histset[105] = 1;          // TrkPatRec BestTrackID tracks events Ecl > 60

  book_track_histset[110] = 1;          // TrkPatRec TrackID[0] - SetC
  book_track_histset[111] = 1;          // TrkPatRec TrackID[1]
  book_track_histset[112] = 1;          // TrkPatRec TrackID[2]
  book_track_histset[113] = 1;          // TrkPatRec TrackID[3]
  book_track_histset[114] = 1;          // TrkPatRec TrackID[1]
  book_track_histset[115] = 1;          // TrkPatRec TrackID[2]
  book_track_histset[116] = 1;          // TrkPatRec TrackID[3]

  book_track_histset[200] = 1;		// CalPatRec all  tracks 
  book_track_histset[201] = 1;		// CalPatRec BestTrackID tracks 
  book_track_histset[202] = 1;		// CalPatRec BestTrackID no fitCons&momErr&t0Err tracks 
  book_track_histset[203] = 1;		// CalPatRec BestTrackID, dpf>1 tracks 
  book_track_histset[204] = 1;          // CalPatRec all  tracks events Ecl > 60
  book_track_histset[205] = 1;          // CalPatRec BestTrackID tracks events Ecl > 60

  book_track_histset[210] = 1;          // CalPatRec TrackID[0] - SetC
  book_track_histset[211] = 1;          // CalPatRec TrackID[1]
  book_track_histset[212] = 1;          // CalPatRec TrackID[2]
  book_track_histset[213] = 1;          // CalPatRec TrackID[3]
  book_track_histset[214] = 1;          // CalPatRec TrackID[1]
  book_track_histset[215] = 1;          // CalPatRec TrackID[2]
  book_track_histset[216] = 1;          // CalPatRec TrackID[3]

  book_track_histset[300] = 1;		// TrkPatRec not CalPatRec all  tracks 
  book_track_histset[301] = 1;		// TrkPatRec not CalPatRec BestTrackID tracks 
  book_track_histset[302] = 1;		// TrkPatRec not CalPatRec BestTrackID no fitCons&momErr&t0Err tracks 
  book_track_histset[303] = 1;		// TrkPatRec not CalPatRec BestTrackID, dpf>1 tracks 
  book_track_histset[304] = 1;          // TrkPatRec not CalPatRec all  tracks events Ecl > 60
  book_track_histset[305] = 1;          // TrkPatRec not CalPatRec BestTrackID tracks events Ecl > 60

  book_track_histset[310] = 1;          // TrkPatRec not CalPatRec TrackID[0] - SetC
  book_track_histset[311] = 1;          // TrkPatRec not CalPatRec TrackID[1]
  book_track_histset[312] = 1;          // TrkPatRec not CalPatRec TrackID[2]
  book_track_histset[313] = 1;          // TrkPatRec not CalPatRec TrackID[3]
  book_track_histset[314] = 1;          // TrkPatRec not CalPatRec TrackID[1]
  book_track_histset[315] = 1;          // TrkPatRec not CalPatRec TrackID[2]
  book_track_histset[316] = 1;          // TrkPatRec not CalPatRec TrackID[3]

  book_track_histset[400] = 1;		// CalPatRec not TrkPatRec all  tracks 
  book_track_histset[401] = 1;		// CalPatRec not TrkPatRec BestTrackID tracks 
  book_track_histset[402] = 1;		// CalPatRec not TrkPatRec BestTrackID no fitCons&momErr&t0Err tracks 
  book_track_histset[403] = 1;		// CalPatRec not TrkPatRec BestTrackID, dpf>1 tracks 
  book_track_histset[404] = 1;          // CalPatRec not TrkPatRec all  tracks events Ecl > 60
  book_track_histset[405] = 1;          // CalPatRec not TrkPatRec BestTrackID tracks events Ecl > 60

  book_track_histset[410] = 1;          // CalPatRec not TrkPatRec TrackID[0] - SetC
  book_track_histset[411] = 1;          // CalPatRec not TrkPatRec TrackID[1]
  book_track_histset[412] = 1;          // CalPatRec not TrkPatRec TrackID[2]
  book_track_histset[413] = 1;          // CalPatRec not TrkPatRec TrackID[3]
  book_track_histset[414] = 1;          // CalPatRec not TrkPatRec TrackID[1]
  book_track_histset[415] = 1;          // CalPatRec not TrkPatRec TrackID[2]
  book_track_histset[416] = 1;          // CalPatRec not TrkPatRec TrackID[3]

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
// need MC truth branch
//-----------------------------------------------------------------------------
void TTrackCompModule::FillEventHistograms(EventHist_t* Hist) {
  double            cos_th, xv(-1.e6), yv(-1.e6), rv(-1.e6), zv(-1.e6), p;
  TLorentzVector    mom (1.,0.,0.,0);

  if (fParticle) { 
    fParticle->Momentum(mom);
    xv = fParticle->Vx()+3904.;
    yv = fParticle->Vy();
    rv = sqrt(xv*xv+yv*yv);
    zv = fParticle->Vz();
  }

  p      = mom.P();
  cos_th = mom.Pz()/p;

  xv = fParticle->Vx()+3904.;
  yv = fParticle->Vy();
  rv = sqrt(xv*xv+yv*yv);
  zv = fParticle->Vz();

  Hist->fLumWt->Fill(fLumWt);
  Hist->fRv->Fill(rv);
  Hist->fZv->Fill(zv);

  TSimParticle* simp = fSimPar.fParticle;

  Hist->fPdgCode->Fill(simp->fPdgCode);
  Hist->fMomTargetEnd->Fill(simp->fMomTargetEnd);
  Hist->fMomTrackerFront->Fill(simp->fMomTrackerFront);	// 
  Hist->fNshCE->Fill(simp->fNStrawHits);

  Hist->fEleMom->Fill(p);
  Hist->fEleCosTh->Fill(cos_th);

  Hist->fNTracks[0]->Fill(fNTracks[0]);
  Hist->fNTracks[1]->Fill(fNTracks[1]);

  int nsh_tot = GetHeaderBlock()->fNStrawHits;

  Hist->fNshTot[0]->Fill(nsh_tot);
  Hist->fNshTot[1]->Fill(nsh_tot);

  Hist->fNClusters->Fill(fNClusters);
  Hist->fEClMax->Fill(fEClMax);
  Hist->fTClMax->Fill(fTClMax);

  double dp(1.e6);

  if ((fNTracks[0] == 1) && (fNTracks[1] == 1)) {
    TrackPar_t* tp = &fTrackPar[0][0];
    TrackPar_t* cp = &fTrackPar[1][0];

    dp = tp->fP-cp->fP;
  }
					// momentum difference
  Hist->fDp->Fill(dp);
}

//-----------------------------------------------------------------------------
// fill efficiency histograms : need 10 histogram sets
// pitch = 1./tan(dip)
//-----------------------------------------------------------------------------
void TTrackCompModule::FillEfficiencyHistograms(TStnTrackBlock*  TrackBlock, 
						TStnTrackID*     TrackID   , 
						TrackPar_t*      TPar      ,
						int              HistSet   ) {
  if (fSimPar.fParticle->NStrawHits() >= 20) {
    FillEventHistograms(fHist.fEvent[HistSet]);

    if (fSimp->fMomTrackerFront > 100.) {
      FillEventHistograms(fHist.fEvent[HistSet+1]);

      TLorentzVector vdmom;
      vdmom.SetXYZM(fSimPar.fTFront->McMomentumX(),
		    fSimPar.fTFront->McMomentumY(),		      
		    fSimPar.fTFront->McMomentumZ(),
		    fSimPar.fTFront->Mass());

      float ce_pitch  = vdmom.Pt()/vdmom.Pz();
      float min_pitch = 1./TrackID->MaxTanDip();
      float max_pitch = 1./TrackID->MinTanDip();

      if ((min_pitch < ce_pitch) && (ce_pitch < max_pitch)) {
	FillEventHistograms(fHist.fEvent[HistSet+2]);
	  
	if (TrackBlock->NTracks() > 0) {
	  TStnTrack* track = TrackBlock->Track(0);
	  int id_word      = TrackID->IDWord(track);

	  FillEventHistograms(fHist.fEvent[HistSet+3]);
	  
	  if ((id_word & TStnTrackID::kTrkQualBit) == 0) {
	    FillEventHistograms(fHist.fEvent[HistSet+4]);
	    
	    if ((id_word & TStnTrackID::kT0Bit) == 0) {
	      FillEventHistograms(fHist.fEvent[HistSet+5]);
	      
	      if ((id_word & TStnTrackID::kTanDipBit) == 0) {
		FillEventHistograms(fHist.fEvent[HistSet+6]);
		
		if (((id_word & TStnTrackID::kD1Bit) == 0) && 
		    ((id_word & TStnTrackID::kD1Bit) == 0)    ) {
		  
		  FillEventHistograms(fHist.fEvent[HistSet+7]);
		  
		  if ((id_word & TStnTrackID::kTanDipBit) == 0) {
		    FillEventHistograms(fHist.fEvent[HistSet+8]);

		    if ((103.5 < TPar->fP) && (TPar->fP < 105)) {
		      FillEventHistograms(fHist.fEvent[HistSet+9]);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

//-----------------------------------------------------------------------------
// for DIO : ultimately, one would need to renormalize the distribution
//-----------------------------------------------------------------------------
void TTrackCompModule::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track, TrackPar_t* Tp) {

  TLorentzVector  mom;
					// pointer to local track parameters
  //  int itrk = Track->Number();

  //  TrackPar_t* tp = fTrackPar+itrk;

					// fP - corrected momentum, fP0 and fP2 - not corrected
  Hist->fP[0]->Fill (Tp->fP);
  Hist->fP[1]->Fill (Tp->fP);
  Hist->fP[2]->Fill (Tp->fP);
					// track fP0 and track fP2 are supposed to be the same...
  Hist->fP0->  Fill (Track->fP0);
  Hist->fP2->  Fill (Track->fP2);

  Hist->fPDio->Fill(Tp->fP,Tp->fDioWt);

  Hist->fPlw->Fill   (Tp->fP, Tp->fLumWt);
  Hist->fPDiolw->Fill(Tp->fP, Tp->fTotWt);

  Hist->fFitMomErr->Fill(Track->fFitMomErr);

  Hist->fPt    ->Fill(Track->fPt    );
  Hist->fPFront->Fill(Track->fPFront);
  Hist->fPStOut->Fill(Track->fPStOut);
					// dp: Tracker-only resolution

  Hist->fDpFront ->Fill(Tp->fDpF);
  Hist->fXDpF    ->Fill(Tp->fDpF/Track->fFitMomErr);
  Hist->fDpFDio  ->Fill(Tp->fDpF,Tp->fDioWt);
  Hist->fDpFront0->Fill(Tp->fDp0);
  Hist->fDpFront2->Fill(Tp->fDp2);
  Hist->fDpFSt   ->Fill(Tp->fDpFSt);
  Hist->fDpFVsZ1 ->Fill(Track->fZ1,Tp->fDpF);

  Hist->fCosTh->Fill(Track->Momentum()->CosTheta());
  Hist->fChi2->Fill (Track->fChi2);
  Hist->fNDof->Fill(Track->NActive()-5.);
  Hist->fChi2Dof->Fill(Track->fChi2/(Track->NActive()-5.));
  Hist->fNActive->Fill(Track->NActive());
  Hist->fNWrong->Fill(Track->NWrong());
  Hist->fNDoublets->Fill(Track->NOSDoublets()+Track->NSSDoublets());
  Hist->fNOSDoublets->Fill(Track->NOSDoublets());
  Hist->fNSSDoublets->Fill(Track->NSSDoublets());
  Hist->fNAmb0->Fill(Track->fNDoublets >> 16);
  Hist->fT0->Fill(Track->fT0);
  Hist->fT0Err->Fill(Track->fT0Err);
  Hist->fQ->Fill(-1);
  Hist->fFitCons[0]->Fill(Track->fFitCons);
  Hist->fFitCons[1]->Fill(Track->fFitCons);

  Hist->fD0->Fill(Track->fD0);
  Hist->fZ0->Fill(Track->fZ0);
  Hist->fTanDip->Fill(Track->fTanDip);
  Hist->fDtZ0->Fill(Tp->fDtZ0);
  
  Hist->fAlgMask->Fill(Track->AlgMask());

  Hist->fChi2Match->Fill(Tp->fChi2Match);
  Hist->fChi2XY->Fill(Tp->fChi2XY);
  Hist->fChi2T->Fill (Tp->fChi2T);

  Hist->fDt->Fill(Tp->fDt);
  Hist->fDx->Fill(Tp->fDx);
  Hist->fDy->Fill(Tp->fDy);
  Hist->fDz->Fill(Tp->fDz);
  Hist->fDu->Fill(Tp->fDu);
  Hist->fDv->Fill(Tp->fDv);
  Hist->fPath->Fill(Tp->fPath);

  Hist->fFConsVsNActive->Fill(Track->NActive(),Track->fFitCons);
  Hist->fDaveTrkQual->Fill(Track->DaveTrkQual());

}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TTrackCompModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("TrkPatRec"     ,"TStnTrackBlock"     ,&fTrackBlock[0]);
  RegisterDataBlock("CalPatRec"     ,"TStnTrackBlock"     ,&fTrackBlock[1]);
  RegisterDataBlock("ClusterBlock"  ,"TStnClusterBlock"   ,&fClusterBlock );
  RegisterDataBlock("SimpBlock"     ,"TSimpBlock"         ,&fSimpBlock    );
  RegisterDataBlock("GenpBlock"     ,"TGenpBlock"         ,&fGenpBlock    );
  RegisterDataBlock("VdetBlock"     ,"TVdetDataBlock"     ,&fVdetBlock    );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// initialize likelihood histograms
//-----------------------------------------------------------------------------
  return 0;
}


//_____________________________________________________________________________
int TTrackCompModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}

//_____________________________________________________________________________

void TTrackCompModule::FillHistograms() {

  TStnTrack*   trk;
  TrackPar_t*  tp;
  int          ihist, n_setc_tracks[2];
//-----------------------------------------------------------------------------
// event histograms
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// EVT_1: events with 50 MeV+ cluster and T0 > 550
//-----------------------------------------------------------------------------
  if ((fEClMax > 50.) && (fTClMax > 550)) {
    FillEventHistograms(fHist.fEvent[1]);
  }
//-----------------------------------------------------------------------------
// EVT_2: events with 50 MeV+ cluster and both tracks T0 > 550
//-----------------------------------------------------------------------------
  if ((fEClMax > 50.) && (fTClMax > 550)) {
    FillEventHistograms(fHist.fEvent[1]);
  }

//-----------------------------------------------------------------------------
// what does CalPatRec add ?
// TRK_0 : TrkPatRec tracks, BEST_ID
//-----------------------------------------------------------------------------
  if ((fTrackBlock[0]->NTracks() > 0) && (fTrackPar[0][0].fIDWord[fBestID] == 0)) {
    TStnTrack* trk = fTrackBlock[0]->Track(0);
    TrackPar_t* tp = &fTrackPar[0][0];
    FillTrackHistograms(fHist.fTrack[0],trk,tp);
  }
  else {
    if (fTrackBlock[1]->NTracks() > 0) {
      if (fTrackPar[1][0].fIDWord[3] == 0) {
//-----------------------------------------------------------------------------
// TRK_1: have CalPatRec track TrkQual > 0.3 and no TrkPatRec TrkQual > 0.4 track
//-----------------------------------------------------------------------------
	TStnTrack* trk = fTrackBlock[1]->Track(0);
	TrackPar_t* tp = &fTrackPar[1][0];
	FillTrackHistograms(fHist.fTrack[1],trk,tp);
      }
      if (fTrackPar[1][0].fIDWord[2] == 0) {
//-----------------------------------------------------------------------------
// TRK_2: have CalPatRec track TrkQual > 0.2 and no TrkPatRec TrkQual > 0.4 track
//-----------------------------------------------------------------------------
	TStnTrack* trk = fTrackBlock[1]->Track(0);
	TrackPar_t* tp = &fTrackPar[1][0];
	FillTrackHistograms(fHist.fTrack[2],trk,tp);
      }
    }
  }
//-----------------------------------------------------------------------------
// TRK_3 : an attempt to define best track
//-----------------------------------------------------------------------------
    TStnTrack*  best_track(0);
    TrackPar_t* best_tp;

    if (fTrackBlock[1]->NTracks() > 0) {
      if (fTrackPar[1][0].fIDWord[fBestID] == 0) {
	best_track = fTrackBlock[1]->Track(0);
	best_tp    = &fTrackPar[1][0];
      }
      else if ((fTrackBlock[0]->NTracks() > 0) && (fTrackPar[0][0].fIDWord[fBestID] == 0)) {
	best_track = fTrackBlock[0]->Track(0);
	best_tp    = &fTrackPar[0][0];
      }
      else if (fTrackPar[1][0].fIDWord[2] == 0) {
	best_track = fTrackBlock[1]->Track(0);
	best_tp    = &fTrackPar[1][0];
      }
    }
    else {
      if ((fTrackBlock[0]->NTracks() > 0) && (fTrackPar[0][0].fIDWord[fBestID] == 0)) {
	best_track = fTrackBlock[0]->Track(0);
	best_tp    = &fTrackPar[0][0];
      }
    }

    if (best_track != 0) {
      FillTrackHistograms(fHist.fTrack[3],best_track,best_tp);
    }

//-----------------------------------------------------------------------------
// TrkPatRec and CalPatRec histograms, inclusive, ihist defines the offset
// i=0:TrkPatRec, i=1:CalPatRec
//-----------------------------------------------------------------------------
  for (int i=0; i<2; i++) {
    n_setc_tracks[i] = 0;
    ihist            = 100*(i+1);
    for (int itrk=0; itrk<fNTracks[i]; itrk++) {
      trk = fTrackBlock[i]->Track(itrk);
      tp  = fTrackPar[i]+itrk;
//-----------------------------------------------------------------------------
// set IHIST+0: all tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[ihist+0],trk,tp);
//-----------------------------------------------------------------------------
// TRK_101, TRK_201: BestID 
//-----------------------------------------------------------------------------
      if (tp->fIDWord[fBestID] == 0) {
	FillTrackHistograms(fHist.fTrack[ihist+1],trk,tp);
	n_setc_tracks[i] += 1;
      }
//-----------------------------------------------------------------------------
// TRK_102, TRK_202: (BestID - FitConsBit - T0ErrBit - MomErrBit) tracks 
//-----------------------------------------------------------------------------
      int mask = TStnTrackID::kFitConsBit || TStnTrackID::kT0ErrBit || TStnTrackID::kFitMomErrBit;
      if ((tp->fIDWord[fBestID] & ~mask) == 0) {
	FillTrackHistograms(fHist.fTrack[ihist+2],trk,tp);
      }
//-----------------------------------------------------------------------------
// IHIST+3: (SetC + (dpf > 1)tracks 
//-----------------------------------------------------------------------------
      if ((tp->fIDWord[fBestID] == 0) & (tp->fDpF > 1)) {
	FillTrackHistograms(fHist.fTrack[ihist+3],trk,tp);
      }
//-----------------------------------------------------------------------------
// IHIST+4: add   Ecl > 60 requirement
//-----------------------------------------------------------------------------
      if (fEClMax > 60.) {
	FillTrackHistograms(fHist.fTrack[ihist+4],trk,tp);

	if (tp->fIDWord[0] == 0) {
	  FillTrackHistograms(fHist.fTrack[ihist+5],trk,tp);
	}
      }

      for (int idd=0; idd<fNID; idd++) {
	if (tp->fIDWord[idd] == 0) FillTrackHistograms(fHist.fTrack[ihist+10+idd],trk,tp);
      }
    }
//-----------------------------------------------------------------------------
// either ID=300 (TrkPatRec not CalPatRec) or 400(CalPatRec not TrkPatRec)
//-----------------------------------------------------------------------------
    if ((fNTracks[i] > 0) && (fNTracks[1-i] == 0)) {

      ihist = 100*(i+1)+200;

      for (int itrk=0; itrk<fNTracks[i]; itrk++) {
	trk = fTrackBlock[i]->Track(itrk);
	tp  = fTrackPar[i]+itrk;
	//-----------------------------------------------------------------------------
	// set IHIST+0: all tracks
	//-----------------------------------------------------------------------------
	FillTrackHistograms(fHist.fTrack[ihist+0],trk,tp);
	//-----------------------------------------------------------------------------
	// IHIST+1: Set C selection
	//-----------------------------------------------------------------------------
	if (tp->fIDWord[fBestID] == 0) {
	  FillTrackHistograms(fHist.fTrack[ihist+1],trk,tp);
	  n_setc_tracks[i] += 1;
	}
	//-----------------------------------------------------------------------------
	// IHIST+2: (SetC - FitConsBit - T0ErrBit - MomErrBit) tracks 
	//-----------------------------------------------------------------------------
	int mask = TStnTrackID::kFitConsBit || TStnTrackID::kT0ErrBit || TStnTrackID::kFitMomErrBit;
	if ((tp->fIDWord[fBestID] & ~mask) == 0) {
	  FillTrackHistograms(fHist.fTrack[ihist+2],trk,tp);
	}
	//-----------------------------------------------------------------------------
	// IHIST+3: (SetC + (dpf > 1)tracks 
	//-----------------------------------------------------------------------------
	if ((tp->fIDWord[fBestID] == 0) & (tp->fDpF > 1)) {
	  FillTrackHistograms(fHist.fTrack[ihist+3],trk,tp);
	}
	//-----------------------------------------------------------------------------
	// IHIST+4: add   Ecl > 60 requirement
	//-----------------------------------------------------------------------------
	if (fEClMax > 60.) {
	  FillTrackHistograms(fHist.fTrack[ihist+4],trk,tp);

	  if (tp->fIDWord[0] == 0) {
	    FillTrackHistograms(fHist.fTrack[ihist+5],trk,tp);
	  }
	}

	for (int idd=0; idd<fNID; idd++) {
	  if (tp->fIDWord[idd] == 0) FillTrackHistograms(fHist.fTrack[ihist+10+idd],trk,tp);
	}
      }

    }
  }
//-----------------------------------------------------------------------------
// efficiency histograms, use fDaveTrkQual > 0.4 for the cuts
//-----------------------------------------------------------------------------
  FillEfficiencyHistograms(fTrackBlock[0],fTrackID[4],&fTrackPar[0][0],10);
  FillEfficiencyHistograms(fTrackBlock[1],fTrackID[4],&fTrackPar[1][0],20);
}



//-----------------------------------------------------------------------------
// assume less than 20 tracks 
//-----------------------------------------------------------------------------
int TTrackCompModule::InitTrackPar(TStnTrackBlock*   TrackBlock  , 
				   TStnClusterBlock* ClusterBlock, 
				   TrackPar_t*       TrackPar    ) {
  TrackPar_t*           tp;
  TStnTrack*            track;
  int                   id_word, icorr;
  double                xs;
  TEmuLogLH::PidData_t  dat;
//-----------------------------------------------------------------------------
// momentum corrections for TrkPatRec and CalPatRec
//-----------------------------------------------------------------------------
  const double kMomentumCorr[2] = { 0.049, 0.020 };

  const char* block_name = TrackBlock->GetNode()->GetName();

  if      (strcmp(block_name,"TrkPatRec" ) == 0) icorr = 0;
  else if (strcmp(block_name,"CalPatRec" ) == 0) icorr = 1;
  else if (strcmp(block_name,"TrackBlock") == 0) icorr = 2;
  else {
    icorr = -999;
    Error("TTrackCompModule::InitTrackPar","IN TROUBLE");
    return -1;
  }
//-----------------------------------------------------------------------------
// loop over tracks
//-----------------------------------------------------------------------------
  int ntrk = TrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    tp             = TrackPar+itrk;
    track          = TrackBlock->Track(itrk);

    for (int idd=0; idd<fNID; idd++) {
      tp->fIDWord[idd] = fTrackID[idd]->IDWord(track);
    }

    id_word = tp->fIDWord[fBestID];
    track->fIDWord = id_word;
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
//-----------------------------------------------------------------------------
// in this scheme correction is set right before the call
// in case of MergePatRec use BestAlg - 
// hopefully, TTrackComp will never use MergePatRec branch
//-----------------------------------------------------------------------------
    if (icorr == 2) icorr = track->BestAlg();

    tp->fP     = track->fP     +kMomentumCorr[icorr];		// correcting
    tp->fDpF   = tp->fP        -track->fPFront;
    tp->fDp0   = track->fP0    -track->fPFront;
    tp->fDp2   = track->fP2    -track->fPFront;
    tp->fDpFSt = track->fPFront-track->fPStOut;

    if (fFillDioHist == 0) tp->fDioWt = 1.;
    else                   tp->fDioWt = TStntuple::DioWeightAl(fEleE);

    tp->fLumWt   = GetHeaderBlock()->LumWeight();
    tp->fTotWt   = tp->fLumWt*tp->fDioWt;
    tp->fDioWtRC = tp->fDioWt;
    tp->fTotWtRC = tp->fLumWt*tp->fDioWtRC;

    tp->fDtZ0 = -1.e6;
    if (fSimPar.fTMid) tp->fDtZ0 = track->T0()-fSimPar.fTMid->Time();
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
    tp->fChi2XY    = -1.e6;
    tp->fChi2T     = -1.e6;
    tp->fPath      = -1.e6;
    tp->fSinTC     = -1.e6;
    tp->fDrTC      = -1.e6;
    tp->fSInt      = -1.e6;

    if (vr) {
      tp->fEcl = vr->fEnergy;
      tp->fEp  = tp->fEcl/tp->fP;

      tp->fDx  = vr->fDx;
      tp->fDy  = vr->fDy;
      tp->fDz  = vr->fDz;
//-----------------------------------------------------------------------------
// v4_2_4: correct by additional 0.22 ns - track propagation by 6 cm
//-----------------------------------------------------------------------------
      tp->fDt  = vr->fDt ; // v4_2_4: - 0.22; // - 1.;

      nx  = vr->fNxTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);
      ny  = vr->fNyTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);

      tp->fDu        = vr->fDx*nx+vr->fDy*ny;
      tp->fDv        = vr->fDx*ny-vr->fDy*nx;
      tp->fChi2Match = vr->fChi2Match;
					// from now on the matching chi2 has XY part only
      tp->fChi2XY    = vr->fChi2Match;
      tp->fChi2T     = vr->fChi2Time;
      tp->fPath      = vr->fPath;
//-----------------------------------------------------------------------------
// angle
//-----------------------------------------------------------------------------
      TStnCluster* cl = fClusterBlock->Cluster(vr->fClusterIndex);
      tp->fSinTC = nx*cl->fNy-ny*cl->fNx;
      tp->fDrTC  = vr->fDr;
      tp->fSInt  = vr->fSInt;
    }

    if ((tp->fEp > 0) && (track->fEp > 0) && (fabs(tp->fEp-track->fEp) > 1.e-6)) {
      GetHeaderBlock()->Print(Form(" TTrackAnaModule ERROR: tp->fEp = %10.5f  track->fEp = %10.5f\n ",tp->fEp,track->fEp));
    }
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

  return 0;
}

//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TTrackCompModule::Event(int ientry) {

  fTrackBlock[0]->GetEntry(ientry);
  fTrackBlock[1]->GetEntry(ientry);
  fClusterBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fVdetBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// luminosity weight
//-----------------------------------------------------------------------------
  fLumWt = GetHeaderBlock()->LumWeight();
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fNGenp     = fGenpBlock->NParticles();
  fNClusters = fClusterBlock->NClusters();

  fCluster = NULL;
  fEClMax  = -1;
  fTClMax  = -1;
  if (fNClusters > 0) {
    fCluster = fClusterBlock->Cluster(0);
    fEClMax  = fCluster->Energy();
    fTClMax  = fCluster->Time  ();
  }
//-----------------------------------------------------------------------------
// MC generator info
//-----------------------------------------------------------------------------
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

  if (fParticle) fEleE = fParticle->Energy();
  else           fEleE = -1.;
//-----------------------------------------------------------------------------
// may want to revisit the definition of fSimp in the future
//-----------------------------------------------------------------------------
  fSimp             = fSimpBlock->Particle(0);
  fSimPar.fParticle = fSimp;
//-----------------------------------------------------------------------------
// virtual detectors - for fSimp need parameters at the tracker front
//-----------------------------------------------------------------------------
  int nvdhits = fVdetBlock->NHits();
  for (int i=0; i<nvdhits; i++) {
    TVdetHitData* vdhit = fVdetBlock->Hit(i);
    if (vdhit->PdgCode() == fSimp->fPdgCode) {
      if ((vdhit->Index() == 13) || (vdhit->Index() == 14)) {
	fSimPar.fTFront = vdhit;
      }
      else if ((vdhit->Index() == 11) || (vdhit->Index() == 12)) {
	fSimPar.fTMid = vdhit;
      }
    }
  }
//-----------------------------------------------------------------------------
// initialize additional track parameters
//-----------------------------------------------------------------------------
  for (int i=0; i<2; i++) {
    fNTracks    [i] = fTrackBlock[i]->NTracks();
    fNGoodTracks[i] = 0;
    InitTrackPar(fTrackBlock[i],fClusterBlock,fTrackPar[i]);
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
// looking mostly at the CalPatRec tracks
//-----------------------------------------------------------------------------
// diagnostics is related to CalPatRec
//-----------------------------------------------------------------------------
void TTrackCompModule::Debug() {

  TStnTrack* trk;
  TrackPar_t* tp(NULL);
  char        text[500];
  int         calpatrec(1);
//-----------------------------------------------------------------------------
// bit 0: All Events
//-----------------------------------------------------------------------------
  if (GetDebugBit(0) == 1) {
    GetHeaderBlock()->Print(Form("TTrackCompModule :bit000:"));
  }

  TStnTrackBlock* cprb = fTrackBlock[calpatrec];

  int ntrk = cprb->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = cprb->Track(itrk);
    tp  = &fTrackPar[calpatrec][itrk];
//-----------------------------------------------------------------------------
// bit 3: Set C CALPATREC tracks with large DPF > 5 MeV
//-----------------------------------------------------------------------------
    if (GetDebugBit(3) == 1) {
      if ((trk->fIDWord == 0) && (tp->fDpF > 5.)) {
	if (tp->fDpF > 5) {
	  GetHeaderBlock()->Print(Form("TTrackCompModule bit003: tp->DpF = %10.3f trk->fP = %10.3f trk->fPFront = %10.3f",
				       tp->fDpF, trk->fP,trk->fPFront));
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// bit 4: events with TrkPatRec track, a 60 MeV+ cluster and no CalPatRec track
//-----------------------------------------------------------------------------
  if (GetDebugBit(4) == 1) {
    if ((fNTracks[0] > 0) && (fNClusters > 0) && (fNTracks[1] == 0)) {
      if (fEClMax > 60.) {
	sprintf(text,"TTrackCompModule bit004: N(TPR) = %i N(CPR) = %i E(cl) = %10.3f",
		fNTracks[0],fNTracks[1],fEClMax);
	GetHeaderBlock()->Print(text);
      }
    }
  }
//-----------------------------------------------------------------------------
// bit 5: Set C CALPATREC tracks with P > 106
//-----------------------------------------------------------------------------
  if ((GetDebugBit(5) == 1) && (ntrk > 0)) {
    trk = cprb->Track(0);

    if ((trk->fIDWord == 0) && (trk->fP > fDebugCut[5].fXMin) && (trk->fP < fDebugCut[5].fXMax)) {
      GetHeaderBlock()->Print(Form("TTrackCompModule bit005: tp->DpF = %10.3f trk->fP = %10.3f trk->fPFront = %10.3f",
				   tp->fDpF, trk->fP,trk->fPFront));
    }
  }
//-----------------------------------------------------------------------------
// bit 6: Set C CALPATREC tracks with 1.5 < tp->fDFp < 5
//-----------------------------------------------------------------------------
  if ((GetDebugBit(6) == 1) && (ntrk > 0)) {
    trk = cprb->Track(0);

    if ((trk->fIDWord == 0) && (tp->fDpF >= fDebugCut[6].fXMin) && (tp->fDpF < fDebugCut[6].fXMax)) {
      GetHeaderBlock()->Print(Form("TTrackCompModule bit006: tp->DpF = %10.3f trk->fP = %10.3f trk->fPFront = %10.3f",
				   tp->fDpF, trk->fP,trk->fPFront));
    }
  }
//-----------------------------------------------------------------------------
// bit 7: events with N(Set C CalPatRec tracks) > 0
//-----------------------------------------------------------------------------
  if ((GetDebugBit(7) == 1) && (fNGoodTracks[calpatrec] > 0)) {
    GetHeaderBlock()->Print(Form("TTrackCompModule :bit007:"));
  }
//-----------------------------------------------------------------------------
// bit 8: Set C CALPATREC tracks with P > 105
//-----------------------------------------------------------------------------
  if ((GetDebugBit(8) == 1) && (ntrk > 0)) {
    trk = cprb->Track(0);

    if ((trk->fIDWord == 0) && (trk->fP > 105.)) {
      GetHeaderBlock()->Print(Form("TTrackCompModule :bit008: tp->DpF = %10.3f trk->fP = %10.3f trk->fPFront = %10.3f",
				   tp->fDpF, trk->fP,trk->fPFront));
    }
  }


}

//_____________________________________________________________________________
int TTrackCompModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TTrackCompModule::Test001() {
}

