//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//  3  : events with CalPatRec tracks with DFF > 5
//  4  : events with TrkPatRec track and a cluster, but with no CalPatRec track
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
#include "ana/TTrackCompModule.hh"

ClassImp(TTrackCompModule)
//-----------------------------------------------------------------------------
TTrackCompModule::TTrackCompModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fTrackNumber.Set(100);

  fMinT0 = 0; // do not cut on time by default

  fTrackID = new TStnTrackID();
//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
//-----------------------------------------------------------------------------
  fPdgCode       = 11;
  fGeneratorCode = 2;			// conversionGun, 28:StoppedParticleReactionGun
}

//-----------------------------------------------------------------------------
TTrackCompModule::~TTrackCompModule() {
}

//-----------------------------------------------------------------------------
void TTrackCompModule::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP[0]       ,"p"        ,Form("%s: Track P(Z1)"       ,Folder), 400,  80  ,120. ,Folder);
  HBook1F(Hist->fP[1]       ,"p_1"      ,Form("%s: Track P(total)[1]" ,Folder), 100, 100  ,105. ,Folder);
  HBook1F(Hist->fP[2]       ,"p_2"      ,Form("%s: Track P(total)[1]" ,Folder),2000,   0  ,200. ,Folder);
  HBook1F(Hist->fP0         ,"p0"       ,Form("%s: Track P(Z0)"       ,Folder),1000,   0  ,200. ,Folder);
  HBook1F(Hist->fP2         ,"p2"       ,Form("%s: Track P(z=-1540)"  ,Folder),1000,   0  ,200. ,Folder);

  //  HBook1D(Hist->fPDio       ,"pdio"     ,Form("%s: Track P(DIO WT)"   ,Folder), 400,  90  ,110. ,Folder);
  //  Hist->fPDio->Sumw2(kTRUE);

  HBook1F(Hist->fFitMomErr  ,"momerr"   ,Form("%s: Track FitMomError" ,Folder), 200,   0  ,  1. ,Folder);
  HBook1F(Hist->fPFront     ,"pf"       ,Form("%s: Track P(front)   " ,Folder), 400,  90  ,110. ,Folder);
  HBook1F(Hist->fDpFront    ,"dpf"      ,Form("%s: Track P-P(front) " ,Folder), 200,  -5. ,  5. ,Folder);
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
  HBook1F(Hist->fT0         ,"t0"       ,Form("%s: track T0"          ,Folder), 200, 0,2000,Folder);
  HBook1F(Hist->fT0Err      ,"t0err"    ,Form("%s: track T0Err"       ,Folder), 100, 0,  10,Folder);
  HBook1F(Hist->fQ          ,"q"        ,Form("%s: track Q"           ,Folder),   4,-2,   2,Folder);
  HBook1F(Hist->fFitCons[0] ,"fcon"     ,Form("%s: track fit cons [0]",Folder), 200, 0,   1,Folder);
  HBook1F(Hist->fFitCons[1] ,"fcon1"    ,Form("%s: track fit cons [1]",Folder), 1000, 0,   0.1,Folder);
  HBook1F(Hist->fD0         ,"d0"       ,Form("%s: track D0      "    ,Folder), 200,-200, 200,Folder);
  HBook1F(Hist->fZ0         ,"z0"       ,Form("%s: track Z0      "    ,Folder), 200,-2000,2000,Folder);
  HBook1F(Hist->fTanDip     ,"tdip"     ,Form("%s: track tan(dip)"    ,Folder), 200, 0.0 ,2.0,Folder);
  HBook1F(Hist->fResid      ,"resid"    ,Form("%s: hit residuals"     ,Folder), 500,-0.5 ,0.5,Folder);
  HBook1F(Hist->fAlgMask    ,"alg"      ,Form("%s: algorithm mask"    ,Folder),  10,  0, 10,Folder);

  //  HBook1F(Hist->fZ1         ,"z1"       ,Form("%s: track Z1      "    ,Folder), 200,-2000,2000,Folder);
  //  HBook2F(Hist->fNHVsStation,"nh_vs_st" ,Form("%s: N(hits) Vs Station",Folder),  50, 0,50,10,-0.5,9.5,Folder);
  //  HBook2F(Hist->fNHVsNSt    ,"nh_vs_nst",Form("%s: N(hits) Vs NSt"    ,Folder),  10,-0.5,9.5,40,-0.5,39.5,Folder);

  //  HBook1F(Hist->fPdgCode    ,"pdg"      ,Form("%s: track PDG code"    ,Folder), 100,-50,50,Folder);
}

//-----------------------------------------------------------------------------
void TTrackCompModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
  HBook1F(Hist->fEleMom    ,"ce_mom"   ,Form("%s: Conversion Electron Momentum"    ,Folder),1000,  0,200,Folder);
  HBook1F(Hist->fRv        ,"rv"       ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 1000,Folder);
  HBook1F(Hist->fZv        ,"zv"       ,Form("%s: Z(Vertex)"                       ,Folder), 300, 0,15000,Folder);
  HBook1F(Hist->fNTracks[0],"ntrk_0"   ,Form("%s: N(Reconstructed Tracks)[0]"      ,Folder),100,0,100,Folder);
  HBook1F(Hist->fNTracks[1],"ntrk_1"   ,Form("%s: N(Reconstructed Tracks)[1]"      ,Folder),100,0,100,Folder);
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
}

//-----------------------------------------------------------------------------
void TTrackCompModule::BookSimpHistograms(SimpHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fPdgCode         ,"pdg"   ,Form("%s: PDG code"                     ,Folder),200,-100,100,Folder);
  HBook1F(Hist->fMomTargetEnd    ,"ptarg" ,Form("%s: CE mom after Stopping Target" ,Folder),400,  90,110,Folder);
  HBook1F(Hist->fMomTrackerFront ,"pfront",Form("%s: CE mom at the Tracker Front"  ,Folder),400,  90,110,Folder);
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
// book crystal histograms
//-----------------------------------------------------------------------------
//   HBook1F(fHist.fCrystalR[0],"rc_0"     ,Form("disk [0] crystal radius"),100,0,1000,"Hist");
//   HBook1F(fHist.fCrystalR[1],"rc_1"     ,Form("disk [1] crystal radius"),100,0,1000,"Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events

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

  book_track_histset[  0] = 1;		// TrkPatRec all tracks 
  book_track_histset[  1] = 1;		// TrkPatRec all tracks 

  book_track_histset[100] = 1;		// CalPatRec all tracks 
  book_track_histset[101] = 1;		// CalPatRec Set "C" tracks

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
  double            cos_th, xv, yv, rv, zv, p;
  TLorentzVector    mom;

  fParticle->Momentum(mom);

  p      = mom.P();

  cos_th = mom.Pz()/p;

  xv = fParticle->Vx()+3904.;
  yv = fParticle->Vy();
  rv = sqrt(xv*xv+yv*yv);
  zv = fParticle->Vz();

  Hist->fEleMom->Fill(p);
  Hist->fEleCosTh->Fill(cos_th);
  Hist->fRv->Fill(rv);
  Hist->fZv->Fill(zv);

  Hist->fNTracks[0]->Fill(fNTracks[0]);
  Hist->fNTracks[1]->Fill(fNTracks[1]);
}

//-----------------------------------------------------------------------------
void TTrackCompModule::FillSimpHistograms(SimpHist_t* Hist, TSimParticle* Simp) {

  Hist->fPdgCode->Fill(Simp->fPdgCode);
  Hist->fMomTargetEnd->Fill(Simp->fMomTargetEnd);
  Hist->fMomTrackerFront->Fill(Simp->fMomTrackerFront);
  Hist->fNStrawHits->Fill(Simp->fNStrawHits);
}

//-----------------------------------------------------------------------------
// for DIO : ultimately, one would need to renormalize the distribution
//-----------------------------------------------------------------------------
void TTrackCompModule::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track, TrackPar_t* Tp) {

  TLorentzVector  mom;
					// pointer to local track parameters
  //  itrk = Track->Number();

  Hist->fP[0]->Fill (Track->fP);
  Hist->fP[1]->Fill (Track->fP);
  Hist->fP[2]->Fill (Track->fP);
  Hist->fP0->  Fill (Track->fP0);
  Hist->fP2->  Fill (Track->fP2);

  Hist->fFitMomErr->Fill(Track->fFitMomErr);

  Hist->fPt    ->Fill(Track->fPt    );
  Hist->fPFront->Fill(Track->fPFront);
  Hist->fPStOut->Fill(Track->fPStOut);
					// dp: Tracker-only resolution

  Hist->fDpFront ->Fill(Tp->fDpF);
  Hist->fDpFront0->Fill(Tp->fDp0);
  Hist->fDpFront2->Fill(Tp->fDp2);
  Hist->fDpFSt   ->Fill(Tp->fDpFSt);
  Hist->fDpFVsZ1 ->Fill(Track->fZ1,Tp->fDpF);

  Hist->fCosTh->Fill(Track->Momentum()->CosTheta());
  Hist->fChi2->Fill (Track->fChi2);
  Hist->fNDof->Fill(Track->fNActive-5.);
  Hist->fChi2Dof->Fill(Track->fChi2/(Track->fNActive-5.));
  Hist->fNActive->Fill(Track->fNActive);
  Hist->fT0->Fill(Track->fT0);
  Hist->fT0Err->Fill(Track->fT0Err);
  Hist->fQ->Fill(-1);
  Hist->fFitCons[0]->Fill(Track->fFitCons);
  Hist->fFitCons[1]->Fill(Track->fFitCons);

  Hist->fD0->Fill(Track->fD0);
  Hist->fZ0->Fill(Track->fZ0);
  Hist->fTanDip->Fill(Track->fTanDip);
  Hist->fAlgMask->Fill(Track->AlgMask());

//   chi2c = Track->fChi2C/(Track->fNActive-5.);
//   Hist->fChi2DofC->Fill(chi2c);
//-----------------------------------------------------------------------------
// there is an inconsistency in the SIMP block filling - in Mu2e offline 
// the particle momentumis is kept in MeV/c, while the PDG mass  -in GeV/c^2..
// thus the energy is screwed up... kludge around
// assign muon mass
//-----------------------------------------------------------------------------
//   double ekin(-1.);
//   if (fSimp) {
//     double p, m;
//     //    p    = fSimp->fStartMom.P();
//     p = Track->fP;
//     m    = 105.658; // in MeV
//     ekin = sqrt(p*p+m*m)-m;
//   }
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
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// initialize likelihood histograms
//-----------------------------------------------------------------------------
  fTrackID->SetMinT0(fMinT0);
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
  int          ihist;
//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// Simp histograms
//-----------------------------------------------------------------------------
//   if (fSimp) {
//     FillSimpHistograms(fHist.fSimp[0],fSimp);
//   }
//-----------------------------------------------------------------------------
// TrkPatRec histograms, ihist defines the offset
//-----------------------------------------------------------------------------
  for (int i=0; i<2; i++) {
    ihist = 100*i;
    for (int itrk=0; itrk<fNTracks[i]; itrk++) {
      trk = fTrackBlock[i]->Track(itrk);
      tp  = fTrackPar[i]+itrk;
//-----------------------------------------------------------------------------
// set IHIST+0: all tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[ihist+0],trk,tp);

      if (trk->fIDWord == 0) {
//-----------------------------------------------------------------------------
// IHIST+1: Set C selection
//-----------------------------------------------------------------------------
	FillTrackHistograms(fHist.fTrack[ihist+1],trk,tp);
      }
    }
  }
}



//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TTrackCompModule::Event(int ientry) {

  double                p;
  TStnTrack*            track;
  int                   id_word; //, alg_mask;
  TLorentzVector        mom;

  fTrackBlock[0]->GetEntry(ientry);
  fTrackBlock[1]->GetEntry(ientry);
  fClusterBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fNGenp     = fGenpBlock->NParticles();
  fNClusters = fClusterBlock->NClusters();

  fCluster = NULL;
  fEclMax  = -1;
  if (fNClusters > 0) {
    fCluster = fClusterBlock->Cluster(0);
    fEclMax  = fCluster->Energy();
  }

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

  fParticle->Momentum(mom);
					// this is a kludge, to be removed at the next 
					// ntupling 
  //  fEleE     = fParticle->Energy();
  p         = mom.P();
  fEleE     = sqrt(p*p+0.511*0.511);
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  TrackPar_t*   tp;
  int           ntrk;

  for (int i=0; i<2; i++) {
    fNTracks    [i] = fTrackBlock[i]->NTracks();
    fNGoodTracks[i] = 0;

    ntrk = fNTracks[i];

    for (int itrk=0; itrk<ntrk; itrk++) {
//-----------------------------------------------------------------------------
// assume less 20 tracks
//-----------------------------------------------------------------------------
      tp             = fTrackPar[i]+itrk;
      track          = fTrackBlock[i]->Track(itrk);

      id_word        = fTrackID->IDWord(track);
      track->fIDWord = id_word;
      if (id_word == 0) {
	fNGoodTracks[i] += 1;
      }

      //      alg_mask = track->AlgMask();
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

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
// looking mostly at the CalPatRec tracks
//-----------------------------------------------------------------------------
void TTrackCompModule::Debug() {

  TStnTrack* trk;
  TrackPar_t* tp;
  char        text[500];

  int calpatrec(1);

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
// bit 4: events with TrkPatRec track and no CalPatRec track
//-----------------------------------------------------------------------------
  if (GetDebugBit(0) == 1) {
    sprintf(text,"TTrackCompModule bit004: N(TPR) = %i N(CPR) = %i E(cl) = %10.3f",
	    fNTracks[0],fNTracks[1],fEclMax);
    GetHeaderBlock()->Print(text);
  }
//-----------------------------------------------------------------------------
// bit 4: events with TrkPatRec track and no CalPatRec track
//-----------------------------------------------------------------------------
  if (GetDebugBit(4) == 1) {
    if ((fNTracks[0] > 0) && (fNClusters > 0) && (fNTracks[1] == 0)) {
      if (fEclMax > 60.) {
	sprintf(text,"TTrackCompModule bit004: N(TPR) = %i N(CPR) = %i E(cl) = %10.3f",
		fNTracks[0],fNTracks[1],fEclMax);
	GetHeaderBlock()->Print(text);
	
      }
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

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}

