///////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
// 0  : all events
// 1  : passed events
// 2  : rejected events
// 
// 3  : events with set tracks passing the quality cuts and N(tracks) > 1
// 4  : unused
// 5  : unused
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
// #include "CalorimeterGeom/inc/HexMap.hh"

#include "ana/TTrackRecoEffAnaModule.hh"

ClassImp(TTrackRecoEffAnaModule)
//-----------------------------------------------------------------------------
TTrackRecoEffAnaModule::TTrackRecoEffAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fPdgCode       = 11;
  fGeneratorCode = 28;
  fTrkPatRecOnly = 0;

  fMinNMCHits    = 20;
  fMinMCMomentum = 100.;
  fMinMCPitch    = 1.;
  fMaxMCPitch    = sqrt(3.);
}

//-----------------------------------------------------------------------------
TTrackRecoEffAnaModule::~TTrackRecoEffAnaModule() {
}


//-----------------------------------------------------------------------------
void TTrackRecoEffAnaModule::BookStrawHitHistograms(StrawHitHist_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  //-----------------------------------------------------------------------------
  //  
  //-----------------------------------------------------------------------------

  HBook1F(Hist->fPdgCode,"pdg_code",Form("%s: PDG code",Folder),  200,-1000, 1000,Folder);
  HBook1F(Hist->fGenCode,"gen_code",Form("%s: generator code",Folder),  100, -10, 90,Folder);
  HBook1F(Hist->fEnergy,"energy",Form("%s: Hit Energy" ,Folder),200, 0, 0.1,Folder);
  HBook1F(Hist->fTime  ,"time"  ,Form("%s: Hit Time  " ,Folder),200, 0,2000,Folder);
  HBook1F(Hist->fDt    ,"dt"    ,Form("%s: Hit Dt"     ,Folder),1000, -10,10,Folder);
  HBook1F(Hist->fMcMomentum,"mom"    ,Form("%s: MC Momentum",Folder),200, 0,200,Folder);
}


//-----------------------------------------------------------------------------
void TTrackRecoEffAnaModule::BookGenpHistograms(GenpHist_t* Hist, const char* Folder) {
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
void TTrackRecoEffAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fNStrawHits[0] ,"nsh_0"   ,Form("%s: N Straw Hits[0]" ,Folder),  200, 0,  200,Folder);
  HBook1F(Hist->fNStrawHits[1] ,"nsh_1"   ,Form("%s: N Straw Hits[1]" ,Folder), 1000, 0,10000,Folder);
  HBook1F(Hist->fNHitsSignal   ,"nhs"     ,Form("%s: N Hits   CE"     ,Folder),  200, 0,  200,Folder);
  HBook1F(Hist->fMomTF         ,"mom_tf"  ,Form("%s: Momentum CE F"   ,Folder), 1100, 0,  110,Folder);
  
  HBook1F(Hist->fNTracks       ,"ntrk"    ,Form("%s: NTracks"         ,Folder),   20, 0,   20,Folder);
  HBook1F(Hist->fFitCons       ,"fcons"   ,Form("%s: Fit Consistency" ,Folder), 1000, 0,    1,Folder);
  HBook1F(Hist->fT0            ,"t0"      ,Form("%s: T0"              ,Folder),  200, 0, 2000,Folder);
  HBook1F(Hist->fT0Err         ,"t0err"   ,Form("%s: T0 Err"          ,Folder),  200, 0,    2,Folder);
  HBook1F(Hist->fNActive       ,"nactive" ,Form("%s: N Active"        ,Folder),  100, 0,  100,Folder);
  HBook1F(Hist->fTanDipMC      ,"tdip_mc" ,Form("%s: Tan(DIP) MC"     ,Folder),  200, 0,    2,Folder);
  HBook1F(Hist->fTanDip        ,"tdip"    ,Form("%s: Tan(DIP)"        ,Folder),  200, 0,    2,Folder);
  HBook1F(Hist->fD0            ,"d0"      ,Form("%s: D0"              ,Folder),  200,-200,200,Folder);
  HBook1F(Hist->fRMax          ,"rmax"    ,Form("%s: RMax"            ,Folder),  200, 0, 2000,Folder);
  HBook1F(Hist->fP             ,"p"       ,Form("%s: Track Momentum"  ,Folder), 1100, 0,  110,Folder);
  HBook1F(Hist->fAlgMask       ,"alg"     ,Form("%s: Algorithm Mask"  ,Folder),   10, 0,   10,Folder);
  HBook1F(Hist->fClusterE      ,"ecl"     ,Form("%s: Cluster E"       ,Folder),  200, 0,  200,Folder);
  
}

//-----------------------------------------------------------------------------
void TTrackRecoEffAnaModule::BookRecoEffHistograms(RecoEffHist_t* Hist, const char* Folder) {
  char name [200];
  char title[200];

  for (int i=0; i<5; i++) {
    sprintf(name ,"ntracks_%i",i);
    sprintf(title,"%s: ntracks[%i]",Folder,i);
    HBook1F(Hist->fNTracks[i],name,title, 10, 0,   10,Folder);

    sprintf(name ,"fcons_%i",i);
    sprintf(title,"%s: fcons[%i]",Folder,i);
    HBook1F(Hist->fFitCons[i],name,title, 1000, 0,   1,Folder);

    sprintf(name ,"nactive_%i",i);
    sprintf(title,"%s: nactive[%i]",Folder,i);
    HBook1F(Hist->fNActive[i],name,title, 200, 0, 200,Folder);

    sprintf(name ,"fitmomerr_%i",i);
    sprintf(title,"%s: fitmomerr[%i]",Folder,i);
    HBook1F(Hist->fFitMomErr[i],name,title, 100, 0, 1,Folder);

    sprintf(name ,"t0err_%i",i);
    sprintf(title,"%s: t0err[%i]",Folder,i);
    HBook1F(Hist->fT0Err[i],name,title, 200, 0,   2,Folder);

    sprintf(name ,"t0_%i",i);
    sprintf(title,"%s: t0[%i]",Folder,i);
    HBook1F(Hist->fT0[i],name,title, 200, 0,   2000,Folder);

    sprintf(name ,"tandip_%i",i);
    sprintf(title,"%s: tandip[%i]",Folder,i);
    HBook1F(Hist->fTanDip[i],name,title, 200, 0,  2,Folder);

    sprintf(name ,"d0_%i",i);
    sprintf(title,"%s: d0[%i]",Folder,i);
    HBook1F(Hist->fD0[i],name,title, 400, -200,  200,Folder);

    sprintf(name ,"rmax_%i",i);
    sprintf(title,"%s: rmax[%i]",Folder,i);
    HBook1F(Hist->fRMax[i],name,title, 200, 0,  1000,Folder);

    sprintf(name ,"p_%i",i);
    sprintf(title,"%s: p[%i]",Folder,i);
    HBook1F(Hist->fP[i],name,title, 1200, 0,  120,Folder);
  }

  HBook1F(Hist->fFailedBits,"failed_bits" ,Form("%s: Failed Bits",Folder),32, 0, 32,Folder);
  HBook1F(Hist->fPassed    ,"passed"      ,Form("%s: passed"     ,Folder),2 , 0, 2 ,Folder);
}

//-----------------------------------------------------------------------------
void TTrackRecoEffAnaModule::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fP,"p"   ,Form("%s: track momentum" ,Folder),  1200, 0,  120,Folder);
}

//_____________________________________________________________________________
void TTrackRecoEffAnaModule::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");

// //-----------------------------------------------------------------------------
// // book crystal histograms
// //-----------------------------------------------------------------------------
//   HBook1F(fHist.fCrystalR[0],"rc_0"     ,Form("disk [0] crystal radius"),100,0,1000,"Hist");
//   HBook1F(fHist.fCrystalR[1],"rc_1"     ,Form("disk [1] crystal radius"),100,0,1000,"Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;		// N(hits CE >= 20)
  book_event_histset[ 2] = 1;		// P(MC) > 100
  book_event_histset[ 3] = 1;		// 1 < pitch < sqrt(3) 
  book_event_histset[ 4] = 1;		// reconstructed track 
  book_event_histset[ 5] = 1;		// quality cuts passed
  book_event_histset[ 6] = 1;		// quality cuts passed
  book_event_histset[ 7] = 1;		// quality cuts passed
  book_event_histset[ 8] = 1;		// quality cuts passed
  book_event_histset[ 9] = 1;		// mom signal window 

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
// book straw hit histograms
//-----------------------------------------------------------------------------
  int book_strawhit_histset[kNStrawHitHistSets];
  for (int i=0; i<kNStrawHitHistSets; i++) book_strawhit_histset[i] = 0;

  book_strawhit_histset[0] = 1;		// all clusters

  for (int i=0; i<kNStrawHitHistSets; i++) {
    if (book_strawhit_histset[i] != 0) {
      sprintf(folder_name,"sh_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fStrawHit[i] = new StrawHitHist_t;
      BookStrawHitHistograms(fHist.fStrawHit[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book Genp histograms
//-----------------------------------------------------------------------------
  int book_genp_histset[kNGenpHistSets];
  for (int i=0; i<kNGenpHistSets; i++) book_genp_histset[i] = 0;

  book_genp_histset[0] = 1;		// all particles

  for (int i=0; i<kNGenpHistSets; i++) {
    if (book_genp_histset[i] != 0) {
      sprintf(folder_name,"gen_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fGenp[i] = new GenpHist_t;
      BookGenpHistograms(fHist.fGenp[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book TRACK histograms
//-----------------------------------------------------------------------------
  int book_track_histset[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

  book_track_histset[0] = 1;		// all tracks
  book_track_histset[1] = 1;		// tracks in events with N(hits CE) >= 20

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
// book RECO EFF histograms
//-----------------------------------------------------------------------------
  int book_recoeff_histset[kNRecoEffHistSets];
  for (int i=0; i<kNRecoEffHistSets; i++) book_recoeff_histset[i] = 0;

  book_recoeff_histset[0] = 1;		// all tracks

  for (int i=0; i<kNRecoEffHistSets; i++) {
    if (book_recoeff_histset[i] != 0) {
      sprintf(folder_name,"eff_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fRecoEff[i] = new RecoEffHist_t;
      BookRecoEffHistograms(fHist.fRecoEff[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TTrackRecoEffAnaModule::FillStrawHitHistograms(StrawHitHist_t* Hist, TStrawHitData* Hit) {

  Hist->fPdgCode->Fill(Hit->PdgCode());
  Hist->fGenCode->Fill(Hit->GeneratorCode());
  Hist->fEnergy->Fill(Hit->Energy());
  Hist->fTime  ->Fill(Hit->Time  ());
  Hist->fDt    ->Fill(Hit->Dt    ());
  Hist->fMcMomentum->Fill(Hit->McMomentum());
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TTrackRecoEffAnaModule::FillEventHistograms(EventHist_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  Hist->fNStrawHits[0]->Fill(fNStrawHits);
  Hist->fNStrawHits[1]->Fill(fNStrawHits);
  Hist->fNHitsSignal->Fill(fNHitsSignal);
  Hist->fMomTF->Fill(fMomTF);
  Hist->fNTracks->Fill(fNTracks);

  Hist->fFitCons->Fill(fFitCons);
  Hist->fT0Err->Fill(fT0Err);
  Hist->fT0->Fill(fT0);
  Hist->fD0->Fill(fD0);
  Hist->fRMax->Fill(fRMax);
  Hist->fNActive->Fill(fNActive);
  Hist->fTanDipMC->Fill(fTanDipMC);
  Hist->fTanDip->Fill(fTanDip);
  Hist->fP->Fill(fP);
  Hist->fAlgMask->Fill(fAlgMask);
  Hist->fClusterE->Fill(fClusterE);
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TTrackRecoEffAnaModule::FillRecoEffHistograms(RecoEffHist_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  fIDWord = 0;

  if (fNTracks   <=   0                      ) fIDWord |= kNTracksBit;
  if (fFitCons   <  2.e-3                    ) fIDWord |= kFitConsBit;
  if (fNActive   <     20                    ) fIDWord |= kNActiveBit;
  if (fFitMomErr >= 0.25                     ) fIDWord |= kFitMomErrBit;
  if (fT0Err     >= 0.9                      ) fIDWord |= kT0ErrBit;
  if (fT0        <  700                      ) fIDWord |= kT0Bit;
  if ( (fTanDip < 1)  || (fTanDip > sqrt(3.))) fIDWord |= kTanDipBit;
  if ( (fD0    < -80) || (fD0 > 105.)        ) fIDWord |= kD0Bit;
  if ( (fRMax < 450)  || (fRMax > 680.)      ) fIDWord |= kRMaxBit;
  if ( (fP < 103.5)   || (fP > 105.0)        ) fIDWord |= kPBit;

  
  Hist->fNTracks[0]->Fill(fNTracks);
  if ((fIDWord & ~kNTracksBit) == 0) Hist->fNTracks[1]->Fill(fNTracks);
  if (fIDWord == 0) Hist->fNTracks[4]->Fill(fNTracks);

  Hist->fFitCons[0]->Fill(fFitCons);
  if ((fIDWord & ~kFitConsBit) == 0) Hist->fFitCons[1]->Fill(fFitCons);
  if (fIDWord == 0) Hist->fFitCons[4]->Fill(fFitCons);

  Hist->fNActive[0]->Fill(fNActive);
  if ((fIDWord & ~kNActiveBit) == 0) Hist->fNActive[1]->Fill(fNActive);
  if (fIDWord == 0) Hist->fNActive[4]->Fill(fNActive);

  Hist->fFitMomErr[0]->Fill(fFitMomErr);
  if ((fIDWord & ~kFitMomErrBit) == 0) Hist->fFitMomErr[1]->Fill(fFitMomErr);
  if (fIDWord == 0) Hist->fFitMomErr[4]->Fill(fFitMomErr);

  Hist->fT0Err[0]->Fill(fT0Err);
  if ((fIDWord & ~kT0ErrBit) == 0) Hist->fT0Err[1]->Fill(fT0Err);
  if (fIDWord == 0) Hist->fT0Err[4]->Fill(fT0Err);

  Hist->fT0[0]->Fill(fT0);
  if ((fIDWord & ~kT0Bit) == 0) Hist->fT0[1]->Fill(fT0);
  if (fIDWord == 0) Hist->fT0[4]->Fill(fT0);

  Hist->fTanDip[0]->Fill(fTanDip);
  if ((fIDWord & ~kTanDipBit) == 0) Hist->fTanDip[1]->Fill(fTanDip);
  if (fIDWord == 0) Hist->fTanDip[4]->Fill(fTanDip);

  Hist->fD0[0]->Fill(fD0);
  if ((fIDWord & ~kD0Bit) == 0) Hist->fD0[1]->Fill(fD0);
  if (fIDWord == 0) Hist->fD0[4]->Fill(fD0);

  Hist->fRMax[0]->Fill(fRMax);
  if ((fIDWord & ~kRMaxBit) == 0) Hist->fRMax[1]->Fill(fRMax);
  if (fIDWord == 0) Hist->fRMax[4]->Fill(fRMax);

  Hist->fP[0]->Fill(fP);
  if ((fIDWord & ~kPBit) == 0) Hist->fP[1]->Fill(fP);
  if (fIDWord == 0) Hist->fP[4]->Fill(fP);

//-----------------------------------------------------------------------------
//  ***** now histogram effect of cuts in the order they are applied
//-----------------------------------------------------------------------------
  Hist->fNTracks[2]->Fill(fNTracks);
  if ((fIDWord & kNTracksBit) != 0)                         goto END;
  Hist->fNTracks[3]->Fill(fNTracks);

  Hist->fFitCons[2]->Fill(fFitCons);
  if ((fIDWord & kFitConsBit) != 0)                         goto END;
  Hist->fFitCons[3]->Fill(fFitCons);

  Hist->fNActive[2]->Fill(fNActive);
  if ((fIDWord & kNActiveBit) != 0)                         goto END;
  Hist->fNActive[3]->Fill(fNActive);

  Hist->fFitMomErr[2]->Fill(fFitMomErr);
  if ((fIDWord & kFitMomErrBit) != 0)                         goto END;
  Hist->fFitMomErr[3]->Fill(fFitMomErr);

  Hist->fT0Err[2]->Fill(fT0Err);
  if ((fIDWord & kT0ErrBit) != 0)                         goto END;
  Hist->fT0Err[3]->Fill(fT0Err);

  Hist->fT0[2]->Fill(fT0);
  if ((fIDWord & kT0Bit) != 0)                         goto END;
  Hist->fT0[3]->Fill(fT0);

  Hist->fTanDip[2]->Fill(fTanDip);
  if ((fIDWord & kTanDipBit) != 0)                         goto END;
  Hist->fTanDip[3]->Fill(fTanDip);

  Hist->fD0[2]->Fill(fD0);
  if ((fIDWord & kD0Bit) != 0)                         goto END;
  Hist->fD0[3]->Fill(fD0);

  Hist->fRMax[2]->Fill(fRMax);
  if ((fIDWord & kRMaxBit) != 0)                         goto END;
  Hist->fRMax[3]->Fill(fRMax);

  Hist->fP[2]->Fill(fP);
  if ((fIDWord & kPBit) != 0)                         goto END;
  Hist->fP[3]->Fill(fP);

 END:
//-----------------------------------------------------------------------------
//  single histogram showing how often every particular cut failed
//-----------------------------------------------------------------------------
  for (int bit=0; bit<32; bit++) {
    if (((fIDWord >> bit) & 0x1) == 1) {
      Hist->fFailedBits->Fill(bit);
    }
  }
  Hist->fPassed->Fill(fIDWord == 0);

}

//-----------------------------------------------------------------------------
void TTrackRecoEffAnaModule::FillGenpHistograms(GenpHist_t* Hist, TGenParticle* Genp) {
//   int    gen_id;
//   float  p, cos_th, z0, t0, r0, x0, y0;

//   TLorentzVector mom, v;

  Hist->fPdgCode[0]->Fill(Genp->GetPdgCode());
  Hist->fPdgCode[1]->Fill(Genp->GetPdgCode());
}

//-----------------------------------------------------------------------------
void TTrackRecoEffAnaModule::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track) {
//   int    gen_id;
//   float  p, cos_th, z0, t0, r0, x0, y0;

//   TLorentzVector mom, v;

//   Hist->fPdgCode[0]->Fill(Genp->GetPdgCode());
//   Hist->fPdgCode[1]->Fill(Genp->GetPdgCode());
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TTrackRecoEffAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("StrawDataBlock"  ,"TStrawDataBlock" ,&fStrawDataBlock);
  RegisterDataBlock("VDetBlock"       ,"TVDetDataBlock"  ,&fVDetDataBlock);
  RegisterDataBlock("GenpBlock"       ,"TGenpBlock"      ,&fGenpBlock);
  RegisterDataBlock("TrackBlock"      ,"TStnTrackBlock"  ,&fTrackBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  printf(" MinNMCHits   : %10i\n"  ,fMinNMCHits);
  printf(" MinMCMomentum: %10.3f\n",fMinMCMomentum);
  printf(" MinMCPitch   : %10.3f\n",fMinMCPitch);
  printf(" MaxMCPitch   : %10.3f\n",fMaxMCPitch);

  return 0;
}



//-----------------------------------------------------------------------------
// event histograms (efficiency, aim to reproduce Dave's numbers)
//-----------------------------------------------------------------------------
void TTrackRecoEffAnaModule::FillHistograms() {
//-----------------------------------------------------------------------------
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);

  if (fNHitsSignal >= fMinNMCHits) {
    FillEventHistograms(fHist.fEvent[1]);

    if (fMomTF > fMinMCMomentum) {
      FillEventHistograms(fHist.fEvent[2]);
//-----------------------------------------------------------------------------
// EVT_3: MC pitch
//-----------------------------------------------------------------------------
     if ((fPitchTF > fMinMCPitch) && (fPitchTF < fMaxMCPitch)) {
	FillEventHistograms(fHist.fEvent[3]);
//-----------------------------------------------------------------------------
// EVT_4: track reconstructed
//-----------------------------------------------------------------------------
	if ((fNTracks > 0) && (fRecoAlgFlag == 1)) { 
	  FillEventHistograms(fHist.fEvent[4]);
//-----------------------------------------------------------------------------
// EVT_5: track passes quality cuts
//-----------------------------------------------------------------------------
	  if ((fFitCons      >= 2.e-3) && 
	      (fNActive      >= 20   ) && 
	      (fTrack->FitMomErr() < 0.25) && 
	      (fT0Err     < 0.9 )    ) {

	    FillEventHistograms(fHist.fEvent[5]);
	    
	    if (fT0 >= 700) {
	      FillEventHistograms(fHist.fEvent[6]);

	      if ((fTanDip > 1) && (fTanDip < sqrt(3.))) {
		FillEventHistograms(fHist.fEvent[7]);

		int rmin_cut_fails = (fD0   <  -80.) || (fD0 >= 105.);
		int rmax_cut_fails = (fRMax <  450.) || (fRMax >= 680.);

		int cosmics_flag (0);

		if (rmin_cut_fails || rmax_cut_fails) cosmics_flag = 1;
		
		if (cosmics_flag == 0) {
		  FillEventHistograms(fHist.fEvent[8]);
		  
		  if ((fP > 103.5) && (fP < 105.0)) {
		    FillEventHistograms(fHist.fEvent[9]);
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
// now the track reco cuts only, the ones which one can apply
//-----------------------------------------------------------------------------
  FillRecoEffHistograms(fHist.fRecoEff[0]);

//-----------------------------------------------------------------------------
// straw hit histograms
//-----------------------------------------------------------------------------
//  int            nh;
  TStrawHitData* hit;

  for (int i=0; i<fNStrawHits; i++) {
    hit = fStrawDataBlock->Hit(i);
    FillStrawHitHistograms(fHist.fStrawHit[0],hit);
  }
//-----------------------------------------------------------------------------
// fill GENP histograms
// GEN_0: all particles
//-----------------------------------------------------------------------------
//   TGenParticle* genp;
//   for (int i=0; i<fNGenp; i++) {
//     genp = fGenpBlock->Particle(i);
//     FillGenpHistograms(fHist.fGenp[0],genp);
//  }
}



//_____________________________________________________________________________
int TTrackRecoEffAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TTrackRecoEffAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;
  TStrawHitData*        hit(NULL);

  fStrawDataBlock->GetEntry(ientry);
  fVDetDataBlock->GetEntry(ientry);
  fTrackBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fNGenp      = fGenpBlock->NParticles();
  fNStrawHits = fStrawDataBlock->NHits();

  fNHitsSignal = 0;

  for (int i=0; i<fNStrawHits; i++) {
    hit = fStrawDataBlock->Hit(i);
    if ((hit->PdgCode() == fPdgCode) && (hit->GeneratorCode() == fGeneratorCode)) {
      fNHitsSignal += 1;
    }
  }
//-----------------------------------------------------------------------------
// define MC particle parameters at the tracker entrance: absolute momentum and
// the pitch angle
// use the first TF hit - in principle, teh same particle may have more than 
// one hit in VD (it can turn back)
//-----------------------------------------------------------------------------
  TVDetHitData* vdhit;

  fNHitsTF   = 0;
  fNHitsTB   = 0;
  fMomTF     = -1.;
  fMomTB     = -1.;
  fPitchTF   = -1.;

  fNVDetHits = fVDetDataBlock->NHits();

  float  px, py, pz;

  for (int i=0; i<fNVDetHits; i++) {
    vdhit = fVDetDataBlock->Hit(i);
    if (vdhit->Index() == 13) {
					// tracker FRONT
      fNHitsTF += 1;
      if (fMomTF < 0) { 
	fMomTF   = vdhit->McMomentum();
	px       = vdhit->McMomentumX();
	py       = vdhit->McMomentumY();
	pz       = vdhit->McMomentumZ();
	fPitchTF = sqrt(px*px+py*py)/pz;
      }
    }
    else if (hit->Index() == 15) {
					// tracker END
      fNHitsTB += 1;
      if (fMomTB < 0) fMomTB = vdhit->McMomentum();
    }
  }

  fTanDipMC = fPitchTF;

  fNTracks   = fTrackBlock->NTracks();
  fTrack     = NULL;
  fNActive   = -1.;
  fFitCons   = -1.;
  fT0        = -1.;
  fT0Err     = -1.;
  fTanDip    = -1.;
  fT0        = -1.;
  fD0        = -1.e6;
  fRMax      = +1.e6;
  fP         = -1.;
  fT0Err     = 1.e6;
  fAlgMask   = -1;
  fClusterE  = -1.;
  fFitMomErr = 1.e6;

  if (fNTracks > 0) {
    fTrack     = fTrackBlock->Track(0);
    fNActive   = fTrack->NActive();
    fFitCons   = fTrack->FitCons();
    fT0        = fTrack->T0();
    fT0Err     = fTrack->T0Err();
    fD0        = fTrack->D0();
    fRMax      = fD0 + 2./fTrack->fC0;
    fP         = fTrack->P();
    fTanDip    = 1./fTrack->TanDip();
    fAlgMask   = fTrack->AlgMask();
    fClusterE  = fTrack->ClusterE();
    fFitMomErr = fTrack->FitMomErr();
  }

//-----------------------------------------------------------------------------
// check if TrkPatRec-only efficiency has been requested
//-----------------------------------------------------------------------------
  fRecoAlgFlag   = 1;
  if (fTrkPatRecOnly) {
    fRecoAlgFlag = (fAlgMask != 2);
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TTrackRecoEffAnaModule::Debug() {

//-----------------------------------------------------------------------------
// bit 3: events with N(tracks) > 1 and first track passing quality cuts
//-----------------------------------------------------------------------------
  if (GetDebugBit(3) == 1) {
    int ntrk = fTrackBlock->NTracks();
    if ((ntrk > 1) && (fFitCons > 2.e-3)) {
      GetHeaderBlock()->Print(Form("NTracks = %5i fFitCons = %10.3f",ntrk,fFitCons));
    }
  }
}

//_____________________________________________________________________________
int TTrackRecoEffAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TTrackRecoEffAnaModule::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}

