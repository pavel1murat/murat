//////////////////////////////////////////////////////////////////////////////
// use of TStnTrack::fTmp:
//
// Tmp(0) : corrected momentum at the tracker front - not yet
// 
// use of debug bits: 
//  0  : one line per track
//  1  : passed events
//  2  : rejected events
// 
//  3  : events with good tracks P > 70 and T0 > 300
//  5  : events with N(tracks) > 1
//
// 3 different ID : 
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TDirectory.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------

#include "murat/ana/TMomscaleAnaModule.hh"

ClassImp(murat::TMomscaleAnaModule)

namespace murat {
//-----------------------------------------------------------------------------
TMomscaleAnaModule::TMomscaleAnaModule(const char* name, const char* title): 
  TAnaModule(name,title) { 
  fTrackBlockName         = "TrackBlock";
  fTrackStrawHitBlockName = "TrackStrawHitBlock";
  fTrackNumber.Set(100);

  //  fDiskCalorimeter = new TDiskCalorimeter();
  fCalorimeterType = 2;

  fMinT0           = 700; 
  //  fMinETrig        = 50.;               // MeV
					// track-cluster matching timing cut
  // fMinDtTcm        = -5.;
  // fMaxDtTcm        =  8.;
//-----------------------------------------------------------------------------
// initialize Track ID
// 0: SetC  1: TrkQual>0.1 2:TrkQual>0.4
// what about number of hits ? - 3: no cuts on the number of hits
//-----------------------------------------------------------------------------
  fNID             = 2;

  fTrackID[0]      = new TStnTrackID();

  int mask0 = TStnTrackID::kNActiveBit | TStnTrackID::kChi2DofBit | TStnTrackID::kMomErrBit ;
  fTrackID[0]->SetMaxChi2Dof(3.0 );
  fTrackID[0]->SetMinNActive(25  );
  fTrackID[0]->SetMaxMomErr (0.25);

  mask0    |= TStnTrackID::kTanDipBit;
  fTrackID[0]->SetMinTanDip (0.6);
  fTrackID[0]->SetMaxTanDip (1.0);

  mask0    |= TStnTrackID::kD0Bit;
  fTrackID[0]->SetMinD0 (-100);         // mm
  fTrackID[0]->SetMaxD0 ( 100);

  fTrackID[0]->SetUseMask(mask0);
//-----------------------------------------------------------------------------
// tight ID
//-----------------------------------------------------------------------------
  fTrackID[1]      = new TStnTrackID();

  int mask1 = TStnTrackID::kNActiveBit | TStnTrackID::kChi2DofBit | TStnTrackID::kMomErrBit ;
  fTrackID[1]->SetMaxChi2Dof(3.0 );
  fTrackID[1]->SetMinNActive(30  );     // cut tighter on the number of hits
  fTrackID[1]->SetMaxMomErr (0.25);

  mask1    |= TStnTrackID::kTanDipBit;
  fTrackID[1]->SetMinTanDip (0.65);
  fTrackID[1]->SetMaxTanDip (1.20);

  mask1    |= TStnTrackID::kD0Bit;
  fTrackID[1]->SetMinD0 ( -50);         // mm
  fTrackID[1]->SetMaxD0 ( 100);

  fTrackID[1]->SetUseMask(mask1);

  fBestID           = 0;

  for (int i=0; i<20; i++) {
    TrackPar_t* tp  = fTrackPar+i;
    tp->fTrackID[0] = TAnaModule::fTrackID[0];
    tp->fTrackID[1] = TAnaModule::fTrackID[1];
  }

//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
// 2:conversionGun, 28:StoppedParticleReactionGun - see 
//-----------------------------------------------------------------------------
  fDirection        = 1;

  fDnMax            = 15;
//-----------------------------------------------------------------------------
// this is redefiend by the dataset catalog anyway
//-----------------------------------------------------------------------------
  fPDGCode          = -11;
  fMCProcessCode    = 181;   // mu2e::ProcessCode::mu2ePienu;
}

//-----------------------------------------------------------------------------
TMomscaleAnaModule::~TMomscaleAnaModule() {
  for (int i=0; i<fNID; i++) delete fTrackID[i];
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TMomscaleAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockName.Data()        ,"TStnTrackBlock"     ,&fTrackBlock  );
  // RegisterDataBlock(fTrackStrawHitBlockName.Data(),"TTrackStrawHitBlock",&fTrackStrawHitBlock);

  RegisterDataBlock("KSFBlock"        ,"TStnTrackSeedBlock"  ,&fTrackSeedBlock  );
  RegisterDataBlock("HelixBlock"      ,"TStnHelixBlock"      ,&fHelixBlock      );
  RegisterDataBlock("TimeClusterBlock","TStnTimeClusterBlock",&fTimeClusterBlock);

  RegisterDataBlock("ClusterBlock" ,"TStnClusterBlock" ,&fClusterBlock );
  RegisterDataBlock("CalDataBlock" ,"TCalDataBlock"    ,&fCalDataBlock );
  // RegisterDataBlock("StrawHitBlock","TStrawHitBlock"   ,&fStrawHitBlock);
  RegisterDataBlock("GenpBlock"    ,"TGenpBlock"       ,&fGenpBlock    );
  RegisterDataBlock("SimpBlock"    ,"TSimpBlock"       ,&fSimpBlock    );
  RegisterDataBlock("SpmcBlockVDet","TStepPointMCBlock",&fVDetBlock    );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}

//-----------------------------------------------------------------------------
void TMomscaleAnaModule::BookMomscaleHistograms(MomscaleHist_t* Hist, const char* Folder) {
  HBook2F(Hist->fEpVsP ,"ep_vs_p",Form("%s: E/p vs P "       ,Folder), 120, 0,  120, 110, 0, 1.1,Folder);
}

//-----------------------------------------------------------------------------
void TMomscaleAnaModule::BookTimeClusterHistograms(TimeClusterHist_t* Hist, const char* Folder) {
  HBook1F(Hist->fNsh  ,"nsh"  ,Form("%s: N(single straw hits)",Folder), 200, 0,  200,Folder);
  HBook1F(Hist->fNch  ,"nch"  ,Form("%s: N(combo hits)"       ,Folder), 200, 0,  200,Folder);
  HBook1F(Hist->fT0   ,"t0"   ,Form("%s: T0"                  ,Folder), 200, 0, 2000,Folder);
  HBook1F(Hist->fT0Err,"t0err",Form("%s: T0Err"               ,Folder), 200, 0, 2000,Folder);
}

//_____________________________________________________________________________
void TMomscaleAnaModule::BookHistograms() {

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
  // HBook1F(fHist.fCrystalR[0],"rc_0"     ,Form("disk [0] crystal radius"),100,0,1000,"Hist");
  // HBook1F(fHist.fCrystalR[1],"rc_1"     ,Form("disk [1] crystal radius"),100,0,1000,"Hist");
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
// book helix histograms
//-----------------------------------------------------------------------------
  int book_helix_histset[kNHelixHistSets];
  for (int i=0; i<kNHelixHistSets; i++) book_helix_histset[i] = 0;

  book_helix_histset[ 0] = 1;		// all events
  book_helix_histset[ 1] = 1;		// events with N(reconstructed tracks) > 0
  book_helix_histset[ 2] = 1;		// helices with a TrackSeed

  for (int i=0; i<kNHelixHistSets; i++) {
    if (book_helix_histset[i] != 0) {
      sprintf(folder_name,"hel_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fHelix[i] = new HelixHist_t;
      BookHelixHistograms(fHist.fHelix[i],Form("Hist/%s",folder_name));
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

  book_track_histset[  0] = 1;		// all tracks
  book_track_histset[  1] = 1;		// fIDWord[0]=0 tracks
  book_track_histset[  2] = 1;		// fIDWord[0]=0 tracks tdip<0.75
  book_track_histset[  3] = 1;		// fIDWord[0]=0 tracks tdip>0.75
  book_track_histset[  4] = 1;		// fIDWord[0]=0 tracks p>50
  book_track_histset[  5] = 1;		// fIDWord[0]=0 tracks 45<p<50
  book_track_histset[  6] = 1;		// fIDWord[0]=0 tracks 40<p<45

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
// book track histograms
//-----------------------------------------------------------------------------
  int book_track_id_histset[kNTrackIDHistSets];
  for (int i=0; i<kNTrackIDHistSets; i++) book_track_id_histset[i] = 0;

  book_track_id_histset[  0] = 1;		// all tracks

  for (int i=0; i<kNTrackIDHistSets; i++) {
    if (book_track_id_histset[i] == 0)                      continue;
    sprintf(folder_name,"tid_%i",i);
    fol = (TFolder*) hist_folder->FindObject(folder_name);
    if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
    fHist.fTrackID[i] = new TStnTrackID::Hist_t;
    BookTrackIDHistograms(fHist.fTrackID[i],Form("Hist/%s",folder_name));
  }
//-----------------------------------------------------------------------------
// book momscale histograms
//-----------------------------------------------------------------------------
  int book_momscale_histset[kNMomscaleHistSets];
  for (int i=0; i<kNMomscaleHistSets; i++) book_momscale_histset[i] = 0;

  book_momscale_histset[  0] = 1;		// all          tracks
  book_momscale_histset[  1] = 1;		// fIDWord[0]=0 tracks
  book_momscale_histset[  2] = 1;		// fIDWord[0]=0 tracks tdip>0.75
  book_momscale_histset[  3] = 1;		// fIDWord[0]=0 tracks tdip<0.75
  book_momscale_histset[  4] = 1;		// fIDWord[0]=0 tracks p>50
  book_momscale_histset[  5] = 1;		// fIDWord[0]=0 tracks 45<p<50
  book_momscale_histset[  6] = 1;		// fIDWord[0]=0 tracks 40<p<45

  for (int i=0; i<kNMomscaleHistSets; i++) {
    if (book_momscale_histset[i] == 0)                        continue;
    sprintf(folder_name,"pip_%i",i);
    fol = (TFolder*) hist_folder->FindObject(folder_name);
    if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
    fHist.fMomscale[i] = new MomscaleHist_t;
    BookMomscaleHistograms(fHist.fMomscale[i],Form("Hist/%s",folder_name));
  }
//-----------------------------------------------------------------------------
// book time cluster histograms
//-----------------------------------------------------------------------------
  int book_tc_histset[kNTimeClusterHistSets];
  for (int i=0; i<kNTimeClusterHistSets; i++) book_tc_histset[i] = 0;

  book_tc_histset[0] = 1;		// all clusters

  for (int i=0; i<kNTimeClusterHistSets; i++) {
    if (book_tc_histset[i] == 0)                            continue;
    sprintf(folder_name,"tc_%i",i);
    fol = (TFolder*) hist_folder->FindObject(folder_name);
    if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
    fHist.fTimeCluster[i] = new TimeClusterHist_t;
    BookTimeClusterHistograms(fHist.fTimeCluster[i],Form("Hist/%s",folder_name));
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
}


//-----------------------------------------------------------------------------
// fill efficiency histograms : need 10 histogram sets
// pitch = 1./tan(dip)
//-----------------------------------------------------------------------------
void TMomscaleAnaModule::FillEfficiencyHistograms(TStnTrackBlock*  TrackBlock, 
                                                TStnTrackID*     TrackID   , 
                                                int              HistSet   ) {
  // having MC truth is a must for calculating efficiency!
  if (fEvtPar.fSimp == nullptr) return;

  if (fEvtPar.fSimp->NStrawHits() >= 20) {
    FillEventHistograms(fHist.fEvent[HistSet],&fEvtPar);

    if (fEvtPar.fSimp->fMomTrackerFront > 100.) {
      FillEventHistograms(fHist.fEvent[HistSet+1],&fEvtPar);

      TVector3 vdmom;

      if (fSimPar.fTFront != NULL) {
//-----------------------------------------------------------------------------
// regular case
//-----------------------------------------------------------------------------
	vdmom.SetXYZ(fSimPar.fTFront->Mom()->X(),
		     fSimPar.fTFront->Mom()->Y(),		      
		     fSimPar.fTFront->Mom()->Z());
      }
      else {
//-----------------------------------------------------------------------------
// pathological case when an upstream simulated MC particle is reconstructed 
// as the downstream one and there is no hit... efficiency doens't make sense
// just make sure the code doesnt' crash
//-----------------------------------------------------------------------------
	vdmom.SetXYZ(fSimPar.fParticle->fStartMom.X(),
		     fSimPar.fParticle->fStartMom.Y(),
		     fSimPar.fParticle->fStartMom.Z());
      }

      float ce_pitch  = vdmom.Pt()/vdmom.Pz();
      float min_pitch = 1./TrackID->MaxTanDip();
      float max_pitch = 1./TrackID->MinTanDip();

      if ((min_pitch < ce_pitch) && (ce_pitch < max_pitch)) {
	FillEventHistograms(fHist.fEvent[HistSet+2],&fEvtPar);
	  
	if (TrackBlock->NTracks() > 0) {
	  TStnTrack* track = TrackBlock->Track(0);
	  int id_word      = TrackID->IDWord(track);

	  FillEventHistograms(fHist.fEvent[HistSet+3],&fEvtPar);
	  
	  if ((id_word & TStnTrackID::kTrkQualBit) == 0) {
	    FillEventHistograms(fHist.fEvent[HistSet+4],&fEvtPar);
	    
	    if ((id_word & TStnTrackID::kT0Bit) == 0) {
	      FillEventHistograms(fHist.fEvent[HistSet+5],&fEvtPar);
	      
	      if ((id_word & TStnTrackID::kTanDipBit) == 0) {
		FillEventHistograms(fHist.fEvent[HistSet+6],&fEvtPar);
		
		if (((id_word & TStnTrackID::kD0Bit  ) == 0) && 
		    ((id_word & TStnTrackID::kRMaxBit) == 0)    ) {
		  
		  FillEventHistograms(fHist.fEvent[HistSet+7],&fEvtPar);
		  
		  if ((id_word & TStnTrackID::kTanDipBit) == 0) {
		    FillEventHistograms(fHist.fEvent[HistSet+8],&fEvtPar);

		    if ((103.5 < track->fP) && (track->fP < 105)) {
		      FillEventHistograms(fHist.fEvent[HistSet+9],&fEvtPar);
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
void    TMomscaleAnaModule::FillMomscaleHistograms    (MomscaleHist_t*  Hist, 
                                                       TStnTrack*       Trk, 
                                                       TrackPar_t*      Tp, 
                                                       SimPar_t*        SimPar,
                                                       double           Weight) {
  Hist->fEpVsP->Fill(Tp->fP,Tp->fEp,Weight);
}

//-----------------------------------------------------------------------------
void    TMomscaleAnaModule::FillTimeClusterHistograms(TimeClusterHist_t* Hist  , 
                                                      TStnTimeCluster*   Tc    ,
                                                      double             Weight) {
  Hist->fNsh->Fill  (Tc->NHits()     ,Weight);
  Hist->fNch->Fill  (Tc->NComboHits(),Weight);
  Hist->fT0->Fill   (Tc->T0()        ,Weight);
  Hist->fT0Err->Fill(Tc->T0Err()     ,Weight);
}

//_____________________________________________________________________________
void TMomscaleAnaModule::FillHistograms() {

  // double       cos_th (-2.); //,  cl_e(-1.);
  // int          disk_id(-1);
//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);
  if (fEvtPar.fNTracks[0]> 0) FillEventHistograms(fHist.fEvent[1],&fEvtPar);
//-----------------------------------------------------------------------------
// helix histograms
// HEL_0: all
// HEL_1: helices in events with reconstructed tracks
//-----------------------------------------------------------------------------
  int nhel = fHelixBlock->NHelices();
  for (int i=0; i<nhel; i++) {
    TStnHelix*  hel = fHelixBlock->Helix(i);
    HelixPar_t* hp  = fHelixPar+i;
    FillHelixHistograms(fHist.fHelix[0],hel,hp);
    if (fEvtPar.fNTracks[0]  >  0) FillHelixHistograms(fHist.fHelix[1],hel,hp);
    if (hel->fTrackSeedIndex >= 0) FillHelixHistograms(fHist.fHelix[2],hel,hp);
  }
//-----------------------------------------------------------------------------
// SIMP histograms
//-----------------------------------------------------------------------------
  if (fEvtPar.fSimp) {
    SimpData_t sd;

    // TSimParticle* mom = fSimpBlock->Particle(fEvtPar.fSimp->ParentID());

    FillSimpHistograms(fHist.fSimp[0],fEvtPar.fSimp,&sd);
  }
//-----------------------------------------------------------------------------
// track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
  TStnTrack*   trk;
  TrackPar_t*  tp;

  for (int i=0; i<fEvtPar.fNTracks[0]; ++i ) {
    trk       = fTrackBlock->Track(i);
    tp        = fTrackPar+i;

    fTrackID[0]->FillHistograms(fHist.fTrackID[0],trk);

    tp->fDnTrackTc = tp->fTimeCluster->NHits()-trk->NActive();
    tp->fDtTrackTc = tp->fTimeCluster->NHits()-trk->NActive();
//-----------------------------------------------------------------------------
// find closest time cluster 
//-----------------------------------------------------------------------------
    float dt_tc_min = 1.e6;   // signed !!!

    for (int itc=0; itc<fEvtPar.fNTimeClusters; itc++) {
      TStnTimeCluster* tc = fTimeClusterBlock->TimeCluster(itc);
      if (tc == tp->fTimeCluster)                 continue;
      float dt = tp->fTimeCluster->T0()-tc->T0();
      if (fabs(dt) < fabs(dt_tc_min)) {
        dt_tc_min = dt;
      }
    }

    tp->fDtTcTc = dt_tc_min;

    FillTrackHistograms   (fHist.fTrack   [  0],trk,tp,&fSimPar);
    FillMomscaleHistograms(fHist.fMomscale[  0],trk,tp,&fSimPar);

    if (tp->fIDWord[fBestID] == 0) {
					// track passes the quality selection

      FillTrackHistograms   (fHist.fTrack   [  1],trk,tp,&fSimPar);
      FillMomscaleHistograms(fHist.fMomscale[  1],trk,tp,&fSimPar);

      if (trk->TanDip() < 0.75) {
        FillTrackHistograms   (fHist.fTrack   [  2],trk,tp,&fSimPar);
        FillMomscaleHistograms(fHist.fMomscale[  2],trk,tp,&fSimPar);
      }
      else {
        FillTrackHistograms   (fHist.fTrack   [  3],trk,tp,&fSimPar);
        FillMomscaleHistograms(fHist.fMomscale[  3],trk,tp,&fSimPar);
      }

      if (tp->fP > 50) {
        FillTrackHistograms   (fHist.fTrack   [  4],trk,tp,&fSimPar);
        FillMomscaleHistograms(fHist.fMomscale[  4],trk,tp,&fSimPar);
      }
      else if ((tp->fP > 45) and (tp->fP < 50)) {
        FillTrackHistograms   (fHist.fTrack   [  5],trk,tp,&fSimPar);
        FillMomscaleHistograms(fHist.fMomscale[  5],trk,tp,&fSimPar);
      }
      else if ((tp->fP > 40) and (tp->fP < 45)) {
        FillTrackHistograms   (fHist.fTrack   [  6],trk,tp,&fSimPar);
        FillMomscaleHistograms(fHist.fMomscale[  6],trk,tp,&fSimPar);
      }
    }
  }
//-----------------------------------------------------------------------------
// track reco efficiency - filles Event histograms
//-----------------------------------------------------------------------------
  FillEfficiencyHistograms(fTrackBlock,TAnaModule::fTrackID[fBestID],11);
//-----------------------------------------------------------------------------
// fill GENP histograms
// GEN_0: all particles
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  for (int i=0; i<fNGenp; i++) {
    genp = fGenpBlock->Particle(i);
    FillGenpHistograms(fHist.fGenp[0],genp);
  }
}

//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TMomscaleAnaModule::Event(int ientry) {

  //  TDiskCalorimeter::GeomData_t disk_geom;

  fTrackBlock->GetEntry(ientry);
  fTrackSeedBlock->GetEntry(ientry);
  fHelixBlock->GetEntry(ientry);
  fTimeClusterBlock->GetEntry(ientry);
  // fClusterBlock->GetEntry(ientry);
  // fCalDataBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fVDetBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// at some point, the TVdetBlock class got obsolete, and now the virtual detector 
// hits are stored in TStepPointMCBlock 
// dont' try to read it as that would fail
//-----------------------------------------------------------------------------
  fEvtPar.fDioLOWt          =  1.;
  fEvtPar.fDioLLWt          =  1.;
  fEvtPar.fNCrvClusters     = -1;
  fEvtPar.fNCrvPulses       = -1;
  fEvtPar.fNCrvCoincidences = -1;
  fEvtPar.fNStrawHits       = GetHeaderBlock()->NStrawHits();
  fEvtPar.fNComboHits       = GetHeaderBlock()->NComboHits();
  fEvtPar.fNGenp            = fGenpBlock->NParticles();
  fEvtPar.fNSimp            = fSimpBlock->NParticles();
  fEvtPar.fSimp             = nullptr;
  fEvtPar.fPartE            = -1.;
  fEvtPar.fPionSurvProb     = 1.;
  fEvtPar.fNHelices         = fHelixBlock->NHelices();
  fEvtPar.fNTimeClusters    = fTimeClusterBlock->NTimeClusters();
//-----------------------------------------------------------------------------
// luminosity weight
//-----------------------------------------------------------------------------
  fLumWt = GetHeaderBlock()->LumWeight();
//-----------------------------------------------------------------------------
// figure the MC truth
// pi+ --> e+ nu case : determine the event weight
//-----------------------------------------------------------------------------
  for (int i=fEvtPar.fNSimp-1; i>=0; i--) {
    TSimParticle* simp = fSimpBlock->Particle(i);
    int pdg_code       = simp->PDGCode();
    int generator_id   = simp->GeneratorID();         // MC process code

    if ((pdg_code == fPDGCode) and (generator_id == fMCProcessCode)) {
      fEvtPar.fSimp  = simp;
      fEvtPar.fPartE = simp->StartMom()->Energy();
    }

    if ((abs(pdg_code) == 211) && (simp->GeneratorID() == 56)) {
//-----------------------------------------------------------------------------
// found the pion, survival probability
//-----------------------------------------------------------------------------
      fEvtPar.fPionSurvProb  = exp(-simp->fEndProperTime);
      break;
    }
  }
//-----------------------------------------------------------------------------
// calculate DIO weights once per event
//-----------------------------------------------------------------------------
  if (fEvtPar.fPartE > 0) {
    fEvtPar.fDioLOWt = TStntuple::DioWeightAl   (fEvtPar.fPartE);
    fEvtPar.fDioLLWt = TStntuple::DioWeightAl_LL(fEvtPar.fPartE);
  }
//-----------------------------------------------------------------------------
// may want to revisit the definition of fSimp in the future
//-----------------------------------------------------------------------------
  fSimPar.fParticle = fEvtPar.fSimp;
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
//-----------------------------------------------------------------------------
// process virtual detectors - for fSimp need parameters at tracker entrance
// use the first fit
//-----------------------------------------------------------------------------
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  int nvdhits = fVDetBlock->NStepPoints();
  for (int i=0; i<nvdhits; i++) {
    TStepPointMC* vdhit = fVDetBlock->StepPointMC(i);
    if (vdhit->PDGCode() == fEvtPar.fSimp->fPdgCode) {
      if ((vdhit->VolumeID() == 13) || (vdhit->VolumeID() == 14)) {
	if (fDirection*vdhit->Mom()->Z() > 0) {
	  if (fSimPar.fTFront == 0) fSimPar.fTFront = vdhit;
	}
      }
      else if ((vdhit->VolumeID() == 11) || (vdhit->VolumeID() == 12)) {
	if (fDirection*vdhit->Mom()->Z() > 0) {
	  if (fSimPar.fTMid == 0) fSimPar.fTMid = vdhit;
	}
      }
      else if (vdhit->VolumeID() == mu2e::VirtualDetectorId::TT_Back) {
	if (fDirection*vdhit->Mom()->Z() > 0) {
	  if (fSimPar.fTBack == NULL) fSimPar.fTBack = vdhit;
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// init calorimeter
//-----------------------------------------------------------------------------
  fNClusters  = fClusterBlock->NClusters();
  fNCalHits   = fCalDataBlock->NHits();
//-----------------------------------------------------------------------------
// init calorimeter clusters - remember, the first one not necessarily is the 
// most energetic
//-----------------------------------------------------------------------------
  fNClusters = fClusterBlock->NClusters();
  if (fNClusters == 0) fCluster = 0;
  else                 fCluster = fClusterBlock->Cluster(0);
//-----------------------------------------------------------------------------
// loop over tracks and calculate needed parameters
//-----------------------------------------------------------------------------
  fNStrawHits          = GetHeaderBlock()->fNStrawHits;

  fEvtPar.fNTracks[0]  = fTrackBlock->NTracks();
  if (fEvtPar.fNTracks[0] == 0) fTrack = 0;
  else                          fTrack = fTrackBlock->Track(0);

  fEvtPar.fNGoodTracks[0] = 0;
  fEvtPar.fNGoodTracks[1] = 0;
  fNMatchedTracks         = 0;
//-----------------------------------------------------------------------------
// determine the number of CalPatRec tracks - this assumes that the list of 
// tracks has been created by MergePatRec
//-----------------------------------------------------------------------------
  int ntrk = fTrackBlock->NTracks();

  InitTrackPar(fTrackBlock,fClusterBlock,fTrackPar,&fSimPar);
//-----------------------------------------------------------------------------
// additional initializations - helices and time clusters
//-----------------------------------------------------------------------------
  fN300 = 0;
  for (int itrk=0; itrk<ntrk; itrk++) {
    TStnTrack*   trk = fTrackBlock->Track(itrk);
    TrackPar_t*  tp  = fTrackPar+itrk;

    tp->fHelix       = nullptr;
    tp->fTimeCluster = nullptr;

    int ih = trk->fHelixIndex;
    if (ih < fEvtPar.fNHelices) {
      tp->fHelix = fHelixBlock->Helix(ih);
    }
//-----------------------------------------------------------------------------
// fN300 : the number of late (T0>300) good tracks above 60 MeV/c
//-----------------------------------------------------------------------------
    if ((tp->fIDWord[fBestID] == 0) and (tp->fP > 60) and (trk->T0() > 300)) { 
      fN300++;
    }
//-----------------------------------------------------------------------------
// currently, at least in the KK case, an TStnHelix doesn't have a helix index stored
// find a time cluster closest to teh track in time
//-----------------------------------------------------------------------------
    TStnTimeCluster* closest_tc(nullptr);
    float min_dt = 1.e6;
    for (int itc=0; itc<fEvtPar.fNTimeClusters; itc++) {
      TStnTimeCluster* tc = fTimeClusterBlock->TimeCluster(itc);
      float dt = trk->T0()-tc->T0();
      if (fabs(dt) < fabs(min_dt)) {
        closest_tc = tc;
        min_dt     = dt;
      }
    }

    tp->fTimeCluster = closest_tc;
    tp->fDtTrackTc   = min_dt;
  }

  fEventPassedSelections = 0;

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TMomscaleAnaModule::Debug() {

  TStnTrack* trk;
  TrackPar_t* tp;
  int ntrk = fTrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = fTrackBlock->Track(itrk);
    tp  = &fTrackPar[itrk];
//-----------------------------------------------------------------------------
// bit 0: all tracks
//-----------------------------------------------------------------------------
    if (GetDebugBit(0) == 1) {
      GetHeaderBlock()->Print(Form("bit_000: All p = %10.3f",trk->Momentum()->P()));
    }
//-----------------------------------------------------------------------------
// bit 1: IDWord =0 0 tracks
//-----------------------------------------------------------------------------
    if (GetDebugBit(1) == 1) {
      if (tp->fIDWord[fBestID] == 0) {
	GetHeaderBlock()->Print(Form("bit_001: IDWord=0 p = %10.3f",trk->Momentum()->P()));
      }
    }
//-----------------------------------------------------------------------------
// bit 3: good tracks with P > 70 T0 > 300
//-----------------------------------------------------------------------------
    if (GetDebugBit(3) == 1) {
      if (fEventPassedSelections == 1) {
        if ((tp->fP > 70) and (trk->fT0 > 300)) {
	  GetHeaderBlock()->Print(Form("large P: %f",tp->fP));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 4: good events with 65 < P < 70 , T0 > 300
//-----------------------------------------------------------------------------
    if (GetDebugBit(4) == 1) {
      if (tp->fIDWord[fBestID] == 0) {
        if ((tp->fDnTrackTc <= fDnMax) and (fabs(tp->fDtTcTc) > 100)) {
          if ((tp->fP > 65) and (tp->fP < 700) and (trk->fT0 > 300)) {
            GetHeaderBlock()->Print(Form(":bit_04: candidate: i=%i trk->P0=%10.3f tp->fP=%10.3f",itrk,trk->P0(),tp->fP));
          }
        }
      }
    }
  }
//-----------------------------------------------------------------------------
// bit 5: events with N(tracks) > 1
//-----------------------------------------------------------------------------
  if (GetDebugBit(5) == 1) {
    int ntrk = fTrackBlock->NTracks();
    if (ntrk > 1) {
      GetHeaderBlock()->Print(Form("NTracks = %5i",ntrk));
    }
  }
}

//_____________________________________________________________________________
int TMomscaleAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TMomscaleAnaModule::Test001() {
}

}
