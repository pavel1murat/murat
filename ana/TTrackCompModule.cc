//////////////////////////////////////////////////////////////////////////////
// compare track fitting with two different ambiguity resolvers
//
// use of tmp:
//
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//  3  : events with KDAR tracks with DPF > 3
//  4  : events with KPAR track and a 50+ MeV cluster, but with no KDAR track
//  5  : events with KDAR tracks with P > 106
//  6  : events with KDAR tracks with 1.5 < DPF < 5
//  7  : events with N(set "C" KDAR tracks)  > 0
//  8  : events with KDAR tracks with P > 105
//  9  : all events - print DMVA
// 10  : events with KPAR DMVA > 0.3 - print DMVA
// 11  : validate MergePatRec
// 12  : print events with dtZ0 < -10
// 13  : events with KDAR tracks with DPF > 10
// 14  : misreconstructed events in dsid=e11s721z
//
// call: "track_comp(28,4)
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/base/TStnDataset.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
#include "murat/ana/mva_data.hh"

// framework
#include "fhiclcpp/ParameterSet.h"

// Xerces XML Parser
#include <xercesc/dom/DOM.hpp>

#include "Mu2eUtilities/inc/MVATools.hh"
#include <string>
#include "math.h"

#include "murat/ana/TTrackCompModule.hh"

using std::string;
using std::vector;

ClassImp(murat::TTrackCompModule)

namespace murat {
//-----------------------------------------------------------------------------
TTrackCompModule::TTrackCompModule(const char* name, const char* title): TAnaModule(name,title) {

  fPdgCode           = 11;		// electron
  fGeneratorCode     =  2; // 2:ConversionGun 28:StoppedParticleReactionGun

  fTrackBlockName[0] = "TrackBlockPar";
  fTrackBlockName[1] = "TrackBlockDar";
//-----------------------------------------------------------------------------
// TrackID[0] : "SetC"
// i = 1..6 : cut on DaveTrkQual > 0.1*i instead
// fTrackID[1-19] - check MVA
//-----------------------------------------------------------------------------
  fNID  = 20;  // this is the limit ....
  for (int i=0; i<fNID; i++) {
    fTrackID[i] = new TStnTrackID();
    if (i > 0) {
      fTrackID[i]->SetMaxMomErr (100);
      fTrackID[i]->SetMaxT0Err  (100);
      fTrackID[i]->SetMinNActive( -1);
      fTrackID[i]->SetMinFitCons(-1.);
      fTrackID[i]->SetMinTrkQual(0.05*i);
    }
  }
//-----------------------------------------------------------------------------
// TStntuple 
//-----------------------------------------------------------------------------
  fStnt = TStntuple::Instance();

//-----------------------------------------------------------------------------
// track ID for RMC background estimates
//-----------------------------------------------------------------------------
  fTrackID_RMC = new TStnTrackID();

  fTrackID_RMC->SetMaxChi2Dof(4. );
  // fTrackID_RMC->SetMaxT0Err  (1.5);
  // fTrackID_RMC->SetMaxMomErr (0.4);
  fTrackID_RMC->SetMaxT0Err  (1.3);
  fTrackID_RMC->SetMaxMomErr (0.3);
  fTrackID_RMC->SetMinNActive(20 );
  fTrackID_RMC->SetMaxDNa    ( 5 );
  fTrackID_RMC->SetMinTanDip (1./sqrt(3.));
  fTrackID_RMC->SetMaxTanDip (1.5);
  fTrackID_RMC->SetMinD0     (-100.);
  fTrackID_RMC->SetMaxD0     ( 100.);

  int mask = TStnTrackID::kNActiveBit | TStnTrackID::kChi2DofBit | TStnTrackID::kT0Bit     | 
             TStnTrackID::kT0ErrBit   | TStnTrackID::kMomErrBit  | TStnTrackID::kTanDipBit | 
             TStnTrackID::kD0Bit      | TStnTrackID::kRMaxBit    | TStnTrackID::kDNaBit    ;

  fTrackID_RMC->SetUseMask(mask);

  fKMaxRMC        = 90.;
//-----------------------------------------------------------------------------
// 000 : no Z, use chi2d 
//-----------------------------------------------------------------------------
  fMinETrig       = 50.;
  fLogLH          = new TEmuLogLH();
//-----------------------------------------------------------------------------
// debugging information
//-----------------------------------------------------------------------------
  fDebugCut[5].fXMin = 106.;
  fDebugCut[5].fXMax = 200.;

  fDebugCut[6].fXMin = 1.5;
  fDebugCut[6].fXMax = 10.0;
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  fMbTime            = 1695.;
  fFillHistograms    = 1;
}

//-----------------------------------------------------------------------------
TTrackCompModule::~TTrackCompModule() {
  delete fLogLH;
  for (int i=0; i<fNID; i++) delete fTrackID[i];

  delete fTprMVA;
  delete fCprMVA;
}

//-----------------------------------------------------------------------------
// TrkRecAlgorithm : "cpr" or "tpr"
// TrainingDataset : just 'fele2s51b1' - goes into a file name
// MVAType : 
// ---------
// 000-004 : log(fcons) with different weights 
// 100-104 : chi2d      with different weights
// 202     : chi2       training by Arpan
//-----------------------------------------------------------------------------
void TTrackCompModule::SetMVA(const char* TrkRecAlgorithm, const char* TrainingDataset, int MvaType) {

  printf(" [TTrackCompModule::SetMVA] TrkRecAlgorithm:%s  TrainingDataset:%s MvaType:%i\n",
	 TrkRecAlgorithm,TrainingDataset,MvaType);

  fUseMVA         = 1;

  TString trk_alg = TrkRecAlgorithm;
  trk_alg.ToUpper();

  if (trk_alg == "CPR") {
    if (fCprMVA) delete fCprMVA;
    fCprMVA           = new mva_data("CPR",TrainingDataset,MvaType);
    fTmvaAlgorithmCpr = MvaType;
  }
  else if (trk_alg == "TPR") {
    if (fTprMVA) delete fTprMVA;
    fTprMVA           = new mva_data("TPR",TrainingDataset,MvaType);
    fTmvaAlgorithmTpr = MvaType;
  }
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TTrackCompModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockName[0].Data(), "TStnTrackBlock"   , &fTrackBlock[0] );
  RegisterDataBlock(fTrackBlockName[1].Data(), "TStnTrackBlock"   , &fTrackBlock[1] );
  RegisterDataBlock("ClusterBlock"           , "TStnClusterBlock" , &fClusterBlock  );
  RegisterDataBlock("SimpBlock"              , "TSimpBlock"       , &fSimpBlock     );
  RegisterDataBlock("GenpBlock"              , "TGenpBlock"       , &fGenpBlock     );
  RegisterDataBlock("TrackSeedBlock"         ,"TStnTrackSeedBlock", &fTrackSeedBlock);
  RegisterDataBlock("HelixBlock"             , "TStnHelixBlock"   , &fHelixBlock    );
  RegisterDataBlock("SpmcBlockVDet"          , "TStepPointMCBlock", &fSpmcBlockVDet );
//-----------------------------------------------------------------------------
// for validation purposes
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// initialize likelihood histograms
//-----------------------------------------------------------------------------
  if (fWriteTmvaTree >= 0) {
    TDirectory* dir = gDirectory;

    const char* dsname = GetAna()->GetInputModule()->GetDataset(0)->GetName();

    int algo = fWriteTmvaTree;

    if      (algo == 0) {
      fTmvaFile  = new TFile(Form("%s.tmva_training_trkpatrec.root",dsname),"recreate");
    }
    else if (algo == 1) {
      fTmvaFile  = new TFile(Form("%s.tmva_training_calpatrec.root",dsname),"recreate");
    }

    fSigTree  = new TTree("tmva_training_tree","TMVA Signal Tree");

    //    fBgrTree  = new TTree("background","TMVA Background Tree");

					// signal tree

    fSigBranch.fP          = fSigTree->Branch("p"       ,&fTmvaData.fP         ,"F");
    fSigBranch.fPMC        = fSigTree->Branch("pmc"     ,&fTmvaData.fPMC       ,"F");
    fSigBranch.fTanDip     = fSigTree->Branch("tdip"    ,&fTmvaData.fTanDip    ,"F");
    fSigBranch.fNActive    = fSigTree->Branch("nactive" ,&fTmvaData.fNActive   ,"F");
    fSigBranch.fNaFract    = fSigTree->Branch("nafract" ,&fTmvaData.fNaFract   ,"F");
    fSigBranch.fChi2Dof    = fSigTree->Branch("chi2d"   ,&fTmvaData.fChi2Dof   ,"F");
    fSigBranch.fFitCons    = fSigTree->Branch("fcons"   ,&fTmvaData.fFitCons   ,"F");
    fSigBranch.fMomErr     = fSigTree->Branch("momerr"  ,&fTmvaData.fMomErr    ,"F");
    fSigBranch.fT0Err      = fSigTree->Branch("t0err"   ,&fTmvaData.fT0Err     ,"F");
    fSigBranch.fD0         = fSigTree->Branch("d0"      ,&fTmvaData.fD0        ,"F");
    fSigBranch.fRMax       = fSigTree->Branch("rmax"    ,&fTmvaData.fRMax      ,"F");
    fSigBranch.fNdaOverNa  = fSigTree->Branch("nda_o_na",&fTmvaData.fNdaOverNa ,"F");
    fSigBranch.fNzaOverNa  = fSigTree->Branch("nza_o_na",&fTmvaData.fNzaOverNa ,"F");
    fSigBranch.fNmaOverNa  = fSigTree->Branch("nma_o_na",&fTmvaData.fNmaOverNa ,"F");
    fSigBranch.fZ1         = fSigTree->Branch("z1"      ,&fTmvaData.fZ1        ,"F");
    fSigBranch.fWeight     = fSigTree->Branch("wt"      ,&fTmvaData.fWeight    ,"F");

    dir->cd();
  }
//-----------------------------------------------------------------------------
// init two MVA-based classifiers - KPAR (TPR) and KDAR (CPR) 
//-----------------------------------------------------------------------------
  if (fUseMVA) {
    fhicl::ParameterSet pset_tpr, pset_cpr;

    fNMVA = 1;

    string              s1(fTprMVA->XmlWeightsFile());
    string              s2(fCprMVA->XmlWeightsFile());

    printf(">>> [TTrackCompModule::BeginJob] Init TrkPatRec MVA from %s\n",s1.data());
    pset_tpr.put<string>("MVAWeights",s1);
    fTprQualMva = new mu2e::MVATools(pset_tpr);
    fTprQualMva->initMVA();
    //    fTrkQualMva->showMVA();

    printf(">>> [TTrackCompModule::BeginJob] Init CalPatRec MVA from %s\n",s2.data());
    pset_cpr.put<string>("MVAWeights",s2);
    fCprQualMva = new mu2e::MVATools(pset_cpr);
    fCprQualMva->initMVA();
  }

  // string hist_dir      = gEnv->GetValue("mu2e.TrkQual.HistDir","_none_");
  // string trk_qual_dsid = gEnv->GetValue("mu2e.TrkQual.Dsid"   ,"_none_");

  // fTrkQualFile    = Form("%s/%s.track_comp_use_mva_%03i.hist",
  // 			 hist_dir.data(),trk_qual_dsid.data(),fUseMVA);

  // fTrackProb[0]   = new prob_dist(fTrkQualFile.Data(),"TrackComp","trk_100/mvaout");
  // fTrackProb[1]   = new prob_dist(fTrkQualFile.Data(),"TrackComp","trk_200/mvaout");

  fBestID[0]      = fTprMVA->BestID();		  // Dave's default: DaveTrkQual > 0.4
  fBestID[1]      = fCprMVA->BestID();		  // KDAR     : CprQual     > 0.85

  return 0;
}


//_____________________________________________________________________________
int TTrackCompModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//-----------------------------------------------------------------------------
void TTrackCompModule::BookDTrackHistograms(HistBase_t* HistR, const char* Folder) {
  DTrackHist_t* Hist =  (DTrackHist_t*) HistR;

  HBook1F(Hist->fDp       ,"dp"      ,Form("%s: P(1)-P(2)"               ,Folder),200, -1,1,Folder);
  HBook1F(Hist->fRMomErr10,"rmomerr" ,Form("%s: MomEff(1)/MomErr(0)"     ,Folder),200,  0,2,Folder);
  HBook1F(Hist->fDT0      ,"dt0"     ,Form("%s: T0(1)-T0(2)"             ,Folder),200, -5,5,Folder);


}

//_____________________________________________________________________________
void TTrackCompModule::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();

  TH1::SetDefaultSumw2(kTRUE);	

  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;		// events with EclMax > fMinETrig and TClMax > 550
  book_event_histset[ 2] = 1;           // *** fill in
  book_event_histset[ 3] = 1;           // events with at least one helix
  book_event_histset[ 4] = 1;           // events with at least one reconstructed track seed
  book_event_histset[ 5] = 1;           // events with at least one reconstructed track seed nhits >= 15 chi2(ZPHI)<4

					// 10:19: KPAR efficiency, 20:29: KDAR efficiency

  for (int i=10; i<20; i++) book_event_histset[i] = 1;
  for (int i=20; i<30; i++) book_event_histset[i] = 1;

  book_event_histset[41] = 1;           // 41: events with at least one positron P > 85 MeV/c

  book_event_histset[51] = 1;           // 51: events with at least one reconstructed PAR track
  book_event_histset[52] = 1;           // 52: events with at least one reconstructed PAR track passing all cuts 

  book_event_histset[61] = 1;           // 61: events with at least one reconstructed DAR track
  book_event_histset[62] = 1;           // 62: events with at least one reconstructed DAR track passing all cuts 

  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
    }
  }
//--------------------------------------------------------------------------------
// book trackSeed histograms
//--------------------------------------------------------------------------------
  int book_trackSeed_histset[kNTrackSeedHistSets];
  for (int i=0; i<kNTrackSeedHistSets; ++i)  book_trackSeed_histset[i] = 0;

  book_trackSeed_histset[0] = 1;   // events with at least one trackSeed
  book_trackSeed_histset[1] = 1;   // events with at least one trackSeed with p > 80 MeV/c
  book_trackSeed_histset[2] = 1;   // events with at least one trackSeed with p > 90 MeV/c
  book_trackSeed_histset[3] = 1;   // events with at least one trackSeed with p > 100 MeV/c
  book_trackSeed_histset[4] = 1;   // events with at least one trackSeed with 10 < nhits < 15
  book_trackSeed_histset[5] = 1;   // events with at least one trackSeed with nhits >= 15
  book_trackSeed_histset[6] = 1;   // events with at least one trackSeed with nhits >= 15 and chi2(ZPhi)<4

  for (int i=0; i<kNTrackSeedHistSets; i++) {
    if (book_trackSeed_histset[i] != 0) {
      sprintf(folder_name,"trkseed_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrackSeed[i] = new TrackSeedHist_t;
      BookTrackSeedHistograms(fHist.fTrackSeed[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track histograms
//-----------------------------------------------------------------------------
  TString*  track_selection   [kNTrackHistSets];

  for (int i=0; i<kNTrackHistSets; i++) { track_selection[i] = NULL; }

  track_selection[  0] = new TString("good PAR tracks");
  track_selection[  1] = new TString("good DAR tracks in events with no good kPAR tracks TrkQual>0.3");
  track_selection[  2] = new TString("good DAR tracks in events with no good kPAR tracks TrkQual>0.2");
  //  track_selection[  3] = new TString("best track");

  track_selection[100] = new TString("PAR all tracks");
  track_selection[101] = new TString("PAR BestTrackID");
  track_selection[102] = new TString("PAR BestTrackID no fitCons&momErr&t0Err tracks");
  track_selection[103] = new TString("PAR BestTrackID, dpf>1 tracks");
  track_selection[104] = new TString("PAR all  tracks events Ecl > 60");
  track_selection[105] = new TString("PAR BestTrackID tracks events Ecl > 60");
  track_selection[106] = new TString("PAR tracks with dtz0 < -10");
  track_selection[109] = new TString("PAR tracks with XDpF > +10 MeV");

  track_selection[110] = new TString("PAR TrackID[0] - SetC");
  track_selection[111] = new TString("PAR TrackID[1]");
  track_selection[112] = new TString("PAR TrackID[2]");
  track_selection[113] = new TString("PAR TrackID[3]");
  track_selection[114] = new TString("PAR TrackID[4]");
  track_selection[115] = new TString("PAR TrackID[5]");
  track_selection[116] = new TString("PAR TrackID[6]");
  track_selection[117] = new TString("PAR TrackID[7]");
  track_selection[118] = new TString("PAR TrackID[8]");
  track_selection[119] = new TString("PAR TrackID[9]");
  track_selection[120] = new TString("PAR TrackID[10]");
  track_selection[121] = new TString("PAR TrackID[11]");
  track_selection[122] = new TString("PAR TrackID[12]");
  track_selection[123] = new TString("PAR TrackID[13]");
  track_selection[124] = new TString("PAR TrackID[14]");
  track_selection[125] = new TString("PAR TrackID[15]");
  track_selection[126] = new TString("PAR TrackID[16]");
  track_selection[127] = new TString("PAR TrackID[17]");
  track_selection[128] = new TString("PAR TrackID[18]");
  track_selection[129] = new TString("PAR TrackID[19]");

  track_selection[141] = new TString("PAR tracks N(active) > 20");
  track_selection[142] = new TString("PAR tracks N(active) > 20, |D0| < 100");
  track_selection[143] = new TString("PAR tracks N(active) > 20, |D0| < 100, DN(active) < 6");
  track_selection[144] = new TString("PAR tracks N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");

  track_selection[150] = new TString("PAR tracks |DpF| < 1");
  track_selection[151] = new TString("PAR tracks |DpF| < 1 ; N(active) > 20");
  track_selection[152] = new TString("PAR tracks |DpF| < 1 ; N(active) > 20, |D0| < 100");
  track_selection[153] = new TString("PAR tracks |DpF| < 1 ; N(active) > 20, |D0| < 100, DN(active) < 6");
  track_selection[154] = new TString("PAR tracks |DpF| < 1 ; N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");

  track_selection[160] = new TString("PAR tracks DpF   > 2");
  track_selection[161] = new TString("PAR tracks DpF   > 2 ; N(active) > 20");
  track_selection[162] = new TString("PAR tracks DpF   > 2 ; N(active) > 20, |D0| < 100");
  track_selection[163] = new TString("PAR tracks DpF   > 2 ; N(active) > 20, |D0| < 100, DN(active) < 6");
  track_selection[164] = new TString("PAR tracks DpF   > 2 ; N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");

  track_selection[171] = new TString("PAR+ tracks with final selections");
  track_selection[172] = new TString("PAR- tracks with final selections");
  track_selection[173] = new TString("PAR+ tracks with final selections and process-defined weight");
  track_selection[174] = new TString("PAR- tracks with final selections and process-defined weight");

  track_selection[177] = new TString("cosmics+: #171 + 90 < p < 93");
  track_selection[178] = new TString("cosmics-: #172 + 90 < p < 93");

  track_selection[179] = new TString("PAR- tracks with final selections and DIO weight");

  track_selection[180] = new TString("PAR+ tracks all");                                       
  track_selection[181] = new TString("PAR+ tracks N(active) > 20");                                       
  track_selection[182] = new TString("PAR+ tracks N(active) > 20 and |D0| < 100");			       
  track_selection[183] = new TString("PAR+ tracks N(active) > 20, |D0| < 100, DN(active) < 6");	       
  track_selection[184] = new TString("PAR+ tracks N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");
  track_selection[185] = new TString("PAR+ tracks with final selections, E/P < 0.7 (mu+/pi+)");
  track_selection[186] = new TString("PAR+ tracks with final selections, E/P > 0.7 (e+)");

  track_selection[190] = new TString("PAR- tracks all");                                       
  track_selection[191] = new TString("PAR- tracks N(active) > 20");                                       
  track_selection[192] = new TString("PAR- tracks N(active) > 20 and |D0| < 100");
  track_selection[193] = new TString("PAR- tracks N(active) > 20, |D0| < 100, DN(active) < 6");
  track_selection[194] = new TString("PAR- tracks N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");
  track_selection[195] = new TString("PAR- tracks with final selections, E/P < 0.7 (mu-/pi-)");
  track_selection[196] = new TString("PAR- tracks with final selections, E/P > 0.7 (e-)");

  track_selection[200] = new TString("DAR all tracks");
  track_selection[201] = new TString("DAR BestTrackID");
  track_selection[202] = new TString("DAR BestTrackID no fitCons&momErr&t0Err tracks");
  track_selection[203] = new TString("DAR BestTrackID, dpf>1 tracks");
  track_selection[204] = new TString("DAR all  tracks events Ecl > 60");
  track_selection[205] = new TString("DAR BestTrackID tracks events Ecl > 60");
  track_selection[206] = new TString("DAR tracks with dtz0 < -10");
  track_selection[209] = new TString("DAR tracks with XDpF > +10 MeV");

  track_selection[210] = new TString("DAR TrackID[0] - SetC");
  track_selection[211] = new TString("DAR TrackID[1]");
  track_selection[212] = new TString("DAR TrackID[2]");
  track_selection[213] = new TString("DAR TrackID[3]");
  track_selection[214] = new TString("DAR TrackID[4]");
  track_selection[215] = new TString("DAR TrackID[5]");
  track_selection[216] = new TString("DAR TrackID[6]");
  track_selection[217] = new TString("DAR TrackID[7]");
  track_selection[218] = new TString("DAR TrackID[8]");
  track_selection[219] = new TString("DAR TrackID[9]");
  track_selection[220] = new TString("DAR TrackID[10]");
  track_selection[221] = new TString("DAR TrackID[11]");
  track_selection[222] = new TString("DAR TrackID[12]");
  track_selection[223] = new TString("DAR TrackID[13]");
  track_selection[224] = new TString("DAR TrackID[14]");
  track_selection[225] = new TString("DAR TrackID[15]");
  track_selection[226] = new TString("DAR TrackID[16]");
  track_selection[227] = new TString("DAR TrackID[17]");
  track_selection[228] = new TString("DAR TrackID[18]");
  track_selection[229] = new TString("DAR TrackID[19");

  track_selection[241] = new TString("DAR tracks N(active) > 20");
  track_selection[242] = new TString("DAR tracks N(active) > 20, |D0| < 100");
  track_selection[243] = new TString("DAR tracks N(active) > 20, |D0| < 100, DN(active) < 6");
  track_selection[244] = new TString("DAR tracks N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");

  track_selection[250] = new TString("DAR tracks |DpF| < 1");
  track_selection[251] = new TString("DAR tracks |DpF| < 1 ; N(active) > 20");
  track_selection[252] = new TString("DAR tracks |DpF| < 1 ; N(active) > 20, |D0| < 100");
  track_selection[253] = new TString("DAR tracks |DpF| < 1 ; N(active) > 20, |D0| < 100, DN(active) < 6");
  track_selection[254] = new TString("DAR tracks |DpF| < 1 ; N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");

  track_selection[260] = new TString("DAR tracks DpF   > 2");
  track_selection[261] = new TString("DAR tracks DpF   > 2 ; N(active) > 20");
  track_selection[262] = new TString("DAR tracks DpF   > 2 ; N(active) > 20, |D0| < 100");
  track_selection[263] = new TString("DAR tracks DpF   > 2 ; N(active) > 20, |D0| < 100, DN(active) < 6");
  track_selection[264] = new TString("DAR tracks DpF   > 2 ; N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");

  track_selection[271] = new TString("DAR+ tracks with final selections");
  track_selection[272] = new TString("DAR- tracks with final selections");
  track_selection[273] = new TString("DAR+ tracks with final selections and process-defined weight");
  track_selection[274] = new TString("DAR- tracks with final selections and process-defined weight");

  track_selection[277] = new TString("cosmics+: #271 + 90 < p < 93");
  track_selection[278] = new TString("cosmics-: #272 + 90 < p < 93");

  track_selection[279] = new TString("DAR- tracks with final selections and DIO weight");

  track_selection[280] = new TString("DAR+ tracks all");                                       
  track_selection[281] = new TString("DAR+ tracks N(active) > 20");                                       
  track_selection[282] = new TString("DAR+ tracks N(active) > 20 and |D0| < 100");			       
  track_selection[283] = new TString("DAR+ tracks N(active) > 20, |D0| < 100, DN(active) < 6");	       
  track_selection[284] = new TString("DAR+ tracks N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");
  track_selection[285] = new TString("DAR+ tracks with final selections, E/P < 0.7 (mu+/pi+)");
  track_selection[286] = new TString("DAR+ tracks with final selections, E/P > 0.7 (e+)");

  track_selection[290] = new TString("DAR- tracks all");                                       
  track_selection[291] = new TString("DAR- tracks N(active) > 20");                                       
  track_selection[292] = new TString("DAR- tracks N(active) > 20 and |D0| < 100");
  track_selection[293] = new TString("DAR- tracks N(active) > 20, |D0| < 100, DN(active) < 6");
  track_selection[294] = new TString("DAR- tracks N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");
  track_selection[295] = new TString("DAR- tracks with final selections, E/P < 0.7 (mu-/pi-)");
  track_selection[296] = new TString("DAR- tracks with final selections, E/P > 0.7 (e-)");

  track_selection[301] = new TString("PAR- tracks final selections, dr/dz(cal) <  0");
  track_selection[302] = new TString("PAR- tracks final selections, dr/dz(cal) >= 0");
  track_selection[305] = new TString("PAR- tracks final selections, cluster in the 1st disk");
  track_selection[306] = new TString("PAR- tracks final selections, cluster in the 2nd disk");
  track_selection[307] = new TString("PAR- tracks final selections, dt <  -4 ns");
  track_selection[308] = new TString("PAR- tracks final selections, dt >= -4 ns");

  track_selection[401] = new TString("DAR- tracks final selections, dr/dz(cal) <  0");
  track_selection[402] = new TString("DAR- tracks final selections, dr/dz(cal) >= 0");
  track_selection[405] = new TString("DAR- tracks final selections, cluster in the 1st disk");
  track_selection[406] = new TString("DAR- tracks final selections, cluster in the 2nd disk");
  track_selection[407] = new TString("DAR- tracks final selections, dt <  -4 ns");
  track_selection[408] = new TString("DAR- tracks final selections, dt >= -4 ns");

  const char* folder_title;
  for (int i=0; i<kNTrackHistSets; i++) {
    if (track_selection[i] != 0) {
      sprintf(folder_name,"trk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title = folder_name;
      folder_title = track_selection[i]->Data();
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fTrack[i] = new TrackHist_t;
      BookTrackHistograms(fHist.fTrack[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book delta track histograms
//-----------------------------------------------------------------------------
  TString*  dtrack_selection   [kNDTrackHistSets];

  for (int i=0; i<kNDTrackHistSets; i++) { dtrack_selection[i] = NULL; }

  dtrack_selection[  0] = new TString("all PAR-DAR tracks");

  //  const char* folder_title;
  for (int i=0; i<kNDTrackHistSets; i++) {
    if (dtrack_selection[i] != NULL) {
      sprintf(folder_name,"dtrk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title = folder_name;
      folder_title = dtrack_selection[i]->Data();
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fDTrack[i] = new DTrackHist_t;
      BookDTrackHistograms(fHist.fDTrack[i],Form("Hist/%s",folder_name));
    }
  }

}


//-----------------------------------------------------------------------------
// TmvaAlgorithm: 100*usez + algorigthm
//-----------------------------------------------------------------------------
int TTrackCompModule::FillTmvaTree() {
  int rc(0), loc(-1);

  loc = fWriteTmvaTree;			// 0:kPAR , 1:kDAR

  if (fTrackBlock[loc]->NTracks() != 1) return rc;

					// first KDAR track

  TStnTrack* trk = fTrackBlock[loc]->Track(0);
  TrackPar_t* tp = &fTrackPar[loc][0];

  float na       = trk->NActive();

  fTmvaData.fP          = tp->fP;
  fTmvaData.fPMC        = trk->fPFront;
  fTmvaData.fTanDip     = trk->TanDip();
  fTmvaData.fNActive    = na;
  fTmvaData.fNaFract    = na/trk->NHits();
  fTmvaData.fChi2Dof    = trk->Chi2Dof();
  fTmvaData.fFitCons    = trk->FitCons();
  fTmvaData.fMomErr     = trk->FitMomErr();
  fTmvaData.fT0Err      = trk->T0Err();
  fTmvaData.fD0         = trk->D0();
  fTmvaData.fRMax       = trk->RMax();
  fTmvaData.fNdaOverNa  = trk->NDoubletsAct()/na;
  fTmvaData.fNzaOverNa  = trk->NHitsAmbZero()/na;
  fTmvaData.fNmaOverNa  = trk->NMatActive()/na;
  fTmvaData.fZ1         = trk->fZ1;
  fTmvaData.fWeight     = 1.;

// try one! // there should be two trees - signal and background
// 
//   if (tp->fDpF > 0.7) {
//     fBgrTree->Fill();
//   }
//   else if (fabs(tp->fDpF) < 0.25) {
//     fSigTree->Fill();
//   }

  fSigTree->Fill();

  return rc;
}

//-----------------------------------------------------------------------------
// for DIO : ultimately, one would need to renormalize the distribution
//-----------------------------------------------------------------------------
void TTrackCompModule::FillDTrackHistograms(HistBase_t* HistR, TStnTrack* Trk1, TrackPar_t* Tp1, TStnTrack* Trk2, TrackPar_t* Tp2) {

  //  TLorentzVector  mom;

  DTrackHist_t* Hist = (DTrackHist_t*) HistR;

					// Tp->fP - corrected momentum, fP0 and fP2 - not corrected
  float dp = Tp1->fP-Tp2->fP;
  Hist->fDp->Fill (dp);

  float rmomerr10 = Trk2->FitMomErr()/Trk1->FitMomErr();
  Hist->fRMomErr10->Fill (rmomerr10);

  float dt0 = Trk1->T0()-Trk2->T0();
  Hist->fDT0->Fill (dt0);
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
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);
//-----------------------------------------------------------------------------
// EVT_1: events with 50 MeV+ cluster and T0 > 550
//-----------------------------------------------------------------------------
  if ((fEClMax > 50.) && (fTClMax > 550)) {
    FillEventHistograms(fHist.fEvent[1],&fEvtPar);
  }
//-----------------------------------------------------------------------------
// EVT_2: events with 50 MeV+ cluster and both tracks T0 > 550
//-----------------------------------------------------------------------------
  if ((fEClMax > fMinETrig) && (fTClMax > 550)) {
    FillEventHistograms(fHist.fEvent[2],&fEvtPar);
  }
//-----------------------------------------------------------------------------
// EVT_3: events with at least one helix
//-----------------------------------------------------------------------------
  if (fNHelices > 0) FillEventHistograms(fHist.fEvent[3],&fEvtPar);
//-----------------------------------------------------------------------------
// EVT_4: events with at least one track seed
// EVT_5: events with at least one GOOD track seed
//-----------------------------------------------------------------------------
  if (fNTrackSeeds     > 0) FillEventHistograms(fHist.fEvent[4],&fEvtPar);
  if (fNGoodTrackSeeds > 0) FillEventHistograms(fHist.fEvent[5],&fEvtPar);
//-----------------------------------------------------------------------------
// EVT_41: events with good DAR positron track 
//-----------------------------------------------------------------------------
  if ((fNTracks[kDAR] > 0) && (fWeight > 0)) {
    TStnTrack*  trk = fTrackBlock[kDAR]->Track(0);
    TrackPar_t* tp  = &fTrackPar[kDAR][0];
    if ((trk->P0() > 85) && (tp->fIDWord_RMC == 0) && (trk->Charge() > 0)) {
      FillEventHistograms(fHist.fEvent[41],&fEvtPar);
    }
  }
//-----------------------------------------------------------------------------
// EVT_51 : events with PAR tracks 
// EVT_52 : events with good PAR tracks 
//-----------------------------------------------------------------------------
  // TLorentzVector    mom (1.,0.,0.,0);
  // if (fParticle) fParticle->Momentum(mom);
  // double mc_mom = mom.P();

  if (fNTracks[kPAR] > 0) {
    FillEventHistograms(fHist.fEvent[51],&fEvtPar);
    TrackPar_t* tp  = &fTrackPar[kPAR][0];
    if (tp->fIDWord_RMC == 0) FillEventHistograms(fHist.fEvent[52],&fEvtPar);
  }
//-----------------------------------------------------------------------------
// EVT_61 : events with DAR tracks 
// EVT_62 : events with good DAR tracks 
//-----------------------------------------------------------------------------
  if (fNTracks[kDAR] > 0) {
    FillEventHistograms(fHist.fEvent[61],&fEvtPar);
    TrackPar_t* tp  = &fTrackPar[kDAR][0];
    if (tp->fIDWord_RMC == 0) FillEventHistograms(fHist.fEvent[62],&fEvtPar);
  }

//--------------------------------------------------------------------------------
// track seed histograms
//--------------------------------------------------------------------------------
  TStnTrackSeed* trkSeed;
  for (int i=0; i<fNTrackSeeds; ++i) {
    trkSeed = fTrackSeedBlock->TrackSeed(i);
    
    FillTrackSeedHistograms(fHist.fTrackSeed[0], trkSeed);
    
    int     nhits    = trkSeed->NHits();
    double  p        = trkSeed->P();
    double  chi2     = trkSeed->Chi2();
   
    if (p >  80.) FillTrackSeedHistograms(fHist.fTrackSeed[1], trkSeed);
    if (p >  90.) FillTrackSeedHistograms(fHist.fTrackSeed[2], trkSeed);
    if (p > 100.) FillTrackSeedHistograms(fHist.fTrackSeed[3], trkSeed);

    if ( (nhits>10) && (nhits <  15)) FillTrackSeedHistograms(fHist.fTrackSeed[4], trkSeed);
    if ( nhits>=15                  ) FillTrackSeedHistograms(fHist.fTrackSeed[5], trkSeed);
    if ( (chi2 < 4) && (nhits >= 15)) FillTrackSeedHistograms(fHist.fTrackSeed[6], trkSeed);
  }

//-----------------------------------------------------------------------------
// what does kDAR add ?
// TRK_0 : kPAR tracks, BEST_ID
//-----------------------------------------------------------------------------
  if ((fTrackBlock[kPAR]->NTracks() > 0) && (fTrackPar[kPAR][0].fIDWord[fBestID[0]] == 0)) {
    TStnTrack* trk = fTrackBlock[kPAR]->Track(0);
    TrackPar_t* tp = &fTrackPar[kPAR][0];
    FillTrackHistograms(fHist.fTrack[0],trk,tp,&fSimPar);
  }
  else {
    if (fTrackBlock[kDAR]->NTracks() > 0) {
      if (fTrackPar[kDAR][0].fIDWord[3] == 0) {
//-----------------------------------------------------------------------------
// TRK_1: have KDAR track TrkQual > 0.3 and no KPAR TrkQual > 0.4 track
//-----------------------------------------------------------------------------
	TStnTrack* trk = fTrackBlock[kDAR]->Track(0);
	TrackPar_t* tp = &fTrackPar[kDAR][0];
	FillTrackHistograms(fHist.fTrack[1],trk,tp,&fSimPar);
      }
      if (fTrackPar[1][0].fIDWord[2] == 0) {
//-----------------------------------------------------------------------------
// TRK_2: have KDAR track TrkQual > 0.2 and no KPAR TrkQual > 0.4 track
//-----------------------------------------------------------------------------
	TStnTrack* trk = fTrackBlock[1]->Track(0);
	TrackPar_t* tp = &fTrackPar[1][0];
	FillTrackHistograms(fHist.fTrack[2],trk,tp,&fSimPar);
      }
    }
  }
//-----------------------------------------------------------------------------
// KPAR and KDAR histograms, inclusive, ihist defines the offset
// i=0:KPAR, i=1:KDAR
//-----------------------------------------------------------------------------
  for (int i=0; i<2; i++) {
    n_setc_tracks[i] = 0;
    ihist            = 100*(i+1);
    for (int itrk=0; itrk<fNTracks[i]; itrk++) {
      trk = fTrackBlock[i]->Track(itrk);
      tp  = fTrackPar[i]+itrk;
//-----------------------------------------------------------------------------
// TRK_100, TRK_200: all reconstructed tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[ihist+0],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// TRK_101, TRK_201: BestID 
//-----------------------------------------------------------------------------
      if (tp->fIDWord[fBestID[i]] == 0) {
//-----------------------------------------------------------------------------
// cosmics
//-----------------------------------------------------------------------------
	FillTrackHistograms(fHist.fTrack[ihist+1],trk,tp,&fSimPar);
	n_setc_tracks[i] += 1;
      }
//-----------------------------------------------------------------------------
// TRK_102, TRK_202: (BestID - FitConsBit - T0ErrBit - MomErrBit) tracks 
//-----------------------------------------------------------------------------
      int mask = (TStnTrackID::kFitConsBit | TStnTrackID::kT0ErrBit | TStnTrackID::kMomErrBit);
      if ((tp->fIDWord[fBestID[i]] & ~mask) == 0) {
	FillTrackHistograms(fHist.fTrack[ihist+2],trk,tp,&fSimPar);
      }
//-----------------------------------------------------------------------------
// IHIST+3: (SetC + (dpf > 1)tracks 
//-----------------------------------------------------------------------------
      if ((tp->fIDWord[fBestID[i]] == 0) & (tp->fDpF > 1)) {
	FillTrackHistograms(fHist.fTrack[ihist+3],trk,tp,&fSimPar);
      }
//-----------------------------------------------------------------------------
// IHIST+4: add   Ecl > 60 requirement
//-----------------------------------------------------------------------------
      if (fEClMax > 60.) {
	FillTrackHistograms(fHist.fTrack[ihist+4],trk,tp,&fSimPar);

	if (tp->fIDWord[0] == 0) {
	  FillTrackHistograms(fHist.fTrack[ihist+5],trk,tp,&fSimPar);
	}
      }
//-----------------------------------------------------------------------------
// IHIST+6: oddly, tracks with seemingly wrong reconstruced T0
//-----------------------------------------------------------------------------
      if (tp->fDtZ0 < -10) {
	FillTrackHistograms(fHist.fTrack[ihist+6],trk,tp,&fSimPar);
      }
//-----------------------------------------------------------------------------
// IHIST+9: tracks with XDpF > 10
//-----------------------------------------------------------------------------
      if (tp->fXDpF   > 10) FillTrackHistograms(fHist.fTrack[ihist+9],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// different cuts on track quality variable
//-----------------------------------------------------------------------------
      for (int idd=0; idd<fNID; idd++) {
	if (tp->fIDWord[idd] == 0) FillTrackHistograms(fHist.fTrack[ihist+10+idd],trk,tp,&fSimPar);
      }
//-----------------------------------------------------------------------------
// IHIST+41: N(a) > 20
// IHIST+42: N(a) > 20, |D0| < 100
// IHIST+43: N(a) > 20, |D0| < 100, DNa < 6
// IHIST+44: N(a) > 20, |D0| < 100, DNa < 6, chi2d < 4
//-----------------------------------------------------------------------------
      if (trk->NActive() > 20) { 
	FillTrackHistograms(fHist.fTrack[ihist+41],trk,tp,&fSimPar);
	if (fabs(trk->D0()) < 100) {
	  FillTrackHistograms(fHist.fTrack[ihist+42],trk,tp,&fSimPar);
	  int dna = trk->NHits()-trk->NActive();
	  if (dna < 6) {
	    FillTrackHistograms(fHist.fTrack[ihist+43],trk,tp,&fSimPar);
	    if (trk->Chi2Dof() < 4.) {
	      FillTrackHistograms(fHist.fTrack[ihist+44],trk,tp,&fSimPar);
	    }
	  }
	}
      }

      if (trk->Charge() > 0) {
//-----------------------------------------------------------------------------
// positive tracks only:
// ---------------------
// IHIST+80: all
// IHIST+81: N(a) > 20
// IHIST+82: N(a) > 20, |D0| < 100
// IHIST+83: N(a) > 20, |D0| < 100, DNa < 6
// IHIST+84: N(a) > 20, |D0| < 100, DNa < 6, chi2d < 4
//-----------------------------------------------------------------------------
	FillTrackHistograms(fHist.fTrack[ihist+80],trk,tp,&fSimPar);
	if (trk->NActive() > 20) { 
	  FillTrackHistograms(fHist.fTrack[ihist+81],trk,tp,&fSimPar);
	  if (fabs(trk->D0()) < 100) {
	    FillTrackHistograms(fHist.fTrack[ihist+82],trk,tp,&fSimPar);
	    int dna = trk->NHits()-trk->NActive();
	    if (dna < 6) {
	      FillTrackHistograms(fHist.fTrack[ihist+83],trk,tp,&fSimPar);
	      if (trk->Chi2Dof() < 4.) {
		FillTrackHistograms(fHist.fTrack[ihist+84],trk,tp,&fSimPar);
	      }
	    }
	  }
	}
      }
      else {
//-----------------------------------------------------------------------------
// negative tracks only:
// ---------------------
// IHIST+90: all
// IHIST+91: N(a) > 20
// IHIST+92: N(a) > 20, |D0| < 100
// IHIST+93: N(a) > 20, |D0| < 100, DNa < 6
// IHIST+94: N(a) > 20, |D0| < 100, DNa < 6, chi2d < 4
//-----------------------------------------------------------------------------
	FillTrackHistograms(fHist.fTrack[ihist+90],trk,tp,&fSimPar);
	if (trk->NActive() > 20) { 
	  FillTrackHistograms(fHist.fTrack[ihist+91],trk,tp,&fSimPar);
	  if (fabs(trk->D0()) < 100) {
	    FillTrackHistograms(fHist.fTrack[ihist+92],trk,tp,&fSimPar);
	    int dna = trk->NHits()-trk->NActive();
	    if (dna < 6) {
	      FillTrackHistograms(fHist.fTrack[ihist+93],trk,tp,&fSimPar);
	      if (trk->Chi2Dof() < 4.) {
		FillTrackHistograms(fHist.fTrack[ihist+94],trk,tp,&fSimPar);
	      }
	    }
	  }
	}
      }
//-----------------------------------------------------------------------------
// |DPF| < 1 MeV/c
// IHIST+51: N(a) > 20
// IHIST+52: N(a) > 20, |D0| < 100
// IHIST+53: N(a) > 20, |D0| < 100, DNa < 6
// IHIST+54: N(a) > 20, |D0| < 100, DNa < 6, chi2d < 4
//-----------------------------------------------------------------------------
      if (fabs(tp->fDpF) < 1.) {
	FillTrackHistograms(fHist.fTrack[ihist+50],trk,tp,&fSimPar);
	if (trk->NActive() > 20) { 
	  FillTrackHistograms(fHist.fTrack[ihist+51],trk,tp,&fSimPar);
	  if (fabs(trk->D0()) < 100) {
	    FillTrackHistograms(fHist.fTrack[ihist+52],trk,tp,&fSimPar);
	    int dna = trk->NHits()-trk->NActive();
	    if (dna < 6) {
	      FillTrackHistograms(fHist.fTrack[ihist+53],trk,tp,&fSimPar);
	      if (trk->Chi2Dof() < 4.) {
		FillTrackHistograms(fHist.fTrack[ihist+54],trk,tp,&fSimPar);
	      }
	    }
	  }
	}
      }
//-----------------------------------------------------------------------------
// |DPF| > 2 MeV/c
// IHIST+61: N(a) > 20
// IHIST+62: N(a) > 20 |D0| < 150
// IHIST+63: N(a) > 20, |D0| < 150, DNa < 6
// IHIST+64: N(a) > 20, |D0| < 150, DNa < 6, chi2d < 4
//-----------------------------------------------------------------------------
      if (tp->fDpF > 2.) {
	FillTrackHistograms(fHist.fTrack[ihist+60],trk,tp,&fSimPar);
	if (trk->NActive() > 20) { 
	  FillTrackHistograms(fHist.fTrack[ihist+61],trk,tp,&fSimPar);
	  if (fabs(trk->D0()) < 100) {
	    FillTrackHistograms(fHist.fTrack[ihist+62],trk,tp,&fSimPar);
	    int dna = trk->NHits()-trk->NActive();
	    if (dna < 6) {
	      FillTrackHistograms(fHist.fTrack[ihist+63],trk,tp,&fSimPar);
	      if (trk->Chi2Dof() < 4.) {
		FillTrackHistograms(fHist.fTrack[ihist+64],trk,tp,&fSimPar);
	      }
	    }
	  }
	}
      }
//-----------------------------------------------------------------------------
// final selections
//-----------------------------------------------------------------------------
      if (tp->fIDWord_RMC == 0) {
	if (trk->Charge() == 1) {
//-----------------------------------------------------------------------------
// positive tracks
//-----------------------------------------------------------------------------
	  FillTrackHistograms  (fHist.fTrack[ihist+71],trk,tp,&fSimPar);
	  if (fProcess != -1) FillTrackHistograms(fHist.fTrack[ihist+73],trk,tp,&fSimPar,fWeight); 
	  if ((tp->fP > 90.) && (tp->fP < 93.)) FillTrackHistograms(fHist.fTrack[ihist+77],trk,tp,&fSimPar);            // for cosmics

	  if   (tp->fEp < 0.7) FillTrackHistograms(fHist.fTrack[ihist+85],trk,tp,&fSimPar);   // "e+"
	  else                 FillTrackHistograms(fHist.fTrack[ihist+86],trk,tp,&fSimPar);   // "mu+/pi+"
	}
	else {
//-----------------------------------------------------------------------------
// negative tracks
//-----------------------------------------------------------------------------
	  FillTrackHistograms(fHist.fTrack[ihist+72],trk,tp,&fSimPar);
	  if (fProcess > 0) FillTrackHistograms(fHist.fTrack[ihist+74],trk,tp,&fSimPar,fWeight);                       // weighting
	  
	  if ((tp->fP > 90.) && (tp->fP < 93.)) FillTrackHistograms(fHist.fTrack[ihist+78],trk,tp,&fSimPar);            // for cosmics

	  if   (tp->fEp < 0.7) FillTrackHistograms(fHist.fTrack[ihist+95],trk,tp,&fSimPar);   // "e-"
	  else                 FillTrackHistograms(fHist.fTrack[ihist+96],trk,tp,&fSimPar);   // "mu-/pi-"

	  FillTrackHistograms(fHist.fTrack[ihist+79],trk,tp,&fSimPar,tp->fDioLLWt);                                       // DIO

//-----------------------------------------------------------------------------
// debugging mu-
//-----------------------------------------------------------------------------
	  if (tp->fDrDzCal < 0 ) FillTrackHistograms(fHist.fTrack[ihist+200+ 1            ],trk,tp,&fSimPar);
	  else                   FillTrackHistograms(fHist.fTrack[ihist+200+ 2            ],trk,tp,&fSimPar);

	  if (tp->fEcl     > 0 ) FillTrackHistograms(fHist.fTrack[ihist+200+ 5+tp->fDiskID],trk,tp,&fSimPar);

	  if (tp->fDt      < -4) FillTrackHistograms(fHist.fTrack[ihist+200+ 7            ],trk,tp,&fSimPar);
	  else                   FillTrackHistograms(fHist.fTrack[ihist+200+ 8            ],trk,tp,&fSimPar);
	}
      }
    }
  }

  if (fNTracks[0] == fNTracks[1]) {
    for (int itrk=0; itrk<fNTracks[0]; itrk++) {
      TStnTrack* trk1  = fTrackBlock[0]->Track(itrk);
      TStnTrack* trk2  = fTrackBlock[1]->Track(itrk);
      TrackPar_t* tp1  = &fTrackPar[0][0]+itrk;
      TrackPar_t* tp2  = &fTrackPar[1][0]+itrk;
      FillDTrackHistograms(fHist.fDTrack[0],trk1,tp1,trk2,tp2);
    }
  }
//-----------------------------------------------------------------------------
// efficiency histograms, use fTrkQual > 0.4 for the cuts
//-----------------------------------------------------------------------------
  // FillEfficiencyHistograms(fTrackBlock[0],fTrackID[fBestID[0]],&fTrackPar[0][0],10);
  // FillEfficiencyHistograms(fTrackBlock[1],fTrackID[fBestID[1]],&fTrackPar[1][0],20);

//-----------------------------------------------------------------------------
// fill little tree if requested
//-----------------------------------------------------------------------------
  if (fWriteTmvaTree >= 0) FillTmvaTree();

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
  fTrackSeedBlock->GetEntry(ientry);
  fHelixBlock->GetEntry(ientry);
  fSpmcBlockVDet->GetEntry(ientry);
//-----------------------------------------------------------------------------
// luminosity weight
//-----------------------------------------------------------------------------
  fLumWt = GetHeaderBlock()->LumWeight();
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fNGenp       = fGenpBlock->NParticles();
  fNSimp       = fSimpBlock->NParticles();
  fNClusters   = fClusterBlock->NClusters();
  fNHelices    = fHelixBlock->NHelices();
  fNTrackSeeds = fTrackSeedBlock->NTrackSeeds();

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

  fProcess = -1;
  fWeight  =  1.;
  fPhotonE = -1.;

  if (fNGenp > 0) {
    TGenParticle* p0 = fGenpBlock->Particle(0);
    if ((p0->GetStatusCode() == 41) && (p0->GetPdgCode() == 22)) {
//-----------------------------------------------------------------------------
// RMC
//-----------------------------------------------------------------------------
      fProcess = 41;
      fPhotonE = p0->Energy();
      fWeight  = TStntuple::RMC_ClosureAppxWeight(fPhotonE,fKMaxRMC);
    }
    else if ((p0->GetStatusCode() == 11) && (p0->GetPdgCode() == 22)) {
//-----------------------------------------------------------------------------
// RPC
//-----------------------------------------------------------------------------
      fProcess = 11;
      fPhotonE = p0->Energy();
      double time_wt = fGenpBlock->Weight();
      fWeight  = TStntuple::RPC_PhotonEnergyWeight(fPhotonE)*time_wt;
    }
  }
//-----------------------------------------------------------------------------
// may want to revisit the definition of fSimPar in future
//-----------------------------------------------------------------------------
  fSimPar.fParticle = fSimpBlock->Particle(0);
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  fSimPar.fGenp     = NULL;
//-----------------------------------------------------------------------------
// sometimes, there is no GENP info - as in Bob's pbar dataset
//-----------------------------------------------------------------------------
  if ((fProcess == -1) && (fNSimp > 0)) {
    if (fSimPar.fParticle->PDGCode() == -2212) {
//-----------------------------------------------------------------------------
// pbar production, assume Bob's dataset
// assuming parent particle exists, determine the production cross-section weight
// pbeam, nx, ny: the beam momentum and direction
//-----------------------------------------------------------------------------
      double pbeam(8.9), nx(-0.24192190), nz(-0.97029573);
      const TLorentzVector* sm = fSimPar.fParticle->StartMom();
      double px    = sm->Px();
      double pz    = sm->Pz();
      double costh = (px*nx+pz*nz)/sqrt(px*px+pz*pz);
      double th    = TMath::ACos(costh);
//-----------------------------------------------------------------------------
// convert momentum to GeV/c
//-----------------------------------------------------------------------------
      double plab  = sm->P()/1000.;  

      fWeight      = fStnt->PBar_Striganov_d2N(pbeam,plab,th);
    }
  }
//-----------------------------------------------------------------------------
// virtual detectors - for fSimPar need parameters at the tracker front
//-----------------------------------------------------------------------------
  int nsteps = fSpmcBlockVDet->NStepPoints();
  
  for (int i=0; i<nsteps; i++) {
    TStepPointMC* step = fSpmcBlockVDet->StepPointMC(i);
    if (step->PDGCode() == fSimPar.fParticle->PDGCode()) {
      if ((step->VolumeID() == 13) || (step->VolumeID() == 14)) {
	fSimPar.fTFront = step;
      }
      else if ((step->VolumeID() == 11) || (step->VolumeID() == 12)) {
	fSimPar.fTMid = step;
      }
    }
  }
//-----------------------------------------------------------------------------
// initialize track seed parameters
//-----------------------------------------------------------------------------
  fNGoodTrackSeeds = 0;
  for (int i=0; i<fNTrackSeeds; i++) {
    TStnTrackSeed* ts = fTrackSeedBlock->TrackSeed(i);
    int nhits = ts->NHits();
    if ((nhits > 15) && (ts->Chi2()/(nhits-4.9999) < 4.)) {
      fNGoodTrackSeeds += 1;
    }
  }
//-----------------------------------------------------------------------------
// initialize additional track parameters
//-----------------------------------------------------------------------------
  for (int i=0; i<2; i++) {
    fNTracks    [i] = fTrackBlock[i]->NTracks();
    fNGoodTracks[i] = 0;

    for (int itrk=0; itrk<fNTracks[i]; itrk++) {
      TrackPar_t* tp  = fTrackPar[i]+itrk;

      tp->fTrackID[0] = TAnaModule::fTrackID_BOX;
      tp->fTrackID[1] = TAnaModule::fTrackID_MVA;
      tp->fLogLH      = TAnaModule::fLogLH;
    }

    InitTrackPar(fTrackBlock[i],fClusterBlock,fTrackPar[i],&fSimPar);
  }

  if (fFillHistograms) FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
// looking mostly at the KDAR tracks
//-----------------------------------------------------------------------------
// diagnostics is related to KDAR
//-----------------------------------------------------------------------------
void TTrackCompModule::Debug() {

  TStnTrack* trk;
  TrackPar_t* tp(NULL);
  char        text[500];
//-----------------------------------------------------------------------------
// bit 0: All Events
//-----------------------------------------------------------------------------
  if (GetDebugBit(0) == 1) {
    GetHeaderBlock()->Print(Form("TTrackCompModule :bit000:"));
    printf("KPAR:\n");
    fTrackBlock[kPAR]->Print();

    for (int i=0; i<fNTracks[kPAR]; i++) { 
      PrintTrack(fTrackBlock[kPAR]->Track(i),&fTrackPar[kPAR][i],"data");
    }

    printf("KDAR:\n");
    fTrackBlock[kDAR]->Print();

    for (int i=0; i<fNTracks[kDAR]; i++) { 
      PrintTrack(fTrackBlock[kDAR]->Track(i),&fTrackPar[kDAR][i],"data");
    }
  }
//-----------------------------------------------------------------------------
// bit 3: KDAR tracks with large DPF > 5 MeV
//-----------------------------------------------------------------------------
  TStnTrackBlock* cprb = fTrackBlock[kDAR];
  int ntrk = cprb->NTracks();
  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = cprb->Track(itrk);
    tp  = &fTrackPar[kDAR][itrk];
    //    if ((GetDebugBit(3) == 1) && (tp->fDpF > 5.) && (trk->NActive() > 25) && (trk->Chi2Dof() < 4)) {
    if ((GetDebugBit(3) == 1) && (tp->fIDWord_RMC == 0) && (tp->fDpF > 1.)) {
      GetHeaderBlock()->Print(Form("TTrackCompModule bit003: tp->DpF = %10.3f trk->fP = %10.3f trk->fPFront = %10.3f nactv: %4i chi2d: %10.3f",
				   tp->fDpF, tp->fP,tp->fPFront,trk->NActive(),trk->Chi2Dof()));
    }
  }
//-----------------------------------------------------------------------------
// bit 4: events with KPAR track, a 60 MeV+ cluster and no KDAR track
//-----------------------------------------------------------------------------
  if (GetDebugBit(4) == 1) {
    if ((fNTracks[0] > 0) && (fNClusters > 0) && (fNTracks[1] == 0)) {
      if (fEClMax > fMinETrig) {
	sprintf(text,"TTrackCompModule bit004: N(TPR) = %i N(CPR) = %i E(cl) = %10.3f",
		fNTracks[0],fNTracks[1],fEClMax);
	GetHeaderBlock()->Print(text);
      }
    }
  }
//-----------------------------------------------------------------------------
// bit 5: Set C KDAR tracks with P > 106
//-----------------------------------------------------------------------------
  if ((GetDebugBit(5) == 1) && (ntrk > 0)) {
    trk = cprb->Track(0);

    if ((trk->fIDWord == 0) && (trk->fP > fDebugCut[5].fXMin) && (trk->fP < fDebugCut[5].fXMax)) {
      GetHeaderBlock()->Print(Form("TTrackCompModule bit005: tp->DpF = %10.3f trk->fP = %10.3f trk->fPFront = %10.3f",
				   tp->fDpF, trk->fP,trk->fPFront));
    }
  }
//-----------------------------------------------------------------------------
// bit 6: Set C KDAR tracks with 1.5 < tp->fDFp < 5
//-----------------------------------------------------------------------------
  if ((GetDebugBit(6) == 1) && (ntrk > 0)) {
    trk = cprb->Track(0);
    tp  = &fTrackPar[kDAR][0];

    int best_id = fBestID[kDAR];

    if ((tp->fIDWord[best_id] == 0) && (tp->fDpF >= fDebugCut[6].fXMin) && (tp->fDpF < fDebugCut[6].fXMax)) {
      GetHeaderBlock()->Print(Form("TTrackCompModule bit006: tp->DpF = %10.3f trk->fP = %10.3f trk->fPFront = %10.3f",
				   tp->fDpF, trk->fP,trk->fPFront));
    }
  }
//-----------------------------------------------------------------------------
// bit 7: events with N(Set C KDAR tracks) > 0
//-----------------------------------------------------------------------------
  if ((GetDebugBit(7) == 1) && (fNGoodTracks[kDAR] > 0)) {
    GetHeaderBlock()->Print(Form("TTrackCompModule :bit007:"));
  }
//-----------------------------------------------------------------------------
// bit 8: Set C KDAR tracks with P > 105
//-----------------------------------------------------------------------------
  if ((GetDebugBit(8) == 1) && (ntrk > 0)) {
    trk = cprb->Track(0);

    if ((trk->fIDWord == 0) && (trk->fP > 105.)) {
      GetHeaderBlock()->Print(Form("TTrackCompModule :bit008: tp->DpF = %10.3f trk->fP = %10.3f trk->fPFront = %10.3f",
				   tp->fDpF, trk->fP,trk->fPFront));
    }
  }

  if (GetDebugBit(11) == 1) {
    if ((fNTracks[0] == 1) && (fNTracks[1] == 1)) {
      GetHeaderBlock()->Print(Form("TTrackCompModule :bit011"));
      PrintTrack(fTrackBlock[0]->Track(0),&fTrackPar[0][0],"");
      PrintTrack(fTrackBlock[1]->Track(0),&fTrackPar[1][0],"data");
    }
  }
//-----------------------------------------------------------------------------
// bit 12: DAR tracks with dtZ0 < -10
//-----------------------------------------------------------------------------
  if ((GetDebugBit(12) == 1) && (ntrk > 0)) {
    for (int itrk=0; itrk<ntrk; itrk++) {
      trk = cprb->Track(itrk);
      tp  = &fTrackPar[kDAR][itrk];
      if ((GetDebugBit(12) == 1) && (tp->fDtZ0 < -10)) {
	GetHeaderBlock()->Print(Form("TTrackCompModule :bit012: tp->fDtZ0 = %10.3f", tp->fDtZ0));
      }
    }
  }
//-----------------------------------------------------------------------------
// bit 13: DAR tracks with dpf > 10
//-----------------------------------------------------------------------------
  if (GetDebugBit(13) == 1) {
    for (int itrk=0; itrk<ntrk; itrk++) {
      trk = cprb->Track(itrk);
      tp  = &fTrackPar[kDAR][itrk];
      if ((tp->fDpF > 10) && (trk->fPFront > 90)) {
	GetHeaderBlock()->Print(Form("TTrackCompModule :bit013: tp->fDpf = %10.3f trk->fPFront = %10.3f", 
				     tp->fDpF,trk->fPFront));
      }
    }
  }
}

//-----------------------------------------------------------------------------
int TTrackCompModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());

  if (fWriteTmvaTree >= 0) {
    printf("[TTrackCompModule::EndJob] Writing output TMVA training file %s\n",fTmvaFile->GetName());
    fTmvaFile->Write();
    delete fTmvaFile;
    fTmvaFile  = 0;
    fSigTree   = 0;
    //    fBgrTree   = 0;
  }

  return 0;
}

//_____________________________________________________________________________
void TTrackCompModule::Test001() {
}


}
