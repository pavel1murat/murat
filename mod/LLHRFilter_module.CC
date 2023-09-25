///////////////////////////////////////////////////////////////////////////////
// Class LLHRFilter : fills Stntuple (P.Murat)
// ------------------------------------------
// order of the data blocks is essential - they are filled according to the
// order in which they are declared...
//
///////////////////////////////////////////////////////////////////////////////

#ifdef __GNUG__
#pragma implementation
#endif

#include <string>
#include <cstdio>

#include "Stntuple/obj/AbsEvent.hh"
#include "Stntuple/obj/TStnEvent.hh"

#include <assert.h>
#include <iostream>
#include <iomanip>

#include "TH1.h"
#include "TEnv.h"
#include "TString.h"
#include "TProfile.h"
#include "TFolder.h"
#include "TSystem.h"

#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnErrorLogger.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"

//  #include "Stntuple/obj/TStnTriggerBlock.hh"

#include "murat/mod/LLHRFilter_module.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/mod/StntupleUtilities.hh"

#include "Stntuple/val/stntuple_val_functions.hh"

#include "Stntuple/alg/TStntuple.hh"

// ClassImp(LLHRFilter)

static const char rcsid[] = "$Name:  $";

void stntuple_get_version(char* ver, char* test);
namespace mu2e {
//------------------------------------------------------------------------------
// constructors
//------------------------------------------------------------------------------
LLHRFilter::LLHRFilter(fhicl::ParameterSet const& PSet): 
  StntupleModule        (PSet,"LLHRFilter")
  , fProcessName        (PSet.get<std::string> ("processName"    ,"PROD"          ))

  , fG4ModuleLabel      (PSet.get<std::string> ("g4ModuleLabel"      , "g4run"        ))
  , fStrawHitMaker      (PSet.get<std::string> ("strawHitMaker"      , "makeSH"       ))
  , fTrkPatRecDem       (PSet.get<std::string> ("trkPatRecDem"       , "trkPatRecDem" ))
  , fTrkPatRecUem       (PSet.get<std::string> ("trkPatRecUem"       , "trkPatRecUem" ))
  , fCaloCrystalHitMaker(PSet.get<std::string> ("caloCrystalHitMaker", "CaloCrystalHitsMaker"))
  , fCaloClusterMaker   (PSet.get<std::string> ("caloClusterMaker"   , "makeCaloCluster"))
  , fTrkExtrapol        (PSet.get<std::string> ("trkExtrapol"        , "trkExtrapol"  ))
  , fTrkCalMatch        (PSet.get<std::string> ("trkCalMatch"        , "caloMatching" ))
  , fPidDem             (PSet.get<std::string> ("pidDem"             , "undefined"    ))
  , fFilterEp           (PSet.get<int>         ("filterEp"           , 0    ))
  , fMinEP              (PSet.get<double>      ("minEP"              , 1.1  ))
  , fFilterLLHRCal      (PSet.get<int>         ("filterLLHRCal"      , 1    ))
  , fMaxLLHRCal         (PSet.get<double>      ("maxLLHRCal"         , -20. ))
  , fMinTActive         (PSet.get<double>      ("minTActive"         , 700. ))
{

  char  ver[100], text[100];
  stntuple_get_version(ver,text);

  fVersion      = new TNamed(ver,text);
  TModule::fFolder->Add(fVersion);

  fTrackID      = new TStnTrackID();
  fLogLH        = new TEmuLogLH();

  fPointersInitialized = 0;
}


//------------------------------------------------------------------------------
LLHRFilter::~LLHRFilter() {
  delete fVersion;
}


//------------------------------------------------------------------------------
bool LLHRFilter::beginRun(art::Run& aRun) {

  static int first_begin_run = 1;

  THistModule::beforeBeginRun(aRun);

  if (first_begin_run) {
//-----------------------------------------------------------------------------
// if we runnning stnmaker_prod.exe, save revision of the TCL file in STNTUPLE
//-----------------------------------------------------------------------------
    first_begin_run = 0;
    const char* c = gSystem->Getenv("STNMAKER_PROD_TCL");
    if (c) TModule::fFolder->Add(new TNamed("STNMAKER_PROD_TCL",c));
    else   TModule::fFolder->Add(new TNamed("STNMAKER_PROD_TCL","unknown"));
  }

  THistModule::afterBeginRun(aRun);

  return 1;
}

//------------------------------------------------------------------------------
bool LLHRFilter::endRun(art::Run& aRun ) {
  THistModule::beforeEndRun(aRun);
  THistModule::afterEndRun (aRun);
  return 1;
}


//------------------------------------------------------------------------------
void LLHRFilter::endJob() {

  THistModule::beforeEndJob();
  THistModule::afterEndJob ();

}

//------------------------------------------------------------------------------
void LLHRFilter::beginJob() {

  //  int split_mode, compression_level, buffer_size;

  THistModule::beforeBeginJob();

//-----------------------------------------------------------------------------
// initialize likelihood histograms
//-----------------------------------------------------------------------------
  const char*  mu2e_hist_dir;
  char         fn[200];

  mu2e_hist_dir = gEnv->GetValue("mu2e.HistDir",gSystem->Getenv("MU2E_HIST_DIR"));

  sprintf(fn,"%s/v4_1_8/e00s1202.track_ana.hist",mu2e_hist_dir);
  fLogLH->InitEleDtHist(fn);
  fLogLH->InitEleEpHist(fn);

  sprintf(fn,"%s/v4_1_8/m00s1202.track_ana.hist",mu2e_hist_dir);
  fLogLH->InitMuoDtHist(fn);
  fLogLH->InitMuoEpHist(fn);

  sprintf(fn,"%s/v4_1_8/e00s1202.track_ana.hist",mu2e_hist_dir);
  TH1* he_xs = gh1(fn,"TrackAna","trk_1/xslope");
  fLogLH->SetEleXsHist(he_xs);

  sprintf(fn,"%s/v4_1_8/m00s1202.track_ana.hist",mu2e_hist_dir);
  TH1* hm_xs = gh1(fn,"TrackAna","trk_1/xslope");
  fLogLH->SetMuoXsHist(hm_xs);

  THistModule::afterEndJob();
}

//_____________________________________________________________________________
bool LLHRFilter::filter(AbsEvent& AnEvent) {
  bool rc = false;
  // when execution comes here al the registered data blocks are already
  // initialized with the event data. Left: variables in the data blocks
  // which depend on the variable defined in other blocks, like track number
  // for a muon or an electron - the idea is that these are defined during the
  // 2nd loop in FillStntupleModule, where ResolveLinks methods are called
  // for each data block

//-----------------------------------------------------------------------------
// connect to the error reporting facility
//-----------------------------------------------------------------------------
//  TStnErrorLogger* logger = Event()->GetErrorLogger();
//   logger->Connect("Report(Int_t, const char*)",
// 		  "StntupleModule",
// 		  this,
// 		  "LogError(const char*)");
//-----------------------------------------------------------------------------
// disconnect from the error reporting signal and return back to AC++
//-----------------------------------------------------------------------------
//   logger->Disconnect("Report(Int_t,const char*)",
// 		     this,"LogError(Int_t,const char*)");

  if (fPointersInitialized == 0) {
//-----------------------------------------------------------------------------
// track branches: for ROOT v3 to use streamers one has to specify split=-1
//-----------------------------------------------------------------------------
    fTrackBlock = (TStnTrackBlock*) Event()->GetDataBlock("TrackBlock");

    fTrackBlock->AddCollName("mu2e::KalRepCollection"              ,fTrkPatRecDem.data()    ,"DownstreameMinus");
    fTrackBlock->AddCollName("mu2e::CaloClusterCollection"         ,fCaloClusterMaker.data(),"AlgoCLOSESTSeededByENERGY");
    fTrackBlock->AddCollName("mu2e::TrkToCaloExtrapolCollection"   ,fTrkExtrapol.data()     ,"");
    fTrackBlock->AddCollName("mu2e::TrackClusterLink"              ,fTrkCalMatch.data()     ,"");
    fTrackBlock->AddCollName("mu2e::StrawHitCollection"            ,fStrawHitMaker.data()   ,"");
    fTrackBlock->AddCollName("mu2e::PtrStepPointMCVectorCollection",fStrawHitMaker.data()   ,"StrawHitMCPtr");
    fTrackBlock->AddCollName("mu2e::PIDProductCollection"          ,fPidDem.data()          ,"");
    fTrackBlock->AddCollName("mu2e::StepPointMCCollection"         ,fG4ModuleLabel.data()   ,"");

    fPointersInitialized = 1;
  }

  TStnTrack* track;
  TEmuLogLH::PidData_t  dat;

  double ep, llhr_cal;

  //  int    id_word;

  fTrackBlock = (TStnTrackBlock*) Event()->GetDataBlock("TrackBlock");

  StntupleInitMu2eTrackBlock(fTrackBlock,&AnEvent,0);

  int nt = fTrackBlock->NTracks();
//-----------------------------------------------------------------------------
// calculate likelihoods
//-----------------------------------------------------------------------------
  double xs;
  for (int i=0; i<nt; i++) {
    track   = fTrackBlock->Track(i);

    xs = track->XSlope();
    track->fLogLHRXs    = fLogLH->LogLHRXs(xs);

    dat.fDt   = track->Dt();
    dat.fEp   = track->Ep();
    dat.fPath = -1.;
    if (track->fVMinS) dat.fPath = track->fVMinS->fPath;
      
    track->fEleLogLHCal = fLogLH->LogLHCal(&dat,11);
    track->fMuoLogLHCal = fLogLH->LogLHCal(&dat,13);
  }
//-----------------------------------------------------------------------------
// perform selections
//-----------------------------------------------------------------------------
  if (fFilterEp != 0) {
    for (int i=0; i<nt; i++) {
      track   = fTrackBlock->Track(i);
      //      id_word = fTrackID->IDWord(track);
      ep = track->Ep();

      if (ep > fMinEP) {
	rc = true;
	break;
      }
    }
  }
  else if (fFilterLLHRCal != 0) {
    for (int i=0; i<nt; i++) {
      track = fTrackBlock->Track(i);
      llhr_cal = track->LogLHRCal();
      if (llhr_cal < fMaxLLHRCal) {
	rc = true;
	break;
      }
    }
  }

  return rc;
}

// //_____________________________________________________________________________
// int LLHRFilter::InitCalDataBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eCalDataBlock(Block,event,mode);
// }

// //_____________________________________________________________________________
// int LLHRFilter::InitHeaderBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eHeaderBlock(Block,event,mode);
// }

//_____________________________________________________________________________
// int LLHRFilter::InitTriggerBlock(TStnDataBlock* Block) {
//   int mode = 0;
//   AbsEvent* event = AbsEnv::instance()->theEvent();
//   return StntupleInitMu2eTriggerBlock(Block,event,mode);
// }

} // end namespace mu2e

using mu2e::LLHRFilter;

DEFINE_ART_MODULE(LLHRFilter)
