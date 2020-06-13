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
#include "murat/ana/TTrackCompModule.hh"
#include "murat/ana/mva_data.hh"

// framework
#include "fhiclcpp/ParameterSet.h"

// Xerces XML Parser
#include <xercesc/dom/DOM.hpp>

#include "Mu2eUtilities/inc/MVATools.hh"
#include <string>
#include "math.h"

using std::string;
using std::vector;

ClassImp(TTrackCompModule)
//-----------------------------------------------------------------------------
TTrackCompModule::TTrackCompModule(const char* name, const char* title):
  TStnModule    (name,title),
  fPdgCode      (11),                       // electron
  fGeneratorCode( 2)                        // 2:ConversionGun 28:StoppedParticleReactionGun
{
					    // this is the default for MDC2018 stntuples
  fTrackBlockName[0] = "TrackBlockPar";
  fTrackBlockName[1] = "TrackBlockDar";
//-----------------------------------------------------------------------------
// TrackID[0] : "SetC"
// i = 1..6 : cut on DaveTrkQual > 0.1*i instead
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
// 
//-----------------------------------------------------------------------------
  fMinETrig       = 50.;
  fNMVA           = 0;

  fUseMVA         = 0;

  fTprMVA         = new mva_data("trkpatrec","dave"    ,002);
  fCprMVA         = new mva_data("calpatrec","e11s5731",002);

  fLogLH          = new TEmuLogLH();
//-----------------------------------------------------------------------------
// debugging information
//-----------------------------------------------------------------------------
  fDebugCut[5].fXMin = 106.;
  fDebugCut[5].fXMax = 200.;

  fDebugCut[6].fXMin = 1.5;
  fDebugCut[6].fXMax = 10.0;
//-----------------------------------------------------------------------------
// ntuples for TMVA training
//-----------------------------------------------------------------------------
  fWriteTmvaTree     = -1;
  fTmvaAlgorithmTpr  = -1;
  fTmvaAlgorithmCpr  = -1;
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
void TTrackCompModule::SetMVA(const char* TrkRecAlgorithm, const char* TrainingDataset, int MvaType) {

  fUseMVA         = 1;

  TString trk_alg = TrkRecAlgorithm;
  trk_alg.ToUpper();

  if (trk_alg == "CALPATREC") {
    if (fCprMVA) delete fCprMVA;
    fCprMVA           = new mva_data("CALPATREC",TrainingDataset,MvaType);
    fTmvaAlgorithmCpr = MvaType;
  }
  else if (trk_alg == "TRKPATREC") {
    if (fTprMVA) delete fTprMVA;
    fTprMVA           = new mva_data("TRKPATREC",TrainingDataset,MvaType);
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

  string hist_dir      = gEnv->GetValue("mu2e.TrkQual.HistDir","_none_");
  string trk_qual_dsid = gEnv->GetValue("mu2e.TrkQual.Dsid"   ,"_none_");

  fTrkQualFile    = Form("%s/%s.track_comp_use_mva_%03i.hist",
			 hist_dir.data(),trk_qual_dsid.data(),fUseMVA);

  fTrackProb[0]   = new prob_dist(fTrkQualFile.Data(),"TrackComp","trk_100/mvaout");
  fTrackProb[1]   = new prob_dist(fTrkQualFile.Data(),"TrackComp","trk_200/mvaout");

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
void TTrackCompModule::BookTrackSeedHistograms   (HistBase_t*   HistR, const char* Folder){
  
  TrackSeedHist_t* Hist =  (TrackSeedHist_t*) HistR;

  HBook1F(Hist->fNHits        ,"nhits" ,Form("%s: # of straw hits"              ,Folder),  150,    0,  150,Folder);
  HBook1F(Hist->fClusterTime  ,"clt"   ,Form("%s: cluster time; t_{cluster}[ns]",Folder),  800,  400, 1700,Folder);
  HBook1F(Hist->fClusterEnergy,"cle"   ,Form("%s: cluster energy; E [MeV]      ",Folder),  400,    0,  200,Folder);
  HBook1F(Hist->fRadius       ,"r"     ,Form("%s: curvature radius; r [mm]"     ,Folder),  500,    0,  500,Folder);
  HBook1F(Hist->fMom          ,"p"     ,Form("%s: momentum; p [MeV/c]"          ,Folder),  300,    50, 200,Folder);
  HBook1F(Hist->fPt           ,"pt"    ,Form("%s: pT; pT [MeV/c]"               ,Folder),  600,    0,  150,Folder);
  HBook1F(Hist->fTanDip       ,"tdip"  ,Form("%s: tanDip; tanDip"               ,Folder),  300,    0,    3,Folder);
  HBook1F(Hist->fChi2         ,"chi2"  ,Form("%s: #chi^{2}-XY; #chi^{2}/ndof"   ,Folder),  100,    0,   10,Folder);
  HBook1F(Hist->fFitCons      ,"fcons" ,Form("%s: Fit consistency; Fit-cons"    ,Folder),  100,    0,    1,Folder);
  HBook1F(Hist->fD0           ,"d0"    ,Form("%s: D0; d0 [mm]"                  ,Folder), 1600, -400,  400,Folder);
}

//-----------------------------------------------------------------------------
void TTrackCompModule::BookTrackHistograms(HistBase_t* HistR, const char* Folder) {
//   char name [200];
//   char title[200];

  TrackHist_t* Hist =  (TrackHist_t*) HistR;

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
void TTrackCompModule::BookDTrackHistograms(HistBase_t* HistR, const char* Folder) {
  DTrackHist_t* Hist =  (DTrackHist_t*) HistR;

  HBook1F(Hist->fDp       ,"dp"      ,Form("%s: P(1)-P(2)"               ,Folder),200, -1,1,Folder);
  HBook1F(Hist->fRMomErr10,"rmomerr" ,Form("%s: MomEff(1)/MomErr(0)"     ,Folder),200,  0,2,Folder);
  HBook1F(Hist->fDT0      ,"dt0"     ,Form("%s: T0(1)-T0(2)"             ,Folder),200, -5,5,Folder);


}

//-----------------------------------------------------------------------------
void TTrackCompModule::BookEventHistograms(HistBase_t* HistR, const char* Folder) {
  //  char name [200];
  //  char title[200];

  EventHist_t* Hist =  (EventHist_t*) HistR;

  HBook1D(Hist->fLumWt     ,"lumwt"    ,Form("%s: Luminosity Weight"               ,Folder),200, 0,10,Folder);
  Hist->fLumWt->Sumw2(kTRUE);

  HBook1F(Hist->fRv        ,"rv"       ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 1000,Folder);
  HBook1F(Hist->fZv        ,"zv"       ,Form("%s: Z(Vertex)"                       ,Folder), 300, 0,15000,Folder);

  HBook1F(Hist->fPdgCode         ,"pdg"   ,Form("%s: PDG code"                     ,Folder),200,-100,100,Folder);
  HBook1F(Hist->fMomTargetEnd    ,"ptarg" ,Form("%s: MC mom after Stopping Target" ,Folder),400,  90,110,Folder);
  HBook1F(Hist->fMomTrackerFront ,"pfront",Form("%s: MC mom at the Tracker Front"  ,Folder),400,  90,110,Folder);
  HBook1F(Hist->fNshCE           ,"nsh_ce",Form("%s: CE Number of Straw Hits"      ,Folder),150,0,150,Folder);

  HBook1F(Hist->fMcCosTh   ,"mc_costh" ,Form("%s: MC Particle Cos(Theta) Lab"      ,Folder),100,-1,1,Folder);
  HBook1F(Hist->fMcMom     ,"mc_mom"   ,Form("%s: MC Particle Momentum"            ,Folder),1000,  0,200,Folder);

  HBook1F(Hist->fNHelices  ,"nhel"     ,Form("%s: nhelices"                        ,Folder), 10,0, 10,Folder);
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
  HBook1F(Hist->fInstLumi  ,"inst_lum" ,Form("%s: Inst Luminosity"                 ,Folder),500, 0,1.e8,Folder);
  HBook1F(Hist->fGMom      ,"gmom"     ,Form("%s: Photon Momentum"                 ,Folder),200, 0,200 ,Folder);
  HBook1F(Hist->fGMomRMC   ,"gmom_rmc" ,Form("%s: Photon Momentum, RMC weighted"   ,Folder),200, 0,200 ,Folder);
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
  int       book_track_histset[kNTrackHistSets];
  TString*  track_selection   [kNTrackHistSets];

  for (int i=0; i<kNTrackHistSets; i++) { book_track_histset[i] = 0; track_selection[i] = NULL; }

  book_track_histset[  0] = 1; track_selection[  0] = new TString("good PAR tracks");
  book_track_histset[  1] = 1; track_selection[  1] = new TString("good DAR tracks in events with no good kPAR tracks TrkQual>0.3");
  book_track_histset[  2] = 1; track_selection[  2] = new TString("good DAR tracks in events with no good kPAR tracks TrkQual>0.2");
  book_track_histset[  3] = 1; track_selection[  3] = new TString("best track");

  book_track_histset[100] = 1; track_selection[100] = new TString("PAR all tracks");
  book_track_histset[101] = 1; track_selection[101] = new TString("PAR BestTrackID");
  book_track_histset[102] = 1; track_selection[102] = new TString("PAR BestTrackID no fitCons&momErr&t0Err tracks");
  book_track_histset[103] = 1; track_selection[103] = new TString("PAR BestTrackID, dpf>1 tracks");
  book_track_histset[104] = 1; track_selection[104] = new TString("PAR all  tracks events Ecl > 60");
  book_track_histset[105] = 1; track_selection[105] = new TString("PAR BestTrackID tracks events Ecl > 60");
  book_track_histset[106] = 1; track_selection[106] = new TString("PAR tracks with dtz0 < -10");
  book_track_histset[109] = 1; track_selection[109] = new TString("PAR tracks with XDpF > +10 MeV");

  book_track_histset[110] = 1; track_selection[110] = new TString("PAR TrackID[0] - SetC");
  book_track_histset[111] = 1; track_selection[111] = new TString("PAR TrackID[1]");
  book_track_histset[112] = 1; track_selection[112] = new TString("PAR TrackID[2]");
  book_track_histset[113] = 1; track_selection[113] = new TString("PAR TrackID[3]");
  book_track_histset[114] = 1; track_selection[114] = new TString("PAR TrackID[4]");
  book_track_histset[115] = 1; track_selection[115] = new TString("PAR TrackID[5]");
  book_track_histset[116] = 1; track_selection[116] = new TString("PAR TrackID[6]");
  book_track_histset[117] = 1; track_selection[117] = new TString("PAR TrackID[7]");
  book_track_histset[118] = 1; track_selection[118] = new TString("PAR TrackID[8]");
  book_track_histset[119] = 1; track_selection[119] = new TString("PAR TrackID[9]");
  book_track_histset[120] = 1; track_selection[120] = new TString("PAR TrackID[10]");
  book_track_histset[121] = 1; track_selection[121] = new TString("PAR TrackID[11]");
  book_track_histset[122] = 1; track_selection[122] = new TString("PAR TrackID[12]");
  book_track_histset[123] = 1; track_selection[123] = new TString("PAR TrackID[13]");
  book_track_histset[124] = 1; track_selection[124] = new TString("PAR TrackID[14]");
  book_track_histset[125] = 1; track_selection[125] = new TString("PAR TrackID[15]");
  book_track_histset[126] = 1; track_selection[126] = new TString("PAR TrackID[16]");
  book_track_histset[127] = 1; track_selection[127] = new TString("PAR TrackID[17]");
  book_track_histset[128] = 1; track_selection[128] = new TString("PAR TrackID[18]");
  book_track_histset[129] = 1; track_selection[129] = new TString("PAR TrackID[19]");

  book_track_histset[141] = 1; track_selection[141] = new TString("PAR tracks N(active) > 20");
  book_track_histset[142] = 1; track_selection[142] = new TString("PAR tracks N(active) > 20, |D0| < 100");
  book_track_histset[143] = 1; track_selection[143] = new TString("PAR tracks N(active) > 20, |D0| < 100, DN(active) < 6");
  book_track_histset[144] = 1; track_selection[144] = new TString("PAR tracks N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");

  book_track_histset[150] = 1; track_selection[150] = new TString("PAR tracks |DpF| < 1");
  book_track_histset[151] = 1; track_selection[151] = new TString("PAR tracks |DpF| < 1 ; N(active) > 20");
  book_track_histset[152] = 1; track_selection[152] = new TString("PAR tracks |DpF| < 1 ; N(active) > 20, |D0| < 100");
  book_track_histset[153] = 1; track_selection[153] = new TString("PAR tracks |DpF| < 1 ; N(active) > 20, |D0| < 100, DN(active) < 6");
  book_track_histset[154] = 1; track_selection[154] = new TString("PAR tracks |DpF| < 1 ; N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");

  book_track_histset[160] = 1; track_selection[160] = new TString("PAR tracks DpF   > 2");
  book_track_histset[161] = 1; track_selection[161] = new TString("PAR tracks DpF   > 2 ; N(active) > 20");
  book_track_histset[162] = 1; track_selection[164] = new TString("PAR tracks DpF   > 2 ; N(active) > 20, |D0| < 100");
  book_track_histset[163] = 1; track_selection[163] = new TString("PAR tracks DpF   > 2 ; N(active) > 20, |D0| < 100, DN(active) < 6");
  book_track_histset[164] = 1; track_selection[164] = new TString("PAR tracks DpF   > 2 ; N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");

  book_track_histset[171] = 1; track_selection[171] = new TString("PAR+ tracks with final selections");
  book_track_histset[172] = 1; track_selection[172] = new TString("PAR- tracks with final selections");
  book_track_histset[173] = 1; track_selection[173] = new TString("PAR+ tracks with final selections and process-defined weight");
  book_track_histset[174] = 1; track_selection[174] = new TString("PAR- tracks with final selections and process-defined weight");

  book_track_histset[177] = 1; track_selection[177] = new TString("cosmics+: #171 + 90 < p < 93");
  book_track_histset[178] = 1; track_selection[178] = new TString("cosmics-: #172 + 90 < p < 93");

  book_track_histset[179] = 1; track_selection[179] = new TString("PAR- tracks with final selections and DIO weight");

  book_track_histset[180] = 1; track_selection[180] = new TString("PAR+ tracks all");                                       
  book_track_histset[181] = 1; track_selection[181] = new TString("PAR+ tracks N(active) > 20");                                       
  book_track_histset[182] = 1; track_selection[182] = new TString("PAR+ tracks N(active) > 20 and |D0| < 100");			       
  book_track_histset[183] = 1; track_selection[183] = new TString("PAR+ tracks N(active) > 20, |D0| < 100, DN(active) < 6");	       
  book_track_histset[184] = 1; track_selection[184] = new TString("PAR+ tracks N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");
  book_track_histset[185] = 1; track_selection[185] = new TString("PAR+ tracks with final selections, E/P < 0.7 (mu+/pi+)");
  book_track_histset[186] = 1; track_selection[186] = new TString("PAR+ tracks with final selections, E/P > 0.7 (e+)");

  book_track_histset[190] = 1; track_selection[190] = new TString("PAR- tracks all");                                       
  book_track_histset[191] = 1; track_selection[191] = new TString("PAR- tracks N(active) > 20");                                       
  book_track_histset[192] = 1; track_selection[192] = new TString("PAR- tracks N(active) > 20 and |D0| < 100");
  book_track_histset[193] = 1; track_selection[193] = new TString("PAR- tracks N(active) > 20, |D0| < 100, DN(active) < 6");
  book_track_histset[194] = 1; track_selection[194] = new TString("PAR- tracks N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");
  book_track_histset[195] = 1; track_selection[195] = new TString("PAR- tracks with final selections, E/P < 0.7 (mu-/pi-)");
  book_track_histset[196] = 1; track_selection[196] = new TString("PAR- tracks with final selections, E/P > 0.7 (e-)");

  book_track_histset[200] = 1; track_selection[200] = new TString("DAR all tracks");
  book_track_histset[201] = 1; track_selection[201] = new TString("DAR BestTrackID");
  book_track_histset[202] = 1; track_selection[202] = new TString("DAR BestTrackID no fitCons&momErr&t0Err tracks");
  book_track_histset[203] = 1; track_selection[203] = new TString("DAR BestTrackID, dpf>1 tracks");
  book_track_histset[204] = 1; track_selection[204] = new TString("DAR all  tracks events Ecl > 60");
  book_track_histset[205] = 1; track_selection[205] = new TString("DAR BestTrackID tracks events Ecl > 60");
  book_track_histset[206] = 1; track_selection[206] = new TString("DAR tracks with dtz0 < -10");
  book_track_histset[209] = 1; track_selection[209] = new TString("DAR tracks with XDpF > +10 MeV");

  book_track_histset[210] = 1; track_selection[210] = new TString("DAR TrackID[0] - SetC");
  book_track_histset[211] = 1; track_selection[211] = new TString("DAR TrackID[1]");
  book_track_histset[212] = 1; track_selection[212] = new TString("DAR TrackID[2]");
  book_track_histset[213] = 1; track_selection[213] = new TString("DAR TrackID[3]");
  book_track_histset[214] = 1; track_selection[214] = new TString("DAR TrackID[4]");
  book_track_histset[215] = 1; track_selection[215] = new TString("DAR TrackID[5]");
  book_track_histset[216] = 1; track_selection[216] = new TString("DAR TrackID[6]");
  book_track_histset[217] = 1; track_selection[217] = new TString("DAR TrackID[7]");
  book_track_histset[218] = 1; track_selection[218] = new TString("DAR TrackID[8]");
  book_track_histset[219] = 1; track_selection[219] = new TString("DAR TrackID[9]");
  book_track_histset[220] = 1; track_selection[220] = new TString("DAR TrackID[10]");
  book_track_histset[221] = 1; track_selection[221] = new TString("DAR TrackID[11]");
  book_track_histset[222] = 1; track_selection[222] = new TString("DAR TrackID[12]");
  book_track_histset[223] = 1; track_selection[223] = new TString("DAR TrackID[13]");
  book_track_histset[224] = 1; track_selection[224] = new TString("DAR TrackID[14]");
  book_track_histset[225] = 1; track_selection[225] = new TString("DAR TrackID[15]");
  book_track_histset[226] = 1; track_selection[226] = new TString("DAR TrackID[16]");
  book_track_histset[227] = 1; track_selection[227] = new TString("DAR TrackID[17]");
  book_track_histset[228] = 1; track_selection[228] = new TString("DAR TrackID[18]");
  book_track_histset[229] = 1; track_selection[229] = new TString("DAR TrackID[19");

  book_track_histset[241] = 1; track_selection[241] = new TString("DAR tracks N(active) > 20");
  book_track_histset[242] = 1; track_selection[242] = new TString("DAR tracks N(active) > 20, |D0| < 100");
  book_track_histset[243] = 1; track_selection[243] = new TString("DAR tracks N(active) > 20, |D0| < 100, DN(active) < 6");
  book_track_histset[244] = 1; track_selection[244] = new TString("DAR tracks N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");

  book_track_histset[250] = 1; track_selection[250] = new TString("DAR tracks |DpF| < 1");
  book_track_histset[251] = 1; track_selection[251] = new TString("DAR tracks |DpF| < 1 ; N(active) > 20");
  book_track_histset[252] = 1; track_selection[252] = new TString("DAR tracks |DpF| < 1 ; N(active) > 20, |D0| < 100");
  book_track_histset[253] = 1; track_selection[253] = new TString("DAR tracks |DpF| < 1 ; N(active) > 20, |D0| < 100, DN(active) < 6");
  book_track_histset[254] = 1; track_selection[254] = new TString("DAR tracks |DpF| < 1 ; N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");

  book_track_histset[260] = 1; track_selection[260] = new TString("DAR tracks DpF   > 2");
  book_track_histset[261] = 1; track_selection[261] = new TString("DAR tracks DpF   > 2 ; N(active) > 20");
  book_track_histset[262] = 1; track_selection[264] = new TString("DAR tracks DpF   > 2 ; N(active) > 20, |D0| < 100");
  book_track_histset[263] = 1; track_selection[263] = new TString("DAR tracks DpF   > 2 ; N(active) > 20, |D0| < 100, DN(active) < 6");
  book_track_histset[264] = 1; track_selection[264] = new TString("DAR tracks DpF   > 2 ; N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");

  book_track_histset[271] = 1; track_selection[271] = new TString("DAR+ tracks with final selections");
  book_track_histset[272] = 1; track_selection[272] = new TString("DAR- tracks with final selections");
  book_track_histset[273] = 1; track_selection[273] = new TString("DAR+ tracks with final selections and process-defined weight");
  book_track_histset[274] = 1; track_selection[274] = new TString("DAR- tracks with final selections and process-defined weight");

  book_track_histset[277] = 1; track_selection[277] = new TString("cosmics+: #271 + 90 < p < 93");
  book_track_histset[278] = 1; track_selection[278] = new TString("cosmics-: #272 + 90 < p < 93");

  book_track_histset[279] = 1; track_selection[279] = new TString("DAR- tracks with final selections and DIO weight");

  book_track_histset[280] = 1; track_selection[280] = new TString("DAR+ tracks all");                                       
  book_track_histset[281] = 1; track_selection[281] = new TString("DAR+ tracks N(active) > 20");                                       
  book_track_histset[282] = 1; track_selection[282] = new TString("DAR+ tracks N(active) > 20 and |D0| < 100");			       
  book_track_histset[283] = 1; track_selection[283] = new TString("DAR+ tracks N(active) > 20, |D0| < 100, DN(active) < 6");	       
  book_track_histset[284] = 1; track_selection[284] = new TString("DAR+ tracks N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");
  book_track_histset[285] = 1; track_selection[285] = new TString("DAR+ tracks with final selections, E/P < 0.7 (mu+/pi+)");
  book_track_histset[286] = 1; track_selection[286] = new TString("DAR+ tracks with final selections, E/P > 0.7 (e+)");

  book_track_histset[290] = 1; track_selection[290] = new TString("DAR- tracks all");                                       
  book_track_histset[291] = 1; track_selection[291] = new TString("DAR- tracks N(active) > 20");                                       
  book_track_histset[292] = 1; track_selection[292] = new TString("DAR- tracks N(active) > 20 and |D0| < 100");
  book_track_histset[293] = 1; track_selection[293] = new TString("DAR- tracks N(active) > 20, |D0| < 100, DN(active) < 6");
  book_track_histset[294] = 1; track_selection[294] = new TString("DAR- tracks N(active) > 20, |D0| < 100, DN(active) < 6, chi2d < 4");
  book_track_histset[295] = 1; track_selection[295] = new TString("DAR- tracks with final selections, E/P < 0.7 (mu-/pi-)");
  book_track_histset[296] = 1; track_selection[296] = new TString("DAR- tracks with final selections, E/P > 0.7 (e-)");

  book_track_histset[301] = 1; track_selection[301] = new TString("PAR- tracks final selections, dr/dz(cal) <  0");
  book_track_histset[302] = 1; track_selection[302] = new TString("PAR- tracks final selections, dr/dz(cal) >= 0");
  book_track_histset[305] = 1; track_selection[305] = new TString("PAR- tracks final selections, cluster in the 1st disk");
  book_track_histset[306] = 1; track_selection[306] = new TString("PAR- tracks final selections, cluster in the 2nd disk");
  book_track_histset[307] = 1; track_selection[307] = new TString("PAR- tracks final selections, dt <  -4 ns");
  book_track_histset[308] = 1; track_selection[308] = new TString("PAR- tracks final selections, dt >= -4 ns");

  book_track_histset[401] = 1; track_selection[401] = new TString("DAR- tracks final selections, dr/dz(cal) <  0");
  book_track_histset[402] = 1; track_selection[402] = new TString("DAR- tracks final selections, dr/dz(cal) >= 0");
  book_track_histset[405] = 1; track_selection[405] = new TString("DAR- tracks final selections, cluster in the 1st disk");
  book_track_histset[406] = 1; track_selection[406] = new TString("DAR- tracks final selections, cluster in the 2nd disk");
  book_track_histset[407] = 1; track_selection[407] = new TString("DAR- tracks final selections, dt <  -4 ns");
  book_track_histset[408] = 1; track_selection[408] = new TString("DAR- tracks final selections, dt >= -4 ns");

  const char* folder_title;
  for (int i=0; i<kNTrackHistSets; i++) {
    if (book_track_histset[i] != 0) {
      sprintf(folder_name,"trk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title = folder_name;
      if (track_selection[i] != NULL) folder_title = track_selection[i]->Data();
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fTrack[i] = new TrackHist_t;
      BookTrackHistograms(fHist.fTrack[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book delta track histograms
//-----------------------------------------------------------------------------
//  int       book_dtrack_histset[kNDTrackHistSets];
  TString*  dtrack_selection   [kNDTrackHistSets];

  for (int i=0; i<kNDTrackHistSets; i++) { dtrack_selection[i] = NULL; }

  dtrack_selection[  0] = new TString("all PAR-DAR tracks");

  //  const char* folder_title;
  for (int i=0; i<kNDTrackHistSets; i++) {
    if (dtrack_selection[i] != 0) {
      sprintf(folder_name,"dtrk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title = folder_name;
      if (dtrack_selection[i] != NULL) folder_title = dtrack_selection[i]->Data();
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fDTrack[i] = new DTrackHist_t;
      BookDTrackHistograms(fHist.fDTrack[i],Form("Hist/%s",folder_name));
    }
  }

}

//-----------------------------------------------------------------------------
// need MC truth branch
//-----------------------------------------------------------------------------
void TTrackCompModule::FillEventHistograms(HistBase_t* HistR) {

  double            cos_th, xv(-1.e6), yv(-1.e6), rv(-1.e6), zv(-1.e6), p;
  TLorentzVector    mom (1.,0.,0.,0);

  EventHist_t* Hist = (EventHist_t*) HistR;

  if (fParticle) { 
    fParticle->Momentum(mom);
    xv = fParticle->Vx()+3904.;
    yv = fParticle->Vy();
    rv = sqrt(xv*xv+yv*yv);
    zv = fParticle->Vz();
  }

  p      = mom.P();
  cos_th = mom.Pz()/p;

  // xv = fParticle->Vx()+3904.;
  // yv = fParticle->Vy();
  // rv = sqrt(xv*xv+yv*yv);
  // zv = fParticle->Vz();

  Hist->fLumWt->Fill(fLumWt);
  Hist->fRv->Fill(rv);
  Hist->fZv->Fill(zv);

  TSimParticle* simp = fSimPar.fParticle;

  if (simp) {
    Hist->fPdgCode->Fill(simp->fPdgCode);
    Hist->fMomTargetEnd->Fill(simp->fMomTargetEnd);
    Hist->fMomTrackerFront->Fill(simp->fMomTrackerFront);	// 
    Hist->fNshCE->Fill(simp->fNStrawHits);
  }
  else {
    Hist->fPdgCode->Fill(-1.e6);
    Hist->fMomTargetEnd->Fill(-1.);
    Hist->fMomTrackerFront->Fill(-1.);	// 
    Hist->fNshCE->Fill(-1);
  }

  Hist->fMcMom->Fill(p);
  Hist->fMcCosTh->Fill(cos_th);

  Hist->fNHelices->Fill(fNHelices);

  Hist->fNTracks[0]->Fill(fNTracks[0]);
  Hist->fNTracks[1]->Fill(fNTracks[1]);

  int nsh_tot = GetHeaderBlock()->fNStrawHits;

  Hist->fNshTot[0]->Fill(nsh_tot);
  Hist->fNshTot[1]->Fill(nsh_tot);

  Hist->fNClusters->Fill(fNClusters);
  Hist->fEClMax->Fill(fEClMax);
  Hist->fTClMax->Fill(fTClMax);
//-----------------------------------------------------------------------------
//  difference between the [corrected] PAR and DAR momenta
//-----------------------------------------------------------------------------
  double dp {1.e6};

  if ((fNTracks[0] == 1) && (fNTracks[1] == 1)) {
    TrackPar_t* tp = &fTrackPar[0][0];
    TrackPar_t* cp = &fTrackPar[1][0];

    dp = tp->fP-cp->fP;
  }
  Hist->fDp->Fill(dp);
//-----------------------------------------------------------------------------
// photon momentum plots (for RMC )
//-----------------------------------------------------------------------------
  if (fProcess == 41) {
    Hist->fGMom->Fill   (fPhotonE, 1.);
    Hist->fGMomRMC->Fill(fPhotonE, fWeight);
  }
}

//-----------------------------------------------------------------------------
// fill efficiency histograms : need 10 histogram sets
// pitch = 1./tan(dip)
//-----------------------------------------------------------------------------
void TTrackCompModule::FillEfficiencyHistograms(TStnTrackBlock*  TrackBlock, 
						TStnTrackID*     TrackID   , 
						TrackPar_t*      TPar      ,
						int              HistSet   ) {
  if (fSimPar.fParticle == NULL) return;

  if (fSimPar.fParticle->NStrawHits() >= 20) {
    FillEventHistograms(fHist.fEvent[HistSet]);

    if (fSimp->fMomTrackerFront > 100.) {
      FillEventHistograms(fHist.fEvent[HistSet+1]);

      TVector3 vdmom;
      vdmom.SetXYZ(fSimPar.fTFront->Mom()->X(),
		   fSimPar.fTFront->Mom()->Y(),		      
		   fSimPar.fTFront->Mom()->Z());

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
		
		if (((id_word & TStnTrackID::kD0Bit  ) == 0) && 
		    ((id_word & TStnTrackID::kRMaxBit) == 0)    ) {
		  
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
void TTrackCompModule::FillTrackHistograms(HistBase_t* HistR, TStnTrack* Track, TrackPar_t* Tp, double Weight) {

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


//--------------------------------------------------------------------------------
// function to fill TrasckSeedHit block
//--------------------------------------------------------------------------------
void TTrackCompModule::FillTrackSeedHistograms(HistBase_t*   HistR, TStnTrackSeed* TrkSeed) {
  
  TrackSeedHist_t* Hist = (TrackSeedHist_t*) HistR;
  
  int         nhits    = TrkSeed->NHits      ();
  double      clusterT = TrkSeed->ClusterTime();
  double      clusterE = TrkSeed->ClusterEnergy();
  
  double      mm2MeV   = 3/10.;
  double      pT       = TrkSeed->Pt();
  double      radius   = pT/mm2MeV;

  double      tanDip   = TrkSeed->TanDip();  
  double      p        = pT/std::cos( std::atan(tanDip));
  
  Hist->fNHits      ->Fill(nhits);	 
  Hist->fClusterTime->Fill(clusterT);
  Hist->fClusterEnergy->Fill(clusterE);
  
  Hist->fRadius     ->Fill(radius);    
  Hist->fMom        ->Fill(p);	 
  Hist->fPt         ->Fill(pT);	 
  Hist->fTanDip     ->Fill(tanDip);    
  
  Hist->fChi2       ->Fill(TrkSeed->Chi2()/(nhits-5.));
  Hist->fFitCons    ->Fill(TrkSeed->FitCons());
  Hist->fD0         ->Fill(TrkSeed->D0());
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
  if ((fEClMax > fMinETrig) && (fTClMax > 550)) {
    FillEventHistograms(fHist.fEvent[2]);
  }
//-----------------------------------------------------------------------------
// EVT_3: events with at least one helix
//-----------------------------------------------------------------------------
  if (fNHelices > 0) FillEventHistograms(fHist.fEvent[3]);
//-----------------------------------------------------------------------------
// EVT_4: events with at least one track seed
// EVT_5: events with at least one GOOD track seed
//-----------------------------------------------------------------------------
  if (fNTrackSeeds     > 0) FillEventHistograms(fHist.fEvent[4]);
  if (fNGoodTrackSeeds > 0) FillEventHistograms(fHist.fEvent[5]);
//-----------------------------------------------------------------------------
// EVT_41: events with good DAR positron track 
//-----------------------------------------------------------------------------
  if ((fNTracks[kDAR] > 0) && (fWeight > 0)) {
    TStnTrack*  trk = fTrackBlock[kDAR]->Track(0);
    TrackPar_t* tp  = &fTrackPar[kDAR][0];
    if ((trk->P0() > 85) && (tp->fIDWord_RMC == 0) && (trk->Charge() > 0)) {
      FillEventHistograms(fHist.fEvent[41]);
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
    FillEventHistograms(fHist.fEvent[51]);
    TrackPar_t* tp  = &fTrackPar[kPAR][0];
    if (tp->fIDWord_RMC == 0) FillEventHistograms(fHist.fEvent[52]);
  }
//-----------------------------------------------------------------------------
// EVT_61 : events with DAR tracks 
// EVT_62 : events with good DAR tracks 
//-----------------------------------------------------------------------------
  if (fNTracks[kDAR] > 0) {
    FillEventHistograms(fHist.fEvent[61]);
    TrackPar_t* tp  = &fTrackPar[kDAR][0];
    if (tp->fIDWord_RMC == 0) FillEventHistograms(fHist.fEvent[62]);
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
    FillTrackHistograms(fHist.fTrack[0],trk,tp);
  }
  else {
    if (fTrackBlock[kDAR]->NTracks() > 0) {
      if (fTrackPar[kDAR][0].fIDWord[3] == 0) {
//-----------------------------------------------------------------------------
// TRK_1: have KDAR track TrkQual > 0.3 and no KPAR TrkQual > 0.4 track
//-----------------------------------------------------------------------------
	TStnTrack* trk = fTrackBlock[kDAR]->Track(0);
	TrackPar_t* tp = &fTrackPar[kDAR][0];
	FillTrackHistograms(fHist.fTrack[1],trk,tp);
      }
      if (fTrackPar[1][0].fIDWord[2] == 0) {
//-----------------------------------------------------------------------------
// TRK_2: have KDAR track TrkQual > 0.2 and no KPAR TrkQual > 0.4 track
//-----------------------------------------------------------------------------
	TStnTrack* trk = fTrackBlock[1]->Track(0);
	TrackPar_t* tp = &fTrackPar[1][0];
	FillTrackHistograms(fHist.fTrack[2],trk,tp);
      }
    }
  }
//-----------------------------------------------------------------------------
// TRK_3 : an attempt to define the best track
//-----------------------------------------------------------------------------
  TStnTrack  *tpr(0), *cpr(0), *best_track(0);
  TrackPar_t *best_tp, *tprp, *cprp;

  int ntpr = fTrackBlock[0]->NTracks();
  int ncpr = fTrackBlock[1]->NTracks();

  if ((ntpr > 0) && (fTrackPar[0][0].fIDWord[fBestID[0]] == 0)) {
    tpr  = fTrackBlock[0]->Track(0);
    tprp = &fTrackPar[0][0];
  }

  if ((ncpr > 0) && (fTrackPar[1][0].fIDWord[fBestID[1]] == 0)) {
    cpr  = fTrackBlock[1]->Track(0);
    cprp = &fTrackPar[1][0];
  }

  if      (tpr != NULL) {
    if (cpr == NULL) {
//-----------------------------------------------------------------------------
// only KPAR track is present
//-----------------------------------------------------------------------------
      best_track = tpr;
      best_tp    = tprp;
    }
    else {
//-----------------------------------------------------------------------------
// general case: both tracks present and passsed the ID cuts, figure which one 
// is the best
//-----------------------------------------------------------------------------
      double tpr_prob = fTrackProb[0]->prob(tprp->fMVAOut[0]);
      double cpr_prob = fTrackProb[1]->prob(cprp->fMVAOut[0]);

      if (tpr_prob > cpr_prob) {
	best_track = tpr;
	best_tp    = tprp;
      }
      else {
	best_track = cpr;
	best_tp    = cprp;
      }
    }
  }
  else if (cpr != NULL) {
//-----------------------------------------------------------------------------
// only KDAR track is present
//-----------------------------------------------------------------------------
    best_track = cpr;
    best_tp    = cprp;
  }

  if (best_track != 0) {
    FillTrackHistograms(fHist.fTrack[3],best_track,best_tp);
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
      FillTrackHistograms(fHist.fTrack[ihist+0],trk,tp);
//-----------------------------------------------------------------------------
// TRK_101, TRK_201: BestID 
//-----------------------------------------------------------------------------
      if (tp->fIDWord[fBestID[i]] == 0) {
//-----------------------------------------------------------------------------
// cosmics
//-----------------------------------------------------------------------------
	FillTrackHistograms(fHist.fTrack[ihist+1],trk,tp);
	n_setc_tracks[i] += 1;
      }
//-----------------------------------------------------------------------------
// TRK_102, TRK_202: (BestID - FitConsBit - T0ErrBit - MomErrBit) tracks 
//-----------------------------------------------------------------------------
      int mask = (TStnTrackID::kFitConsBit | TStnTrackID::kT0ErrBit | TStnTrackID::kMomErrBit);
      if ((tp->fIDWord[fBestID[i]] & ~mask) == 0) {
	FillTrackHistograms(fHist.fTrack[ihist+2],trk,tp);
      }
//-----------------------------------------------------------------------------
// IHIST+3: (SetC + (dpf > 1)tracks 
//-----------------------------------------------------------------------------
      if ((tp->fIDWord[fBestID[i]] == 0) & (tp->fDpF > 1)) {
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
//-----------------------------------------------------------------------------
// IHIST+6: oddly, tracks with seemingly wrong reconstruced T0
//-----------------------------------------------------------------------------
      if (tp->fDtZ0 < -10) {
	FillTrackHistograms(fHist.fTrack[ihist+6],trk,tp);
      }
//-----------------------------------------------------------------------------
// IHIST+9: tracks with XDpF > 10
//-----------------------------------------------------------------------------
      if (tp->fXDpF   > 10) FillTrackHistograms(fHist.fTrack[ihist+9],trk,tp);
//-----------------------------------------------------------------------------
// different cuts on track quality variable
//-----------------------------------------------------------------------------
      for (int idd=0; idd<fNID; idd++) {
	if (tp->fIDWord[idd] == 0) FillTrackHistograms(fHist.fTrack[ihist+10+idd],trk,tp);
      }
//-----------------------------------------------------------------------------
// IHIST+41: N(a) > 20
// IHIST+42: N(a) > 20, |D0| < 100
// IHIST+43: N(a) > 20, |D0| < 100, DNa < 6
// IHIST+44: N(a) > 20, |D0| < 100, DNa < 6, chi2d < 4
//-----------------------------------------------------------------------------
      if (trk->NActive() > 20) { 
	FillTrackHistograms(fHist.fTrack[ihist+41],trk,tp);
	if (fabs(trk->D0()) < 100) {
	  FillTrackHistograms(fHist.fTrack[ihist+42],trk,tp);
	  int dna = trk->NHits()-trk->NActive();
	  if (dna < 6) {
	    FillTrackHistograms(fHist.fTrack[ihist+43],trk,tp);
	    if (trk->Chi2Dof() < 4.) {
	      FillTrackHistograms(fHist.fTrack[ihist+44],trk,tp);
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
	FillTrackHistograms(fHist.fTrack[ihist+80],trk,tp);
	if (trk->NActive() > 20) { 
	  FillTrackHistograms(fHist.fTrack[ihist+81],trk,tp);
	  if (fabs(trk->D0()) < 100) {
	    FillTrackHistograms(fHist.fTrack[ihist+82],trk,tp);
	    int dna = trk->NHits()-trk->NActive();
	    if (dna < 6) {
	      FillTrackHistograms(fHist.fTrack[ihist+83],trk,tp);
	      if (trk->Chi2Dof() < 4.) {
		FillTrackHistograms(fHist.fTrack[ihist+84],trk,tp);
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
	FillTrackHistograms(fHist.fTrack[ihist+90],trk,tp);
	if (trk->NActive() > 20) { 
	  FillTrackHistograms(fHist.fTrack[ihist+91],trk,tp);
	  if (fabs(trk->D0()) < 100) {
	    FillTrackHistograms(fHist.fTrack[ihist+92],trk,tp);
	    int dna = trk->NHits()-trk->NActive();
	    if (dna < 6) {
	      FillTrackHistograms(fHist.fTrack[ihist+93],trk,tp);
	      if (trk->Chi2Dof() < 4.) {
		FillTrackHistograms(fHist.fTrack[ihist+94],trk,tp);
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
	FillTrackHistograms(fHist.fTrack[ihist+50],trk,tp);
	if (trk->NActive() > 20) { 
	  FillTrackHistograms(fHist.fTrack[ihist+51],trk,tp);
	  if (fabs(trk->D0()) < 100) {
	    FillTrackHistograms(fHist.fTrack[ihist+52],trk,tp);
	    int dna = trk->NHits()-trk->NActive();
	    if (dna < 6) {
	      FillTrackHistograms(fHist.fTrack[ihist+53],trk,tp);
	      if (trk->Chi2Dof() < 4.) {
		FillTrackHistograms(fHist.fTrack[ihist+54],trk,tp);
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
	FillTrackHistograms(fHist.fTrack[ihist+60],trk,tp);
	if (trk->NActive() > 20) { 
	  FillTrackHistograms(fHist.fTrack[ihist+61],trk,tp);
	  if (fabs(trk->D0()) < 100) {
	    FillTrackHistograms(fHist.fTrack[ihist+62],trk,tp);
	    int dna = trk->NHits()-trk->NActive();
	    if (dna < 6) {
	      FillTrackHistograms(fHist.fTrack[ihist+63],trk,tp);
	      if (trk->Chi2Dof() < 4.) {
		FillTrackHistograms(fHist.fTrack[ihist+64],trk,tp);
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
	  FillTrackHistograms  (fHist.fTrack[ihist+71],trk,tp);
	  if (fProcess != -1) FillTrackHistograms(fHist.fTrack[ihist+73],trk,tp,fWeight); 
	  if ((tp->fP > 90.) && (tp->fP < 93.)) FillTrackHistograms(fHist.fTrack[ihist+77],trk,tp);            // for cosmics

	  if   (tp->fEp < 0.7) FillTrackHistograms(fHist.fTrack[ihist+85],trk,tp);   // "e+"
	  else                 FillTrackHistograms(fHist.fTrack[ihist+86],trk,tp);   // "mu+/pi+"
	}
	else {
//-----------------------------------------------------------------------------
// negative tracks
//-----------------------------------------------------------------------------
	  FillTrackHistograms(fHist.fTrack[ihist+72],trk,tp);
	  if (fProcess > 0) FillTrackHistograms(fHist.fTrack[ihist+74],trk,tp,fWeight);                       // weighting
	  
	  if ((tp->fP > 90.) && (tp->fP < 93.)) FillTrackHistograms(fHist.fTrack[ihist+78],trk,tp);            // for cosmics

	  if   (tp->fEp < 0.7) FillTrackHistograms(fHist.fTrack[ihist+95],trk,tp);   // "e-"
	  else                 FillTrackHistograms(fHist.fTrack[ihist+96],trk,tp);   // "mu-/pi-"

	  FillTrackHistograms(fHist.fTrack[ihist+79],trk,tp,tp->fDioWt);                                       // DIO

//-----------------------------------------------------------------------------
// debugging mu-
//-----------------------------------------------------------------------------
	  if (tp->fDrDzCal < 0 ) FillTrackHistograms(fHist.fTrack[ihist+200+ 1            ],trk,tp);
	  else                   FillTrackHistograms(fHist.fTrack[ihist+200+ 2            ],trk,tp);

	  if (tp->fEcl     > 0 ) FillTrackHistograms(fHist.fTrack[ihist+200+ 5+tp->fDiskID],trk,tp);

	  if (tp->fDt      < -4) FillTrackHistograms(fHist.fTrack[ihist+200+ 7            ],trk,tp);
	  else                   FillTrackHistograms(fHist.fTrack[ihist+200+ 8            ],trk,tp);
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
  FillEfficiencyHistograms(fTrackBlock[0],fTrackID[fBestID[0]],&fTrackPar[0][0],10);
  FillEfficiencyHistograms(fTrackBlock[1],fTrackID[fBestID[1]],&fTrackPar[1][0],20);

//-----------------------------------------------------------------------------
// fill little tree if requested
//-----------------------------------------------------------------------------
  if (fWriteTmvaTree >= 0) FillTmvaTree();

}

//-----------------------------------------------------------------------------
// assume less than 20 tracks 
//-----------------------------------------------------------------------------
int TTrackCompModule::InitTrackPar(TStnTrackBlock*   TrackBlock  , 
				   TStnClusterBlock* ClusterBlock, 
				   TrackPar_t*       TrackPar    ) {
  TrackPar_t*           tp;
  TStnTrack*            track;
  int                   track_type(-999);
  double                xs;
  TEmuLogLH::PidData_t  dat;
//-----------------------------------------------------------------------------
// momentum corrections for KPAR and KDAR
//-----------------------------------------------------------------------------
  const double kMomentumCorr[2] = { 0. , 0. } ; // CD3: { 0.049, 0.020 };
  const double kDtTcmCorr   [2] = { 0. , 0. } ; // { 0.22 , -0.30 }; // ns, sign: fit peak positions

  const char* block_name = TrackBlock->GetNode()->GetName();

  if      (strcmp(block_name,fTrackBlockName[0].Data()) == 0) track_type = 0;
  else if (strcmp(block_name,fTrackBlockName[1].Data()) == 0) track_type = 1;
  else {
    Error("TTrackCompModule::InitTrackPar",Form("IN TROUBLE: unknown track block: %s",block_name));
    return -1;
  }
//-----------------------------------------------------------------------------
// loop over tracks
//-----------------------------------------------------------------------------
  int ntrk = TrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    tp             = TrackPar+itrk;
    track          = TrackBlock->Track(itrk);
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
    if (track_type == 2) track_type = track->BestAlg();

    tp->fP     = track->fP0    +kMomentumCorr[track_type];		// correcting
//-----------------------------------------------------------------------------
// hits on virtual detectors
//-----------------------------------------------------------------------------
     tp->fPStOut = track->fPStOut;
     tp->fPFront = track->fPFront;

//     double t_stout(1.e6), t_front(1.e6);

//     int nvdhits = fVDetBlock->NHits();
//     for (int i=0; i<nvdhits; i++) {
//       TVDetHitData* vdh = fVDetBlock->Hit(i);
//       if (vdh->GeneratorCode() == fGeneratorCode) {
// 	if ((vdh->Index() == 10) && (vdh->Time() < t_stout)) {
// //-----------------------------------------------------------------------------
// // ST exit 
// //-----------------------------------------------------------------------------
// 	  tp->fPStOut = vdh->McMomentum();
// 	  t_stout     = vdh->Time();
// 	}
// 	if ((vdh->Index() == 13) && (vdh->Time() < t_front)) {
// //-----------------------------------------------------------------------------
// // tracker front
// //-----------------------------------------------------------------------------
// 	  tp->fPFront = vdh->McMomentum();
// 	  t_front     = vdh->Time();
// 	}
//       }
//    }
    
    tp->fDpF     = tp->fP     -tp->fPFront;
    tp->fDp0     = track->fP0 -tp->fPFront;
    tp->fDp2     = track->fP2 -tp->fPFront;
    tp->fDpFSt   = tp->fPFront-tp->fPStOut;

    tp->fXDpF    = tp->fDpF/track->fFitMomErr;
    tp->fDioWt   = TStntuple::DioWeightAl   (fEleE)*(105./10);
    tp->fDioWtRC = TStntuple::DioWeightAl_LL(fEleE)*(105./10);
    tp->fLumWt   = GetHeaderBlock()->LumWeight();

    tp->fDtZ0 = -1.e6;
    if (fSimPar.fTMid) {
      double ttrue = fmod(fSimPar.fTMid->Time(),fMbTime);
      tp->fDtZ0 = track->T0()-ttrue;
    }

    tp->fXtZ0 = tp->fDtZ0/track->fT0Err;

    tp->fDtBack = -1.e6;
    if (fSimPar.fTBack) {
      double ttrue = fmod(fSimPar.fTBack->Time(),fMbTime);
      tp->fDtBack = track->T0()-ttrue;
    }
//-----------------------------------------------------------------------------
// track residuals
//-----------------------------------------------------------------------------
    TStnTrack::InterData_t*  vr = track->fVMaxEp; 
    double    nx, ny;

    tp->fEcl       = -1.e6;
    tp->fDiskID    = -1;
    tp->fEp        = -1.e6;
    tp->fDrDzCal   = -1.e6;
    tp->fDtClZ0    = -1.e6;

    tp->fDu        = -1.e6;
    tp->fDv        = -1.e6;
    tp->fDx        = -1.e6;
    tp->fDy        = -1.e6;
    tp->fDz        = -1.e6;
    tp->fDt        = -1.e6;

    tp->fChi2Tcm = -1.e6;
    tp->fChi2XY    = -1.e6;
    tp->fChi2T     = -1.e6;
    tp->fPath      = -1.e6;
    tp->fSinTC     = -1.e6;
    tp->fDrTC      = -1.e6;
    tp->fSInt      = -1.e6;

    if (vr) {
      tp->fDiskID = vr->fID;
      tp->fEcl    = vr->fEnergy;
      tp->fEp     = tp->fEcl/track->fP2;
      tp->fDrDzCal = (vr->fXTrk*vr->fNxTrk+vr->fYTrk+vr->fNyTrk)/sqrt(vr->fXTrk*vr->fXTrk+vr->fYTrk*vr->fYTrk)/vr->fNzTrk;

      tp->fDx     = vr->fDx;
      tp->fDy     = vr->fDy;
      tp->fDz     = vr->fDz;
//-----------------------------------------------------------------------------
// v4_2_4: correct by additional 0.22 ns - track propagation by 6 cm
//-----------------------------------------------------------------------------
      tp->fDt  = vr->fDt ; // v4_2_4: - 0.22; // - 1.;
      if (track_type >= 0) tp->fDt  -= kDtTcmCorr[track_type];

      nx  = vr->fNxTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);
      ny  = vr->fNyTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);

      tp->fDu        = vr->fDx*nx+vr->fDy*ny;
      tp->fDv        = vr->fDx*ny-vr->fDy*nx;
      tp->fChi2Tcm   = vr->fChi2Match;
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

      if (fSimPar.fTMid) {
	tp->fDtClZ0 = tp->fDt-tp->fDtZ0;
      }
    }

    if ((tp->fEp > 0) && (track->fEp > 0) && (fabs(tp->fEp-track->fEp) > 1.e-6)) {
      GetHeaderBlock()->Print(Form(" TTrackAnaModule ERROR: tp->fEp = %10.5f  track->fEp = %10.5f\n ",tp->fEp,track->fEp));
    }
//-----------------------------------------------------------------------------
// on-the-fly MVA calculation
//-----------------------------------------------------------------------------
    tp->fMVAOut[0] = track->DaveTrkQual();

    if (fUseMVA != 0) {
      vector<double>  pmva;
      pmva.resize(10);

      float na = track->NActive();

      pmva[ 0] = na;
      pmva[ 1] = na/track->NHits();
      pmva[ 2] = -1.e6; 			// defined below
      pmva[ 3] = track->FitMomErr();
      pmva[ 4] = track->T0Err();
      pmva[ 5] = track->D0();
      pmva[ 6] = track->RMax();
      pmva[ 7] = track->NDoubletsAct()/na;
      pmva[ 8] = track->NHitsAmbZero()/na;
      pmva[ 9] = track->NMatActive()/na;
      //      pmva[10] = track->fZ1;

      int alg = track->BestAlg();

      if      (alg == 0) {
//-----------------------------------------------------------------------------
// KPAR track - log(fitcons) used for training
//-----------------------------------------------------------------------------
	if (fTmvaAlgorithmTpr < 0) {
	  tp->fMVAOut[0] = track->DaveTrkQual();
	}
	else {
	  int use_chi2d = fTmvaAlgorithmTpr / 100;
	  if      (use_chi2d == 0) pmva[2] = log10(track->FitCons());
	  else                     pmva[2] = track->Chi2Dof();

	  tp->fMVAOut[0] = fTprQualMva->evalMVA(pmva);
	}
	tp->fMVAOut[1] = fTprQualMva->evalMVA(pmva);

	tp->fProb = fTrackProb[0]->prob(tp->fMVAOut[0]);
      }
      else if (alg == 1) {
//-----------------------------------------------------------------------------
// KDAR track - for fTmvaAlgorithm=101 chi2/N(dof) used (our initial training)
//-----------------------------------------------------------------------------
	int use_chi2d = fTmvaAlgorithmCpr / 100;
	if      (use_chi2d == 0) pmva[2] = log10(track->FitCons());
	else                     pmva[2] = track->Chi2Dof();

	tp->fMVAOut[0] = fCprQualMva->evalMVA(pmva);
	tp->fMVAOut[1] = tp->fMVAOut[0];

	tp->fProb      = fTrackProb[1]->prob(tp->fMVAOut[0]);
      }

      if (GetDebugBit(9)) {
	if (alg == 0) {
	  GetHeaderBlock()->Print(Form("KPAR TrkQual, MVAOut: %10.5f %10.5f",
				       track->DaveTrkQual(),tp->fMVAOut[0]));
	}
      }
      if (GetDebugBit(10)) {
	if ((alg == 0) && (tp->fMVAOut[0] - track->DaveTrkQual() > 0.3)) {
	  GetHeaderBlock()->Print(Form("KPAR TrkQual, MVAOut: %10.5f %10.5f",
				       track->DaveTrkQual(),tp->fMVAOut[0]));
	}
      }
    }
//-----------------------------------------------------------------------------
// finally, the track ID
//-----------------------------------------------------------------------------
    for (int idd=0; idd<fNID; idd++) {
      int idw = fTrackID[idd]->IDWord(track);
//-----------------------------------------------------------------------------
// redefine IDWord to use TQ ANN for KPAR and CQ ANN for KDAR tracks
//-----------------------------------------------------------------------------
      idw &= (~TStnTrackID::kTrkQualBit);
      if (tp->fMVAOut[0] < fTrackID[idd]->MinTrkQual()) idw |= TStnTrackID::kTrkQualBit;

      tp->fIDWord[idd] = idw;
    }
//-----------------------------------------------------------------------------
// RMC ID
//-----------------------------------------------------------------------------
    tp->fIDWord_RMC = fTrackID_RMC->IDWord(track);
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

    int id_word = tp->fIDWord[fBestID[track_type]];

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
  // fWtRMC   =  1.;
  // fWtRPC   =  1.;
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
// may want to revisit the definition of fSimp in future
//-----------------------------------------------------------------------------
  fSimp             = fSimpBlock->Particle(0);
  fSimPar.fParticle = fSimp;
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  fSimPar.fGenp     = NULL;
//-----------------------------------------------------------------------------
// sometimes, there is no GENP info - as in Bob's pbar dataset
//-----------------------------------------------------------------------------
  if ((fProcess == -1) && (fNSimp > 0)) {
    if (fSimp->PDGCode() == -2212) {
//-----------------------------------------------------------------------------
// pbar production, assume Bob's dataset
// assuming parent particle exists, determine the production cross-section weight
// pbeam, nx, ny: the beam momentum and direction
//-----------------------------------------------------------------------------
      double pbeam(8.9), nx(-0.24192190), nz(-0.97029573);
      const TLorentzVector* sm = fSimp->StartMom();
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

  if (fParticle) fEleE = fParticle->Energy();
  else           fEleE = -1.;
//-----------------------------------------------------------------------------
// virtual detectors - for fSimp need parameters at the tracker front
//-----------------------------------------------------------------------------
  int nsteps = fSpmcBlockVDet->NStepPoints();
  
  for (int i=0; i<nsteps; i++) {
    TStepPointMC* step = fSpmcBlockVDet->StepPointMC(i);
    if (step->PDGCode() == fSimp->fPdgCode) {
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
    InitTrackPar(fTrackBlock[i],fClusterBlock,fTrackPar[i]);
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
//-----------------------------------------------------------------------------
// bit 14: rare misreconstructed events from e11s721z
//-----------------------------------------------------------------------------
  if (GetDebugBit(14) == 1) {
    for (int itrk=0; itrk<ntrk; itrk++) {
      trk = cprb->Track(itrk);
      tp  = &fTrackPar[kDAR][itrk];
      if ((tp->fIDWord_RMC == 0) && (trk->Charge() < 0)) {
	if (((tp->fP > 103.5) && (tp->fDioWt > 9.e-11)) || 
	    ((tp->fP > 104  ) && (tp->fDioWt > 1.e-12)) ||
	    ((tp->fP > 106  ) && (tp->fDioWt > 1.e-14))    ) {
	  GetHeaderBlock()->Print(Form("TTrackCompModule :bit014: tp->fP = %10.3f tp->fDioWt = %12.3e", 
				       tp->fP,tp->fDioWt));
	}
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


//_____________________________________________________________________________
void TTrackCompModule::PrintTrack(TStnTrack* Track, TrackPar_t* Tp, Option_t* Option) const {

  TString opt(Option);

  opt.ToLower();
					// "non-const *this" for printing purposes
  TStnTrack* t = (TStnTrack*) Track;

  if ((opt == "") || (opt.Index("banner") >= 0)) {
//-----------------------------------------------------------------------------
// print banner
//-----------------------------------------------------------------------------
    printf("------------------------------------------------------------------------------------------------");
    printf("----------------------------------------------------------------------------------\n");
    printf(" i  nh  na nw nosd nssd na0 ncl  alg_mask    id_word   q     p     p(corr) momerr    T0     T0Err     D0");
    printf("      Z0    TanDip   TBack   chi2/dof   fcon  TrkQual MvaOut[0]  MVAOut[1]   Prob \n");
    printf("------------------------------------------------------------------------------------------------");
    printf("----------------------------------------------------------------------------------\n");
  }

  if ((opt == "") || (opt.Index("data") >= 0)) {
    printf("%2i %3i %3i %2i %4i %4i %3i %3i 0x%08x",
	   t->fNumber,t->NHits(), t->NActive(),t->NWrong(), 
	   t->NOSDoublets(), t->NSSDoublets(), t->NHitsAmbZero(),
	   t->NClusters(),
	   t->AlgorithmID());

    printf(" 0x%08x %1.0f %8.3f %8.3f %7.3f %8.3f %6.3f %7.3f %8.3f %7.4f %8.3f %8.2f %8.2e %7.3f %9.4f %9.4f %9.4f",
	   t->fIDWord,
	   t->fCharge, 
	   t->fP*t->fCharge, Tp->fP*t->fCharge, t->fFitMomErr, t->fT0, t->fT0Err, t->fD0, t->fZ0, t->fTanDip, t->TBack(),
	   t->Chi2Dof(),t->FitCons(),t->DaveTrkQual(),Tp->fMVAOut[0], Tp->fMVAOut[1], Tp->fProb);
    printf("\n");
  }
}
