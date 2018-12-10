//////////////////////////////////////////////////////////////////////////////
// CalPatRec validation: compare CalPatRec tracks to TrkPatRec ones
// assume running on stntuple with 3 tracking branches
//
// use of tmp:
//
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//  3  : events with CalPatRec tracks with DPF > 5
//  4  : events with TrkPatRec track and a 50+ MeV cluster, but with no CalPatRec track
//  5  : events with CalPatRec tracks with P > 106
//  6  : events with CalPatRec tracks with 1.5 < DPF < 5
//  7  : events with N(set "C" CalPatRec tracks)  > 0
//  8  : events with CalPatRec tracks with P > 105
//  9  : all events - print DMVA
// 10  : events with TRKPATREC DMVA > 0.3 - print DMVA
// 11  : validate MergePatRec
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
#include<string>

using std::string;
using std::vector;

ClassImp(TTrackCompModule)
//-----------------------------------------------------------------------------
TTrackCompModule::TTrackCompModule(const char* name, const char* title):
  TStnModule    (name,title),
  fPdgCode      (11),                       // electron
  fGeneratorCode( 2)                        // 2:ConversionGun 28:StoppedParticleReactionGun
{
  fFillDioHist = 1;

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

  fMinETrig       = 50.;
  fNMVA           = 0;

  fUseMVA         = 0;

  fTprMVA         = new mva_data("trkpatrec","dave"    ,002);
  fCprMVA         = new mva_data("calpatrec","e11s5731",002);

  fLogLH       = new TEmuLogLH();
//-----------------------------------------------------------------------------
// debugging information
//-----------------------------------------------------------------------------
  fDebugCut[5].fXMin   = 106.;
  fDebugCut[5].fXMax   = 200.;

  fDebugCut[6].fXMin   = 1.5;
  fDebugCut[6].fXMax   = 10.0;
//-----------------------------------------------------------------------------
// ntuples for TMVA training
//-----------------------------------------------------------------------------
  fWriteTmvaTree     = -1;
  fTmvaAlgorithmTpr  = -1;
  fTmvaAlgorithmCpr  = -1;
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
  RegisterDataBlock(fTrackBlockName[0].Data(), "TStnTrackBlock"     ,&fTrackBlock[0]);
  RegisterDataBlock(fTrackBlockName[1].Data(), "TStnTrackBlock"     ,&fTrackBlock[1]);
  RegisterDataBlock("ClusterBlock"           , "TStnClusterBlock"   ,&fClusterBlock );
  RegisterDataBlock("SimpBlock"              , "TSimpBlock"         ,&fSimpBlock    );
  RegisterDataBlock("GenpBlock"              , "TGenpBlock"         ,&fGenpBlock    );
  RegisterDataBlock("VDetBlock"              , "TVDetDataBlock"     ,&fVDetBlock    );
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
// init two MVA-based classifiers - TrkPatRec (TPR) and CalPatRec (CPR) 
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
  fBestID[1]      = fCprMVA->BestID();		  // CalPatRec     : CprQual     > 0.85

  return 0;
}


//_____________________________________________________________________________
int TTrackCompModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
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
  HBook1F(Hist->fChi2Dof    ,"chi2d"    ,Form("%s: track chi2/N(dof)" ,Folder), 500, 0, 10,Folder);

  HBook1F(Hist->fNActive    ,"nactv"    ,Form("%s: N(active)"         ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fNaFract    ,"nafr"     ,Form("%s: N(active fraction)",Folder), 110, 0.5,1.05,Folder);
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
  HBook1F(Hist->fDtZ0       ,"dtz0"     ,Form("%s: DT(Z0), MC"        ,Folder), 200, -10.0 ,10.0,Folder);

  HBook1F(Hist->fResid      ,"resid"    ,Form("%s: hit residuals"     ,Folder), 500,-0.5 ,0.5,Folder);
  HBook1F(Hist->fAlgMask    ,"alg"      ,Form("%s: algorithm mask"    ,Folder),  10,  0, 10,Folder);

  HBook1F(Hist->fChi2Tcm  ,"chi2tcm"  ,Form("%s: chi2(t-c match)"   ,Folder), 250,  0  ,250 ,Folder);
  HBook1F(Hist->fChi2XY     ,"chi2xy"   ,Form("%s: chi2(t-c match) XY",Folder), 300,-50  ,250 ,Folder);
  HBook1F(Hist->fChi2T      ,"chi2t"    ,Form("%s: chi2(t-c match) T" ,Folder), 250,  0  ,250 ,Folder);

  HBook1F(Hist->fDt         ,"dt"       ,Form("%s: track delta(T)"    ,Folder), 400,-20  ,20 ,Folder);
  HBook1F(Hist->fDx         ,"dx"       ,Form("%s: track delta(X)"    ,Folder), 200,-500 ,500,Folder);
  HBook1F(Hist->fDy         ,"dy"       ,Form("%s: track delta(Y)"    ,Folder), 200,-500 ,500,Folder);
  HBook1F(Hist->fDz         ,"dz"       ,Form("%s: track delta(Z)"    ,Folder), 200,-250 ,250,Folder);
  HBook1F(Hist->fDu         ,"du"       ,Form("%s: track-cluster DU)" ,Folder), 250,-250 ,250,Folder);
  HBook1F(Hist->fDv         ,"dv"       ,Form("%s: track-cluster DV)" ,Folder), 200,-100 ,100,Folder);
  HBook1F(Hist->fPath       ,"path"     ,Form("%s: track sdisk"       ,Folder),  50,   0 ,500,Folder);

  HBook1F(Hist->fECl        ,"ecl"      ,Form("%s: cluster E"         ,Folder), 300, 0   ,150,Folder);
  HBook1F(Hist->fEClEKin    ,"ecl_ekin" ,Form("%s: cluster E/Ekin(mu)",Folder), 200, 0   ,2,Folder);
  HBook1F(Hist->fEp         ,"ep"       ,Form("%s: track E/P"         ,Folder), 300, 0   ,1.5,Folder);

  HBook2F(Hist->fFConsVsNActive,"fc_vs_na" ,Form("%s: FitCons vs NActive",Folder),  150, 0, 150, 200,0,1,Folder);
  HBook1F(Hist->fDaveTrkQual,"dtqual"   ,Form("%s:DaveTrkQual"        ,Folder), 200, -0.5, 1.5,Folder);
  HBook1F(Hist->fMVAOut     ,"mvaout"   ,Form("%s:MVAOut[0]"          ,Folder), 200, -0.5, 1.5,Folder);
  HBook1F(Hist->fDeltaMVA   ,"dmva"     ,Form("%s:MVAOut[0]-TrkQual"  ,Folder), 200, -1.0, 1.0,Folder);
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
  HBook1F(Hist->fInstLumi  ,"inst_lum" ,Form("%s: Inst Luminosity"                 ,Folder),500, 0,1.e8,Folder);
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
  book_event_histset[ 1] = 1;		// events with EclMax > fMinETrig and TClMax > 550

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
  book_track_histset[114] = 1;          // TrkPatRec TrackID[4]
  book_track_histset[115] = 1;          // TrkPatRec TrackID[5]
  book_track_histset[116] = 1;          // TrkPatRec TrackID[6]
  book_track_histset[117] = 1;          // TrkPatRec TrackID[7]
  book_track_histset[118] = 1;          // TrkPatRec TrackID[8]
  book_track_histset[119] = 1;          // TrkPatRec TrackID[9]
  book_track_histset[120] = 1;          // TrkPatRec TrackID[10]
  book_track_histset[121] = 1;          // TrkPatRec TrackID[11]
  book_track_histset[122] = 1;          // TrkPatRec TrackID[12]
  book_track_histset[123] = 1;          // TrkPatRec TrackID[13]
  book_track_histset[124] = 1;          // TrkPatRec TrackID[14]
  book_track_histset[125] = 1;          // TrkPatRec TrackID[15]
  book_track_histset[126] = 1;          // TrkPatRec TrackID[16]
  book_track_histset[127] = 1;          // TrkPatRec TrackID[17]
  book_track_histset[128] = 1;          // TrkPatRec TrackID[18]
  book_track_histset[129] = 1;          // TrkPatRec TrackID[19]

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
  book_track_histset[214] = 1;          // CalPatRec TrackID[4]
  book_track_histset[215] = 1;          // CalPatRec TrackID[5]
  book_track_histset[216] = 1;          // CalPatRec TrackID[6]
  book_track_histset[217] = 1;          // CalPatRec TrackID[7]
  book_track_histset[218] = 1;          // CalPatRec TrackID[8]
  book_track_histset[219] = 1;          // CalPatRec TrackID[9]
  book_track_histset[220] = 1;          // CalPatRec TrackID[10]
  book_track_histset[221] = 1;          // CalPatRec TrackID[11]
  book_track_histset[222] = 1;          // CalPatRec TrackID[12]
  book_track_histset[223] = 1;          // CalPatRec TrackID[13]
  book_track_histset[224] = 1;          // CalPatRec TrackID[14]
  book_track_histset[225] = 1;          // CalPatRec TrackID[15]
  book_track_histset[226] = 1;          // CalPatRec TrackID[16]
  book_track_histset[227] = 1;          // CalPatRec TrackID[17]
  book_track_histset[228] = 1;          // CalPatRec TrackID[18]
  book_track_histset[229] = 1;          // CalPatRec TrackID[19]

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
  book_track_histset[314] = 1;          // TrkPatRec not CalPatRec TrackID[4]
  book_track_histset[315] = 1;          // TrkPatRec not CalPatRec TrackID[5]
  book_track_histset[316] = 1;          // TrkPatRec not CalPatRec TrackID[6]
  book_track_histset[317] = 1;          // TrkPatRec not CalPatRec TrackID[7]
  book_track_histset[318] = 1;          // TrkPatRec not CalPatRec TrackID[8]
  book_track_histset[319] = 1;          // TrkPatRec not CalPatRec TrackID[9]
  book_track_histset[320] = 1;          // TrkPatRec not CalPatRec TrackID[10]
  book_track_histset[321] = 1;          // TrkPatRec not CalPatRec TrackID[11]
  book_track_histset[322] = 1;          // TrkPatRec not CalPatRec TrackID[12]
  book_track_histset[323] = 1;          // TrkPatRec not CalPatRec TrackID[13]
  book_track_histset[324] = 1;          // TrkPatRec not CalPatRec TrackID[14]
  book_track_histset[325] = 1;          // TrkPatRec not CalPatRec TrackID[15]
  book_track_histset[326] = 1;          // TrkPatRec not CalPatRec TrackID[16]
  book_track_histset[327] = 1;          // TrkPatRec not CalPatRec TrackID[17]
  book_track_histset[328] = 1;          // TrkPatRec not CalPatRec TrackID[18]
  book_track_histset[329] = 1;          // TrkPatRec not CalPatRec TrackID[19]

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
  book_track_histset[414] = 1;          // CalPatRec not TrkPatRec TrackID[4]
  book_track_histset[415] = 1;          // CalPatRec not TrkPatRec TrackID[5]
  book_track_histset[416] = 1;          // CalPatRec not TrkPatRec TrackID[6]
  book_track_histset[417] = 1;          // CalPatRec not TrkPatRec TrackID[7]
  book_track_histset[418] = 1;          // CalPatRec not TrkPatRec TrackID[8]
  book_track_histset[419] = 1;          // CalPatRec not TrkPatRec TrackID[9]
  book_track_histset[420] = 1;          // TrkPatRec not CalPatRec TrackID[10]
  book_track_histset[421] = 1;          // CalPatRec not TrkPatRec TrackID[11]
  book_track_histset[422] = 1;          // CalPatRec not TrkPatRec TrackID[12]
  book_track_histset[423] = 1;          // CalPatRec not TrkPatRec TrackID[13]
  book_track_histset[424] = 1;          // CalPatRec not TrkPatRec TrackID[14]
  book_track_histset[425] = 1;          // CalPatRec not TrkPatRec TrackID[15]
  book_track_histset[426] = 1;          // CalPatRec not TrkPatRec TrackID[16]
  book_track_histset[427] = 1;          // CalPatRec not TrkPatRec TrackID[17]
  book_track_histset[428] = 1;          // CalPatRec not TrkPatRec TrackID[18]
  book_track_histset[429] = 1;          // CalPatRec not TrkPatRec TrackID[19]

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
  if (fSimPar.fParticle == NULL) return;

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
void TTrackCompModule::FillTrackHistograms(HistBase_t* HistR, TStnTrack* Track, TrackPar_t* Tp) {

  TLorentzVector  mom;

  TrackHist_t* Hist = (TrackHist_t*) HistR;

					// pointer to local track parameters
  //  int itrk = Track->Number();

  //  TrackPar_t* tp = fTrackPar+itrk;

					// Tp->fP - corrected momentum, fP0 and fP2 - not corrected
  Hist->fP[0]->Fill (Tp->fP);
  Hist->fP[1]->Fill (Tp->fP);
  Hist->fP[2]->Fill (Tp->fP);
					// fP0: momentum in the first point,  fP2 - in the last
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
  Hist->fChi2Dof->Fill(Track->fChi2/(Track->NActive()-5.));

  float na = Track->NActive();

  Hist->fNActive->Fill(na);
  Hist->fNaFract->Fill(na/(Track->NHits()+0.));
  Hist->fNWrong->Fill(Track->NWrong());

  float nd = Track->NDoublets();

  float nad = Track->NDoubletsAct();
  Hist->fNDoublets->Fill(nd);
  Hist->fNadOverNd->Fill(nad/nd);
  Hist->fNOSD->Fill(Track->NOSDoublets());
  Hist->fNSSD->Fill(Track->NSSDoublets());
  Hist->fNdOverNa->Fill(nd/na);
  Hist->fNosdOverNa->Fill(Track->NOSDoublets()/na);
  Hist->fNssdOverNa->Fill(Track->NSSDoublets()/na);
  Hist->fNZeroAmb->Fill(Track->NHitsAmbZero());
  Hist->fNzaOverNa->Fill(Track->NHitsAmbZero()/na);

  int nma = Track->NMatActive();

  Hist->fNMatActive->Fill(nma);
  Hist->fNmaOverNa->Fill(nma/na);

  Hist->fNBend->Fill(Track->NBend());

  Hist->fT0->Fill(Track->fT0);
  Hist->fT0Err->Fill(Track->fT0Err);
  Hist->fQ->Fill(-1);
  Hist->fFitCons[0]->Fill(Track->fFitCons);
  Hist->fFitCons[1]->Fill(Track->fFitCons);

  Hist->fD0->Fill(Track->fD0);
  Hist->fZ0->Fill(Track->fZ0);
  Hist->fTanDip->Fill(Track->fTanDip);
  Hist->fDtZ0->Fill(Tp->fDtZ0);
  Hist->fRMax->Fill(Track->RMax());
  
  Hist->fAlgMask->Fill(Track->AlgMask());

  Hist->fChi2Tcm->Fill(Tp->fChi2Tcm);
  Hist->fChi2XY->Fill(Tp->fChi2XY);
  Hist->fChi2T->Fill (Tp->fChi2T);

  Hist->fDt->Fill(Tp->fDt);
  Hist->fDx->Fill(Tp->fDx);
  Hist->fDy->Fill(Tp->fDy);
  Hist->fDz->Fill(Tp->fDz);
  Hist->fDu->Fill(Tp->fDu);
  Hist->fDv->Fill(Tp->fDv);
  Hist->fPath->Fill(Tp->fPath);

  double    ekin(-1.);
  if (fSimp) {
    double p, m;
    p    = Tp->fP;
    m    = 105.658; // in MeV
    ekin = sqrt(p*p+m*m)-m;
  }

  Hist->fECl->Fill(Tp->fEcl);
  Hist->fEClEKin->Fill(Tp->fEcl/ekin);
  Hist->fEp->Fill(Tp->fEp);

  Hist->fFConsVsNActive->Fill(Track->NActive(),Track->fFitCons);
//-----------------------------------------------------------------------------
// MVA variables
//-----------------------------------------------------------------------------
  Hist->fDaveTrkQual->Fill(Track->DaveTrkQual());

  Hist->fMVAOut->Fill(Tp->fMVAOut[0]);
  float dmva=Tp->fMVAOut[0]-Track->DaveTrkQual();
  Hist->fDeltaMVA->Fill(dmva);
}


//-----------------------------------------------------------------------------
// TmvaAlgorithm: 100*usez + algorigthm
//-----------------------------------------------------------------------------
int TTrackCompModule::FillTmvaTree() {
  int rc(0), loc(-1);

  loc = fWriteTmvaTree;			// 0:trkpatrec , 1:calpatrec

  if (fTrackBlock[loc]->NTracks() != 1) return rc;

					// first CalPatRec track

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
    FillEventHistograms(fHist.fEvent[1]);
  }

//-----------------------------------------------------------------------------
// what does CalPatRec add ?
// TRK_0 : TrkPatRec tracks, BEST_ID
//-----------------------------------------------------------------------------
  if ((fTrackBlock[0]->NTracks() > 0) && (fTrackPar[0][0].fIDWord[fBestID[0]] == 0)) {
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
// only TrkPatRec track is present
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
// only CalPatRec track is present
//-----------------------------------------------------------------------------
    best_track = cpr;
    best_tp    = cprp;
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
// TRK_100, TRK_200: all reconstructed tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[ihist+0],trk,tp);
//-----------------------------------------------------------------------------
// TRK_101, TRK_201: BestID 
//-----------------------------------------------------------------------------
      if (tp->fIDWord[fBestID[i]] == 0) {
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
	if (tp->fIDWord[fBestID[i]] == 0) {
	  FillTrackHistograms(fHist.fTrack[ihist+1],trk,tp);
	  n_setc_tracks[i] += 1;
	}
	//-----------------------------------------------------------------------------
	// IHIST+2: (SetC - FitConsBit - T0ErrBit - MomErrBit) tracks 
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

	for (int idd=0; idd<fNID; idd++) {
	  if (tp->fIDWord[idd] == 0) FillTrackHistograms(fHist.fTrack[ihist+10+idd],trk,tp);
	}
      }

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
// momentum corrections for TrkPatRec and CalPatRec
//-----------------------------------------------------------------------------
  const double kMomentumCorr[2] = { 0.049, 0.020 };
  const double kDtTcmCorr   [2] = { 0.22 , -0.30 }; // ns, sign: fit peak positions

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

    tp->fChi2Tcm = -1.e6;
    tp->fChi2XY    = -1.e6;
    tp->fChi2T     = -1.e6;
    tp->fPath      = -1.e6;
    tp->fSinTC     = -1.e6;
    tp->fDrTC      = -1.e6;
    tp->fSInt      = -1.e6;

    if (vr) {
      tp->fEcl = vr->fEnergy;
      tp->fEp  = tp->fEcl/track->fP2;

      tp->fDx  = vr->fDx;
      tp->fDy  = vr->fDy;
      tp->fDz  = vr->fDz;
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
// TrkPatRec track - log(fitcons) used for training
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
// CalPatRec track - for fTmvaAlgorithm=101 chi2/N(dof) used (our initial training)
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
	  GetHeaderBlock()->Print(Form("TRKPATREC TrkQual, MVAOut: %10.5f %10.5f",
				       track->DaveTrkQual(),tp->fMVAOut[0]));
	}
      }
      if (GetDebugBit(10)) {
	if ((alg == 0) && (tp->fMVAOut[0] - track->DaveTrkQual() > 0.3)) {
	  GetHeaderBlock()->Print(Form("TRKPATREC TrkQual, MVAOut: %10.5f %10.5f",
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
// redefine IDWord to use TQ ANN for TrkPatRec and CQ ANN for CalPatRec tracks
//-----------------------------------------------------------------------------
      idw &= (~TStnTrackID::kTrkQualBit);
      if (tp->fMVAOut[0] < fTrackID[idd]->MinTrkQual()) idw |= TStnTrackID::kTrkQualBit;

      tp->fIDWord[idd] = idw;
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
  fVDetBlock->GetEntry(ientry);
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
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  fSimPar.fGenp     = NULL;
//-----------------------------------------------------------------------------
// virtual detectors - for fSimp need parameters at the tracker front
//-----------------------------------------------------------------------------
  int nvdhits = fVDetBlock->NHits();
  for (int i=0; i<nvdhits; i++) {
    TVDetHitData* vdhit = fVDetBlock->Hit(i);
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
  int         trkpatrec(0), calpatrec(1);
//-----------------------------------------------------------------------------
// bit 0: All Events
//-----------------------------------------------------------------------------
  if (GetDebugBit(0) == 1) {
    GetHeaderBlock()->Print(Form("TTrackCompModule :bit000:"));
    printf("TrkPatRec:\n");
    fTrackBlock[trkpatrec]->Print();

    for (int i=0; i<fNTracks[trkpatrec]; i++) { 
      PrintTrack(fTrackBlock[trkpatrec]->Track(i),&fTrackPar[trkpatrec][i],"data");
    }

    printf("CalPatRec:\n");
    fTrackBlock[calpatrec]->Print();

    for (int i=0; i<fNTracks[calpatrec]; i++) { 
      PrintTrack(fTrackBlock[calpatrec]->Track(i),&fTrackPar[calpatrec][i],"data");
    }
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
      if (fEClMax > fMinETrig) {
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
    tp  = &fTrackPar[calpatrec][0];

    int best_id = fBestID[calpatrec];

    if ((tp->fIDWord[best_id] == 0) && (tp->fDpF >= fDebugCut[6].fXMin) && (tp->fDpF < fDebugCut[6].fXMax)) {
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


  if (GetDebugBit(11) == 1) {
    if ((fNTracks[0] == 1) && (fNTracks[1] == 1)) {
      GetHeaderBlock()->Print(Form("TTrackCompModule :bit011"));
      PrintTrack(fTrackBlock[0]->Track(0),&fTrackPar[0][0],"");
      PrintTrack(fTrackBlock[1]->Track(0),&fTrackPar[1][0],"data");
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
