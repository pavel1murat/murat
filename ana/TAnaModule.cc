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

#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"

#include "Stntuple/base/TStnDataset.hh"
#include "Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/loop/TStnAna.hh"

#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
#include "murat/ana/TAnaModule.hh"

ClassImp(murat::TAnaModule)

namespace murat {

//-----------------------------------------------------------------------------
TAnaModule::TAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fMbTime      = 1695.;
  fMinT0       = 0.;                   // analysis cuts
  fApplyCorr   = 1;                    // by default, corfect momentum and delta(T)     
  fEventWeight = 1.;
  fBatchMode   = 1;
//-----------------------------------------------------------------------------
// track quality : box cuts
//-----------------------------------------------------------------------------
  fTrackID_BOX = new TStnTrackID();

  fTrackID_BOX->SetMaxChi2Dof(3.5 );
  fTrackID_BOX->SetMaxT0Err  (1.2);
  fTrackID_BOX->SetMaxMomErr (0.27);
  fTrackID_BOX->SetMinNActive(18 );
  fTrackID_BOX->SetMaxDNa    ( 6 );
  fTrackID_BOX->SetMinFNa    (0.87);
  fTrackID_BOX->SetMinTanDip (0.5);     // (1./sqrt(3.));
  fTrackID_BOX->SetMaxTanDip (1.0);
  fTrackID_BOX->SetMinD0     (-100.);
  fTrackID_BOX->SetMaxD0     ( 100.);
//-----------------------------------------------------------------------------
// track quality
//-----------------------------------------------------------------------------
  fTrackID_MVA = new TStnTrackID();

  fTrackID_MVA->SetMinTrkQual(0.8);
  fTrackID_MVA->SetMinTanDip (0.5);
  fTrackID_MVA->SetMaxTanDip (1.0);
  fTrackID_MVA->SetMinD0     (-100.);
  fTrackID_MVA->SetMaxD0     ( 100.);

  int mask = TStnTrackID::kTrkQualBit | TStnTrackID::kD0Bit | TStnTrackID::kTanDipBit | TStnTrackID::kT0Bit ;
  fTrackID_MVA->SetUseMask(mask);

  fNID           = 2;  			// BOX and MVA
//-----------------------------------------------------------------------------
// particle ID
//-----------------------------------------------------------------------------
//  fLogLH   = new TEmuLogLH();
//-----------------------------------------------------------------------------
// e/mu PID MVA
//-----------------------------------------------------------------------------
  fUsePidMVA = 0;
  fPidMVA[0] = nullptr;
  fPidMVA[1] = nullptr;

  fPidReader = nullptr;
//-----------------------------------------------------------------------------
// multivariate ways of figuring out the track quality
//-----------------------------------------------------------------------------
  fUseTrqMVA = 0;
  fTrqMVA[0] = nullptr;
  fTrqMVA[1] = nullptr;
//-----------------------------------------------------------------------------
// TStntuple 
//-----------------------------------------------------------------------------
  fStnt      = TStntuple::Instance();
}

//-----------------------------------------------------------------------------
TAnaModule::~TAnaModule() {
  if (fTrqMVA[0]) delete fTrqMVA[0];
  if (fTrqMVA[1]) delete fTrqMVA[1];
}

//-----------------------------------------------------------------------------
// TrkRecAlgorithm : "PAR" or "DAR"
// TrainingDataset : just 'fele2s51b1' - goes into the file name
// MVA TrainingCode : 
// ---------------------
// 0060 : PAR dPf > 0.60
// 0070 : PAR dPf > 0.70
// 1060 : DAR dPf > 0.60
// 1070 : DAR dPf > 0.70
//-----------------------------------------------------------------------------
void TAnaModule::SetTrqMVA(int Block, const char* TrainingDataset, int TrainingCode) {

  printf(" [%s::SetTrqMVA] TrainingDataset:%s TrainingCode:%i\n",GetName(),TrainingDataset,TrainingCode);

  if (TrainingCode > 0) { 
    fUseTrqMVA = 1;

    if (fTrqMVA[Block]) delete fTrqMVA[Block];
    fTrqMVA[Block] = new mva_data(TrainingDataset,TrainingCode);
  }
  else { 
    fUseTrqMVA = 0;
  }
}


//-----------------------------------------------------------------------------
void TAnaModule::SetPidMVA(const char* TrainingDataset, int MVATrainingCode) {

  printf(" [%s::SetPidMVA] TrainingDataset:%s TrainingCode:%i\n",GetName(),TrainingDataset,MVATrainingCode);

  if (MVATrainingCode > 0) { 
    fUsePidMVA = 1;

    if (fPidMVA[0]) delete fPidMVA[0];

    fPidMVA[0] = new mva_data(TrainingDataset,MVATrainingCode);  // TMVA::Reader
  }
  else { 
    fUsePidMVA = 0;
  }
}

//-----------------------------------------------------------------------------
// common initializations
//-----------------------------------------------------------------------------
int TAnaModule::BeginRun() {

  TStnDataset* ds = GetAna()->GetInputModule()->GetDataset(0);

  if (ds->GetMcFlag() != 0) {
    int pdg_code     = ds->GetPDGCode();
    int process_code = ds->GetMCProcessCode();

    if (pdg_code     != 0) fPDGCode       = pdg_code;
    if (process_code >= 0) fMCProcessCode = process_code;
  }

  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);

  printf("%s::BeginRun: run: %6i, MCProcessCode: %5i PDGCode: %5i\n",GetName(),rn,fMCProcessCode,fPDGCode);
  
  return 0;
}

//-----------------------------------------------------------------------------
// common initializations
//-----------------------------------------------------------------------------
int TAnaModule::BeginJob() {
					// make sure T0 is the same - it could've been redefined
  fTrackID_BOX->SetMinT0(fMinT0);
  fTrackID_MVA->SetMinT0(fMinT0);
					// PID initialization: read the likelihood templates
  //  fLogLH->Init("v5_7_0");
  return 0;
}

//-----------------------------------------------------------------------------
int TAnaModule::EndJob() {
  return 0;
}


//-----------------------------------------------------------------------------
double TAnaModule::BatchModeWeight(float lumi, int mode) {
  if(mode <= 0) return 1.;
  if(mode > 2 ) return 1.;
//-----------------------------------------------------------------------------
// Batch mode 1/2 log normal initialization
//-----------------------------------------------------------------------------
  const static double mean_b1 = 1.6e7;
  const static double mean_b2 = 3.9e7;
  const static double sigma = 0.7147;
  const static double mub1 = log(mean_b1) - 0.5*sigma*sigma;
  const static double mub2 = log(mean_b2) - 0.5*sigma*sigma;
  const static double cut_off_norm_b1 = ROOT::Math::lognormal_cdf(1.2e8, mub1, sigma); //Due to max cutoff in generation 
  const static double cut_off_norm_b2 = ROOT::Math::lognormal_cdf(1.2e8, mub2, sigma); //Due to max cutoff in generation
  if(mode == 1) {
    const double p1 = ROOT::Math::lognormal_pdf(lumi, mub1, sigma)/cut_off_norm_b1;
    return p1;
  }
  const double p2 = ROOT::Math::lognormal_pdf(lumi, mub2, sigma)/cut_off_norm_b2;
  return p2;
}

//-----------------------------------------------------------------------------
void TAnaModule::BookClusterHistograms(HistBase_t* Hist, const char* Folder) {

  ClusterHist_t* hist = (ClusterHist_t*) Hist;

  HBook1F(hist->fDiskID ,"disk_id",Form("%s: Disk ID"       ,Folder), 10, 0,  10,Folder);
  HBook1F(hist->fEnergy ,"energy" ,Form("%s: Cluster Energy",Folder),500, 0, 250,Folder);
  HBook1F(hist->fEnergyDiff ,"energydiff" ,Form("%s: Cluster Energy - Gen Energy",Folder),500, -50, 50,Folder);
  HBook1F(hist->fT0     ,"t0"     ,Form("%s: cluster T0"    ,Folder),200, 0,2000,Folder);
  HBook1F(hist->fRow    ,"row"    ,Form("%s: cluster Row"   ,Folder),200, 0, 200,Folder);
  HBook1F(hist->fCol    ,"col"    ,Form("%s: cluster column",Folder),200, 0, 200,Folder);
  HBook1F(hist->fX      ,"x"      ,Form("%s: cluster X"     ,Folder),200, -5000,5000,Folder);
  HBook1F(hist->fY      ,"y"      ,Form("%s: cluster Y"     ,Folder),200,-1000,1000,Folder);
  HBook1F(hist->fZ      ,"z"      ,Form("%s: cluster Z"     ,Folder),200, 11500,13500,Folder);
  HBook1F(hist->fR      ,"r"      ,Form("%s: cluster Radius",Folder),100, 0,  1000,Folder);
  HBook1F(hist->fYMean  ,"ymean"  ,Form("%s: cluster YMean" ,Folder),400,-200,200,Folder);
  HBook1F(hist->fZMean  ,"zmean"  ,Form("%s: cluster ZMean" ,Folder),400,-200,200,Folder);
  HBook1F(hist->fSigY   ,"sigy"   ,Form("%s: cluster SigY"  ,Folder),100, 0,100,Folder);
  HBook1F(hist->fSigZ   ,"sigz"   ,Form("%s: cluster SigZ"  ,Folder),100, 0,100,Folder);
  HBook1F(hist->fSigR   ,"sigr"   ,Form("%s: cluster SigR"  ,Folder),100, 0,100,Folder);
  HBook1F(hist->fNCr0   ,"ncr0"   ,Form("%s: cluster NCR[0]",Folder),100, 0,100,Folder);
  HBook1F(hist->fNCr1   ,"ncr1"   ,Form("%s: cluster NCR[1]",Folder),100, 0,100,Folder);
  HBook1F(hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
  HBook1F(hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);
  HBook1F(hist->fSigE1  ,"sige1"   ,Form("%s: SigmaE/Etot"  ,Folder),200, 0, 10,Folder);
  HBook1F(hist->fSigE2  ,"sige2"   ,Form("%s: SigmaE/Emean" ,Folder),200, 0, 10,Folder);
}

//-----------------------------------------------------------------------------
void TAnaModule::BookCrvClusterHistograms   (CrvClusterHist_t*   Hist, const char* Folder){
 
  HBook1F(Hist->fSectorType    ,"sector"     ,Form("%s: CRV sector type" ,Folder),  15,   0,    15,Folder);
  HBook1F(Hist->fNPulses       ,"npulses"    ,Form("%s: N(pulses)"       ,Folder), 100,   0,   100,Folder);
  HBook1F(Hist->fNPe           ,"npe"        ,Form("%s: N(PE)"           ,Folder), 500,   0,  5000,Folder);
  HBook1F(Hist->fStartTime     ,"tstart"     ,Form("%s: start time, ns"  ,Folder), 400,   0,  2000,Folder);  
  HBook1F(Hist->fEndTime       ,"tend"       ,Form("%s: end time, ns"    ,Folder), 400,   0,  2000,Folder);  
  HBook1F(Hist->fWidth         ,"wwidth"     ,Form("%s: width, ns"       ,Folder), 200,   0,  200,Folder);  
  HBook2F(Hist->fXVsZ          ,"x_vs_z"     ,Form("%s: X vs Z"          ,Folder), 250,   -5000,20000,200,-10000,10000,Folder);  
  HBook2F(Hist->fYVsZ          ,"y_vs_z"     ,Form("%s: Y vs Z"          ,Folder), 250,   -5000,20000,200,     0,4000,Folder);  
}

//-----------------------------------------------------------------------------
void TAnaModule::BookCrvPulseHistograms   (CrvPulseHist_t*   Hist, const char* Folder){

  HBook1F(Hist->fNPe       ,"npe"        ,Form("%s: N(Pe)"     ,Folder), 250,   0,   500,Folder);
  HBook1F(Hist->fNPeHeight ,"npe_height" ,Form("%s: NPE_HEIGHT",Folder), 250,   0,   500,Folder);
  HBook1F(Hist->fNDigis    ,"ndigis"     ,Form("%s: N(digis)"  ,Folder),  10,   0,    10,Folder);
  HBook1F(Hist->fBar       ,"bar"        ,Form("%s: bar"       ,Folder), 600,   0,  6000,Folder);
  HBook1F(Hist->fSipm      ,"sipm"       ,Form("%s: sipm"      ,Folder),  10,   0,    10,Folder);

  HBook1F(Hist->fTime      ,"time"       ,Form("%s: time"      ,Folder), 500,   0,  2000,Folder);
  HBook1F(Hist->fHeight    ,"height"     ,Form("%s: height"    ,Folder), 500,   0,  2000,Folder);
  HBook1F(Hist->fWidth     ,"width"      ,Form("%s: width"     ,Folder), 500,   0,   500,Folder);
  HBook1F(Hist->fChi2      ,"chi2"       ,Form("%s: chi2"      ,Folder), 500,   0,   500,Folder);
  HBook1F(Hist->fLeTime    ,"le_time"    ,Form("%s: LE time"   ,Folder), 500,   0,  2000,Folder);
  HBook1F(Hist->fDt        ,"dt"         ,Form("%s: LEt-t"     ,Folder), 100, -50,    50,Folder);
}

//-----------------------------------------------------------------------------
void TAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  EventHist_t* hist = (EventHist_t*) Hist;
  //  char name [200];
  //  char title[200];

  HBook1F(hist->fEventWeight[0],"eventweight_0" ,Form("%s: Event Weight"            ,Folder), 100, -1.,   2.,Folder);
  HBook1F(hist->fEventWeight[1],"eventweight_1" ,Form("%s: Log10(Event Weight)"     ,Folder), 100,-15.,   5.,Folder);
  HBook1F(hist->fEventE        ,"event_energy"  ,Form("%s: Relevant Event Energy"   ,Folder), 400,  0., 200.,Folder);
  HBook1F(hist->fInstLumi[0]   ,"inst_lumi_0"   ,Form("%s: POT"                     ,Folder), 300,  0.,1.5e8,Folder);
  HBook1F(hist->fInstLumi[1]   ,"inst_lumi_1"   ,Form("%s: POT"                     ,Folder), 300,  0.,1.5e8,Folder);
  HBook1F(hist->fInstLumi[2]   ,"inst_lumi_2"   ,Form("%s: POT"                     ,Folder), 300,  0.,1.5e8,Folder);
  HBook1F(hist->fBatchWeight[0],"batchweight_0" ,Form("%s: Log10(One Batch Weight)" ,Folder), 100, -20., 2.,Folder);
  HBook1F(hist->fBatchWeight[1],"batchweight_1" ,Form("%s: Log10(Two Batch Weight)" ,Folder), 100, -20., 2.,Folder);

  HBook1F(hist->fMcCosTh   ,"mc_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
  HBook1F(hist->fMcMom     ,"mc_mom"   ,Form("%s: Conversion Electron Momentum"    ,Folder),1000,  0,200,Folder);
  HBook1D(hist->fDioMom    ,"dio_mom"  ,Form("%s: DIO momentum"                    ,Folder),1000, 50,150,Folder);
  HBook1F(hist->fRv        ,"rv"      ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 1000,Folder);
  HBook1F(hist->fZv        ,"zv"      ,Form("%s: Z(Vertex)"                       ,Folder), 300, 0,15000,Folder);
  //  HBook1F(hist->fNClusters ,"ncl"      ,Form("%s: Number of Reconstructed Clusters",Folder),200,0,200,Folder);
  HBook1F(hist->fNTracks[0] ,"ntrk_0"     ,Form("%s: Number of Reconstructed Tracks[0]"  ,Folder),100,0,100,Folder);
  HBook1F(hist->fNTracks[1] ,"ntrk_1"     ,Form("%s: Number of Reconstructed Tracks[1]"  ,Folder),100,0,100,Folder);
  HBook1F(hist->fNShTot[0],"nsh_0" ,Form("%s: Number of Straw Hits [0]"        ,Folder),250,0,250,Folder);
  HBook1F(hist->fNShTot[1],"nsh_1" ,Form("%s: Number of Straw Hits [1]"        ,Folder),500,0,10000,Folder);
  HBook1F(hist->fNGoodSH   ,"nsh50"    ,Form("%s: N(SH) +/-50"                     ,Folder),300,0,1500,Folder);
  HBook1F(hist->fNChTot[0],"nch_0" ,Form("%s: Number of Combo Hits [0]"        ,Folder),250,0,250,Folder);
  HBook1F(hist->fNChTot[1],"nch_1" ,Form("%s: Number of Combo Hits [1]"        ,Folder),500,0,10000,Folder);
  HBook1F(hist->fDtClT     ,"dt_clt"   ,Form("%s: DT(cluster-track)"               ,Folder),100,-100,100,Folder);
  HBook1F(hist->fDtClS     ,"dt_cls"   ,Form("%s: DT(cluster-straw hit)"           ,Folder),200,-200,200,Folder);
  HBook1F(hist->fSHTime    ,"shtime"   ,Form("%s: Straw Hit Time"                  ,Folder),400,0,2000,Folder);
  HBook1F(hist->fEClMax    ,"eclmax"   ,Form("%s: Max cluster energy"              ,Folder),150,0,150,Folder);
  HBook1F(hist->fNHyp      ,"nhyp"     ,Form("%s: N(fit hypotheses)"               ,Folder),5,0,5,Folder);
  HBook1F(hist->fBestHyp[0],"bfh0"     ,Form("%s: Best Fit Hyp[0](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(hist->fBestHyp[1],"bfh1"     ,Form("%s: Best Fit Hyp[1](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(hist->fNGenp     ,"ngenp"    ,Form("%s: N(Gen Particles)"                ,Folder),500,0,500,Folder);

  HBook1F(hist->fNCrvClusters    ,"ncrvcl"  ,Form("%s: N(CRV Clusters)"                 ,Folder),100,0,  100,Folder);
  HBook1F(hist->fNCrvCoincidences[0],"ncrvco_0"  ,Form("%s: N(CRV coincidences)[0]"     ,Folder),200,0, 1000,Folder);
  HBook1F(hist->fNCrvCoincidences[1],"ncrvco_1"  ,Form("%s: N(CRV coincidences)[1]"     ,Folder),200,0,  200,Folder);
  HBook1F(hist->fNCrvPulses[0]   ,"ncrvp_0" ,Form("%s: N(CRV pulses)[0]"                ,Folder),500,0,10000,Folder);
  HBook1F(hist->fNCrvPulses[1]   ,"ncrvp_1" ,Form("%s: N(CRV pulses)[1]"                ,Folder),500,0,  500,Folder);
}

//-----------------------------------------------------------------------------
void TAnaModule::BookGenpHistograms(GenpHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP         ,"p"     ,Form("%s: Momentum"     ,Folder),1000,     0, 200,Folder);
  HBook1F(Hist->fPdgCode[0],"pdg_0" ,Form("%s: PDG Code[0]"     ,Folder),200, -100, 100,Folder);
  HBook1F(Hist->fPdgCode[1],"pdg_1" ,Form("%s: PDG Code[1]"     ,Folder),500, -2500, 2500,Folder);
  HBook1F(Hist->fGenID     ,"gen_id",Form("%s: Generator ID" ,Folder), 100,     0, 100,Folder);
  HBook1F(Hist->fZ0        ,"z0"    ,Form("%s: Z0"           ,Folder), 500,  5400, 6400,Folder);
  HBook1F(Hist->fT0        ,"t0"    ,Form("%s: T0"           ,Folder), 200,     0, 2000,Folder);
  HBook1F(Hist->fR0        ,"r"     ,Form("%s: R0"           ,Folder), 100,     0,  100,Folder);
  HBook1F(Hist->fCosTh     ,"cos_th",Form("%s: Cos(Theta)"   ,Folder), 200,   -1.,   1.,Folder);
}

//-----------------------------------------------------------------------------
void TAnaModule::BookHelixHistograms(HistBase_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HelixHist_t* hist = (HelixHist_t*) Hist;

  HBook1F(hist->fBestAlg ,"best_alg",Form("%s: best algorithm"       ,Folder),  10,     0  ,  10  ,Folder);
  HBook1F(hist->fChi2XY  ,"chi2xy"  ,Form("%s: Chi2(XY)/DOF"         ,Folder), 200,     0  ,  20  ,Folder);
  HBook1F(hist->fChi2ZPhi,"chi2zp"  ,Form("%s: Chi2(ZP)/DOF"         ,Folder), 200,     0  ,  20  ,Folder);
  HBook1F(hist->fSinTh   ,"sin_th"  ,Form("%s: sin_th = Pt/P"        ,Folder), 200,    -1  ,   1  ,Folder);
  HBook1F(hist->fD0      ,"d0"      ,Form("%s: D0"                   ,Folder), 250,   250  , 250  ,Folder);
  HBook1F(hist->fHelicity,"hel"     ,Form("%s: helicity"             ,Folder),   3,    -1.5,  1.5 ,Folder);
  HBook1F(hist->fLambda  ,"lam"     ,Form("%s: lambda=dphi/dz"       ,Folder), 200,   -500 , 500  ,Folder);
  HBook1F(hist->fNCh     ,"nch"     ,Form("%s: N(combo hits)"        ,Folder), 200,    -0.5, 199.5,Folder);
  HBook1F(hist->fNSh     ,"nsh"     ,Form("%s: N(straw hits)"        ,Folder), 200,    -0.5, 199.5,Folder);
  HBook1F(hist->fP       ,"p"       ,Form("%s: P"                    ,Folder), 250,     0  , 250  ,Folder);
  HBook1F(hist->fRadius  ,"r"       ,Form("%s: Radius"               ,Folder), 200,     0  , 1000 ,Folder);
}

//-----------------------------------------------------------------------------
void TAnaModule::BookSimpHistograms(SimpHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fPdgCode[0]     ,"pdg_0"       ,Form("%s: PDG code[0]"             ,Folder), 200,-100, 100,Folder);
  HBook1F(Hist->fPdgCode[1]     ,"pdg_1"       ,Form("%s: PDG code[1]"             ,Folder),1000,-500, 500,Folder);
  HBook1F(Hist->fNStrawHits     ,"nsth"        ,Form("%s: n straw hits"            ,Folder), 200,   0, 200,Folder);
  HBook1F(Hist->fMomTargetEnd   ,"ptarg"       ,Form("%s: mom after ST"            ,Folder), 400,   0, 200,Folder);
  HBook1F(Hist->fMomTrackerFront,"pfront"      ,Form("%s: mom at the Tracker Front",Folder), 400,   0, 200,Folder);
  HBook1F(Hist->fTime           ,"time"        ,Form("%s: production time"         ,Folder), 200,   0,2000,Folder);
  HBook1F(Hist->fCosTh          ,"costh"       ,Form("%s: cos(theta)"              ,Folder), 200,  -1,   1,Folder);
  HBook1F(Hist->fStartVolumeID  ,"vid0"        ,Form("%s: start volume ID"         ,Folder), 200,2900,3100,Folder);
  HBook1F(Hist->fEndVolumeID    ,"vid1"        ,Form("%s: end volume ID"           ,Folder), 200,2900,3100,Folder);

  HBook2F(Hist->fNshVsCosTh     ,"nsh_vs_costh",Form("%s: nsh vs costh"            ,Folder), 20 ,-1,1,40,0,200,Folder);
}

//-----------------------------------------------------------------------------
void TAnaModule::BookTrackIDHistograms(TStnTrackID::Hist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200]

  for (int i=0; i<5; i++) {
    HBook1F(Hist->fNActive[i],Form("nactive_%i",i) ,Form("%s: Nactive [%i]"        ,Folder,i), 150,   0  , 150. ,Folder);
    HBook1F(Hist->fFitCons[i],Form("fcons_%i"  ,i) ,Form("%s: FitCons [%i]"        ,Folder,i), 200,   0  ,   1. ,Folder);
    HBook1F(Hist->fChi2Dof[i],Form("chi2d_%i"  ,i) ,Form("%s: Chi2/Dof[%i]"        ,Folder,i), 200,   0  ,  20. ,Folder);
    HBook1F(Hist->fT0     [i],Form("t0_%i"     ,i) ,Form("%s: T0      [%i]"        ,Folder,i), 200,   0  ,2000. ,Folder);
    HBook1F(Hist->fT0Err  [i],Form("t0err_%i"  ,i) ,Form("%s: T0Err   [%i]"        ,Folder,i), 200,   0  ,   2. ,Folder);
    HBook1F(Hist->fMomErr [i],Form("momerr_%i" ,i) ,Form("%s: MomErr  [%i]"        ,Folder,i), 200,   0  ,   1. ,Folder);
    HBook1F(Hist->fDNa    [i],Form("dna_%i"    ,i) ,Form("%s: DNa     [%i]"        ,Folder,i), 100,   0  , 100. ,Folder);
    HBook1F(Hist->fTanDip [i],Form("tandip_%i" ,i) ,Form("%s: TanDip  [%i]"        ,Folder,i), 400,   0  ,   4. ,Folder);
    HBook1F(Hist->fD0     [i],Form("d0_%i"     ,i) ,Form("%s: D0      [%i]"        ,Folder,i), 400, -200., 200. ,Folder);
    HBook1F(Hist->fRMax   [i],Form("rmax_%i"   ,i) ,Form("%s: RMax    [%i]"        ,Folder,i), 400,    0., 800. ,Folder);
    HBook1F(Hist->fTrkQual[i],Form("trkqual_%i",i) ,Form("%s: TrkQual [%i]"        ,Folder,i), 200,  -0.5,   1.5,Folder);
  }

  HBook1F(Hist->fPassed    ,"passed"     ,Form("%s: Passed     "         ,Folder),  5,  0,   5,Folder);
  HBook1F(Hist->fFailedBits,"failed_bits",Form("%s: Failed Bits"         ,Folder), 40,  0,  40,Folder);
}

//-----------------------------------------------------------------------------
void TAnaModule::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP[0]       ,"p"        ,Form("%s: Track P(Z1)"       ,Folder), 800,  80  ,120. ,Folder);
  HBook1F(Hist->fP[1]       ,"p_1"      ,Form("%s: Track P(total)[1]" ,Folder), 100, 100  ,105. ,Folder);
  HBook1F(Hist->fP[2]       ,"p_2"      ,Form("%s: Track P(total)[1]" ,Folder),2000,   0  ,200. ,Folder);
  HBook1F(Hist->fP0         ,"p0"       ,Form("%s: Track P(Z0)"       ,Folder),1000,   0  ,200. ,Folder);
  HBook1F(Hist->fP2         ,"p2"       ,Form("%s: Track P(z=-1540)"  ,Folder),1000,   0  ,200. ,Folder);
  HBook1F(Hist->fPDio       ,"pdio"     ,Form("%s: Track P Wt=LL(DIO)",Folder), 100, 100  ,105. ,Folder);
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
  HBook2F(Hist->fPVsGenE    ,"p_vs_gene",Form("%s: Track P vs Gen E"  ,Folder), 300, 0., 150., 300, 0., 150., Folder);

  HBook1F(Hist->fPt         ,"pt"       ,Form("%s: Track Pt"          ,Folder), 600, 75,95,Folder);
  HBook1F(Hist->fCosTh      ,"costh"    ,Form("%s: Track cos(theta)"  ,Folder), 100,-1,1  ,Folder);
  HBook1F(Hist->fChi2       ,"chi2"     ,Form("%s: Track chi2 total"  ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fChi2Dof    ,"chi2d"    ,Form("%s: track chi2/N(dof)" ,Folder), 500, 0, 10,Folder);

  HBook1F(Hist->fNActive    ,"nactv"    ,Form("%s: N(active)"         ,Folder), 200,  0  , 200 ,Folder);
  HBook1F(Hist->fNaFract    ,"nafr"     ,Form("%s: N(active fraction)",Folder), 110,  0.5,1.05 ,Folder);
  HBook1F(Hist->fDNa        ,"dna"      ,Form("%s: Nhits-Nactive"     ,Folder), 100, -0.5 ,99.5,Folder);
  HBook1F(Hist->fNtch       ,"ntch"     ,Form("%s: N timecluster hits",Folder), 200, -0.5 ,199.5,Folder);
  HBook1F(Hist->fNWrong     ,"nwrng"    ,Form("%s: N(wrong drift sgn)",Folder), 100,  0   ,100,Folder);
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

  HBook1F(Hist->fT0         ,"t0"       ,Form("%s: track T0"          ,Folder), 200,    0, 2000,Folder);
  HBook1F(Hist->fT0Err      ,"t0err"    ,Form("%s: track T0Err"       ,Folder), 100,    0,   10,Folder);
  HBook1F(Hist->fDtTc       ,"dttc"     ,Form("%s: T(trk)-T(TC)"      ,Folder), 100, -100,  100,Folder);

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
  HBook1F(Hist->fSeedFr     ,"seed_fr"  ,Form("%s: Eseed/Etot"        ,Folder), 100, 0   ,1  ,Folder);
  HBook1F(Hist->fNCrystals  ,"ncr"      ,Form("%s: N(crystals)"       ,Folder), 100, 0   ,100,Folder);
  HBook1F(Hist->fEClEKin    ,"ecl_ekin" ,Form("%s: cluster E/Ekin(mu)",Folder), 500, 0   ,5,Folder);
  HBook1F(Hist->fEp         ,"ep"       ,Form("%s: track E/P"         ,Folder), 300, 0   ,1.5,Folder);
  HBook1F(Hist->fDrDzCal    ,"drdzcal"  ,Form("%s: track dr/dz cal"   ,Folder), 200, -5  ,5  ,Folder);
  HBook1F(Hist->fDtClZ0     ,"dtclz0"   ,Form("%s: T(cl_z0)-T(Z0)"    ,Folder), 250, -5 , 5,Folder);
  HBook2F(Hist->fDtClZ0VsECl,"dtclz0_vs_ecl",Form("%s: DtClZ0 vs ECl" ,Folder), 100, 0 , 200, 250, -5 , 5,Folder);
  HBook2F(Hist->fDtClZ0VsP  ,"dtclz0_vs_p"  ,Form("%s: DtClZ0 vs p"   ,Folder), 100, 0 , 200, 250, -5 , 5,Folder);

  HBook2F(Hist->fFConsVsNActive,"fc_vs_na" ,Form("%s: FitCons vs NActive",Folder),  150, 0, 150, 200,0,1,Folder);
  HBook1F(Hist->fDaveTrkQual,"dtqual"   ,Form("%s:DaveTrkQual"        ,Folder), 200, -0.5, 1.5,Folder);
  HBook1F(Hist->fTrqMvaOut  ,"trqmvaout",Form("%s:TRQ MVA output"     ,Folder), 200, -0.5, 1.5,Folder);

  HBook2F(Hist->fPVsTime    ,"p_vs_time",Form("%s: P(Z1) vs T(0)"     ,Folder), 800,  80  ,120. , 200, 0, 2000,Folder);

  HBook1F(Hist->fDtCRV      ,"dtcrv"     ,Form("%s:dT(CRV)"            ,Folder), 400,-1000, 1000,Folder);
  HBook1F(Hist->fDtCRV2     ,"dtcrv2"    ,Form("%s:dT(CRV2)"           ,Folder), 400,-1000, 1000,Folder);
  HBook2F(Hist->fZVsDtCRV   ,"dtcrv_vs_z",Form("%s: dT(CRV) vs Z(CRV)" ,Folder), 250,-5000,20000,500,-500, 500,Folder);
//-----------------------------------------------------------------------------
// TrkCaloHit histograms
//-----------------------------------------------------------------------------
  HBook1F(Hist->fTchTime    ,"tch_time"  ,Form("%s:TCH time"           ,Folder), 200,    0, 2000,Folder);
  HBook1F(Hist->fTchPath    ,"tch_path"  ,Form("%s:TCH path"           ,Folder), 200, -500,  500,Folder);
  HBook1F(Hist->fTchDx      ,"tch_dx"    ,Form("%s:TCH Dx"             ,Folder), 200, -500,  500,Folder);
  HBook1F(Hist->fTchDy      ,"tch_dy"    ,Form("%s:TCH Dy"             ,Folder), 200, -500,  500,Folder);
  HBook1F(Hist->fTchDz      ,"tch_dz"    ,Form("%s:TCH DZ"             ,Folder), 200, -500,  500,Folder);
  HBook1F(Hist->fTchDr      ,"tch_dr"    ,Form("%s:TCH DR"             ,Folder), 200, -100,  100,Folder);
  HBook1F(Hist->fTchDt      ,"tch_dt"    ,Form("%s:TCH Dt"             ,Folder), 400,  -10,   10,Folder);

  HBook2F(Hist->fEpVsPath   ,"ep_vs_path",Form("%s: E/P vs Path"       ,Folder), 100, -100,  400, 120,  0,  1.2,Folder);
  HBook2F(Hist->fEpVsTchDz  ,"ep_vs_tchdz",Form("%s: E/P vs Tch_DZ"    ,Folder), 100, -100,  400, 120,  0,  1.2,Folder);
  HBook2F(Hist->fTchDtVsDz  ,"tchdt_vs_dz",Form("%s: TCH_DT vs Tch_DZ" ,Folder), 100, -100,  400, 200, -10, 10 ,Folder);
//-----------------------------------------------------------------------------
// PID histograms
//-----------------------------------------------------------------------------
  HBook1F(Hist->fPidMvaOut  ,"pidmvaout",Form("%s:PID MVA output"     ,Folder), 200, -0.5, 1.5,Folder);
}

//-----------------------------------------------------------------------------
void TAnaModule::BookTrackSeedHistograms   (HistBase_t*   HistR, const char* Folder){
  
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
void TAnaModule::FillClusterHistograms(HistBase_t* Hist, TStnCluster* Cluster, double Weight) {
  // if(Weight == 0.) return; //ignore 0 weight events

  ClusterHist_t* hist = (ClusterHist_t*) Hist;
  int   row, col;
  float  x, y, z, r;

  row = Cluster->Ix1();
  col = Cluster->Ix2();

  x   = Cluster->fX;
  y   = Cluster->fY;
  z   = Cluster->fZ;
  r   = sqrt(x*x+y*y);

  if ((row < 0) || (row > 9999)) row = -9999;
  if ((col < 0) || (col > 9999)) col = -9999;

  hist->fDiskID->Fill(Cluster->DiskID(), Weight);
  hist->fEnergy->Fill(Cluster->Energy(), Weight);
  hist->fEnergyDiff->Fill((Cluster->Energy() - fEleE), Weight); //assuming generated energy is ideal cluster energy
  hist->fT0->Fill(Cluster->Time(), Weight);
  hist->fRow->Fill(row, Weight);
  hist->fCol->Fill(col, Weight);
  hist->fX->Fill(x, Weight);
  hist->fY->Fill(y, Weight);
  hist->fZ->Fill(z, Weight);
  hist->fR->Fill(r, Weight);

  hist->fYMean->Fill(Cluster->fYMean, Weight);
  hist->fZMean->Fill(Cluster->fZMean, Weight);
  hist->fSigY->Fill(Cluster->fSigY, Weight);
  hist->fSigZ->Fill(Cluster->fSigZ, Weight);
  hist->fSigR->Fill(Cluster->fSigR, Weight);
  hist->fNCr0->Fill(Cluster->fNCrystals, Weight);
  hist->fNCr1->Fill(Cluster->fNCr1, Weight);
  hist->fFrE1->Fill(Cluster->fFrE1, Weight);
  hist->fFrE2->Fill(Cluster->fFrE2, Weight);
  hist->fSigE1->Fill(Cluster->fSigE1, Weight);
  hist->fSigE2->Fill(Cluster->fSigE2, Weight);
}

//-----------------------------------------------------------------------------
void TAnaModule::FillCrvClusterHistograms(CrvClusterHist_t* Hist, TCrvCoincidenceCluster* CrvCluster) {

  Hist->fSectorType->Fill(CrvCluster->SectorType());
  Hist->fNPulses->Fill(CrvCluster->NPulses());
  Hist->fNPe->Fill(CrvCluster->NPe());
  Hist->fStartTime->Fill(CrvCluster->StartTime());
  Hist->fEndTime->Fill(CrvCluster->EndTime());

  float width  = CrvCluster->EndTime()-CrvCluster->StartTime();
  Hist->fWidth->Fill(width);

  float x = CrvCluster->Position()->X();
  float y = CrvCluster->Position()->Y();
  float z = CrvCluster->Position()->Z();

  Hist->fXVsZ->Fill(z,x);
  Hist->fYVsZ->Fill(z,y);
}

//-----------------------------------------------------------------------------
void TAnaModule::FillCrvPulseHistograms(CrvPulseHist_t* Hist, TCrvRecoPulse* Pulse) {

  Hist->fNPe->Fill(Pulse->NPe());
  Hist->fNPeHeight->Fill(Pulse->NPeHeight());
  Hist->fNDigis->Fill(Pulse->NDigis());
  Hist->fBar->Fill(Pulse->Bar());
  Hist->fSipm->Fill(Pulse->Sipm());
  Hist->fTime->Fill(Pulse->Time());
  Hist->fHeight->Fill(Pulse->Height());
  Hist->fWidth->Fill(Pulse->Width());
  Hist->fChi2->Fill(Pulse->Chi2());
  Hist->fLeTime->Fill(Pulse->LeTime());

  float dt = Pulse->LeTime()-Pulse->Time();
  Hist->fDt->Fill(dt);
}

//-----------------------------------------------------------------------------
// need MC truth branch
//-----------------------------------------------------------------------------
  void TAnaModule::FillEventHistograms(HistBase_t* Hist, EventPar_t* Evp) {

    EventHist_t* hist = (EventHist_t*) Hist;
    double            cos_th(-2), xv(-1.e6), yv(-1.e6), rv(-1.e6), zv(-1.e6), p(-1.), dio_wt(-1.);
  //  double            e, m, r;
  TLorentzVector    mom;

  if (Evp->fSimp) {
    mom = *Evp->fSimp->StartMom();
    p      = mom.P();
    cos_th = mom.Pz()/p;
    xv     = Evp->fSimp->StartPos()->X()+3904.;
    yv     = Evp->fSimp->StartPos()->Y();
    rv     = sqrt(xv*xv+yv*yv);
    zv     = Evp->fSimp->StartPos()->Z();
    dio_wt = TStntuple::DioWeightAl(p);
  }

  hist->fEventWeight[0]->Fill(fEventWeight);
  hist->fEventWeight[1]->Fill(log10(fEventWeight));
  hist->fEventE        ->Fill(Evp->fPartE  , fEventWeight);
  hist->fInstLumi[0]   ->Fill(Evp->fInstLum, fEventWeight);

  if (fBatchMode == 1) {
    hist->fInstLumi[1]->Fill(Evp->fInstLum, fEventWeight/Evp->fOneBatchWeight);
    hist->fInstLumi[2]->Fill(Evp->fInstLum, fEventWeight/Evp->fOneBatchWeight*Evp->fTwoBatchWeight);
  } else if(fBatchMode == 2) {
    hist->fInstLumi[1]->Fill(Evp->fInstLum, fEventWeight/Evp->fTwoBatchWeight);
    hist->fInstLumi[2]->Fill(Evp->fInstLum, fEventWeight/Evp->fTwoBatchWeight*Evp->fOneBatchWeight);
  }

  hist->fBatchWeight[0]->Fill(log10(Evp->fOneBatchWeight));
  hist->fBatchWeight[1]->Fill(log10(Evp->fTwoBatchWeight));

  hist->fMcMom->Fill(p);
  hist->fDioMom->Fill(p,dio_wt);
  hist->fMcCosTh->Fill(cos_th);
  hist->fRv->Fill(rv);
  hist->fZv->Fill(zv);

  // hist->fNClusters->Fill(fNClusters);
  hist->fNTracks[0]->Fill  (Evp->fNTracks[0]);
  hist->fNTracks[1]->Fill  (Evp->fNTracks[1]);
  hist->fNShTot[0]->Fill(Evp->fNStrawHits);
  hist->fNShTot[1]->Fill(Evp->fNStrawHits);

  hist->fNChTot[0]->Fill(Evp->fNComboHits);
  hist->fNChTot[1]->Fill(Evp->fNComboHits);

  double emax   = -1;
  double dt     = 9999.;

  hist->fDtClT->Fill(dt);
  hist->fEClMax->Fill(emax);

  // hist->fNHyp->Fill(fNHyp);
  // hist->fBestHyp[0]->Fill(fBestHyp[0]);
  // hist->fBestHyp[1]->Fill(fBestHyp[1]);
  hist->fNGenp->Fill(Evp->fNGenp);

  hist->fNCrvClusters->Fill(Evp->fNCrvClusters);
  hist->fNCrvCoincidences[0]->Fill(Evp->fNCrvCoincidences);
  hist->fNCrvCoincidences[1]->Fill(Evp->fNCrvCoincidences);
  hist->fNCrvPulses[0]->Fill(Evp->fNCrvPulses);
  hist->fNCrvPulses[1]->Fill(Evp->fNCrvPulses);
}

//-----------------------------------------------------------------------------
void TAnaModule::FillGenpHistograms(GenpHist_t* Hist, TGenParticle* Genp) {
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
  void TAnaModule::FillHelixHistograms(HistBase_t* Hist, TStnHelix* Helix, HelixPar_t* Help, double Weight) {
//   char name [200];
//   char title[200];

    HelixHist_t* hist = (HelixHist_t*) Hist;

    hist->fBestAlg->Fill(Helix->BestAlg(),Weight);
    hist->fChi2XY->Fill(Helix->Chi2XY(),Weight);
    hist->fChi2ZPhi->Fill(Helix->Chi2ZPhi(),Weight);

    float sinth = Helix->Pt()/Helix->P();
    hist->fSinTh->Fill(sinth,Weight);

    hist->fD0->Fill(Helix->D0(),Weight);

    hist->fHelicity->Fill(Helix->Helicity(),Weight);
    hist->fLambda->Fill(Helix->Lambda(),Weight);
    hist->fNCh->Fill(Helix->fNComboHits,Weight);
    hist->fNSh->Fill(Helix->fNHits,Weight);
    hist->fP->Fill(Helix->P(),Weight);
    hist->fRadius->Fill(Helix->Radius(),Weight);
}

//-----------------------------------------------------------------------------
void TAnaModule::FillSimpHistograms(SimpHist_t* Hist, TSimParticle* Simp, double Weight) {

  Hist->fPdgCode[0]->Fill(Simp->fPdgCode,Weight);
  Hist->fPdgCode[1]->Fill(Simp->fPdgCode,Weight);
  Hist->fMomTargetEnd->Fill(Simp->fMomTargetEnd,Weight);
  Hist->fMomTrackerFront->Fill(Simp->fMomTrackerFront,Weight);

  int nsh = Simp->NStrawHits();
  Hist->fNStrawHits->Fill(nsh,Weight);

  Hist->fTime->Fill(Simp->StartPos()->T(),Weight);

  double costh = Simp->StartMom()->Z()/Simp->StartMom()->P();
  Hist->fCosTh->Fill(costh,Weight);

  Hist->fStartVolumeID->Fill(Simp->fStartVolumeIndex,Weight);
  Hist->fEndVolumeID->Fill(Simp->fEndVolumeIndex,Weight);

  Hist->fNshVsCosTh->Fill(costh,nsh,Weight);
}


//-----------------------------------------------------------------------------
// for DIO : ultimately, one would need to renormalize the distribution
//-----------------------------------------------------------------------------
void TAnaModule::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track, 
				     TrackPar_t* Tp, SimPar_t* SimPar, double Weight) {

  TLorentzVector  mom;
					// Tp->fP - corrected momentum, fP0 and fP2 - not corrected
  Hist->fP[0]->Fill (Tp->fP,Weight);
  Hist->fP[1]->Fill (Tp->fP,Weight);
  Hist->fP[2]->Fill (Tp->fP,Weight);
					// fP0: uncorrected momentum in the first point,  fP2 - in the last
  Hist->fP0->  Fill (Track->fP0,Weight);
  Hist->fP2->  Fill (Track->fP2,Weight);

  Hist->fPDio->  Fill (Tp->fP,Tp->fDioLLWt);     // fixed for debugging

  Hist->fFitMomErr->Fill(Track->fFitMomErr,Weight);

  Hist->fPt    ->Fill(Track->fPt    , Weight);
  Hist->fPFront->Fill(Track->fPFront, Weight);
  Hist->fPStOut->Fill(Track->fPStOut, Weight);
  Hist->fPVsGenE->Fill(fEleE, Track->fP, Weight);
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
  Hist->fChi2Dof->Fill(Track->Chi2Dof(), Weight);

  float na  = Track->NActive();
  float dna = Track->NHits()-na;

  Hist->fNActive->Fill(na, Weight);
  Hist->fNaFract->Fill(na/(Track->NHits()+0.), Weight);
  Hist->fDNa->Fill(dna, Weight);

  if (Tp->fTimeCluster != nullptr) Hist->fNtch->Fill(Tp->fTimeCluster->NHits(),Weight);
  else                             Hist->fNtch->Fill(-1,Weight);
  Hist->fDtTc->Fill(Tp->fDtTc,Weight);

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
  Hist->fSeedFr->Fill(Tp->fSeedFr, Weight);
  Hist->fNCrystals->Fill(Tp->fNCrystals, Weight);
//-----------------------------------------------------------------------------
// assume muon hypothesis
//-----------------------------------------------------------------------------
  double    ekin(-1.);
  if (SimPar->fParticle) {
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
  Hist->fTrqMvaOut->Fill(Tp->fTrqMvaOut[1], Weight);
//-----------------------------------------------------------------------------
// 2D momentum() vs time (tracker center) histogram for sensitivity calculation
//-----------------------------------------------------------------------------
  Hist->fPVsTime->Fill(Tp->fP, Track->fT0, Weight);

  Hist->fDtCRV   ->Fill(Tp->fDtCRV , Weight);
  Hist->fDtCRV2  ->Fill(Tp->fDtCRV2, Weight);
  Hist->fZVsDtCRV->Fill(Tp->fZCRV  , Tp->fDtCRV, Weight);
//-----------------------------------------------------------------------------
// TrkCaloHit histograms
//-----------------------------------------------------------------------------
  Hist->fTchTime->Fill(Tp->fTchTime, Weight);
  Hist->fTchPath->Fill(Tp->fTchPath, Weight);
  Hist->fTchDx  ->Fill(Tp->fTchDx  , Weight);
  Hist->fTchDy  ->Fill(Tp->fTchDy  , Weight);
  Hist->fTchDz  ->Fill(Tp->fTchDz  , Weight);
  Hist->fTchDr  ->Fill(Tp->fTchDr  , Weight);
  Hist->fTchDt  ->Fill(Tp->fTchDt  , Weight);

  Hist->fEpVsPath ->Fill(Tp->fPath ,Tp->fEp   ,Weight);
  Hist->fEpVsTchDz->Fill(Tp->fTchDz,Tp->fEp   ,Weight);
  Hist->fTchDtVsDz->Fill(Tp->fTchDz,Tp->fTchDt,Weight);
//-----------------------------------------------------------------------------
// MVA-based PID
//-----------------------------------------------------------------------------
  Hist->fPidMvaOut->Fill(Tp->fPidMvaOut[0], Weight);
}


//--------------------------------------------------------------------------------
// function to fill TrasckSeedHit block
//--------------------------------------------------------------------------------
  void TAnaModule::FillTrackSeedHistograms(HistBase_t*   HistR, TStnTrackSeed* TrkSeed, double Weight) {
  
  TrackSeedHist_t* Hist = (TrackSeedHist_t*) HistR;
  
  int         nhits    = TrkSeed->NHits      ();
  double      clusterT = TrkSeed->ClusterTime();
  double      clusterE = TrkSeed->ClusterEnergy();
  
  double      mm2MeV   = 3/10.;
  double      pT       = TrkSeed->Pt();
  double      radius   = pT/mm2MeV;

  double      tanDip   = TrkSeed->TanDip();  
  double      p        = pT/std::cos( std::atan(tanDip));
  
  Hist->fNHits      ->Fill(nhits,Weight);	 
  Hist->fClusterTime->Fill(clusterT,Weight);
  Hist->fClusterEnergy->Fill(clusterE,Weight);
  
  Hist->fRadius     ->Fill(radius,Weight);    
  Hist->fMom        ->Fill(p,Weight);	 
  Hist->fPt         ->Fill(pT,Weight);	 
  Hist->fTanDip     ->Fill(tanDip,Weight);    
  
  Hist->fChi2       ->Fill(TrkSeed->Chi2()/(nhits-5.),Weight);
  Hist->fFitCons    ->Fill(TrkSeed->FitCons(),Weight);
  Hist->fD0         ->Fill(TrkSeed->D0(),Weight);
}

//-----------------------------------------------------------------------------
// on input: TrackPar->fFitType = 0: PAR, =1: DAR
//-----------------------------------------------------------------------------
int TAnaModule::InitTrackPar(TStnTrackBlock*     TrackBlock  , 
			     TStnClusterBlock*   ClusterBlock, 
			     TrackPar_t*         TrackPar    ,
			     SimPar_t*           SimPar      ) {
//-----------------------------------------------------------------------------
// momentum corrections for KPAR and KDAR
//-----------------------------------------------------------------------------
  const double kMomentumCorr[2] = { 0.0 , 0.0 } ; // SU2020 // CD3: { 0.049, 0.020 };
  const double kDtTcmCorr   [2] = { 0.  , 0.  } ; // { 0.22 , -0.30 }; // ns, sign: fit peak positions
//-----------------------------------------------------------------------------
// loop over tracks
//-----------------------------------------------------------------------------
  int ntrk = TrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    TrackPar_t*   tp = TrackPar+itrk;
    TStnTrack* track = TrackBlock->Track(itrk);
    int fit_type = 0;
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
//-----------------------------------------------------------------------------
    if (fApplyCorr != 0) tp->fP = track->fP0 + kMomentumCorr[fit_type]; 
//-----------------------------------------------------------------------------
// hits on virtual detectors
//-----------------------------------------------------------------------------
    tp->fPStOut  = track->fPStOut;
    tp->fPFront  = track->fPFront;

    tp->fDpF     = tp->fP     -tp->fPFront;
    tp->fDp0     = track->fP0 -tp->fPFront;
    tp->fDp2     = track->fP2 -tp->fPFront;
    tp->fDpFSt   = tp->fPFront-tp->fPStOut;

    tp->fXDpF    = tp->fDpF/track->fFitMomErr;

    tp->fLumWt   = GetHeaderBlock()->LumWeight();

    tp->fDtZ0 = -1.e6;
    if (SimPar->fTMid) {
      double ttrue = fmod(SimPar->fTMid->Time(),fMbTime);
      tp->fDtZ0 = track->T0()-ttrue;
    }

    tp->fXtZ0 = tp->fDtZ0/track->fT0Err;

    tp->fDtBack = -1.e6;
    if (SimPar->fTBack) {
      double ttrue = fmod(SimPar->fTBack->Time(),fMbTime);
      tp->fDtBack = track->T0()-ttrue;
    }
//-----------------------------------------------------------------------------
// track residuals
//-----------------------------------------------------------------------------
    TStnTrack::InterData_t*  vr = track->fVMaxEp; 
    double    nx, ny;

    tp->fCluster   = nullptr;

    tp->fEcl       = -1.e6;
    tp->fSeedFr    = -1.e6;
    tp->fNCrystals = -1.e6;
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

    tp->fChi2Tcm   = -1.e6;
    tp->fChi2XY    = -1.e6;
    tp->fChi2T     = -1.e6;
    tp->fPath      = -1.e6;
    tp->fSinTC     = -1.e6;
    tp->fDrTC      = -1.e6;
    tp->fSInt      = -1.e6;

    if (vr) {
      tp->fDiskID  = vr->fID;
      tp->fCluster = ClusterBlock->Cluster(vr->fClusterIndex);
      tp->fSeedFr  = tp->fCluster->SeedFr();
      tp->fNCrystals = tp->fCluster->NCrystals();
      tp->fEcl     = vr->fEnergy;
      tp->fEp      = tp->fEcl/track->fP2;
      tp->fDrDzCal = (vr->fXTrk*vr->fNxTrk+vr->fYTrk+vr->fNyTrk)/sqrt(vr->fXTrk*vr->fXTrk+vr->fYTrk*vr->fYTrk)/vr->fNzTrk;

      tp->fDx      = vr->fDx;
      tp->fDy      = vr->fDy;
      tp->fDz      = vr->fDz;
//-----------------------------------------------------------------------------
// v4_2_4: correct by additional 0.22 ns - track propagation by 6 cm
//-----------------------------------------------------------------------------
      tp->fDt  = vr->fDt ;                                       // v4_2_4: - 0.22; // - 1.;
      if (fApplyCorr != 0) tp->fDt  -= kDtTcmCorr[fit_type];

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
      TStnCluster* cl = ClusterBlock->Cluster(vr->fClusterIndex);
      tp->fSinTC = nx*cl->fNy-ny*cl->fNx;
      tp->fDrTC  = vr->fDr;
      tp->fSInt  = vr->fSInt;

      if (SimPar->fTMid) {
	tp->fDtClZ0 = tp->fDt-tp->fDtZ0;
      }
    }

    if ((tp->fEp > 0) && (track->fEp > 0) && (fabs(tp->fEp-track->fEp) > 1.e-6)) {
      GetHeaderBlock()->Print(Form(" TAnaModule ERROR: tp->fEp = %10.5f  track->fEp = %10.5f\n ",tp->fEp,track->fEp));
    }
//-----------------------------------------------------------------------------
// on-the-fly MVA calculation
// tp->fMVAOut[0] : comes from offline
// tp->fMVAOut[1] : calculated on the fly
//-----------------------------------------------------------------------------
    tp->fTrqMvaOut[0] = track->DaveTrkQual();        // comes from Offline
    tp->fTrqMvaOut[1] = track->DaveTrkQual(); 

    if (fUseTrqMVA != 0) {
//-----------------------------------------------------------------------------
// MVA output calculated on the fly - in principle, should be charge-symmetric
// can have several TRQ MVA's defined, each track block uses only one of them,
// indexed with tp->fTrqMvaIndex
//-----------------------------------------------------------------------------
      float na = track->NActive();
      float nm = track->NMat();

      int imva = tp->fTrqMvaIndex;	                          // 
	
      if (fTrqMVA[imva] != nullptr) {
	fTrqMVA[imva]->fVar[0] = na;
	fTrqMVA[imva]->fVar[1] = na/track->NHits();
	fTrqMVA[imva]->fVar[2] = log10(track->FitCons());
	fTrqMVA[imva]->fVar[3] = track->FitMomErr();
	fTrqMVA[imva]->fVar[4] = track->T0Err();
	fTrqMVA[imva]->fVar[5] = track->NDoubletsAct()/na;
	fTrqMVA[imva]->fVar[6] = track->NHitsAmbZero()/na;
	fTrqMVA[imva]->fVar[7] = track->NMatActive()/nm;
      // pmva[ 8] = track->D0();                          // low-rank, do not use
      // pmva[ 9] = track->RMax();                        // low-rank, do not use

	tp->fTrqMvaOut[1] = fTrqMVA[imva]->Eval();

	if (tp->fTrqMvaOut[1] < -990) {
	  GetHeaderBlock()->Print(Form("BAD TRQ: ntrk, itrk= %2i %2i",ntrk,itrk));
	  PrintTrack(track,tp,"");
	}

	track->SetTmp(0,tp->fTrqMvaOut[1]);
      }
    }
//-----------------------------------------------------------------------------
// track ID 
//-----------------------------------------------------------------------------
    for (int k=0; k<fNID; k++) {
      tp->fIDWord[k]  = tp->fTrackID[k]->IDWord(track);
    }
//-----------------------------------------------------------------------------
// track-CRV residuals
//-----------------------------------------------------------------------------
    tp->fDtCRV       = -1.e6;
    tp->fZCRV        = -1.e6;
    tp->fDtCRV2      = -1.e6;
//-----------------------------------------------------------------------------
// TrackCaloHit
//-----------------------------------------------------------------------------
    TStnTrack::InterData_t*  tch = track->fVTCH;
  
    tp->fTchTime   = tch->fTime;
    tp->fTchPath   = tch->fPath;
    tp->fTchDr     = tch->fDr;
    tp->fTchDx     = tch->fDx;
    tp->fTchDy     = tch->fDy;
    tp->fTchDz     = tch->fDz;
    tp->fTchDt     = tch->fDt;
//-----------------------------------------------------------------------------
// MVA-based PID - can't calculate here, as need two reco outputs
//-----------------------------------------------------------------------------
    tp->fPidMvaOut[0] = -1.e6;
//-----------------------------------------------------------------------------
// both electron and muon tracks are there, need to apply preselections
//-----------------------------------------------------------------------------
    if (fUsePidMVA != 0) {
      if ((tp->fTchDt > -100) && (tp->fTchDz > -50) && (tp->fTchDz < 250)) {

	fPidMVA[0]->fVar[0] = tp->fEp;
	fPidMVA[0]->fVar[1] = tp->fNCrystals;
	fPidMVA[0]->fVar[2] = tp->fSeedFr;
	fPidMVA[0]->fVar[3] = tp->fTchDt;
	fPidMVA[0]->fVar[4] = tp->fTchDz;
	fPidMVA[0]->fVar[5] = tp->fTchDr;
	fPidMVA[0]->fVar[6] = tp->fPath;

	tp->fPidMvaOut[0] = fPidMVA[0]->Eval();
      }
    }
  }

  return 0;
}


//_____________________________________________________________________________
void TAnaModule::PrintTrack(TStnTrack* Track, TrackPar_t* Tp, Option_t* Option) const {

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
    printf("      Z0    TanDip   TBack   chi2/dof   fcon  TrkQ    MvaOut[0]  MVAOut[1] \n");
    printf("------------------------------------------------------------------------------------------------");
    printf("----------------------------------------------------------------------------------\n");
  }

  if ((opt == "") || (opt.Index("data") >= 0)) {
    printf("%2i %3i %3i %2i %4i %4i %3i %3i 0x%08x",
	   t->fNumber,t->NHits(), t->NActive(),t->NWrong(), 
	   t->NOSDoublets(), t->NSSDoublets(), t->NHitsAmbZero(),
	   t->NClusters(),
	   t->AlgorithmID());

    printf(" 0x%08x %1.0f %8.3f %8.3f %7.3f %8.3f %6.3f %7.3f %8.3f %7.4f %8.3f %8.2f %8.2e %7.3f %9.4f %9.4f",
	   t->fIDWord,
	   t->fCharge, 
	   t->fP*t->fCharge, Tp->fP*t->fCharge, t->fFitMomErr, t->fT0, t->fT0Err, t->fD0, t->fZ0, t->fTanDip, t->TBack(),
	   t->Chi2Dof(),t->FitCons(),t->DaveTrkQual(),Tp->fTrqMvaOut[0], Tp->fTrqMvaOut[1]);
    printf("\n");
  }
}

}
