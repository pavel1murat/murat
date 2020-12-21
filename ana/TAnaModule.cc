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
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
#include "fhiclcpp/ParameterSet.h"
#include <xercesc/dom/DOM.hpp>
#include "Mu2eUtilities/inc/MVATools.hh"
//------------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
#include "murat/ana/TAnaModule.hh"

using std::string;
using std::vector;

using namespace murat;

namespace murat {

//-----------------------------------------------------------------------------
TAnaModule::TAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fMbTime      = 1695.;

  fMinT0       = 0.; // analysis cuts
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

  int mask = TStnTrackID::kTrkQualBit | TStnTrackID::kD0Bit | TStnTrackID::kTanDipBit | TStnTrackID::kRMaxBit | TStnTrackID::kT0Bit ;
  fTrackID_MVA->SetUseMask(mask);

  fNID     = 2;  			// BOX and MVA
//-----------------------------------------------------------------------------
// particle ID
//-----------------------------------------------------------------------------
  fLogLH   = new TEmuLogLH();
}

//-----------------------------------------------------------------------------
TAnaModule::~TAnaModule() {
}

//-----------------------------------------------------------------------------
// common initializations
//-----------------------------------------------------------------------------
int TAnaModule::BeginJob() {
					// make sure T0 is the same - it could've been redefined
  fTrackID_BOX->SetMinT0(fMinT0);
  fTrackID_MVA->SetMinT0(fMinT0);
					// PID initialization: read the likelihood templates
  fLogLH->Init("v5_7_0");

  return 0;
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
void TAnaModule::BookGenpHistograms(GenpHist_t* Hist, const char* Folder) {
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
  HBook1F(Hist->fTchDz      ,"tch_dz"    ,Form("%s:TCH DZ"             ,Folder), 100, -250,  250,Folder);
  HBook1F(Hist->fTchDr      ,"tch_dr"    ,Form("%s:TCH DR"             ,Folder), 200, -100,  100,Folder);
  HBook1F(Hist->fTchDt      ,"tch_dt"    ,Form("%s:TCH Dt"             ,Folder), 400,  -10,   10,Folder);
  
}

//-----------------------------------------------------------------------------
void TAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fMcCosTh    ,"mc_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
  HBook1F(Hist->fMcMom      ,"mc_mom"   ,Form("%s: Conversion Electron Momentum"    ,Folder),1000,  0,200,Folder);
  //  HBook1D(Hist->fDioMom     ,"dio_mom"  ,Form("%s: DIO momentum"                    ,Folder),1000, 50,150,Folder);
  HBook1F(Hist->fRv         ,"rv"      ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 1000,Folder);
  HBook1F(Hist->fZv         ,"zv"      ,Form("%s: Z(Vertex)"                       ,Folder), 300, 0,15000,Folder);
  //  HBook1F(Hist->fNClusters ,"ncl"      ,Form("%s: Number of Reconstructed Clusters",Folder),200,0,200,Folder);
  HBook1F(Hist->fNTracks[0] ,"ntrk_0"     ,Form("%s: Number of Reconstructed Tracks[0]"  ,Folder),100,0,100,Folder);
  HBook1F(Hist->fNTracks[1] ,"ntrk_1"     ,Form("%s: Number of Reconstructed Tracks[1]"  ,Folder),100,0,100,Folder);
  HBook1F(Hist->fNShTot[0],"nsh_0" ,Form("%s: Number of Straw Hits [0]"        ,Folder),250,0,250,Folder);
  HBook1F(Hist->fNShTot[1],"nsh_1" ,Form("%s: Number of Straw Hits [1]"        ,Folder),250,0,5000,Folder);
  HBook1F(Hist->fNGoodSH   ,"nsh50"    ,Form("%s: N(SH) +/-50"                     ,Folder),300,0,1500,Folder);
  HBook1F(Hist->fDtClT     ,"dt_clt"   ,Form("%s: DT(cluster-track)"               ,Folder),100,-100,100,Folder);
  HBook1F(Hist->fDtClS     ,"dt_cls"   ,Form("%s: DT(cluster-straw hit)"           ,Folder),200,-200,200,Folder);
  HBook1F(Hist->fSHTime    ,"shtime"   ,Form("%s: Straw Hit Time"                  ,Folder),400,0,2000,Folder);
  HBook1F(Hist->fEClMax    ,"eclmax"   ,Form("%s: Max cluster energy"              ,Folder),150,0,150,Folder);
  HBook1F(Hist->fNHyp      ,"nhyp"     ,Form("%s: N(fit hypotheses)"               ,Folder),5,0,5,Folder);
  HBook1F(Hist->fBestHyp[0],"bfh0"     ,Form("%s: Best Fit Hyp[0](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(Hist->fBestHyp[1],"bfh1"     ,Form("%s: Best Fit Hyp[1](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(Hist->fNGenp     ,"ngenp"    ,Form("%s: N(Gen Particles)"                ,Folder),500,0,500,Folder);

  HBook1F(Hist->fNCrvClusters    ,"ncrvcl"  ,Form("%s: N(CRV Clusters)"                 ,Folder),100,0,  100,Folder);
  HBook1F(Hist->fNCrvCoincidences[0],"ncrvco_0"  ,Form("%s: N(CRV coincidences)[0]"     ,Folder),200,0, 1000,Folder);
  HBook1F(Hist->fNCrvCoincidences[1],"ncrvco_1"  ,Form("%s: N(CRV coincidences)[1]"     ,Folder),200,0,  200,Folder);
  HBook1F(Hist->fNCrvPulses[0]   ,"ncrvp_0" ,Form("%s: N(CRV pulses)[0]"                ,Folder),500,0,10000,Folder);
  HBook1F(Hist->fNCrvPulses[1]   ,"ncrvp_1" ,Form("%s: N(CRV pulses)[1]"                ,Folder),500,0,  500,Folder);
}

//-----------------------------------------------------------------------------
void TAnaModule::BookSimpHistograms(SimpHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fPdgCode   ,"pdg"         ,Form("%s: PDG code"                     ,Folder),200,-100,100,Folder);
  HBook1F(Hist->fNStrawHits,"nsth"        ,Form("%s: n straw hits"                 ,Folder),200,   0,200,Folder);
  HBook1F(Hist->fMomTargetEnd    ,"ptarg" ,Form("%s: CE mom after Stopping Target" ,Folder),400,  90,110,Folder);
  HBook1F(Hist->fMomTrackerFront ,"pfront",Form("%s: CE mom at the Tracker Front"  ,Folder),400,  90,110,Folder);
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

//_____________________________________________________________________________
// void TAnaModule::BookHistograms() {
// }


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
  void TAnaModule::FillEventHistograms(EventHist_t* Hist, EventPar_t* Evp) {
  double            cos_th(-2), xv(-1.e6), yv(-1.e6), rv(-1.e6), zv(-1.e6), p(-1.);
  //  double            e, m, r;
  TLorentzVector    mom;

  if (Evp->fParticle) {
    Evp->fParticle->Momentum(mom);
    p      = mom.P();
    cos_th = mom.Pz()/p;
    xv     = Evp->fParticle->Vx()+3904.;
    yv     = Evp->fParticle->Vy();
    rv     = sqrt(xv*xv+yv*yv);
    zv     = Evp->fParticle->Vz();
    //    dio_wt = TStntuple::DioWeightAl(p);
  }

  Hist->fMcMom->Fill(p);
  Hist->fMcCosTh->Fill(cos_th);
  Hist->fRv->Fill(rv);
  Hist->fZv->Fill(zv);

  // Hist->fNClusters->Fill(fNClusters);
  Hist->fNTracks[0]->Fill  (Evp->fNTracks[0]);
  Hist->fNTracks[1]->Fill  (Evp->fNTracks[1]);
  Hist->fNShTot[0]->Fill(Evp->fNStrawHits);
  Hist->fNShTot[1]->Fill(Evp->fNStrawHits);

  double emax   = -1;
  double dt     = 9999.;

  Hist->fDtClT->Fill(dt);
  Hist->fEClMax->Fill(emax);

  // Hist->fNHyp->Fill(fNHyp);
  // Hist->fBestHyp[0]->Fill(fBestHyp[0]);
  // Hist->fBestHyp[1]->Fill(fBestHyp[1]);
  Hist->fNGenp->Fill(Evp->fNGenp);

  Hist->fNCrvClusters->Fill(Evp->fNCrvClusters);
  Hist->fNCrvCoincidences[0]->Fill(Evp->fNCrvCoincidences);
  Hist->fNCrvCoincidences[1]->Fill(Evp->fNCrvCoincidences);
  Hist->fNCrvPulses[0]->Fill(Evp->fNCrvPulses);
  Hist->fNCrvPulses[1]->Fill(Evp->fNCrvPulses);
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
void TAnaModule::FillSimpHistograms(SimpHist_t* Hist, TSimParticle* Simp) {

  Hist->fPdgCode->Fill(Simp->fPdgCode);
  Hist->fMomTargetEnd->Fill(Simp->fMomTargetEnd);
  Hist->fMomTrackerFront->Fill(Simp->fMomTrackerFront);
  Hist->fNStrawHits->Fill(Simp->fNStrawHits);
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
//-----------------------------------------------------------------------------
// 2D momentum() vs time (tracker center) histogram for sensitivity calculation
//-----------------------------------------------------------------------------
  Hist->fPVsTime->Fill(Tp->fP, Track->fT0, Weight);

  Hist->fDtCRV   ->Fill(Tp->fDtCRV, Weight);
  Hist->fZVsDtCRV->Fill(Tp->fZCRV , Tp->fDtCRV, Weight);
  Hist->fDtCRV2  ->Fill(Tp->fDtCRV2, Weight);
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
}


//--------------------------------------------------------------------------------
// function to fill TrasckSeedHit block
//--------------------------------------------------------------------------------
void TAnaModule::FillTrackSeedHistograms(HistBase_t*   HistR, TStnTrackSeed* TrkSeed) {
  
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
int TAnaModule::InitTrackPar(TStnTrackBlock*           TrackBlock  , 
			     TStnClusterBlock*         ClusterBlock, 
			     TrackPar_t*       TrackPar    ,
			     SimPar_t*         SimPar      ) {
  TrackPar_t*           tp;
  TStnTrack*            track;
  int                   track_type(-999);
  double                xs;
  TEmuLogLH::PidData_t  dat;
  static int            printed(0);
//-----------------------------------------------------------------------------
// momentum corrections for KPAR and KDAR
//-----------------------------------------------------------------------------
  const double kMomentumCorr[2] = { 0. , 0. } ; // CD3: { 0.049, 0.020 };
  const double kDtTcmCorr   [2] = { 0. , 0. } ; // { 0.22 , -0.30 }; // ns, sign: fit peak positions

  const char* block_name = TrackBlock->GetNode()->GetName();

  track_type = 0;
  if      (strcmp(block_name,"TrackBlockPar") == 0) track_type = 0;
  else if (strcmp(block_name,"TrackBlockDar") == 0) track_type = 1;
  else {
    if (printed <= 10) {
      Error("murat::TAnaModule::InitTrackPar",Form("COULD BE IN TROUBLE: unknown track block: %s",block_name));
      printed++;
    }
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
//-----------------------------------------------------------------------------
    if (track_type == 2) track_type = track->BestAlg();

    tp->fP     = track->fP0    +kMomentumCorr[track_type];		// correcting
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
    if (SimPar->fParticle) {
      double e     = SimPar->fParticle->StartMom()->E();
      tp->fDioLOWt = TStntuple::DioWeightAl   (e);
      tp->fDioLLWt = TStntuple::DioWeightAl_LL(e);
    }
    else {
      tp->fDioLOWt = 1. ; // TStntuple::DioWeightAl   (SimPar->fEleE);
      tp->fDioLLWt = 1. ; // TStntuple::DioWeightAl_LL(SimPar->fEleE);
    }
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
//-----------------------------------------------------------------------------
    tp->fMVAOut[0] = -1.e6;
    tp->fMVAOut[1] = -1.e6;

    if (fUseMVA != 0) {
      vector<double>  pmva;
//-----------------------------------------------------------------------------
// alg=0: TrkPatRec track - log(fitcons) used for training
// tp->fMVAOut[0] : comes from offline
// tp->fMVAOut[1] : calculated on the fly
//-----------------------------------------------------------------------------
      int alg = track->BestAlg();
      
      if      (alg == 0) {
	tp->fMVAOut[0] = track->DaveTrkQual();        // comes from Offline

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

	if (fTmvaAlgorithmTpr >= 0) {
	  int use_chi2d = fTmvaAlgorithmTpr / 100;
	  if      (use_chi2d == 0) pmva[2] = log10(track->FitCons());
	  else                     pmva[2] = track->Chi2Dof();

	  tp->fMVAOut[1] = fTprQualMva->evalMVA(pmva);
	}
      }

//-----------------------------------------------------------------------------
// CalPatRec track - for fTmvaAlgorithm=101 chi2/N(dof) used (our initial training)
//-----------------------------------------------------------------------------
      else if (alg == 1) {

	pmva.resize(11);

	float na = track->NActive();

	pmva[ 0] = na;
	pmva[ 1] = na/track->NHits();
	pmva[ 2] = -1.e6; 			// defined below
	pmva[ 3] = track->FitMomErr();
	pmva[ 4] = track->T0Err();
	pmva[ 5] = track->D0();
	pmva[ 6] = track->RMax();
	pmva[ 7] = track->TanDip();
	pmva[ 8] = track->NDoubletsAct()/na;
	pmva[ 9] = track->NHitsAmbZero()/na;
	pmva[10] = track->NMatActive()/na;


	int use_chi2d = fTmvaAlgorithmCpr / 100;
	if      (use_chi2d == 0) pmva[2] = log10(track->FitCons());
	else                     pmva[2] = track->Chi2Dof();

	tp->fMVAOut[1] = fCprQualMva->evalMVA(pmva);
      }
    }
//-----------------------------------------------------------------------------
// track ID 
//-----------------------------------------------------------------------------
    for (int k=0; k<fNID; k++) {
      tp->fIDWord[k]  = tp->fTrackID[k]->IDWord(track);
    }
//-----------------------------------------------------------------------------
// PID likelihoods
//-----------------------------------------------------------------------------
    dat.fDt   = tp->fDt;
    dat.fEp   = tp->fEp;
    dat.fPath = tp->fPath;
      
    xs = track->XSlope();

    track->fEleLogLHCal = tp->fLogLH->LogLHCal(&dat,11);
    track->fMuoLogLHCal = tp->fLogLH->LogLHCal(&dat,13);

    double llhr_cal = track->fEleLogLHCal-track->fMuoLogLHCal;

    int id_word = tp->fIDWord[0];

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

    track->fLogLHRXs = tp->fLogLH->LogLHRXs(xs); 
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

}
