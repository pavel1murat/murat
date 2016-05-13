//////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/obj/TDisk.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------

#include "ana/TTrackAnaModuleA.hh"

//-----------------------------------------------------------------------------
void TTrackAnaModuleA::BookCaloHistograms(HistBase_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];

  CaloHist_t* hist = (CaloHist_t*) Hist;

  HBook1F(hist->fDiskID ,"disk_id",Form("%s: Disk ID"       ,Folder), 10, 0,  10,Folder);

  for (int i=0; i<kNDisks; i++) {
    HBook1F(hist->fEnergy  [i],Form("energy_%i",i),Form("%s: Hit Energy[%i]",Folder,i),200, 0, 100,Folder);
    HBook1F(hist->fTime    [i],Form("time_%i"  ,i),Form("%s: Hit time  [%i]",Folder,i),200, 0,2000,Folder);
    HBook1F(hist->fNHits   [i],Form("nhits_%i" ,i),Form("%s: NHits     [%i]",Folder,i), 50, 0,  50,Folder);
    HBook1F(hist->fRadius  [i],Form("r_%i"     ,i),Form("%s: Radius    [%i]",Folder,i),100, 0,1000,Folder);
    HBook1F(hist->fRadiusWE[i],Form("rwe_%i"   ,i),Form("%s: RadiusWE  [%i]",Folder,i),100, 0,1000,Folder);

    HBook1F(hist->fE700    [i],Form("e700_%i",i),Form("%s: Hit Energy[%i] (T > 700ns)",Folder,i),200, 0, 100,Folder);
    HBook1F(hist->fT700    [i],Form("t700_%i",i),Form("%s: Hit time  [%i] (T > 700ns)",Folder,i),200, 0,2000,Folder);
    HBook1F(hist->fN700    [i],Form("n700_%i",i),Form("%s: NHits     [%i] (T > 700ns)",Folder,i), 50, 0,  50,Folder);

    HBook1F(hist->fR700  [i],Form("r700_%i"  ,i),Form("%s: Radius (T>700) [%i]",Folder,i),100, 0,1000,Folder);
    HBook1F(hist->fRWE700[i],Form("rwe700_%i",i),Form("%s: Radius*E(T>700)[%i]",Folder,i),100, 0,1000,Folder);
  }
}


//-----------------------------------------------------------------------------
void TTrackAnaModuleA::BookClusterHistograms(HistBase_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  ClusterHist_t* hist = (ClusterHist_t*) Hist;

  HBook1F(hist->fDiskID ,"disk_id",Form("%s: Disk ID"       ,Folder), 10, 0,  10,Folder);
  HBook1F(hist->fEnergy ,"energy" ,Form("%s: Cluster Energy",Folder),500, 0, 250,Folder);
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
void TTrackAnaModuleA::BookGenpHistograms(HistBase_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  GenpHist_t* hist = (GenpHist_t*) Hist;

  HBook1F(hist->fP      ,"p"       ,Form("%s: Momentum"     ,Folder),1000,     0, 200,Folder);
  HBook1F(hist->fPdgCode[0],"pdg_code_0",Form("%s: PDG Code[0]"     ,Folder),200, -100, 100,Folder);
  HBook1F(hist->fPdgCode[1],"pdg_code_1",Form("%s: PDG Code[1]"     ,Folder),500, -2500, 2500,Folder);
  HBook1F(hist->fGenID  ,"gen_id"  ,Form("%s: Generator ID" ,Folder), 100,     0, 100,Folder);
  HBook1F(hist->fZ0     ,"z0"      ,Form("%s: Z0"           ,Folder), 500,  5400, 6400,Folder);
  HBook1F(hist->fT0     ,"t0"      ,Form("%s: T0"           ,Folder), 200,     0, 2000,Folder);
  HBook1F(hist->fR0     ,"r"       ,Form("%s: R0"           ,Folder), 100,     0,  100,Folder);
  HBook1F(hist->fCosTh  ,"cos_th"  ,Form("%s: Cos(Theta)"   ,Folder), 200,   -1.,   1.,Folder);
}

//-----------------------------------------------------------------------------
void TTrackAnaModuleA::BookTrackHistograms(HistBase_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  TrackHist_t* hist = (TrackHist_t*) Hist;

  HBook1F(hist->fP[0]       ,"p"        ,Form("%s: Track P(Z1)"       ,Folder), 800,  80  ,120. ,Folder);
  HBook1F(hist->fP[1]       ,"p_1"      ,Form("%s: Track P(total)[1]" ,Folder), 100, 100  ,105. ,Folder);
  HBook1F(hist->fP[2]       ,"p_2"      ,Form("%s: Track P(total)[1]" ,Folder),2000,   0  ,200. ,Folder);
  HBook1F(hist->fP0         ,"p0"       ,Form("%s: Track P(Z0)"       ,Folder),1000,   0  ,200. ,Folder);
  HBook1F(hist->fP2         ,"p2"       ,Form("%s: Track P(z=-1540)"  ,Folder),1000,   0  ,200. ,Folder);
  HBook1D(hist->fPDio       ,"pdio"     ,Form("%s: Track P(DIO WT)"   ,Folder), 800,  80  ,120. ,Folder);
  hist->fPDio->Sumw2(kTRUE);
//-----------------------------------------------------------------------------
// luminosity-weighted distributions for signal and background
//-----------------------------------------------------------------------------
  HBook1D(hist->fPlw        ,"plw"     ,Form("%s: Track P(Lumi-WT)"   ,Folder), 800,  80  ,120. ,Folder);
  hist->fPlw->Sumw2(kTRUE);
  HBook1D(hist->fPDiolw     ,"pdiolw"  ,Form("%s: Trk P WT(Lumi+DIO)" ,Folder), 800,  80  ,120. ,Folder);
  hist->fPDiolw->Sumw2(kTRUE);
//-----------------------------------------------------------------------------
  HBook1F(hist->fFitMomErr  ,"momerr"   ,Form("%s: Track FitMomError" ,Folder), 200,   0  ,  1. ,Folder);
  HBook1F(hist->fPFront     ,"pf"       ,Form("%s: Track P(front)   " ,Folder), 400,  90  ,110. ,Folder);
  HBook1F(hist->fDpFront    ,"dpf"      ,Form("%s: Track P-P(front) " ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(hist->fXDpF       ,"xdpf"     ,Form("%s: DpF/momErr"        ,Folder),1000, -50. , 50. ,Folder);
  HBook1F(hist->fDpFDio     ,"dpfdio"   ,Form("%s: Track DpF(DIO Wt)" ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(hist->fXDpFDio    ,"xdpfdio"  ,Form("%s: DpF/momErr DioWt"  ,Folder),1000, -50. , 50. ,Folder);
  HBook1F(hist->fDpFront0   ,"dp0f"     ,Form("%s: Track P0-P(front)" ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(hist->fDpFront2   ,"dp2f"     ,Form("%s: Track P2-P(front)" ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(hist->fPStOut     ,"pstout"   ,Form("%s: Track P(ST_Out)  " ,Folder), 400,  90. ,110. ,Folder);
  HBook1F(hist->fDpFSt      ,"dpfst"    ,Form("%s: Track Pf-Psto"     ,Folder), 200,  -5  ,  5. ,Folder);
  HBook2F(hist->fDpFVsZ1    ,"dpf_vs_z1",Form("%s: Track DPF Vs Z1"   ,Folder), 200, -2000.,0,200,-5.,5,Folder);

  HBook1F(hist->fPt         ,"pt"       ,Form("%s: Track Pt"          ,Folder), 600, 75,95,Folder);
  HBook1F(hist->fCosTh      ,"costh"    ,Form("%s: Track cos(theta)"  ,Folder), 100,-1,1,Folder);
  HBook1F(hist->fChi2       ,"chi2"     ,Form("%s: Track chi2 total"  ,Folder), 200, 0,200,Folder);
  HBook1F(hist->fChi2Dof    ,"chi2d"    ,Form("%s: track chi2/N(dof)" ,Folder), 500, 0, 10,Folder);

  HBook1F(hist->fNActive    ,"nactv"    ,Form("%s: N(active)"         ,Folder), 200, 0,200,Folder);
  HBook1F(hist->fNaFract    ,"nafr"     ,Form("%s: N(active fraction)",Folder), 110, 0.5,1.05,Folder);
  HBook1F(hist->fNWrong     ,"nwrng"    ,Form("%s: N(wrong drift sgn)",Folder), 100, 0,100,Folder);
  HBook1F(hist->fNDoublets  ,"nd"       ,Form("%s: N(doublets)"       ,Folder),  50, 0, 50,Folder);
  HBook1F(hist->fNadOverNd  ,"nad_nd"   ,Form("%s: Nad/N(doublets)"   ,Folder), 110, 0,1.1,Folder);
  HBook1F(hist->fNSSD       ,"nssd"     ,Form("%s: N(SS doublets)"    ,Folder),  50, 0, 50,Folder);
  HBook1F(hist->fNOSD       ,"nosd"     ,Form("%s: N(OS doublets)"    ,Folder),  50, 0, 50,Folder);
  HBook1F(hist->fNdOverNa   ,"nd_na"    ,Form("%s: NDoublets/Nactive" ,Folder), 100, 0,0.5,Folder);
  HBook1F(hist->fNssdOverNa ,"nssd_na"  ,Form("%s: NSSD/Nactive"      ,Folder), 100, 0,0.5,Folder);
  HBook1F(hist->fNosdOverNa ,"nosd_na"  ,Form("%s: NOSD/Nactive"      ,Folder), 100, 0,0.5,Folder);
  HBook1F(hist->fNZeroAmb   ,"nza"      ,Form("%s: N (Iamb = 0) hits" ,Folder), 100, 0,100,Folder);
  HBook1F(hist->fNzaOverNa  ,"nza_na"   ,Form("%s: NZeroAmb/Nactive"  ,Folder), 100, 0,  1,Folder);
  HBook1F(hist->fNMatActive ,"nma"      ,Form("%s: N (Mat Active"     ,Folder), 100, 0,100,Folder);
  HBook1F(hist->fNmaOverNa  ,"nma_na"   ,Form("%s: NMatActive/Nactive",Folder), 200, 0,   2,Folder);
  HBook1F(hist->fNBend      ,"nbend"    ,Form("%s: Nbend"             ,Folder), 100, 0,1000,Folder);

  HBook1F(hist->fT0         ,"t0"       ,Form("%s: track T0"          ,Folder), 200, 0,2000,Folder);
  HBook1F(hist->fT0Err      ,"t0err"    ,Form("%s: track T0Err"       ,Folder), 100, 0,  10,Folder);
  HBook1F(hist->fQ          ,"q"        ,Form("%s: track Q"           ,Folder),   4,-2,   2,Folder);
  HBook1F(hist->fFitCons[0] ,"fcon"     ,Form("%s: track fit cons [0]",Folder), 200, 0,   1,Folder);
  HBook1F(hist->fFitCons[1] ,"fcon1"    ,Form("%s: track fit cons [1]",Folder), 1000, 0,   0.1,Folder);
  HBook1F(hist->fD0         ,"d0"       ,Form("%s: track D0      "    ,Folder), 200,-200, 200,Folder);
  HBook1F(hist->fZ0         ,"z0"       ,Form("%s: track Z0      "    ,Folder), 200,-2000,2000,Folder);
  HBook1F(hist->fTanDip     ,"tdip"     ,Form("%s: track tan(dip)"    ,Folder), 200, 0.0 ,2.0,Folder);
  HBook1F(hist->fRMax       ,"rmax"     ,Form("%s: track R(max)  "    ,Folder), 200, 0., 1000,Folder);
  HBook1F(hist->fDtZ0       ,"dtz0"     ,Form("%s: DT(Z0), MC"        ,Folder), 200, -10.0 ,10.0,Folder);

  HBook1F(hist->fResid      ,"resid"    ,Form("%s: hit residuals"     ,Folder), 500,-0.5 ,0.5,Folder);
  HBook1F(hist->fAlgMask    ,"alg"      ,Form("%s: algorithm mask"    ,Folder),  10,  0, 10,Folder);

  HBook1F(hist->fChi2Tcm  ,"chi2tcm"  ,Form("%s: chi2(t-c match)"   ,Folder), 250,  0  ,250 ,Folder);
  HBook1F(hist->fChi2XY     ,"chi2xy"   ,Form("%s: chi2(t-c match) XY",Folder), 300,-50  ,250 ,Folder);
  HBook1F(hist->fChi2T      ,"chi2t"    ,Form("%s: chi2(t-c match) T" ,Folder), 250,  0  ,250 ,Folder);

  HBook1F(hist->fDt         ,"dt"       ,Form("%s: track delta(T)"    ,Folder), 400,-20  ,20 ,Folder);
  HBook1F(hist->fDx         ,"dx"       ,Form("%s: track delta(X)"    ,Folder), 200,-500 ,500,Folder);
  HBook1F(hist->fDy         ,"dy"       ,Form("%s: track delta(Y)"    ,Folder), 200,-500 ,500,Folder);
  HBook1F(hist->fDz         ,"dz"       ,Form("%s: track delta(Z)"    ,Folder), 200,-250 ,250,Folder);
  HBook1F(hist->fDu         ,"du"       ,Form("%s: track-cluster DU)" ,Folder), 250,-250 ,250,Folder);
  HBook1F(hist->fDv         ,"dv"       ,Form("%s: track-cluster DV)" ,Folder), 200,-100 ,100,Folder);
  HBook2F(hist->fDvVsDu     ,"dv_vs_du" ,Form("%s: Track Dv Vs Du"    ,Folder), 100, -250,250,100,-100.,100,Folder);
  HBook1F(hist->fPath       ,"path"     ,Form("%s: track sdisk"       ,Folder),  50,   0 ,500,Folder);
  HBook2F(hist->fDuVsPath   ,"du_vs_path",Form("%s: Track Du Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(hist->fDvVsPath   ,"dv_vs_path",Form("%s: T-C  Dv Vs Path"  ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(hist->fDtVsPath   ,"dt_vs_path",Form("%s: T-C DT Vs Path"   ,Folder),  50,   0 ,500,100,  -5.,  5.,Folder);
  HBook2F(hist->fDuVsTDip   ,"du_vs_tdip",Form("%s: Track Du Vs TDip" ,Folder), 100, 0.5 ,1.5,200,-200.,200.,Folder);
  HBook2F(hist->fDvVsTDip   ,"dv_vs_tdip",Form("%s: Track Dv Vs TDip" ,Folder), 100, 0.5 ,1.5,200,-200.,200.,Folder);

  HBook1F(hist->fZ1         ,"z1"       ,Form("%s: track Z1      "    ,Folder), 200,-2000,2000,Folder);
  HBook1F(hist->fNClusters  ,"ncl"      ,Form("%s: track N(clusters)" ,Folder),  10, 0   , 10,Folder);
  HBook1F(hist->fDiskID     ,"disk_id"  ,Form("%s: track disk ID"     ,Folder),  10,-5   ,  5,Folder);
  HBook1F(hist->fXCal       ,"xcal"     ,Form("%s: track XCal"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(hist->fYCal       ,"ycal"     ,Form("%s: track YCal"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(hist->fZCal       ,"zcal"     ,Form("%s: track ZCal"        ,Folder), 200, 1500,3500,Folder);
  HBook1F(hist->fXTrk       ,"xtrk"     ,Form("%s: track XTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(hist->fYTrk       ,"ytrk"     ,Form("%s: track YTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(hist->fRTrk       ,"rtrk"     ,Form("%s: track RTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(hist->fZTrk       ,"ztrk"     ,Form("%s: track ZTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(hist->fECl        ,"ecl"      ,Form("%s: cluster E"         ,Folder), 300, 0   ,150,Folder);
  HBook1F(hist->fEClEKin    ,"ecl_ekin" ,Form("%s: cluster E/Ekin(mu)",Folder), 200, 0   ,2,Folder);
  HBook1F(hist->fEp         ,"ep"       ,Form("%s: track E/P"         ,Folder), 300, 0   ,1.5,Folder);
  HBook2F(hist->fEpVsPath   ,"ep_vs_path",Form("%s: E/P Vs Path"      ,Folder),  50,   0 ,500,150,  0.,  1.5,Folder);
  HBook2F(hist->fNHVsStation,"nh_vs_st" ,Form("%s: N(hits) Vs Station",Folder),  40, 0,40,10,-0.5,9.5,Folder);
  HBook2F(hist->fNHVsNSt    ,"nh_vs_nst",Form("%s: N(hits) Vs NSt"    ,Folder),  10,-0.5,9.5,40,-0.5,39.5,Folder);

  HBook1F(hist->fRSlope     ,"rslope"   ,Form("%s: Res Slope"         ,Folder), 200,-20 , 20,Folder);
  HBook1F(hist->fXSlope     ,"xslope"   ,Form("%s: Res/Sig Slope"     ,Folder), 200,-20 , 20,Folder);

  HBook2F(hist->fEpVsDt     ,"ep_vs_dt" ,Form("%s: E/P vs Dt"         ,Folder), 200, -10, 10,150,0.,1.5,Folder);
  HBook1F(hist->fEleLogLHCal,"ele_llh_c",Form("%s: ELE Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
  HBook1F(hist->fMuoLogLHCal,"muo_llh_c",Form("%s: MUO Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
  HBook1F(hist->fLogLHRCal  ,"llhr_cal" ,Form("%s: LogLH(e/m) Cal"    ,Folder), 200,-100,100,Folder);
  HBook1F(hist->fLogLHRDeDx ,"llhr_dedx",Form("%s: LogLH(e/m) De/Dx"  ,Folder), 200,-20 , 20,Folder);
  HBook1F(hist->fLogLHRXs   ,"llhr_xs"  ,Form("%s: LogLH(e/m) XSlope" ,Folder), 200,-20 , 20,Folder);
  HBook1F(hist->fLogLHRTrk  ,"llhr_trk" ,Form("%s: LogLH(e/m) Trk"    ,Folder), 200,-20 , 20,Folder);
  HBook1F(hist->fLogLHR     ,"llhr"     ,Form("%s: LogLH(e/m)"        ,Folder), 200,-100 ,100,Folder);

  HBook1F(hist->fPdgCode    ,"pdg"      ,Form("%s: track PDG code"    ,Folder), 100,-50,50,Folder);
  HBook1F(hist->fFrGH       ,"fgh"      ,Form("%s: Fraction Goog Hits",Folder), 100, 0,1,Folder);

  HBook2F(hist->fNEPlVsNHPl ,"nep_vs_nhp",Form("%s: Track NEXP vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
  HBook2F(hist->fNDPlVsNHPl ,"ndp_vs_nhp",Form("%s: Track NDIF vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
  HBook2F(hist->fChi2dVsNDPl,"chi2d_vs_ndp",Form("%s: Track Chi2/Dof vs NDP",Folder), 30, 0,30,100,0.,10,Folder);
  HBook2F(hist->fDpFVsNDPl  ,"dpf_vs_ndp"  ,Form("%s: Track DpF vs NDP",Folder)     , 30, 0,30,100,-5,5,Folder);

  HBook1F(hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
  HBook1F(hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);

  HBook1F(hist->fSinTC  ,"sintc"  ,Form("%s: sin(T^C)"      ,Folder),200, -1, 1,Folder);
  HBook1F(hist->fDrTC   ,"drtc"   ,Form("%s: DeltaR(T^C)"   ,Folder),200, -200, 200,Folder);
  HBook1F(hist->fSInt   ,"sint"   ,Form("%s: SInt"          ,Folder),200, -500, 500,Folder);
  HBook1F(hist->fDaveTrkQual,"dtqual",Form("%s:DaveTrkQual" ,Folder),500, -2.5, 2.5,Folder);
  HBook1F(hist->fNMcStrawHits,"nmc_hits",Form("%s:N(MC Straw Hits by parent",Folder),150, 0, 150,Folder);
}

//-----------------------------------------------------------------------------
void TTrackAnaModuleA::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  char name [200];
  //  char title[200];

  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1D(hist->fLumWt     ,"lumwt"    ,Form("%s: Luminosity Weight"               ,Folder),200, 0,10,Folder);
  hist->fLumWt->Sumw2(kTRUE);

  HBook1F(hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1, 1,Folder);
  HBook1F(hist->fEleMom    ,"ce_mom"   ,Form("%s: Conversion Electron Momentum"    ,Folder),1000,  0,200,Folder);
  HBook1D(hist->fDioMom    ,"dio_mom"  ,Form("%s: DIO momentum"                    ,Folder),1000, 50,150,Folder);
  HBook1F(hist->fRv         ,"rv"      ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 1000,Folder);
  HBook1F(hist->fZv         ,"zv"      ,Form("%s: Z(Vertex)"                       ,Folder), 300, 0,15000,Folder);
  HBook1F(hist->fNClusters ,"ncl"      ,Form("%s: Number of Reconstructed Clusters",Folder),200,0,200,Folder);
  HBook1F(hist->fNTracks   ,"ntrk"     ,Form("%s: Number of Reconstructed Tracks"  ,Folder),100,0,100,Folder);
  HBook1F(hist->fNStrawHits[0],"nsh_0" ,Form("%s: Number of Straw Hits [0]"        ,Folder),250,0,250,Folder);
  HBook1F(hist->fNStrawHits[1],"nsh_1" ,Form("%s: Number of Straw Hits [1]"        ,Folder),400,0,8000,Folder);
  HBook1F(hist->fNGoodSH   ,"nsh50"    ,Form("%s: N(SH) +/-50"                     ,Folder),300,0,1500,Folder);
  HBook1F(hist->fDtClT     ,"dt_clt"   ,Form("%s: DT(cluster-track)"               ,Folder),100,-100,100,Folder);
  HBook1F(hist->fDtClS     ,"dt_cls"   ,Form("%s: DT(cluster-straw hit)"           ,Folder),200,-200,200,Folder);
  HBook1F(hist->fSHTime    ,"shtime"   ,Form("%s: Straw Hit Time"                  ,Folder),400,0,2000,Folder);
  HBook1F(hist->fEMax      ,"emax"     ,Form("%s: Max cluster energy"              ,Folder),150,0,150,Folder);
  HBook1F(hist->fNHyp      ,"nhyp"     ,Form("%s: N(fit hypotheses)"               ,Folder),5,0,5,Folder);
  HBook1F(hist->fBestHyp[0],"bfh0"     ,Form("%s: Best Fit Hyp[0](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(hist->fBestHyp[1],"bfh1"     ,Form("%s: Best Fit Hyp[1](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(hist->fNGenp     ,"ngenp"    ,Form("%s: N(Gen Particles)"                ,Folder),500,0,500,Folder);

  for (int i=0; i<2; i++) {
    sprintf(name,"ncch_%i",i);
    HBook1F(hist->fNCaloCrystalHits[i],name,Form("%s: N(calo crystal hits) [%i]",Folder,i),500,0,1000,Folder);
    sprintf(name,"ncch_vs_disk_%i",i);
    HBook2F(hist->fNCaloHitsVsDisk[i],name,Form("%s: N(calo crystal hits) vs disk[%i]",Folder,i),4,0,4,200,0,200,Folder);
    sprintf(name,"ncch_vs_row_%i",i);
    HBook2F(hist->fNCaloHitsVsRow[i],name,Form("%s: N(calo crystal hits) vs row [%i]",Folder,i),20,0,20,200,0,200,Folder);
    sprintf(name,"ncch_vs_col_%i",i);
    HBook2F(hist->fNCaloHitsVsCol[i],name,Form("%s: N(calo crystal hits) vs col [%i]",Folder,i),50,0,50,200,0,200,Folder);
  }

  for (int i=0; i<kNDisks; i++) {
    HBook1F(hist->fETot        [i],Form("etot_%i"    ,i),Form("%s: Etot[%i]",Folder,i), 300, 0,150,Folder);
    HBook2F(hist->fECrVsR      [i],Form("ecr_vs_r_%i",i),Form("%s: E Cr Vs R [%i]"    ,Folder,i), 100, 0,1000,500,0,100,Folder);
    HBook2F(hist->fNCrVsR      [i],Form("ncr_vs_r_%i",i),Form("%s: N Cr Vs R [%i]"    ,Folder,i), 100, 0,1000,100,0,100,Folder);

    HBook2F(hist->fNCrystalHitsVsR[i],Form("ncrh_vs_r_%i",i),Form("%s: N Crystal Hits[%i] vs R",Folder,i), 100, 0, 1000,100,0,100,Folder);
    HBook2F(hist->fNHitCrystalsVsR[i],Form("nhcr_vs_r_%i",i),Form("%s: N Hit Crystals[%i] vs R",Folder,i), 100, 0, 1000,100,0,100,Folder);
  }

  HBook1F(hist->fNHitCrystalsTot,"nhcr_tot",Form("%s: NHit Crystals Tot",Folder), 100, 0,100,Folder);
  HBook1F(hist->fECal,"ecal",Form("%s: E(cal), sum over both disks",Folder), 500, 0,250,Folder);
  HBook1F(hist->fECalOverEKin,"ec_over_ek",Form("%s: E(cal)/E(kin)",Folder), 200, 0,2,Folder);
  HBook1F(hist->fInstLumi    ,"lumi"      ,Form("%s: N(protons/mubunch)",Folder), 2000, 0,1.e8,Folder);
}

//-----------------------------------------------------------------------------
void TTrackAnaModuleA::BookSimpHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  SimpHist_t* hist = (SimpHist_t*) Hist;

  HBook1F(hist->fPdgCode        ,"pdg"   ,Form("%s: PDG code"                     ,Folder),200,-100,100,Folder);
  HBook1F(hist->fNStrawHits     ,"nsth"  ,Form("%s: n straw hits"                 ,Folder),200,   0,200,Folder);
  HBook1F(hist->fMomTargetEnd   ,"ptarg" ,Form("%s: CE mom after Stopping Target" ,Folder),400,  90,110,Folder);
  HBook1F(hist->fMomTrackerFront,"pfront",Form("%s: CE mom at the Tracker Front"  ,Folder),400,  90,110,Folder);
}
//_____________________________________________________________________________
void TTrackAnaModuleA::BookHistograms() {

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
  HBook1F(fHist.fCrystalR[0],"rc_0"     ,Form("disk [0] crystal radius"),100,0,1000,"Hist");
  HBook1F(fHist.fCrystalR[1],"rc_1"     ,Form("disk [1] crystal radius"),100,0,1000,"Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;	        // events with a reconstructed track
  book_event_histset[ 2] = 1;	        // events without reconstructed tracks
  book_event_histset[ 3] = 1;	        // events with a reconstructed cluster
  book_event_histset[ 4] = 1;	        // events without reconstructed clusters
  book_event_histset[ 5] = 1;	        // events w/o reconstructed tracks, |costh|<0.4
  book_event_histset[ 6] = 1;	        // events with tracks passing "Set C" cuts
  book_event_histset[ 7] = 1;	        // events with E(cluster) > 60 MeV
  book_event_histset[ 8] = 1;	        // events with the highest energy cluster on the 1st disk
  book_event_histset[ 9] = 1;	        // events with the highest energy cluster on the 2nd disk
  book_event_histset[10] = 1;	        // events with SetC tracks 103.5 < p < 105.

  book_event_histset[11] = 1;	        // Track ID ladder : 11-20
  book_event_histset[12] = 1;	        // 
  book_event_histset[13] = 1;	        // 
  book_event_histset[14] = 1;	        // 
  book_event_histset[15] = 1;	        // 
  book_event_histset[16] = 1;	        // 
  book_event_histset[17] = 1;	        // 
  book_event_histset[18] = 1;	        // 
  book_event_histset[19] = 1;	        // 
  book_event_histset[20] = 1;	        // 

  book_event_histset[21] = 0;	        // 
  book_event_histset[22] = 0;	        // 
  book_event_histset[23] = 0;	        // 

					// TrkPatRec tracks
  book_event_histset[24] = 1;	        // events with at least one reco track
  book_event_histset[25] = 1;	        // 
  book_event_histset[26] = 1;	        // 
  book_event_histset[27] = 1;	        // 
  book_event_histset[28] = 1;	        // 

  book_event_histset[31] = 1;		// EFFICIENCY - Dave's definition
  book_event_histset[32] = 1;
  book_event_histset[33] = 1;
  

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

  book_track_histset[  0] = 1;		// all tracks e-
  book_track_histset[  1] = 1;		// all tracks e- passing Set C cuts 
  book_track_histset[  2] = 1;		// all tracks e- passing Set C cuts, events with clusters 
  book_track_histset[  3] = 1;		// all tracks e- passing Set C cuts, events w/o  clusters
  book_track_histset[  4] = 1;		// all tracks e- passing Set C cuts, events with clusters, no closest
  book_track_histset[  5] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
  book_track_histset[  6] = 1;		// all tracks e-  E/P > 0.4
  book_track_histset[  7] = 1;		// all tracks e- passing Set C cuts, E/P > 0.4
  book_track_histset[  8] = 1;		// Set C tracks e- , |xslope| < 3.
  book_track_histset[  9] = 1;		// all  tracks in the event when there is no EM clusters E > 60 MeV
  book_track_histset[ 10] = 1;		// Set C tracks in the event when there is no EM clusters E > 60 MeV
  book_track_histset[ 11] = 1;		// all tracks with P > 103.5
  book_track_histset[ 12] = 1;		// tracks with fcons < 1.e-4
  book_track_histset[ 13] = 1;		// SetC tracks with 100 <= P < 110 
  book_track_histset[ 14] = 1;		// tracks with fcons < 1.e-2
  book_track_histset[ 15] = 1;		// tracks intersecting the 1st disk
  book_track_histset[ 16] = 1;		// tracks intersecting the 2nd disk
  book_track_histset[ 17] = 1;		// tracks with no calorimeter intersections
  book_track_histset[ 18] = 1;		// Set C tracks, T0 > 700
  book_track_histset[ 19] = 1;		// Set C tracks with E/P > 0

  book_track_histset[ 20] = 1;		// tracks with Nhits >= 20
  book_track_histset[ 21] = 1;		// tracks with Nhits >= 20 and chi/Ndof < 3
  book_track_histset[ 22] = 1;		// Set C tracks with E/P > 0 and chi2(match) < 100
  book_track_histset[ 23] = 1;		// Set C tracks with E/P > 0 and chi2(match) < 100 and LLHR(cal) > 0 (interesting for muons)
  book_track_histset[ 24] = 1;		// Set C tracks with E/P > 0 and chi2(match) < 100 and LLHR(cal) < 0 (interesting for electrons)

  book_track_histset[ 25] = 1;		// Set C tracks, 100<p<110, 0<E/P<1.15, chi2tcm<100, -5<DT<8, 
  book_track_histset[ 26] = 1;		// [25] + LLHR_CAL > 0 - interesting for muons
  book_track_histset[ 27] = 1;		// [25] + LLHR_CAL < 0 - interesting for electrons
  book_track_histset[ 28] = 1;		// Set C tracks, E/P > 1.1
  book_track_histset[ 29] = 1;		// Set C tracks 100 < P < 110 , 0 < E/P < 1.15 *precursor for TRK_32*

  book_track_histset[ 30] = 1;		// tracks with Nhits >= 25
  book_track_histset[ 31] = 1;		// tracks with Nhits >= 25 and chi/Ndof < 3
  book_track_histset[ 32] = 1;		// Set C tracks 100 < P < 110 , 0 < E/P < 1.15, chi2tcm<100

  book_track_histset[ 33] = 1;		// DaveTrkQual tracks Q > 0.4
  book_track_histset[ 34] = 1;		// DaveTrkQual tracks Q > 0.1
  
  book_track_histset[ 40] = 1;		// all tracks, alg_mask = 1
  book_track_histset[ 41] = 1;		// Set "C" tracks, alg_mask = 1
  book_track_histset[ 42] = 1;		// Set "C" tracks, alg_mask = 1, T > 700

  book_track_histset[ 50] = 1;		// all tracks, alg_mask = 2
  book_track_histset[ 51] = 1;		// Set "C" tracks, alg_mask = 2
  book_track_histset[ 52] = 1;		// Set "C" tracks, alg_mask = 2, T > 700

  book_track_histset[ 60] = 1;		// all tracks, alg_mask = 3
  book_track_histset[ 61] = 1;		// Set "C" tracks, alg_mask = 3
  book_track_histset[ 62] = 1;		// Set "C" tracks, alg_mask = 3, T > 700
  book_track_histset[ 63] = 1;          // best ID, best = 0 (TrkPatRec)
  book_track_histset[ 64] = 1;          // best ID, best = 1 (CalPatRec)
  book_track_histset[ 65] = 0;
  book_track_histset[ 66] = 0;
  book_track_histset[ 67] = 0;
  book_track_histset[ 68] = 0;
  book_track_histset[ 69] = 0;
  book_track_histset[ 70] = 0;

  book_track_histset[ 71] = 1;		// Set "C" tracks   103.5 < P < 105 : interesting for DIO

  book_track_histset[ 81] = 1;
  book_track_histset[ 82] = 1;
  book_track_histset[ 83] = 1;
  book_track_histset[ 84] = 1;
  book_track_histset[ 85] = 1;
  book_track_histset[ 86] = 1;
  

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
// book cluster histograms
//-----------------------------------------------------------------------------
  int book_cluster_histset[kNClusterHistSets];
  for (int i=0; i<kNClusterHistSets; i++) book_cluster_histset[i] = 0;

  book_cluster_histset[0] = 1;		// all clusters
  book_cluster_histset[1] = 1;		// clusters in events with the reconstructed e-
  book_cluster_histset[2] = 1;		// clusters in events with the track passing SetC cuts
  book_cluster_histset[3] = 1;		// clusters in events w/track passing SetC cuts and |dt|<2.5ns 
  book_cluster_histset[4] = 1;		// clusters > 10 MeV
  book_cluster_histset[5] = 1;		// clusters > 60 MeV
  book_cluster_histset[6] = 1;		// clusters disk#0
  book_cluster_histset[7] = 1;		// clusters disk#1

  for (int i=0; i<kNClusterHistSets; i++) {
    if (book_cluster_histset[i] != 0) {
      sprintf(folder_name,"cls_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCluster[i] = new ClusterHist_t;
      BookClusterHistograms(fHist.fCluster[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book calorimeter histograms
//-----------------------------------------------------------------------------
  int book_calo_histset[kNCaloHistSets];
  for (int i=0; i<kNCaloHistSets; i++) book_calo_histset[i] = 0;

  book_calo_histset[0] = 1;		// all crystals
  book_calo_histset[1] = 1;		// all crystals, e > 0
  book_calo_histset[2] = 1;		// all crystals, e > 0.1
  book_calo_histset[3] = 1;		// all crystals, e > 1.0

  for (int i=0; i<kNCaloHistSets; i++) {
    if (book_calo_histset[i] != 0) {
      sprintf(folder_name,"cal_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCalo[i] = new CaloHist_t;
      BookCaloHistograms(fHist.fCalo[i],Form("Hist/%s",folder_name));
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

}

