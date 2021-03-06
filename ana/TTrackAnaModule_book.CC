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
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------

#include "ana/TTrackAnaModule.hh"

// //-----------------------------------------------------------------------------
// void TTrackAnaModule::BookGenpHistograms(GenpHist_t* Hist, const char* Folder) {
// //   char name [200];
// //   char title[200];

//   HBook1F(Hist->fP      ,"p"       ,Form("%s: Momentum"     ,Folder),1000,     0, 200,Folder);
//   HBook1F(Hist->fPdgCode[0],"pdg_code_0",Form("%s: PDG Code[0]"     ,Folder),200, -100, 100,Folder);
//   HBook1F(Hist->fPdgCode[1],"pdg_code_1",Form("%s: PDG Code[1]"     ,Folder),500, -2500, 2500,Folder);
//   HBook1F(Hist->fGenID  ,"gen_id"  ,Form("%s: Generator ID" ,Folder), 100,     0, 100,Folder);
//   HBook1F(Hist->fZ0     ,"z0"      ,Form("%s: Z0"           ,Folder), 500,  5400, 6400,Folder);
//   HBook1F(Hist->fT0     ,"t0"      ,Form("%s: T0"           ,Folder), 200,     0, 2000,Folder);
//   HBook1F(Hist->fR0     ,"r"       ,Form("%s: R0"           ,Folder), 100,     0,  100,Folder);
//   HBook1F(Hist->fCosTh  ,"cos_th"  ,Form("%s: Cos(Theta)"   ,Folder), 200,   -1.,   1.,Folder);
// }

// //-----------------------------------------------------------------------------
// void TTrackAnaModule::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
// //   char name [200];
// //   char title[200];

//   HBook1F(Hist->fP[0]       ,"p"        ,Form("%s: Track P(Z1)"       ,Folder), 800,  80  ,120. ,Folder);
//   HBook1F(Hist->fP[1]       ,"p_1"      ,Form("%s: Track P(total)[1]" ,Folder), 100, 100  ,105. ,Folder);
//   HBook1F(Hist->fP[2]       ,"p_2"      ,Form("%s: Track P(total)[1]" ,Folder),1000,   0  ,500. ,Folder);
//   HBook1F(Hist->fP0         ,"p0"       ,Form("%s: Track P(Z0)"       ,Folder),1000,   0  ,200. ,Folder);
//   HBook1F(Hist->fP2         ,"p2"       ,Form("%s: Track P(z=-1540)"  ,Folder),1000,   0  ,200. ,Folder);
//   HBook1D(Hist->fPDio       ,"pdio"     ,Form("%s: Track P(DIO WT)"   ,Folder), 800,  80  ,120. ,Folder);
//   Hist->fPDio->Sumw2(kTRUE);
// //-----------------------------------------------------------------------------
// // luminosity-weighted distributions for signal and background
// //-----------------------------------------------------------------------------
//   HBook1D(Hist->fPlw        ,"plw"     ,Form("%s: Track P(Lumi-WT)"   ,Folder), 800,  80  ,120. ,Folder);
//   Hist->fPlw->Sumw2(kTRUE);
//   HBook1D(Hist->fPDiolw     ,"pdiolw"  ,Form("%s: Trk P WT(Lumi+DIO)" ,Folder), 800,  80  ,120. ,Folder);
//   Hist->fPDiolw->Sumw2(kTRUE);
// //-----------------------------------------------------------------------------
//   HBook1F(Hist->fFitMomErr  ,"momerr"   ,Form("%s: Track FitMomError" ,Folder), 200,   0  ,  1. ,Folder);
//   HBook1F(Hist->fPFront     ,"pf"       ,Form("%s: Track P(front)   " ,Folder), 400,  90  ,110. ,Folder);
//   HBook1F(Hist->fDpFront    ,"dpf"      ,Form("%s: Track P-P(front) " ,Folder), 200,  -5. ,  5. ,Folder);
//   HBook1F(Hist->fXDpF       ,"xdpf"     ,Form("%s: DpF/momErr"        ,Folder),1000, -50. , 50. ,Folder);
//   HBook1F(Hist->fDpFDio     ,"dpfdio"   ,Form("%s: Track DpF(DIO Wt)" ,Folder), 200,  -5. ,  5. ,Folder);
//   HBook1F(Hist->fXDpFDio    ,"xdpfdio"  ,Form("%s: DpF/momErr DioWt"  ,Folder),1000, -50. , 50. ,Folder);
//   HBook1F(Hist->fDpFront0   ,"dp0f"     ,Form("%s: Track P0-P(front)" ,Folder), 200,  -5. ,  5. ,Folder);
//   HBook1F(Hist->fDpFront2   ,"dp2f"     ,Form("%s: Track P2-P(front)" ,Folder), 200,  -5. ,  5. ,Folder);
//   HBook1F(Hist->fPStOut     ,"pstout"   ,Form("%s: Track P(ST_Out)  " ,Folder), 400,  90. ,110. ,Folder);
//   HBook1F(Hist->fDpFSt      ,"dpfst"    ,Form("%s: Track Pf-Psto"     ,Folder), 200,  -5  ,  5. ,Folder);
//   HBook2F(Hist->fDpFVsZ1    ,"dpf_vs_z1",Form("%s: Track DPF Vs Z1"   ,Folder), 200, -2000.,0,200,-5.,5,Folder);

//   HBook1F(Hist->fPt         ,"pt"       ,Form("%s: Track Pt"          ,Folder), 600, 75,95,Folder);
//   HBook1F(Hist->fCosTh      ,"costh"    ,Form("%s: Track cos(theta)"  ,Folder), 100,-1,1,Folder);
//   HBook1F(Hist->fChi2       ,"chi2"     ,Form("%s: Track chi2 total"  ,Folder), 200, 0,200,Folder);
//   HBook1F(Hist->fChi2Dof    ,"chi2d"    ,Form("%s: track chi2/N(dof)" ,Folder), 500, 0, 10,Folder);

//   HBook1F(Hist->fNActive    ,"nactv"    ,Form("%s: N(active)"         ,Folder), 200, 0,200,Folder);
//   HBook1F(Hist->fNaFract    ,"nafr"     ,Form("%s: N(active fraction)",Folder), 110, 0.5,1.05,Folder);
//   HBook1F(Hist->fNWrong     ,"nwrng"    ,Form("%s: N(wrong drift sgn)",Folder), 100, 0,100,Folder);
//   HBook1F(Hist->fNDoublets  ,"nd"       ,Form("%s: N(doublets)"       ,Folder),  50, 0, 50,Folder);
//   HBook1F(Hist->fNadOverNd  ,"nad_nd"   ,Form("%s: Nad/N(doublets)"   ,Folder), 110, 0,1.1,Folder);
//   HBook1F(Hist->fNSSD       ,"nssd"     ,Form("%s: N(SS doublets)"    ,Folder),  50, 0, 50,Folder);
//   HBook1F(Hist->fNOSD       ,"nosd"     ,Form("%s: N(OS doublets)"    ,Folder),  50, 0, 50,Folder);
//   HBook1F(Hist->fNdOverNa   ,"nd_na"    ,Form("%s: NDoublets/Nactive" ,Folder), 100, 0,0.5,Folder);
//   HBook1F(Hist->fNssdOverNa ,"nssd_na"  ,Form("%s: NSSD/Nactive"      ,Folder), 100, 0,0.5,Folder);
//   HBook1F(Hist->fNosdOverNa ,"nosd_na"  ,Form("%s: NOSD/Nactive"      ,Folder), 100, 0,0.5,Folder);
//   HBook1F(Hist->fNZeroAmb   ,"nza"      ,Form("%s: N (Iamb = 0) hits" ,Folder), 100, 0,100,Folder);
//   HBook1F(Hist->fNzaOverNa  ,"nza_na"   ,Form("%s: NZeroAmb/Nactive"  ,Folder), 100, 0,  1,Folder);
//   HBook1F(Hist->fNMatActive ,"nma"      ,Form("%s: N (Mat Active"     ,Folder), 100, 0,100,Folder);
//   HBook1F(Hist->fNmaOverNa  ,"nma_na"   ,Form("%s: NMatActive/Nactive",Folder), 200, 0,   2,Folder);
//   HBook1F(Hist->fNBend      ,"nbend"    ,Form("%s: Nbend"             ,Folder), 100, 0,1000,Folder);

//   HBook1F(Hist->fT0         ,"t0"       ,Form("%s: track T0"          ,Folder), 200, 0,2000,Folder);
//   HBook1F(Hist->fT0Err      ,"t0err"    ,Form("%s: track T0Err"       ,Folder), 100, 0,  10,Folder);
//   HBook1F(Hist->fQ          ,"q"        ,Form("%s: track Q"           ,Folder),   4,-2,   2,Folder);
//   HBook1F(Hist->fFitCons[0] ,"fcon"     ,Form("%s: track fit cons [0]",Folder), 200, 0,   1,Folder);
//   HBook1F(Hist->fFitCons[1] ,"fcon1"    ,Form("%s: track fit cons [1]",Folder), 1000, 0,   0.1,Folder);
//   HBook1F(Hist->fD0         ,"d0"       ,Form("%s: track D0      "    ,Folder), 200,-200, 200,Folder);
//   HBook1F(Hist->fZ0         ,"z0"       ,Form("%s: track Z0      "    ,Folder), 200,-2000,2000,Folder);
//   HBook1F(Hist->fTanDip     ,"tdip"     ,Form("%s: track tan(dip)"    ,Folder), 400, -2.0 ,2.0,Folder);
//   HBook1F(Hist->fRMax       ,"rmax"     ,Form("%s: track R(max)  "    ,Folder), 200, 0., 1000,Folder);
//   HBook1F(Hist->fDtZ0       ,"dtz0"     ,Form("%s: DT(Z0), MC"        ,Folder), 400, -20.0 ,20.0,Folder);
//   HBook1F(Hist->fDtBack     ,"dtback"   ,Form("%s: DT(TT_back), MC"   ,Folder), 400, -20.0 ,20.0,Folder);

//   HBook1F(Hist->fResid      ,"resid"    ,Form("%s: hit residuals"     ,Folder), 500,-0.5 ,0.5,Folder);
//   HBook1F(Hist->fAlgMask    ,"alg"      ,Form("%s: algorithm mask"    ,Folder),  10,  0, 10,Folder);

//   HBook1F(Hist->fChi2Tcm    ,"chi2tcm"  ,Form("%s: chi2(t-c match)"   ,Folder), 250,  0  ,250 ,Folder);
//   HBook1F(Hist->fChi2XY     ,"chi2xy"   ,Form("%s: chi2(t-c match) XY",Folder), 300,-50  ,250 ,Folder);
//   HBook1F(Hist->fChi2T      ,"chi2t"    ,Form("%s: chi2(t-c match) T" ,Folder), 250,  0  ,250 ,Folder);

//   HBook1F(Hist->fDt         ,"dt"       ,Form("%s: track delta(T)"    ,Folder), 400,-20  ,20 ,Folder);
//   HBook1F(Hist->fDx         ,"dx"       ,Form("%s: track delta(X)"    ,Folder), 200,-500 ,500,Folder);
//   HBook1F(Hist->fDy         ,"dy"       ,Form("%s: track delta(Y)"    ,Folder), 200,-500 ,500,Folder);
//   HBook1F(Hist->fDz         ,"dz"       ,Form("%s: track delta(Z)"    ,Folder), 200,-250 ,250,Folder);
//   HBook1F(Hist->fDu         ,"du"       ,Form("%s: track-cluster DU)" ,Folder), 250,-250 ,250,Folder);
//   HBook1F(Hist->fDv         ,"dv"       ,Form("%s: track-cluster DV)" ,Folder), 200,-100 ,100,Folder);
//   HBook2F(Hist->fDvVsDu     ,"dv_vs_du" ,Form("%s: Track Dv Vs Du"    ,Folder), 100, -250,250,100,-100.,100,Folder);
//   HBook1F(Hist->fPath       ,"path"     ,Form("%s: track sdisk"       ,Folder),  50,   0 ,500,Folder);
//   HBook2F(Hist->fDuVsPath   ,"du_vs_path",Form("%s: Track Du Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
//   HBook2F(Hist->fDvVsPath   ,"dv_vs_path",Form("%s: T-C  Dv Vs Path"  ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
//   HBook2F(Hist->fDtVsPath   ,"dt_vs_path",Form("%s: T-C DT Vs Path"   ,Folder),  50,   0 ,500,100,  -5.,  5.,Folder);
//   HBook2F(Hist->fDuVsTDip   ,"du_vs_tdip",Form("%s: Track Du Vs TDip" ,Folder), 100, 0.5 ,1.5,200,-200.,200.,Folder);
//   HBook2F(Hist->fDvVsTDip   ,"dv_vs_tdip",Form("%s: Track Dv Vs TDip" ,Folder), 100, 0.5 ,1.5,200,-200.,200.,Folder);

//   HBook1F(Hist->fZ1         ,"z1"       ,Form("%s: track Z1      "    ,Folder), 200,-2000,2000,Folder);
//   HBook1F(Hist->fNClusters  ,"ncl"      ,Form("%s: track N(clusters)" ,Folder),  10, 0   , 10,Folder);
//   HBook1F(Hist->fDiskID     ,"disk_id"  ,Form("%s: track disk ID"     ,Folder),  10,-5   ,  5,Folder);
//   HBook1F(Hist->fXCal       ,"xcal"     ,Form("%s: track XCal"        ,Folder), 200,-1000,1000,Folder);
//   HBook1F(Hist->fYCal       ,"ycal"     ,Form("%s: track YCal"        ,Folder), 200,-1000,1000,Folder);
//   HBook1F(Hist->fZCal       ,"zcal"     ,Form("%s: track ZCal"        ,Folder), 200, 1500,3500,Folder);
//   HBook1F(Hist->fXTrk       ,"xtrk"     ,Form("%s: track XTrk"        ,Folder), 200,-1000,1000,Folder);
//   HBook1F(Hist->fYTrk       ,"ytrk"     ,Form("%s: track YTrk"        ,Folder), 200,-1000,1000,Folder);
//   HBook1F(Hist->fRTrk       ,"rtrk"     ,Form("%s: track RTrk"        ,Folder), 200,-1000,1000,Folder);
//   HBook1F(Hist->fZTrk       ,"ztrk"     ,Form("%s: track ZTrk"        ,Folder), 200,-1000,1000,Folder);
//   HBook1F(Hist->fECl        ,"ecl"      ,Form("%s: cluster E"         ,Folder), 300, 0   ,150,Folder);
//   HBook1F(Hist->fEClEKin    ,"ecl_ekin" ,Form("%s: cluster E/Ekin(mu)",Folder), 500, 0   ,5,Folder);
//   HBook1F(Hist->fECalP      ,"ecalp"    ,Form("%s: E(caltot)/trackP"  ,Folder), 300, 0   ,1.5,Folder);
//   HBook1F(Hist->fEDiskP     ,"ediskp"   ,Form("%s: E(disk)/trackP"    ,Folder), 300, 0   ,1.5,Folder);
//   HBook1F(Hist->fEp         ,"ep"       ,Form("%s: track E/P"         ,Folder), 300, 0   ,1.5,Folder);
//   HBook2F(Hist->fEpVsPath   ,"ep_vs_path",Form("%s: E/P Vs Path"      ,Folder),  50,   0 ,500,150,  0.,  1.5,Folder);
//   HBook2F(Hist->fNHVsStation,"nh_vs_st" ,Form("%s: N(hits) Vs Station",Folder),  40, 0,40,10,-0.5,9.5,Folder);
//   HBook2F(Hist->fNHVsNSt    ,"nh_vs_nst",Form("%s: N(hits) Vs NSt"    ,Folder),  10,-0.5,9.5,40,-0.5,39.5,Folder);

//   HBook1F(Hist->fRSlope     ,"rslope"   ,Form("%s: Res Slope"         ,Folder), 200,-20 , 20,Folder);
//   HBook1F(Hist->fXSlope     ,"xslope"   ,Form("%s: Res/Sig Slope"     ,Folder), 200,-20 , 20,Folder);

//   HBook2F(Hist->fEpVsDt     ,"ep_vs_dt" ,Form("%s: E/P vs Dt"         ,Folder), 200, -10, 10,150,0.,1.5,Folder);
//   HBook1F(Hist->fEleLogLHCal,"ele_llh_c",Form("%s: ELE Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
//   HBook1F(Hist->fMuoLogLHCal,"muo_llh_c",Form("%s: MUO Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
//   HBook1F(Hist->fLogLHRCal  ,"llhr_cal" ,Form("%s: LogLH(e/m) Cal"    ,Folder), 200,-100,100,Folder);
//   HBook1F(Hist->fLogLHRDeDx ,"llhr_dedx",Form("%s: LogLH(e/m) De/Dx"  ,Folder), 200,-20 , 20,Folder);
//   HBook1F(Hist->fLogLHRXs   ,"llhr_xs"  ,Form("%s: LogLH(e/m) XSlope" ,Folder), 200,-20 , 20,Folder);
//   HBook1F(Hist->fLogLHRTrk  ,"llhr_trk" ,Form("%s: LogLH(e/m) Trk"    ,Folder), 200,-20 , 20,Folder);
//   HBook1F(Hist->fLogLHR     ,"llhr"     ,Form("%s: LogLH(e/m)"        ,Folder), 200,-100 ,100,Folder);

//   HBook1F(Hist->fPdgCode[0] ,"pdg_0"    ,Form("%s: track PDG code[0]" ,Folder), 100,-50,50,Folder);
//   HBook1F(Hist->fPdgCode[1] ,"pdg_1"    ,Form("%s: track PDG code[1]" ,Folder), 500,-2500,2500,Folder);
//   HBook1F(Hist->fFrGH       ,"fgh"      ,Form("%s: Fraction Goog Hits",Folder), 100, 0,1,Folder);

//   HBook2F(Hist->fNEPlVsNHPl ,"nep_vs_nhp",Form("%s: Track NEXP vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
//   HBook2F(Hist->fNDPlVsNHPl ,"ndp_vs_nhp",Form("%s: Track NDIF vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
//   HBook2F(Hist->fChi2dVsNDPl,"chi2d_vs_ndp",Form("%s: Track Chi2/Dof vs NDP",Folder), 30, 0,30,100,0.,10,Folder);
//   HBook2F(Hist->fDpFVsNDPl  ,"dpf_vs_ndp"  ,Form("%s: Track DpF vs NDP",Folder)     , 30, 0,30,100,-5,5,Folder);

//   HBook1F(Hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
//   HBook1F(Hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);

//   HBook1F(Hist->fSinTC  ,"sintc"  ,Form("%s: sin(T^C)"      ,Folder),200, -1, 1,Folder);
//   HBook1F(Hist->fDrTC   ,"drtc"   ,Form("%s: DeltaR(T^C)"   ,Folder),200, -200, 200,Folder);
//   HBook1F(Hist->fSInt   ,"sint"   ,Form("%s: SInt"          ,Folder),200, -500, 500,Folder);
//   HBook1F(Hist->fDaveTrkQual,"dtqual",Form("%s:DaveTrkQual" ,Folder),500, -2.5, 2.5,Folder);
//   HBook1F(Hist->fNMcStrawHits,"nmc_hits",Form("%s:N(MC Straw Hits by parent",Folder),150, 0, 150,Folder);

//   HBook1F(Hist->fHitEnergy  ,"hit_e"   ,Form("%s: hit energy"   ,Folder),100, 0, 0.05,Folder);
//   HBook1F(Hist->fHitDt      ,"hit_dt"  ,Form("%s: hit Dt=T1-T2" ,Folder),100, -10, 10,Folder);
//   HBook1F(Hist->fHitTRel    ,"hit_trel",Form("%s: hit T-<T>   " ,Folder),100, -50, 50,Folder);

//   HBook1F(Hist->fT0MinusTM  ,"t0_minus_tm",Form("%s: T0-<T>"    ,Folder),100, -50, 50,Folder);

// }

// //-----------------------------------------------------------------------------
// void TTrackAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
//   char name [200];
//   //  char title[200];

//   HBook1D(Hist->fLumWt     ,"lumwt"    ,Form("%s: Luminosity Weight"               ,Folder),200, 0,10,Folder);
//   Hist->fLumWt->Sumw2(kTRUE);

//   HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1, 1,Folder);
//   HBook1F(Hist->fEleMom    ,"ce_mom"   ,Form("%s: Conversion Electron Momentum"    ,Folder),2500,  0,500,Folder);
//   HBook1D(Hist->fDioMom    ,"dio_mom"  ,Form("%s: DIO momentum"                    ,Folder),1000, 50,150,Folder);
//   HBook1F(Hist->fRv         ,"rv"      ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 1000,Folder);
//   HBook1F(Hist->fZv         ,"zv"      ,Form("%s: Z(Vertex)"                       ,Folder), 300, 0,15000,Folder);
//   HBook1F(Hist->fNClusters ,"ncl"      ,Form("%s: Number of Reconstructed Clusters",Folder),200,0,200,Folder);
//   HBook1F(Hist->fNTracks   ,"ntrk"     ,Form("%s: Number of Reconstructed Tracks"  ,Folder),100,0,100,Folder);
//   HBook1F(Hist->fNStrawHits[0],"nsh_0" ,Form("%s: Number of Straw Hits [0]"        ,Folder),250,0,250,Folder);
//   HBook1F(Hist->fNStrawHits[1],"nsh_1" ,Form("%s: Number of Straw Hits [1]"        ,Folder),400,0,8000,Folder);
//   HBook1F(Hist->fNGoodSH   ,"nsh50"    ,Form("%s: N(SH) +/-50"                     ,Folder),300,0,1500,Folder);
//   HBook1F(Hist->fAlgMask   ,"alg_mask" ,Form("%s: algorithm mask"                  ,Folder), 10,0,  10,Folder);
//   HBook1F(Hist->fDtClT     ,"dt_clt"   ,Form("%s: DT(cluster-track)"               ,Folder),100,-100,100,Folder);
//   HBook1F(Hist->fDtClS     ,"dt_cls"   ,Form("%s: DT(cluster-straw hit)"           ,Folder),200,-200,200,Folder);
//   HBook1F(Hist->fSHTime    ,"shtime"   ,Form("%s: Straw Hit Time"                  ,Folder),400,0,2000,Folder);
//   HBook1F(Hist->fEMax      ,"emax"     ,Form("%s: Max cluster energy"              ,Folder),150,0,150,Folder);
//   HBook1F(Hist->fNHyp      ,"nhyp"     ,Form("%s: N(fit hypotheses)"               ,Folder),5,0,5,Folder);
//   HBook1F(Hist->fBestHyp[0],"bfh0"     ,Form("%s: Best Fit Hyp[0](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
//   HBook1F(Hist->fBestHyp[1],"bfh1"     ,Form("%s: Best Fit Hyp[1](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
//   HBook1F(Hist->fNGenp     ,"ngenp"    ,Form("%s: N(Gen Particles)"                ,Folder),500,0,500,Folder);

//   for (int i=0; i<2; i++) {
//     sprintf(name,"ncch_%i",i);
//     HBook1F(Hist->fNCaloCrystalHits[i],name,Form("%s: N(calo crystal hits) [%i]",Folder,i),500,0,1000,Folder);
//     sprintf(name,"ncch_vs_disk_%i",i);
//     HBook2F(Hist->fNCaloHitsVsDisk[i],name,Form("%s: N(calo crystal hits) vs disk[%i]",Folder,i),4,0,4,200,0,200,Folder);
//     sprintf(name,"ncch_vs_row_%i",i);
//     HBook2F(Hist->fNCaloHitsVsRow[i],name,Form("%s: N(calo crystal hits) vs row [%i]",Folder,i),20,0,20,200,0,200,Folder);
//     sprintf(name,"ncch_vs_col_%i",i);
//     HBook2F(Hist->fNCaloHitsVsCol[i],name,Form("%s: N(calo crystal hits) vs col [%i]",Folder,i),50,0,50,200,0,200,Folder);
//   }

//   for (int i=0; i<kNDisks; i++) {
//     HBook1F(Hist->fETot        [i],Form("etot_%i"    ,i),Form("%s: Etot[%i]",Folder,i), 300, 0,150,Folder);
//     HBook2F(Hist->fECrVsR      [i],Form("ecr_vs_r_%i",i),Form("%s: E Cr Vs R [%i]"    ,Folder,i), 100, 0,1000,500,0,100,Folder);
//     HBook2F(Hist->fNCrVsR      [i],Form("ncr_vs_r_%i",i),Form("%s: N Cr Vs R [%i]"    ,Folder,i), 100, 0,1000,100,0,100,Folder);

//     HBook2F(Hist->fNCrystalHitsVsR[i],Form("ncrh_vs_r_%i",i),Form("%s: N Crystal Hits[%i] vs R",Folder,i), 100, 0, 1000,100,0,100,Folder);
//     HBook2F(Hist->fNHitCrystalsVsR[i],Form("nhcr_vs_r_%i",i),Form("%s: N Hit Crystals[%i] vs R",Folder,i), 100, 0, 1000,100,0,100,Folder);
//   }

//   HBook1F(Hist->fNHitCrystalsTot,"nhcr_tot",Form("%s: NHit Crystals Tot",Folder), 100, 0,100,Folder);
//   HBook1F(Hist->fECal,"ecal",Form("%s: E(cal), sum over both disks",Folder), 500, 0,250,Folder);
//   HBook1F(Hist->fECalOverEKin,"ec_over_ek",Form("%s: E(cal)/E(kin)",Folder), 200, 0,2,Folder);
//   HBook1F(Hist->fInstLumi    ,"lumi"      ,Form("%s: N(protons/mubunch)",Folder), 2000, 0,1.e8,Folder);
// }

// //-----------------------------------------------------------------------------
// void TTrackAnaModule::BookSimpHistograms(SimpHist_t* Hist, const char* Folder) {
//   //  char name [200];
//   //  char title[200];

//   HBook1F(Hist->fPdgCode        ,"pdg"   ,Form("%s: PDG code"                     ,Folder),200,-100,100,Folder);
//   HBook1F(Hist->fNStrawHits     ,"nsth"  ,Form("%s: n straw hits"                 ,Folder),200,   0,200,Folder);
//   HBook1F(Hist->fMomTargetEnd   ,"ptarg" ,Form("%s: CE mom after Stopping Target" ,Folder),400,  90,110,Folder);
//   HBook1F(Hist->fMomTrackerFront,"pfront",Form("%s: CE mom at the Tracker Front"  ,Folder),400,  90,110,Folder);
// }

