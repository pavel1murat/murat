//////////////////////////////////////////////////////////////////////////////
// use of TStnTrack::fTmp:
//
// Tmp(0) : corrected momentum at the tracker front - not yet
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
// 34  : EVT_7: events with E_CL > 60 and no CalPatRec tracks 
// 35  : TRK_1: events with P > 106 MeV/c - misreconstruction
// 36  : TRK_23 events with P < 80: odd misidentified muons - turned out to be DIO electrons
// 37  : TRK_26 LLHR_CAL > 5
// 38  : EVT_7: events with E_CL > 60 and no CalPatRec tracks and TrkPatRec
// 39  : trk_1: events with |SIN_TC| > 0.6
// 40  : EVT_7: events with E_CL > 60 and no tracks at all
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
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------

#include "ana/TTrackAnaModuleA.hh"

ClassImp(TTrackAnaModuleA)
//-----------------------------------------------------------------------------
TTrackAnaModuleA::TTrackAnaModuleA(const char* name, const char* title):
  TTrackAnaModuleBase(name,title)
{
  fTrackBlockName = "TrackBlock";
  fTrackNumber.Set(100);
					// track-cluster matching timing cut
  fMinDtTcm        = -5.;
  fMaxDtTcm        =  8.;
//-----------------------------------------------------------------------------
// initialize Track ID
// 0: SetC  1: TrkQual>0.1 2:TrkQual>0.4
//-----------------------------------------------------------------------------
  fNTrkID          = 4;
  for (int i=0; i<fNTrkID; i++) {
    fTrackID[i]    = new TStnTrackID();
  }

  fTrackID[1]->SetMaxMomErr (100);
  fTrackID[1]->SetMaxT0Err  (100);
  fTrackID[1]->SetMinFitCons(-1.);
  fTrackID[1]->SetMinTrkQual(0.1);

  fTrackID[2]->SetMaxMomErr (100);
  fTrackID[2]->SetMaxT0Err  (100);
  fTrackID[2]->SetMinFitCons(-1.);
  fTrackID[2]->SetMinTrkQual(0.4);

  fTrackID[3]->SetMaxMomErr (100);
  fTrackID[3]->SetMaxT0Err  (100);
  fTrackID[3]->SetMinFitCons(-1.);
  fTrackID[3]->SetMinTrkQual(0.4);
  fTrackID[3]->SetMinNActive(-1 );

  fBestID     = 2;			// best: DaveTrkQual > 0.4
}

//-----------------------------------------------------------------------------
TTrackAnaModuleA::~TTrackAnaModuleA() {
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TTrackAnaModuleA::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockName.Data(),"TStnTrackBlock"   ,&fTrackBlock  );

  RegisterDataBlock("ClusterBlock"        ,"TStnClusterBlock" ,&fClusterBlock);
  RegisterDataBlock("CalDataBlock"        ,"TCalDataBlock"    ,&fCalDataBlock);
  RegisterDataBlock("StrawDataBlock"      ,"TStrawDataBlock"  ,&fStrawDataBlock);
  RegisterDataBlock("GenpBlock"           ,"TGenpBlock"       ,&fGenpBlock);
  RegisterDataBlock("SimpBlock"           ,"TSimpBlock"       ,&fSimpBlock);
  RegisterDataBlock("VDetBlock"           ,"TVDetDataBlock"   ,&fVDetBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//-----------------------------------------------------------------------------
// need MC truth branch
//-----------------------------------------------------------------------------
void TTrackAnaModuleA::FillEventHistograms(HistBase_t* Hist) {
  double            cos_th, dio_wt, xv(-1.e6), yv(-1.e6), rv(-1.e6), zv(-1.e6), p;
  double            e, m, r;
  TLorentzVector    mom(1.,0.,0.,0);

  EventHist_t* hist = (EventHist_t*) Hist;

  if (fParticle) { 
    fParticle->Momentum(mom);
    xv = fParticle->Vx()+3904.;
    yv = fParticle->Vy();
    rv = sqrt(xv*xv+yv*yv);
    zv = fParticle->Vz();
  }

  p      = mom.P();
  cos_th = mom.Pz()/p;
  dio_wt = TStntuple::DioWeightAl(p);

  hist->fLumWt->Fill(fLumWt);
  hist->fRv->Fill(rv);
  hist->fZv->Fill(zv);
  hist->fEleMom->Fill(p);
  hist->fDioMom->Fill(p,dio_wt);
  hist->fEleCosTh->Fill(cos_th);

  hist->fNClusters->Fill(fNClusters);
  hist->fNTracks->Fill  (fNTracks[0]);

  int nsh = GetHeaderBlock()->fNStrawHits;
  hist->fNStrawHits[0]->Fill(nsh);
  hist->fNStrawHits[1]->Fill(nsh);

  double emax   = -1;
  double t0_cls = -1;
  double dt     = 9999.;

  TStnCluster* cluster(0);
  if (fNClusters > 0) cluster = fClusterBlock->Cluster(0);

  TStnTrack* track(0);
  if (fNTracks > 0) track = fTrackBlock->Track(0);

  if (cluster) {
    emax   = cluster->Energy();
    t0_cls = cluster->Time();
  }

  double t0_trk = -1;
  if (track) {
    t0_trk = track->fT0;
  }

  if (track && cluster) {
    dt = t0_cls-t0_trk;
  }

  hist->fDtClT->Fill(dt);
  hist->fEMax->Fill(emax);

  TStrawHitData*  sh;
  int n_good_hits = 0;
  int nstraw_hits = fStrawDataBlock->NHits();
  for (int i=0; i<nstraw_hits; i++ ) {
    sh = fStrawDataBlock->Hit(i);
    dt = t0_cls-sh->Time() + 15;
    hist->fDtClS->Fill(dt);
    hist->fSHTime->Fill(sh->Time());

    if (fabs(dt+15.)< 50) n_good_hits += 1;
  }

  hist->fNGoodSH->Fill(n_good_hits);

  hist->fNHyp->Fill(fNHyp);
  hist->fBestHyp[0]->Fill(fBestHyp[0]);
  hist->fBestHyp[1]->Fill(fBestHyp[1]);
  hist->fNGenp->Fill(fNGenp);
//-----------------------------------------------------------------------------
// crystals - count crystals with E > 1MeV
//-----------------------------------------------------------------------------
  TCalHitData* cch;

  int n_cch_1mev = 0;

  if (fCalorimeterType == 1) {
//-----------------------------------------------------------------------------
// vane calorimeter
//-----------------------------------------------------------------------------
    int  nhits_vane[2][kNDisks], nhits_row [2][20], nhits_col[2][50];
    int  crystal_id, vane_id, local_id, vane_row, vane_col;

    for (int i=0; i<kNDisks; i++) {
      nhits_vane[0][i] = 0;
      nhits_vane[1][i] = 0;
    }
      
    for (int i=0; i<20; i++) {
      nhits_row[0][i] = 0;
      nhits_row[1][i] = 0;
    }

    for (int i=0; i<50; i++) {
      nhits_col[0][i] = 0;
      nhits_col[1][i] = 0;
    }
      
    for (int ic=0; ic<fNCalHits; ic++) {
      cch        = fCalDataBlock->CalHitData(ic);
      crystal_id = cch->ID();

      if (cch->Energy() > 1.) {
	n_cch_1mev += 1;
      }
      // for each crystal determine its row and column
      // the following is for vanes
      vane_id  = crystal_id/484.;
      local_id = crystal_id-vane_id*484;
      vane_row = local_id/44;
      vane_col = local_id-vane_row*44;
      
      nhits_vane[0][vane_id ] += 1;
      nhits_row [0][vane_row] += 1;
      nhits_col [0][vane_col] += 1;
      
      if (cch->Energy() > 1.) {
	nhits_row [1][vane_row] += 1;
	nhits_col [1][vane_col] += 1;
	nhits_vane[1][vane_id ] += 1;
      }
    }

    hist->fNCaloCrystalHits[0]->Fill(fNCalHits);
    hist->fNCaloCrystalHits[1]->Fill(n_cch_1mev);

    for (int iv=0; iv<4; iv++) {
      hist->fNCaloHitsVsDisk[0]->Fill(iv,nhits_vane[0][iv]);
      hist->fNCaloHitsVsDisk[1]->Fill(iv,nhits_vane[1][iv]);
    }

    for (int ir=0; ir<20; ir++) {
      hist->fNCaloHitsVsRow[0]->Fill(ir,nhits_row[0][ir]);
      hist->fNCaloHitsVsRow[1]->Fill(ir,nhits_row[1][ir]);
    }

    for (int ic=0; ic<50; ic++) {
      hist->fNCaloHitsVsCol[0]->Fill(ic,nhits_col[0][ic]);
      hist->fNCaloHitsVsCol[1]->Fill(ic,nhits_col[1][ic]);
    }
  }
  else if (fCalorimeterType == 2) {
//-----------------------------------------------------------------------------
// disk calorimeter
//-----------------------------------------------------------------------------
    int      ndisks, n_hit_crystals[4], n_hit_crystals_tot;
    double   etot[4];

    TCalHitData* hit;

    ndisks = fDiskCalorimeter->NDisks();

    int   bin, hit_id, idisk, nhits;
    int   nhits_r[kNDisks][100], n_hit_crystals_r[kNDisks][100];

    for (int id=0; id<kNDisks; id++) {
      n_hit_crystals[id] = 0;
      etot[id]           = 0;

      for (int ib=0; ib<100; ib++) {
	nhits_r         [id][ib] = 0;
	n_hit_crystals_r[id][ib] = 0;
      }
    }

    nhits = fCalDataBlock->NHits();

    for (int i=0; i< nhits; i++) {
      hit    = fCalDataBlock->CalHitData(i);

      hit_id = hit->ID();
      idisk  = fDiskCalorimeter->DiskNumber(hit_id);
      r      = fDiskCalorimeter->CrystalRadius(hit_id);
      e      = hit->Energy(); 

      etot          [idisk] += e;
      n_hit_crystals[idisk] += 1;

      hist->fECrVsR[idisk]->Fill(r,e);
      hist->fNCrVsR[idisk]->Fill(r,1);

      bin  = (int) (r/10.);

      nhits_r         [idisk][bin] += 1;
//-----------------------------------------------------------------------------
// this is not correct, one needs to check whether this crystal has been hit,
// for the moment, to get going, ignore that
//-----------------------------------------------------------------------------
      n_hit_crystals_r[idisk][bin] += 1;
    }

    n_hit_crystals_tot = 0;

    double ecal = 0;
    for (int id=0; id<ndisks; id++) {
      n_hit_crystals_tot += n_hit_crystals[id];
      ecal += etot[id];
//-----------------------------------------------------------------------------
// fill 'per-disk' histograms
//-----------------------------------------------------------------------------
      hist->fETot[id]->Fill(etot[id]);

//-----------------------------------------------------------------------------
// 100 is fixed by the number of bins in the radial distributions
//-----------------------------------------------------------------------------
      for (int ib=0; ib<100; ib++) {
	r = (ib+0.5)*10.;
	hist->fNCrystalHitsVsR[id]->Fill(r,nhits_r         [id][ib]);
	hist->fNHitCrystalsVsR[id]->Fill(r,n_hit_crystals_r[id][ib]);
      }
    }

    hist->fNHitCrystalsTot->Fill(n_hit_crystals_tot);
    hist->fECal->Fill(ecal);

    double ekin(-1.);
//-----------------------------------------------------------------------------
// there is an inconsistency in the SIMP block filling - in Mu2e offline 
// the particle momentumis is kept in MeV/c, while the PDG mass  -in GeV/c^2..
// thus the energy is screwed up... kludge around
// assign muon mass
//-----------------------------------------------------------------------------
    if (fSimp) {
      p    = fSimp->fStartMom.P();
      m    = 105.658; // in MeV
      ekin = sqrt(p*p+m*m)-m;
    }
    hist->fECalOverEKin->Fill(ecal/ekin);

  }

  hist->fInstLumi->Fill(GetHeaderBlock()->InstLum());
}

//-----------------------------------------------------------------------------
void TTrackAnaModuleA::FillCaloHistograms(HistBase_t* HistBase, TStnCrystal* Cr) {

  int                    nhits;
  float                  t, e, r, e700, n700;
  TCalHitData*           hit;

  CaloHist_t* Hist = (CaloHist_t*) HistBase;

					// determine crystal coordinates
  TDisk* disk = Cr->Disk();
  int   idisk = disk->SectionID();
					// time needs to be defiend
  //    t  = Cr->Time();
  e     = Cr->Energy();
  r     = Cr->Radius();
  nhits = Cr->NHits();

  Hist->fDiskID->Fill(idisk);

  Hist->fEnergy  [idisk]->Fill(e);
  Hist->fNHits   [idisk]->Fill(nhits);
  Hist->fRadius  [idisk]->Fill(r);
  Hist->fRadiusWE[idisk]->Fill(r,e);
    
  e700 = 0;
  n700 = 0;
  for (int i=0; i<nhits; i++) {
    hit  = Cr->CalHitData(i);
    t   = hit->Time();
    Hist->fTime[idisk]->Fill(t);
    if (t > 700.) {
      n700 += 1;
      e700 += hit->Energy();
      Hist->fT700[idisk]->Fill(t);
    }
  }

  Hist->fE700   [idisk]->Fill(e700);
  Hist->fN700   [idisk]->Fill(n700);

  if (n700 > 0) {
    Hist->fR700  [idisk]->Fill(r);
    Hist->fRWE700[idisk]->Fill(r,e700);
  }
}


//-----------------------------------------------------------------------------
void TTrackAnaModuleA::FillClusterHistograms(HistBase_t* Hist, TStnCluster* Cluster) {
  int   row, col;
  float  x, y, z, r;

  ClusterHist_t* hist = (ClusterHist_t*) Hist;

  row = Cluster->Ix1();
  col = Cluster->Ix2();

  x   = Cluster->fX+3904.;
  y   = Cluster->fY;
  z   = Cluster->fZ;
  r   = sqrt(x*x+y*y);

  if ((row < 0) || (row > 9999)) row = -9999;
  if ((col < 0) || (col > 9999)) col = -9999;

  hist->fDiskID->Fill(Cluster->DiskID());
  hist->fEnergy->Fill(Cluster->Energy());
  hist->fT0->Fill(Cluster->Time());
  hist->fRow->Fill(row);
  hist->fCol->Fill(col);
  hist->fX->Fill(x);
  hist->fY->Fill(y);
  hist->fZ->Fill(z);
  hist->fR->Fill(r);

  hist->fYMean->Fill(Cluster->fYMean);
  hist->fZMean->Fill(Cluster->fZMean);
  hist->fSigY->Fill(Cluster->fSigY);
  hist->fSigZ->Fill(Cluster->fSigZ);
  hist->fSigR->Fill(Cluster->fSigR);
  hist->fNCr0->Fill(Cluster->fNCrystals);
  hist->fNCr1->Fill(Cluster->fNCr1);
  hist->fFrE1->Fill(Cluster->fFrE1);
  hist->fFrE2->Fill(Cluster->fFrE2);
  hist->fSigE1->Fill(Cluster->fSigE1);
  hist->fSigE2->Fill(Cluster->fSigE2);
}

//-----------------------------------------------------------------------------
void TTrackAnaModuleA::FillGenpHistograms(HistBase_t* Hist, TGenParticle* Genp) {
  int    gen_id;
  float  p, cos_th, z0, t0, r0, x0, y0;

  TLorentzVector mom, v;

  GenpHist_t* hist = (GenpHist_t*) Hist;

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

  hist->fPdgCode[0]->Fill(Genp->GetPdgCode());
  hist->fPdgCode[1]->Fill(Genp->GetPdgCode());
  hist->fGenID->Fill(gen_id);
  hist->fZ0->Fill(z0);
  hist->fT0->Fill(t0);
  hist->fR0->Fill(r0);
  hist->fP->Fill(p);
  hist->fCosTh->Fill(cos_th);
}

//-----------------------------------------------------------------------------
void TTrackAnaModuleA::FillSimpHistograms(HistBase_t* Hist, TSimParticle* Simp) {

  SimpHist_t* hist = (SimpHist_t*) Hist;
  
  hist->fPdgCode->Fill(Simp->fPdgCode);
  hist->fMomTargetEnd->Fill(Simp->fMomTargetEnd);
  hist->fMomTrackerFront->Fill(Simp->fMomTrackerFront);
  hist->fNStrawHits->Fill(Simp->fNStrawHits);
}

//-----------------------------------------------------------------------------
// for DIO : ultimately, one would need to renormalize the distribution
//-----------------------------------------------------------------------------
void TTrackAnaModuleA::FillTrackHistograms(HistBase_t* Hist, TStnTrack* Track, TrackPar_t* TrackPar) {

  TLorentzVector  mom;
  double          r;
  //  int             itrk;
  TrackPar_t*     tp;

  TrackHist_t* hist = (TrackHist_t*) Hist;

  tp   = TrackPar;

  hist->fP[0]->Fill (tp->fP);		// corrected momentum in the first point
  hist->fP[1]->Fill (tp->fP);
  hist->fP[2]->Fill (tp->fP);

  hist->fP0->  Fill (Track->fP0);	// momentum zt Z0
  hist->fP2->  Fill (Track->fP2);	// momentum at z2 = -1540, the same as in the 1st point

  hist->fPDio->Fill  (tp->fP, tp->fDioWt);
  hist->fPlw->Fill   (tp->fP, tp->fLumWt);
  hist->fPDiolw->Fill(tp->fP, tp->fTotWt);


  hist->fFitMomErr->Fill(Track->fFitMomErr);

  hist->fPt    ->Fill(Track->fPt    );
  hist->fPFront->Fill(Track->fPFront);
  hist->fPStOut->Fill(Track->fPStOut);
					// dp: Tracker-only resolution

  hist->fDpFront ->Fill(tp->fDpF);
  hist->fXDpF    ->Fill(tp->fDpF/Track->fFitMomErr);
  hist->fDpFDio  ->Fill(tp->fDpF,tp->fDioWt);
  hist->fXDpFDio ->Fill(tp->fDpF/Track->fFitMomErr,tp->fDioWt);
  hist->fDpFront0->Fill(tp->fDp0);
  hist->fDpFront2->Fill(tp->fDp2);
  hist->fDpFSt   ->Fill(tp->fDpFSt);
  hist->fDpFVsZ1 ->Fill(Track->fZ1,tp->fDpF);

  hist->fCosTh->Fill(Track->Momentum()->CosTheta());
  hist->fChi2->Fill (Track->fChi2);

  float na = Track->NActive();

  hist->fChi2Dof->Fill(Track->fChi2/(na-5.));
  hist->fNActive->Fill(na);
  hist->fNaFract->Fill(na/(Track->NHits()+0.));
  hist->fNWrong->Fill(Track->NWrong());

  float nd  = Track->NDoublets();
  float nad = Track->NDoubletsAct();
  hist->fNDoublets->Fill(nd);
  hist->fNadOverNd->Fill(nad/nd);
  hist->fNOSD->Fill(Track->NOSDoublets());
  hist->fNSSD->Fill(Track->NSSDoublets());
  hist->fNdOverNa->Fill(nd/na);
  hist->fNosdOverNa->Fill(Track->NOSDoublets()/na);
  hist->fNssdOverNa->Fill(Track->NSSDoublets()/na);
  hist->fNZeroAmb->Fill(Track->NHitsAmbZero());
  hist->fNzaOverNa->Fill(Track->NHitsAmbZero()/na);

  int nma = Track->NMatActive();

  hist->fNMatActive->Fill(nma);
  hist->fNmaOverNa->Fill(nma/na);

  hist->fNBend->Fill(Track->NBend());

  hist->fT0->Fill(Track->fT0);
  hist->fT0Err->Fill(Track->fT0Err);
  //  printf("TTrackAnaModuleA::FillTrackHistograms: track charge is not defined yet\n");
  hist->fQ->Fill(-1);
  hist->fFitCons[0]->Fill(Track->fFitCons);
  hist->fFitCons[1]->Fill(Track->fFitCons);

  hist->fD0->Fill(Track->fD0);
  hist->fZ0->Fill(Track->fZ0);
  hist->fTanDip->Fill(Track->fTanDip);
  hist->fDtZ0->Fill(tp->fDtZ0);
  hist->fRMax->Fill(Track->RMax());
  hist->fAlgMask->Fill(Track->AlgMask());
//-----------------------------------------------------------------------------
// track-cluster matching part: 
// - for residuals, determine intersection with the most energetic cluster
// - for track -only parameters use intersection with lowest trajectory length
//-----------------------------------------------------------------------------
  TStnTrack::InterData_t*    vt = Track->fVMinS;  // track-cluster residuals

  if (vt) {
    hist->fDiskID->Fill(vt->fID  );
    hist->fXTrk->Fill  (vt->fXTrk);
    hist->fYTrk->Fill  (vt->fYTrk);

    r = sqrt(vt->fXTrk*vt->fXTrk+vt->fYTrk*vt->fYTrk);
    hist->fRTrk->Fill  (r);

    hist->fZTrk->Fill  (vt->fZTrk);
  }
  else {
//-----------------------------------------------------------------------------
// fill histograms with numbers easy to recognize as dummy
//-----------------------------------------------------------------------------
    hist->fDiskID->Fill(-1.);
    hist->fXTrk->Fill  (999.);
    hist->fYTrk->Fill  (999.);
    hist->fRTrk->Fill  (999.);
    hist->fZTrk->Fill  (-1. );
  }
//-----------------------------------------------------------------------------
// there is an inconsistency in the SIMP block filling - in Mu2e offline 
// the particle momentumis is kept in MeV/c, while the PDG mass  -in GeV/c^2..
// thus the energy is screwed up... kludge around ... assign muon mass
//-----------------------------------------------------------------------------
  double    ekin(-1.);
  if (fSimp) {
    double p, m;
    p    = tp->fP;
    m    = 105.658; // in MeV
    ekin = sqrt(p*p+m*m)-m;
  }

  hist->fECl->Fill(tp->fEcl);
  hist->fEClEKin->Fill(tp->fEcl/ekin);
  hist->fEp->Fill(tp->fEp);
  hist->fEpVsPath->Fill(tp->fPath,tp->fEp);

  hist->fDx->Fill(tp->fDx);
  hist->fDy->Fill(tp->fDy);
  hist->fDz->Fill(tp->fDz);

  hist->fDt->Fill(tp->fDt);
  hist->fChi2Tcm->Fill(tp->fChi2Tcm);
  hist->fChi2XY->Fill(tp->fChi2XY);
  hist->fChi2T->Fill (tp->fChi2T);

  hist->fDu->Fill    (tp->fDu);
  hist->fDv->Fill    (tp->fDv);
  hist->fDvVsDu->Fill(tp->fDu,tp->fDv);

  hist->fPath->Fill(tp->fPath);
  hist->fDuVsPath->Fill(tp->fPath,tp->fDu);

  hist->fDvVsPath->Fill(tp->fPath,tp->fDv);
  hist->fDtVsPath->Fill(tp->fPath,tp->fDt);

  hist->fDuVsTDip->Fill(Track->fTanDip,tp->fDu);
  hist->fDvVsTDip->Fill(Track->fTanDip,tp->fDv);

  hist->fZ1->Fill(Track->fZ1);

  int ncl = Track->NClusters();
  hist->fNClusters->Fill(ncl);

  hist->fRSlope->Fill(Track->RSlope());
  hist->fXSlope->Fill(Track->XSlope());

  double llhr_dedx, llhr_xs, llhr_cal, llhr_trk, llhr;

  hist->fEleLogLHCal->Fill(Track->EleLogLHCal());
  hist->fMuoLogLHCal->Fill(Track->MuoLogLHCal());

  llhr_cal = Track->LogLHRCal();
  hist->fLogLHRCal->Fill(llhr_cal);

  llhr_dedx = Track->LogLHRDeDx();
  llhr_xs   = Track->LogLHRXs();
  llhr_trk  = Track->LogLHRTrk();
  llhr      = llhr_cal+llhr_trk;

  hist->fEpVsDt->Fill(tp->fDt,tp->fEp);
  hist->fLogLHRDeDx->Fill(llhr_dedx);
  hist->fLogLHRXs->Fill(llhr_xs);
  hist->fLogLHRTrk->Fill(llhr_trk);
  hist->fLogLHR->Fill(llhr);

  hist->fPdgCode->Fill(Track->fPdgCode);
  hist->fFrGH->Fill(Track->fNGoodMcHits/(Track->NActive()+1.e-5));

  hist->fNEPlVsNHPl->Fill(tp->fNEPl,tp->fNHPl);
  hist->fNDPlVsNHPl->Fill(tp->fNDPl,tp->fNHPl);
  hist->fChi2dVsNDPl->Fill(tp->fNDPl,Track->Chi2Dof());
  hist->fDpFVsNDPl  ->Fill(tp->fNDPl,tp->fDpF);

  float        fre1(-1), fre2(-1);
  int          icl;
  TStnCluster* cl;

  if (Track->fVMinS) {
    icl = Track->fVMinS->fClusterIndex;
    if (icl >= 0) {
      cl = fClusterBlock->Cluster(icl);
      fre1 = cl->fFrE1;
      fre2 = cl->fFrE2;
    }
  }

  hist->fFrE1->Fill(fre1);
  hist->fFrE2->Fill(fre2);

  hist->fSinTC->Fill(tp->fSinTC);

  hist->fDrTC->Fill(tp->fDrTC);
  hist->fSInt->Fill(tp->fSInt);
  hist->fDaveTrkQual->Fill(Track->DaveTrkQual());
  hist->fNMcStrawHits->Fill(Track->fNMcStrawHits);
}


//_____________________________________________________________________________
void TTrackAnaModuleA::FillHistograms() {

  double       cos_th (-2.),  cl_e(-1.);
  int          disk_id(-1), alg_mask, best_alg;
  TStnCluster  *cl0(NULL);

  //  cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();

  if (fNClusters > 0) {
    cl0     = fClusterBlock->Cluster(0);
    cl_e    = cl0->Energy();
    disk_id = cl0->DiskID();
  }
//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);

  if (fNTracks[0]> 0) FillEventHistograms(fHist.fEvent[1]);
  else                FillEventHistograms(fHist.fEvent[2]);

  if (fNClusters > 0) FillEventHistograms(fHist.fEvent[3]);
  else                FillEventHistograms(fHist.fEvent[4]);

  if ((fNTracks[0] == 0) && (fabs(cos_th) < 0.4)) {
    FillEventHistograms(fHist.fEvent[5]); 
  }

  if (fNGoodTracks > 0) {
    FillEventHistograms(fHist.fEvent[6]); 

    TLorentzVector    mom(1.,0.,0.,0);
    
    if (fParticle) fParticle->Momentum(mom);

    double p, cos_th;

    p      = mom.P();
    cos_th = mom.Pz()/p;

    if (GetDebugBit(31) && (cos_th > 0.8)) {
      GetHeaderBlock()->Print(Form(" bit:031 cos_th = %10.3f p = %10.3f ntrk = %5i",
				   cos_th, p, fNTracks[0]));
    }
  }

  if ((cl_e > 60.) && (cl0->Time() > 500)) {
    FillEventHistograms(fHist.fEvent[7]); 

    if (GetDebugBit(34)) {
      if (fNCalPatRec <= 0) {
	GetHeaderBlock()->Print(Form(" bit:034 cl_e = %10.3f, cl_time = %10.3f fNCalPatRec = 0",cl_e,cl0->Time()));
      }
    }

    if (GetDebugBit(38)) {
      if ((fNCalPatRec <= 0) && (fNTracks[0] > 0)) {
	GetHeaderBlock()->Print(Form(" bit:038 cl_e = %10.3f, cl_time = %10.3f fNCalPatRec = 0",cl_e,cl0->Time()));
      }
    }

    if (GetDebugBit(40)) {
      if (fNTracks[0] <= 0) {
	GetHeaderBlock()->Print(Form(" bit:040 cl_e = %10.3f, cl_time = %10.3f fNTracks[0] = %3i",
				     cl_e,cl0->Time(),fNTracks[0]));
      }
    }
  }

  if      (disk_id == 0) FillEventHistograms(fHist.fEvent[8]);
  else if (disk_id == 1) FillEventHistograms(fHist.fEvent[9]);
//-----------------------------------------------------------------------------
// EVT_10: events with SetC tracks 103.5 < p < 105.0
//-----------------------------------------------------------------------------
  if (fNGoodTracks > 0) {
    TStnTrack* trk = fTrackBlock->Track(0);
    if ((103.5 < trk->fP) && (trk->fP < 105.)) {
      FillEventHistograms(fHist.fEvent[10]);
    }
  }
//-----------------------------------------------------------------------------
// Simp histograms
//-----------------------------------------------------------------------------
  if (fSimp) {
    FillSimpHistograms(fHist.fSimp[0],fSimp);
  }
//-----------------------------------------------------------------------------
// track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
  TStnTrack*   trk;
  TrackPar_t*  tp;

  for (int i=0; i<fNTracks[0]; ++i ) {
    trk = fTrackBlock->Track(i);
    tp  = fTrackPar+i;

    FillTrackHistograms(fHist.fTrack[0],trk,tp);

    if (trk->fIDWord == 0) {
					// track passes selection "C" 

      FillTrackHistograms(fHist.fTrack[1],trk,tp);

      if (GetDebugBit(32) && (trk->fVMinS != NULL)) {
	if (trk->fVMinS->fChi2Match > 100) {
	  GetHeaderBlock()->Print(Form("bit032: chi2(match) = %10.3lf",trk->fVMinS->fChi2Match));
	}
      }

      if (GetDebugBit(35) && (trk->fP > 106.)) {
	GetHeaderBlock()->Print(Form("bit035: P = %10.3lf",trk->fP));
      }

      if (GetDebugBit(39) && (fabs(tp->fSinTC) < 100) && (fabs(tp->fSinTC) > 0.6)) {
	GetHeaderBlock()->Print(Form("bit039: tp->fSinTC =  %10.3lf",tp->fSinTC));
      }

					// events with at least one  cluster
      if (fNClusters > 0) {
	FillTrackHistograms(fHist.fTrack[2],trk,tp);
      }
      else {
					// events without a cluster
	FillTrackHistograms(fHist.fTrack[3],trk,tp);
      }
//-----------------------------------------------------------------------------
// events with a good track, reconstructed clusters but without a match
//-----------------------------------------------------------------------------
      if ((fNClusters > 0) && (trk->NClusters() == 0)) {
	FillTrackHistograms(fHist.fTrack[4],trk,tp);
      }
//-----------------------------------------------------------------------------
// TRK 5 : events with a good track, reconstructed clusters but without a match
//-----------------------------------------------------------------------------
      if ((trk->fVMaxEp) && (fabs(trk->fVMaxEp->fDt) < 2.5)) {
	FillTrackHistograms(fHist.fTrack[5],trk,tp);
      }
//-----------------------------------------------------------------------------
// TRK 8: good track, |xslope| < 3
//-----------------------------------------------------------------------------
      if (fabs(trk->XSlope()) < 3.) {
	FillTrackHistograms(fHist.fTrack[8],trk,tp);
      }
    }
//-----------------------------------------------------------------------------
// TRK 6 : events with a track and a cluster E/P > 0.4
// TRK 7 : events with a "Set C" track and a cluster E/P > 0.4
//-----------------------------------------------------------------------------
    if ((trk->Ep() > 0.4) && ( trk->Ep() < 1.2)) {
      FillTrackHistograms(fHist.fTrack[6],trk,tp);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[7],trk,tp);
      }
    }
//----------------------------------------------------------------------------
//TRK  9: events with track and with no EM cluster      with E < 60 MeV
//TRK 10: events with good track and with no EM cluster with E < 60 MeV
//TRK 10: events with good track and with no EM cluster with E < 60 MeV
//----------------------------------------------------------------------------
    if (cl_e > 60) {
      FillTrackHistograms(fHist.fTrack[9],trk,tp);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[10],trk,tp);
      }
    }
//-----------------------------------------------------------------------------
// TRK 11: tracks with P > 103.5 MeV
//-----------------------------------------------------------------------------
    if (trk->P() > 103.5) FillTrackHistograms(fHist.fTrack[11],trk,tp);
//-----------------------------------------------------------------------------
// TRK 12: tracks with fcon < 1e-4
// TRK 13: "Set C" tracks with 100 <= P < 110 
// TRK 14: tracks with fcon < 1e-2
//-----------------------------------------------------------------------------
    if (trk->fFitCons < 1.e-4) FillTrackHistograms(fHist.fTrack[12],trk,tp);

    if ((trk->fIDWord == 0) && (trk->P() >= 100.) && (trk->P() < 110.)) {
      FillTrackHistograms(fHist.fTrack[13],trk,tp);
    }

    if (trk->fFitCons < 1.e-2) FillTrackHistograms(fHist.fTrack[14],trk,tp);

    TStnTrack::InterData_t*    vt = trk->fVMinS;  // track-only
//-----------------------------------------------------------------------------
// TRK_15: tracks which have intersection with the 1st disk
// TRK_16: tracks which have intersection with the 2nd disk
// TRK_17: tracks which do not have intersections with the calorimeter
//-----------------------------------------------------------------------------
    if (vt) {
      if      (vt->fID == 0) {
	FillTrackHistograms(fHist.fTrack[15],trk,tp);
      }
      else if (vt->fID == 1) {
	FillTrackHistograms(fHist.fTrack[16],trk,tp);
      }
    }
    else {
      FillTrackHistograms(fHist.fTrack[17],trk,tp);
    }
//-----------------------------------------------------------------------------
// TRK_18: Set "C" tracks with T0 > 700
//-----------------------------------------------------------------------------
    if ((trk->fIDWord == 0) && (trk->T0() > 700)) {
      FillTrackHistograms(fHist.fTrack[18],trk,tp);
    }
//-----------------------------------------------------------------------------
// TRK_19: Set "C" tracks with an associated cluster
//-----------------------------------------------------------------------------
    if ((trk->fIDWord == 0) && (trk->Ep() > 0)) {
      FillTrackHistograms(fHist.fTrack[19],trk,tp);

      if (trk->LogLHRCal() < 0) {
	if (GetDebugBit(29)) {
	  GetHeaderBlock()->Print(Form(" bit:029 LLHR(CAL) = %10.3f ep = %10.3f dt = %10.3f chi2_tcm = %10.3f",
				       trk->LogLHRCal(), trk->Ep(), trk->Dt(), trk->fVMinS->fChi2Match));
	}
      }

      if (tp->fDu < -80.) {
	if (GetDebugBit(33)) {
	  GetHeaderBlock()->Print(Form(" bit:033 DU = %10.3f dv = %10.3f ep = %10.3f dt = %10.3f",
				       tp->fDu, tp->fDv, tp->fEp, trk->Dt()));
	}
      }
    }
//-----------------------------------------------------------------------------
// TRK_20: tracks with >= 20 hits
// TRK_21: tracks with >= 20 hits and Chi2/Ndof < 3
//-----------------------------------------------------------------------------
    if (trk->NActive() >= 20) {
      FillTrackHistograms(fHist.fTrack[20],trk,tp);
      if (trk->Chi2Dof() < 3) {
	FillTrackHistograms(fHist.fTrack[21],trk,tp);
      }
    }
//-----------------------------------------------------------------------------
// TRK 22: Set "C" tracks with an associated cluster and chi2(match) < 100
// TRK 23: Set "C" tracks with an associated cluster and chi2(match) < 100 and LLHR(cal) > 0
//         this is interesting to see which muons are getting misidentified
//-----------------------------------------------------------------------------
    if ((trk->fIDWord == 0) && (tp->fEp > 0) && (tp->fChi2Tcm < 100.)) {
      FillTrackHistograms(fHist.fTrack[22],trk,tp);
      if    (trk->LogLHRCal() > 0) {
	FillTrackHistograms(fHist.fTrack[23],trk,tp);
	
	if (trk->fP < 80.) {
	  if (GetDebugBit(36)) {
	    GetHeaderBlock()->Print(Form(" bit:036 trk p = %10.3f E/P = %10.3f", trk->fP,tp->fEp));
	  }
	}
      }
      else {
//-----------------------------------------------------------------------------
// TRK_24: Set "C" tracks with an associated cluster and chi2(match) < 100 and LLHR(cal) < 0
//         this set allows to see which electrons are getting misidentified
//-----------------------------------------------------------------------------
	FillTrackHistograms(fHist.fTrack[24],trk,tp);
      }
//-----------------------------------------------------------------------------
// TRK_25: Set "C" tracks, 100 < P < 110, 0 < E/p < 1.15,  |dt_corr| < 3, chi2(match) < 100
//-----------------------------------------------------------------------------
//      double dt_corr = trk->Dt(); // -1.;
      if ( (tp->fDt > fMinDtTcm) && (tp->fDt < fMaxDtTcm) && 
	   (tp->fEp <      1.15) && 
	   (trk->fP >      100.) && (trk->fP <       110)   ) {
	FillTrackHistograms(fHist.fTrack[25],trk,tp);
//-----------------------------------------------------------------------------
// more details on the calorimeter-based likelihood 
// TRK_26 : TRK 25 events with LLHR_CAL > 0 ( interesting for muons)
// TRK_27 : TRK 25 events with LLHR_CAL < 0 ( interesting for electrons)
//-----------------------------------------------------------------------------
	double llhr_cal = trk->LogLHRCal();
	if (llhr_cal > 0) {
	  FillTrackHistograms(fHist.fTrack[26],trk,tp);

	  if (llhr_cal > 5.) {
	    if (GetDebugBit(37)) {
	      GetHeaderBlock()->Print(Form("bit:037 trk p = %10.3f E/P = %10.3f", trk->fP,tp->fEp));
	    }
	  }
	}
	else {
	  FillTrackHistograms(fHist.fTrack[27],trk,tp);
	}
      }

    }
//-----------------------------------------------------------------------------
// TRK_28 : events with a "Set C" track and a cluster E/P > 1.1
//-----------------------------------------------------------------------------
    if (trk->fIDWord == 0) {
      if (tp->fEp > 1.1) {
	FillTrackHistograms(fHist.fTrack[28],trk,tp);

	if (GetDebugBit(28)) {
	  GetHeaderBlock()->Print(Form(" bit:028 LLHR(CAL) = %10.3f ep = %10.3f dt = %10.3f",
				       trk->LogLHRCal(), trk->Ep(), trk->Dt()));
	}
      }
    }
//-----------------------------------------------------------------------------
// TRK_29: "Set C" track 100 < P < 110, 0 < E/P < 1.15              : next to TRK_13
// TRK_32: "Set C" track 100 < P < 110, 0 < E/P < 1.15, chi2tcm<100 : next to TRK_29
//-----------------------------------------------------------------------------
    if (trk->fIDWord == 0) {
      if ((tp->fEp > 0) && (trk->P() > 100.) && (trk->P() < 110.) && (tp->fEp < 1.15)) {
	FillTrackHistograms(fHist.fTrack[29],trk,tp);

	if (tp->fChi2Tcm < 100.) {
	  FillTrackHistograms(fHist.fTrack[32],trk,tp);
	}
      }
    }
//-----------------------------------------------------------------------------
// TRK_30: tracks with >= 25 hits
// TRK_31: tracks with >= 25 hits and Chi2/Ndof < 3
//-----------------------------------------------------------------------------
    if (trk->NActive() >= 25) {
      FillTrackHistograms(fHist.fTrack[30],trk,tp);
      if (trk->Chi2Dof() < 3) {
	FillTrackHistograms(fHist.fTrack[31],trk,tp);
      }
    }
//-----------------------------------------------------------------------------
// TRK_33: "DaveTrkQual" Q > 0.4 tracks
//-----------------------------------------------------------------------------
    if (tp->fIDWord[fBestID] == 0) {
      FillTrackHistograms(fHist.fTrack[33],trk,tp);
    }
//-----------------------------------------------------------------------------
// TRK_34: "DaveTrkQual" Q > 0.1 tracks
//-----------------------------------------------------------------------------
    if (tp->fIDWord[1] == 0) {
      FillTrackHistograms(fHist.fTrack[34],trk,tp);
    }
//-----------------------------------------------------------------------------
// split tracks by the algorithm mask: 1 , 2 , or 3
//-----------------------------------------------------------------------------
    alg_mask = trk->AlgMask();
    best_alg = trk->BestAlg();
    if      (alg_mask == 1) {
//-----------------------------------------------------------------------------
// TrkPatRec-only tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[40],trk,tp);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[41],trk,tp);
	// print run/event numbers :
	if (GetDebugBit(6)) {
	  double ep = trk->Ep();
	  if ((ep > 0.8) && (ep < 1.1)) {
	    GetHeaderBlock()->Print(Form(" bit:006 trk_41: track E/P = %8.3f",trk->Ep()));
	  }
	}
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[42],trk,tp);
      }
    }
    else if (alg_mask == 2) {
//-----------------------------------------------------------------------------
// CalPatRec-only tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[50],trk,tp);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[51],trk,tp);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[52],trk,tp);
      }
    }
    else if (alg_mask == 3) {
//-----------------------------------------------------------------------------
// TrkPatRec+CalPatRec tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[60],trk,tp);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[61],trk,tp);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[62],trk,tp);
      }
    }
//-----------------------------------------------------------------------------
// TRK_63 : fIDWORD[bestID] == 0, TRkPatRec
// TRK_64 : fIDWORD[bestID] == 0, CalPatRec
//-----------------------------------------------------------------------------
    if (tp->fIDWord[fBestID] == 0) {
      if      (best_alg == 0) FillTrackHistograms(fHist.fTrack[63],trk,tp);
      else if (best_alg == 1) FillTrackHistograms(fHist.fTrack[64],trk,tp);
    }
//-----------------------------------------------------------------------------
// TRK_71: SetC tracks  103.5 < p < 105 : DIO studies
//-----------------------------------------------------------------------------
    if ((trk->fIDWord == 0) && (trk->fP > 103.5) && (trk->fP < 105.)) {
      FillTrackHistograms(fHist.fTrack[71],trk,tp);
    }
  }

//-----------------------------------------------------------------------------
// Track Reco efficiency, Dave style
//-----------------------------------------------------------------------------
  FillEfficiencyHistograms(fTrackBlock,fTrackID[fBestID],11);
//-----------------------------------------------------------------------------
// cluster histograms 
//-----------------------------------------------------------------------------
  TStnCluster*  cl;
  int           id;
  for (int i=0; i<fNClusters; ++i ) {
    cl = fClusterBlock->Cluster(i);
    id = cl->DiskID();
    FillClusterHistograms(fHist.fCluster[0],cl);

    if (fNTracks[0]     >  0 ) FillClusterHistograms(fHist.fCluster[1],cl);
    if (fNGoodTracks    >  0 ) FillClusterHistograms(fHist.fCluster[2],cl);
    if (fNMatchedTracks >  0 ) FillClusterHistograms(fHist.fCluster[3],cl);
    if (cl->Energy()    > 10.) FillClusterHistograms(fHist.fCluster[4],cl);
    if (cl->Energy()    > 60.) FillClusterHistograms(fHist.fCluster[5],cl);

    if      (id == 0         ) FillClusterHistograms(fHist.fCluster[6],cl);
    else if (id == 1         ) FillClusterHistograms(fHist.fCluster[7],cl);
  }
//-----------------------------------------------------------------------------
// calorimeter histograms
//-----------------------------------------------------------------------------
  TDisk*         disk;
  TStnCrystal*   cr;

  if (fCalorimeterType == 2) {
    int nd = fDiskCalorimeter->NDisks();

    for (int i=0; i<nd; i++) {
      disk = fDiskCalorimeter->Disk(i);
      for (int ic=0; ic<disk->NCrystals(); ic++) {
	cr = disk->Crystal(ic);
	FillCaloHistograms(fHist.fCalo[0],cr);

	if (cr->Energy() > 0) {
	  FillCaloHistograms(fHist.fCalo[1],cr);
	}
	if (cr->Energy() > 0.1) {
	  FillCaloHistograms(fHist.fCalo[2],cr);
	}
	if (cr->Energy() > 1.0) {
	  FillCaloHistograms(fHist.fCalo[3],cr);
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// radial distributions for crystals
//-----------------------------------------------------------------------------
  static int first_entry(1);

  if (first_entry == 1) {
    first_entry = 0;

    if (fCalorimeterType == 2) {
      int nd = fDiskCalorimeter->NDisks();
	
      for (int i=0; i<nd; i++) {
	disk = fDiskCalorimeter->Disk(i);
	for (int ic=0; ic<disk->NCrystals(); ic++) {
	  cr = disk->Crystal(ic);

	  fHist.fCrystalR[i]->Fill(cr->Radius());
	}
      }
    }
  }
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
int TTrackAnaModuleA::Event(int ientry) {

  TDiskCalorimeter::GeomData_t disk_geom;

  fTrackBlock  ->GetEntry(ientry);
  fClusterBlock->GetEntry(ientry);
  fStrawDataBlock->GetEntry(ientry);
  fCalDataBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fVDetBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// luminosity weight
//-----------------------------------------------------------------------------
  fLumWt = GetHeaderBlock()->LumWeight();
//-----------------------------------------------------------------------------
// look for signal particle defined by the PDG code and the generator code
//-----------------------------------------------------------------------------
  fNGenp    = fGenpBlock->NParticles();

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
  fSimPar.fGenp     = fParticle;
//-----------------------------------------------------------------------------
// process virtual detectors - for fSimp need parameters at tracker entrance
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

  if (fDiskCalorimeter->Initialized() == 0) {
    disk_geom.fNDisks = fCalDataBlock->NDisks();

    for (int i=0; i<disk_geom.fNDisks; i++) {
      //      disk_geom.fNCrystals[i] = fCalDataBlock->fNCrystals[i];
      disk_geom.fRMin[i]      = fCalDataBlock->fRMin[i];
      disk_geom.fRMax[i]      = fCalDataBlock->fRMax[i];
      disk_geom.fZ0  [i]      = fCalDataBlock->fZ0  [i];
    }

    disk_geom.fHexSize          = fCalDataBlock->CrystalSize()*2;
    // kludge , so far
    disk_geom.fMinFraction      = 1.; // fCalDataBlock->MinFraction();
    disk_geom.fWrapperThickness = fCalDataBlock->WrapperThickness();
    disk_geom.fShellThickness   = fCalDataBlock->ShellThickness();

    fDiskCalorimeter->Init(&disk_geom);
  }

  fNTracks[0] = fTrackBlock->NTracks();
  if (fNTracks[0] == 0) fTrack = 0;
  else                  fTrack = fTrackBlock->Track(0);

  fNClusters  = fClusterBlock->NClusters();
  fNCalHits   = fCalDataBlock->NHits();
  fNStrawHits = GetHeaderBlock()->fNStrawHits;

  fDiskCalorimeter->InitEvent(fCalDataBlock);
//-----------------------------------------------------------------------------
// loop over tracks and calculate needed parameters
//-----------------------------------------------------------------------------
  fNHyp           = -1;
  fBestHyp[0]     = -1;
  fBestHyp[1]     = -1;
//-----------------------------------------------------------------------------
// determine the number of CalPatRec tracks - this assumes that the list of 
// tracks has been created by MergePatRec
//-----------------------------------------------------------------------------
  TStnTrack*   track;
  int ntrk = fTrackBlock->NTracks();
  fNCalPatRec = 0;
  for (int itrk=0; itrk<ntrk; itrk++) {
    track        = fTrackBlock->Track(itrk);
    int alg_mask = track->AlgMask();
    if (alg_mask & 0x2) fNCalPatRec += 1;
  }
  
  InitTrackPar(fTrackBlock,fClusterBlock,fTrackPar,fNGoodTracks,fNMatchedTracks);
//-----------------------------------------------------------------------------
// init calorimeter clusters - remember, the first one not necessarily is the 
// most energetic
//-----------------------------------------------------------------------------
  fNClusters = fClusterBlock->NClusters();
  if (fNClusters == 0) fCluster = 0;
  else                 fCluster = fClusterBlock->Cluster(0);

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TTrackAnaModuleA::Debug() {

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
	GetHeaderBlock()->Print(Form("bit_000: All p = %10.3f",
				     trk->Momentum()->P()));
    }
//-----------------------------------------------------------------------------
// bit 3: Set C tracks with large DX : 70mm < |DX| < 90mm
//-----------------------------------------------------------------------------
    if (GetDebugBit(3) == 1) {
      if (trk->fIDWord == 0) {
	TStnTrack::InterData_t*    vr = trk->fVMaxEp; // residuals
	if ((vr && (fabs(vr->fDx) > 70) && (fabs(vr->fDx) < 90))) {
	  GetHeaderBlock()->Print(Form("large DX: %f",vr->fDx));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 4: tracks with DpF > 1MeV - positive tail...
//-----------------------------------------------------------------------------
    if (GetDebugBit(4) == 1) {
      if (tp->fDpF > 1.) {
	GetHeaderBlock()->Print(Form("pF pRec, fDpf = %10.3f  %10.3f  %10.3f",
				     trk->fPFront, trk->Momentum()->P(),tp->fDpF));
      }
    }
//-----------------------------------------------------------------------------
// bit 9: Set C tracks with DpF > 1MeV - positive tail...
//-----------------------------------------------------------------------------
    if (GetDebugBit(9) == 1) {
      double ep = trk->Ep();
      if (trk->fIDWord == 0) { 
	if (((ep > 0.42) && (ep < 0.46)) || ((ep > 0.35) && (ep < 0.39))) {
	  GetHeaderBlock()->Print(Form("bit:009 ep = %10.3f e = %10.3f p = %10.3f",
				       trk->fEp,trk->fEp*trk->fP,trk->fP));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 10: Set C tracks with Ecl > 80
//-----------------------------------------------------------------------------
    if (GetDebugBit(10) == 1) {
      double ecl = trk->ClusterE();
      if (trk->fIDWord == 0) { 
	if (ecl > 60) {
	  GetHeaderBlock()->Print(Form("bit:010 e = %10.3f p = %10.3f",
				       ecl,trk->fP));
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
int TTrackAnaModuleA::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TTrackAnaModuleA::Test001() {
}

