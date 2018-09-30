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
#include "DataProducts/inc/VirtualDetectorId.hh"

//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------

#include "ana/TTrackAnaModule.hh"

ClassImp(TTrackAnaModule)
//-----------------------------------------------------------------------------
TTrackAnaModule::TTrackAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fTrackBlockName         = "TrackBlock";
  fTrackStrawHitBlockName = "TrackStrawHitBlock";
  fTrackNumber.Set(100);

  fDiskCalorimeter = new TDiskCalorimeter();
  fCalorimeterType = 2;
  fFillDioHist     = 1;

  fMinT0           = 700; 
  fMinETrig        = 50.;               // MeV
					// track-cluster matching timing cut
  fMinDtTcm        = -5.;
  fMaxDtTcm        =  8.;
//-----------------------------------------------------------------------------
// initialize Track ID
// 0: SetC  1: TrkQual>0.1 2:TrkQual>0.4
// what about number of hits ? - 3: no cuts on the number of hits
//-----------------------------------------------------------------------------
  fNID             = 4;
  for (int i=0; i<fNID; i++) {
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

  fBestID     = 3;			// best: DaveTrkQual > 0.4, no NActive cut
//-----------------------------------------------------------------------------
// PID: initialize likelihood histograms
// TRK 19: "Set C" plus reconstructed and matched cluster
//-----------------------------------------------------------------------------
  fLogLH      = new TEmuLogLH  ();
//   const char   *pid_version;
//   pid_version = gEnv->GetValue("mu2e.PidVersion","_none_");
//   fLogLH->Init(pid_version);
//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
// 2:conversionGun, 28:StoppedParticleReactionGun - see 
//-----------------------------------------------------------------------------
  fPdgCode       = 11;
  fGeneratorCode =  2;
  fDirection     =  1;

  fApplyCorrections = 0;
}

//-----------------------------------------------------------------------------
TTrackAnaModule::~TTrackAnaModule() {
  for (int i=0; i<fNID; i++) delete fTrackID[i];
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TTrackAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockName.Data()        ,"TStnTrackBlock"     ,&fTrackBlock  );
  RegisterDataBlock(fTrackStrawHitBlockName.Data(),"TTrackStrawHitBlock",&fTrackStrawHitBlock);

  RegisterDataBlock("ClusterBlock"        ,"TStnClusterBlock" ,&fClusterBlock);
  RegisterDataBlock("CalDataBlock"        ,"TCalDataBlock"    ,&fCalDataBlock);
  RegisterDataBlock("StrawDataBlock"      ,"TStrawDataBlock"  ,&fStrawDataBlock);
  RegisterDataBlock("GenpBlock"           ,"TGenpBlock"       ,&fGenpBlock);
  RegisterDataBlock("SimpBlock"           ,"TSimpBlock"       ,&fSimpBlock);
  RegisterDataBlock("VdetBlock"           ,"TVdetDataBlock"   ,&fVdetBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//_____________________________________________________________________________
int TTrackAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//-----------------------------------------------------------------------------
// need MC truth branch
//-----------------------------------------------------------------------------
void TTrackAnaModule::FillEventHistograms(EventHist_t* Hist) {
  double            cos_th, dio_wt, xv(-1.e6), yv(-1.e6), rv(-1.e6), zv(-1.e6), p;
  double            e, m, r;
  TLorentzVector    mom(1.,0.,0.,0);

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

  Hist->fLumWt->Fill(fLumWt);
  Hist->fRv->Fill(rv);
  Hist->fZv->Fill(zv);
  Hist->fEleMom->Fill(p);
  Hist->fDioMom->Fill(p,dio_wt);
  Hist->fEleCosTh->Fill(cos_th);

  Hist->fNClusters->Fill(fNClusters);
  Hist->fNTracks->Fill  (fNTracks[0]);

  int nsh = GetHeaderBlock()->fNStrawHits;
  Hist->fNStrawHits[0]->Fill(nsh);
  Hist->fNStrawHits[1]->Fill(nsh);

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
  int    alg_mask(-1);
  if (track) {
    t0_trk   = track->fT0;
    alg_mask = track->AlgMask();
  }

  if (track && cluster) {
    dt = t0_cls-t0_trk;
  }

  Hist->fAlgMask->Fill(alg_mask);
  Hist->fDtClT->Fill(dt);
  Hist->fEMax->Fill(emax);

  TStrawHitData*  sh;
  int n_good_hits = 0;
  int nstraw_hits = fStrawDataBlock->NHits();
  for (int i=0; i<nstraw_hits; i++ ) {
    sh = fStrawDataBlock->Hit(i);
    dt = t0_cls-sh->Time() + 15;
    Hist->fDtClS->Fill(dt);
    Hist->fSHTime->Fill(sh->Time());

    if (fabs(dt+15.)< 50) n_good_hits += 1;
  }

  Hist->fNGoodSH->Fill(n_good_hits);

  Hist->fNHyp->Fill(fNHyp);
  Hist->fBestHyp[0]->Fill(fBestHyp[0]);
  Hist->fBestHyp[1]->Fill(fBestHyp[1]);
  Hist->fNGenp->Fill(fNGenp);
//-----------------------------------------------------------------------------
// crystals - count crystals with E > 1MeV
//-----------------------------------------------------------------------------
  int          ndisks, n_hit_crystals[kNDisks], n_hit_crystals_tot;
  double       etot[kNDisks];
  //  TCalHitData* cch;

  //  int n_cch_1mev = 0;

  //  int   nhits_vane[2][kNDisks], nhits_row [2][20], nhits_col[2][50];
  //  int   crystal_id, vane_id, local_id, vane_row, vane_col;
  int   bin, hit_id, idisk, nhits;
  int   nhits_r[kNDisks][100], n_hit_crystals_r[kNDisks][100];

  if (fCalorimeterType == 2) {
//-----------------------------------------------------------------------------
// disk calorimeter
//-----------------------------------------------------------------------------

    TCalHitData* hit;

    ndisks = fDiskCalorimeter->NDisks();

    if (ndisks > 0) {

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
	
	Hist->fECrVsR[idisk]->Fill(r,e);
	Hist->fNCrVsR[idisk]->Fill(r,1);
	
	bin  = (int) (r/10.);
	
	nhits_r         [idisk][bin] += 1;
//-----------------------------------------------------------------------------
// this is not correct, one needs to check whether this crystal has been hit,
// for the moment, to get going, ignore that
//-----------------------------------------------------------------------------
	n_hit_crystals_r[idisk][bin] += 1;
      }
    }

    n_hit_crystals_tot = 0;

    double ecal = 0;
    for (int id=0; id<ndisks; id++) {
      n_hit_crystals_tot += n_hit_crystals[id];
      ecal += etot[id];
//-----------------------------------------------------------------------------
// fill 'per-disk' histograms
//-----------------------------------------------------------------------------
      Hist->fETot[id]->Fill(etot[id]);

//-----------------------------------------------------------------------------
// 100 is fixed by the number of bins in the radial distributions
//-----------------------------------------------------------------------------
      for (int ib=0; ib<100; ib++) {
	r = (ib+0.5)*10.;
	Hist->fNCrystalHitsVsR[id]->Fill(r,nhits_r         [id][ib]);
	Hist->fNHitCrystalsVsR[id]->Fill(r,n_hit_crystals_r[id][ib]);
      }
    }

    Hist->fNHitCrystalsTot->Fill(n_hit_crystals_tot);
    Hist->fECal->Fill(ecal);

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
    Hist->fECalOverEKin->Fill(ecal/ekin);
  }

  Hist->fInstLumi->Fill(GetHeaderBlock()->InstLum());
}

//-----------------------------------------------------------------------------
void TTrackAnaModule::FillCaloHistograms(CaloHist_t* Hist, TStnCrystal* Cr) {

  int                    nhits;
  float                  t, e, r, e700, n700;
  TCalHitData*           hit;

  // determine crystal coordinates

  TDisk* disk = Cr->Disk();

  int idisk = disk->SectionID();
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
void TTrackAnaModule::FillClusterHistograms(ClusterHist_t* Hist, TStnCluster* Cluster) {
  int   row, col;
  float  x, y, z, r;

  row = Cluster->Ix1();
  col = Cluster->Ix2();

  x   = Cluster->fX+3904.;
  y   = Cluster->fY;
  z   = Cluster->fZ;
  r   = sqrt(x*x+y*y);

  if ((row < 0) || (row > 9999)) row = -9999;
  if ((col < 0) || (col > 9999)) col = -9999;

  Hist->fDiskID->Fill(Cluster->DiskID());
  Hist->fEnergy->Fill(Cluster->Energy());
  Hist->fT0->Fill(Cluster->Time());
  Hist->fRow->Fill(row);
  Hist->fCol->Fill(col);
  Hist->fX->Fill(x);
  Hist->fY->Fill(y);
  Hist->fZ->Fill(z);
  Hist->fR->Fill(r);

  Hist->fYMean->Fill(Cluster->fYMean);
  Hist->fZMean->Fill(Cluster->fZMean);
  Hist->fSigY->Fill(Cluster->fSigY);
  Hist->fSigZ->Fill(Cluster->fSigZ);
  Hist->fSigR->Fill(Cluster->fSigR);
  Hist->fNCr0->Fill(Cluster->fNCrystals);
  Hist->fNCr1->Fill(Cluster->fNCr1);
  Hist->fFrE1->Fill(Cluster->fFrE1);
  Hist->fFrE2->Fill(Cluster->fFrE2);
  Hist->fSigE1->Fill(Cluster->fSigE1);
  Hist->fSigE2->Fill(Cluster->fSigE2);
}

//-----------------------------------------------------------------------------
void TTrackAnaModule::FillGenpHistograms(GenpHist_t* Hist, TGenParticle* Genp) {
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
void TTrackAnaModule::FillSimpHistograms(SimpHist_t* Hist, TSimParticle* Simp) {

  Hist->fPdgCode->Fill(Simp->fPdgCode);
  Hist->fMomTargetEnd->Fill(Simp->fMomTargetEnd);
  Hist->fMomTrackerFront->Fill(Simp->fMomTrackerFront);
  Hist->fNStrawHits->Fill(Simp->fNStrawHits);
}

//-----------------------------------------------------------------------------
// for DIO : ultimately, one would need to renormalize the distribution
//-----------------------------------------------------------------------------
void TTrackAnaModule::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track) {

  TLorentzVector  mom;
  double          r;
  int             itrk;
  TrackPar_t*     tp;
					// pointer to local track parameters
  itrk = Track->Number();
  tp   = fTrackPar+itrk;

  Hist->fP[0]->Fill (tp->fP);		// corrected momentum in the first point
  Hist->fP[1]->Fill (tp->fP);
  Hist->fP[2]->Fill (tp->fP);

  Hist->fP0->  Fill (Track->fP0);	// momentum zt Z0
  Hist->fP2->  Fill (Track->fP2);	// momentum at z2 = -1540, the same as in the 1st point

  Hist->fPDio->Fill  (tp->fP, tp->fDioWt);
  Hist->fPlw->Fill   (tp->fP, tp->fLumWt);
  Hist->fPDiolw->Fill(tp->fP, tp->fTotWt);


  Hist->fFitMomErr->Fill(Track->fFitMomErr);

  Hist->fPt    ->Fill(Track->fPt    );
  Hist->fPFront->Fill(Track->fPFront);
  Hist->fPStOut->Fill(Track->fPStOut);
					// dp: Tracker-only resolution

  Hist->fDpFront ->Fill(tp->fDpF);
  Hist->fXDpF    ->Fill(tp->fDpF/Track->fFitMomErr);
  Hist->fDpFDio  ->Fill(tp->fDpF,tp->fDioWt);
  Hist->fXDpFDio ->Fill(tp->fDpF/Track->fFitMomErr,tp->fDioWt);
  Hist->fDpFront0->Fill(tp->fDp0);
  Hist->fDpFront2->Fill(tp->fDp2);
  Hist->fDpFSt   ->Fill(tp->fDpFSt);
  Hist->fDpFVsZ1 ->Fill(Track->fZ1,tp->fDpF);

  Hist->fCosTh->Fill(Track->Momentum()->CosTheta());
  Hist->fChi2->Fill (Track->fChi2);

  float na = Track->NActive();

  Hist->fChi2Dof->Fill(Track->fChi2/(na-5.));
  Hist->fNActive->Fill(na);
  Hist->fNaFract->Fill(na/(Track->NHits()+0.));
  Hist->fNWrong->Fill(Track->NWrong());

  float nd  = Track->NDoublets();
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
  //  printf("TTrackAnaModule::FillTrackHistograms: track charge is not defined yet\n");
  Hist->fQ->Fill(-1);
  Hist->fFitCons[0]->Fill(Track->fFitCons);
  Hist->fFitCons[1]->Fill(Track->fFitCons);

  Hist->fD0->Fill(Track->fD0);
  Hist->fZ0->Fill(Track->fZ0);
  Hist->fTanDip->Fill(Track->fTanDip);
  Hist->fDtZ0->Fill(tp->fDtZ0);
  Hist->fDtBack->Fill(tp->fDtBack);
  Hist->fRMax->Fill(Track->RMax());
  Hist->fAlgMask->Fill(Track->AlgMask());
//-----------------------------------------------------------------------------
// track-cluster matching part: 
// - for residuals, determine intersection with the most energetic cluster
// - for track -only parameters use intersection with lowest trajectory length
//-----------------------------------------------------------------------------
  TStnTrack::InterData_t*    vt = Track->fVMinS;  // track-cluster residuals

  if (vt) {
    Hist->fDiskID->Fill(vt->fID  );
    Hist->fXTrk->Fill  (vt->fXTrk);
    Hist->fYTrk->Fill  (vt->fYTrk);

    r = sqrt(vt->fXTrk*vt->fXTrk+vt->fYTrk*vt->fYTrk);
    Hist->fRTrk->Fill  (r);

    Hist->fZTrk->Fill  (vt->fZTrk);
  }
  else {
//-----------------------------------------------------------------------------
// fill histograms with numbers easy to recognize as dummy
//-----------------------------------------------------------------------------
    Hist->fDiskID->Fill(-1.);
    Hist->fXTrk->Fill  (999.);
    Hist->fYTrk->Fill  (999.);
    Hist->fRTrk->Fill  (999.);
    Hist->fZTrk->Fill  (-1. );
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

  Hist->fECl->Fill(tp->fEcl);
  Hist->fEClEKin->Fill(tp->fEcl/ekin);

  Hist->fECalP->Fill(tp->fECalP);
  Hist->fEDiskP->Fill(tp->fEDiskP);
  Hist->fEp->Fill(tp->fEp);
  Hist->fEpVsPath->Fill(tp->fPath,tp->fEp);

  Hist->fDx->Fill(tp->fDx);
  Hist->fDy->Fill(tp->fDy);
  Hist->fDz->Fill(tp->fDz);

  Hist->fDt->Fill(tp->fDt);
  Hist->fChi2Tcm->Fill(tp->fChi2Tcm);
  Hist->fChi2XY->Fill(tp->fChi2XY);
  Hist->fChi2T->Fill (tp->fChi2T);

  Hist->fDu->Fill    (tp->fDu);
  Hist->fDv->Fill    (tp->fDv);
  Hist->fDvVsDu->Fill(tp->fDu,tp->fDv);

  Hist->fPath->Fill(tp->fPath);
  Hist->fDuVsPath->Fill(tp->fPath,tp->fDu);

  Hist->fDvVsPath->Fill(tp->fPath,tp->fDv);
  Hist->fDtVsPath->Fill(tp->fPath,tp->fDt);

  Hist->fDuVsTDip->Fill(Track->fTanDip,tp->fDu);
  Hist->fDvVsTDip->Fill(Track->fTanDip,tp->fDv);

  Hist->fZ1->Fill(Track->fZ1);

  int ncl = Track->NClusters();
  Hist->fNClusters->Fill(ncl);

  Hist->fRSlope->Fill(Track->RSlope());
  Hist->fXSlope->Fill(Track->XSlope());

  double llhr_dedx, llhr_xs, llhr_cal, llhr_trk, llhr;

  Hist->fEleLogLHCal->Fill(Track->EleLogLHCal());
  Hist->fMuoLogLHCal->Fill(Track->MuoLogLHCal());

  llhr_cal = Track->LogLHRCal();
  Hist->fLogLHRCal->Fill(llhr_cal);

  llhr_dedx = Track->LogLHRDeDx();
  llhr_xs   = Track->LogLHRXs();
  llhr_trk  = Track->LogLHRTrk();
  llhr      = llhr_cal+llhr_trk;

  Hist->fEpVsDt->Fill(tp->fDt,tp->fEp);
  Hist->fLogLHRDeDx->Fill(llhr_dedx);
  Hist->fLogLHRXs->Fill(llhr_xs);
  Hist->fLogLHRTrk->Fill(llhr_trk);
  Hist->fLogLHR->Fill(llhr);
					// two histograms with different limits and binning
  Hist->fPdgCode[0]->Fill(Track->fPdgCode);
  Hist->fPdgCode[1]->Fill(Track->fPdgCode);
  Hist->fFrGH->Fill(Track->fNGoodMcHits/(Track->NActive()+1.e-5));

  Hist->fNEPlVsNHPl->Fill(tp->fNEPl,tp->fNHPl);
  Hist->fNDPlVsNHPl->Fill(tp->fNDPl,tp->fNHPl);
  Hist->fChi2dVsNDPl->Fill(tp->fNDPl,Track->Chi2Dof());
  Hist->fDpFVsNDPl  ->Fill(tp->fNDPl,tp->fDpF);

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

  Hist->fFrE1->Fill(fre1);
  Hist->fFrE2->Fill(fre2);

  Hist->fSinTC->Fill(tp->fSinTC);

  Hist->fDrTC->Fill(tp->fDrTC);
  Hist->fSInt->Fill(tp->fSInt);
  Hist->fDaveTrkQual->Fill(Track->DaveTrkQual());
  Hist->fNMcStrawHits->Fill(Track->fNMcStrawHits);

  if (fTrackStrawHitBlock->NTracks() > 0) {
    
    int nhits = fTrackStrawHitBlock->NTrackHits(itrk);

    for (int i=0; i<nhits; i++) {
      TTrackStrawHitData* hit = fTrackStrawHitBlock->Hit(itrk,i);
      Hist->fHitEnergy->Fill(hit->Energy());
      Hist->fHitDt->Fill(hit->Dt());
      Hist->fHitTRel->Fill(hit->Time()-tp->fTMean);
    }
  }
  
  Hist->fT0MinusTM->Fill(Track->T0()-tp->fTMean);
}


//-----------------------------------------------------------------------------
// fill efficiency histograms : need 10 histogram sets
// pitch = 1./tan(dip)
//-----------------------------------------------------------------------------
void TTrackAnaModule::FillEfficiencyHistograms(TStnTrackBlock*  TrackBlock, 
					       TStnTrackID*     TrackID   , 
					       int              HistSet   ) {
  // having MC truth is a must for calculating efficiency!
  if (fSimp == nullptr) return;

  if (fSimp->NStrawHits() >= 20) {
    FillEventHistograms(fHist.fEvent[HistSet]);

    if (fSimp->fMomTrackerFront > 100.) {
      FillEventHistograms(fHist.fEvent[HistSet+1]);

      TLorentzVector vdmom;

      if (fSimPar.fTFront != NULL) {
//-----------------------------------------------------------------------------
// regular case
//-----------------------------------------------------------------------------
	vdmom.SetXYZM(fSimPar.fTFront->McMomentumX(),
		      fSimPar.fTFront->McMomentumY(),		      
		      fSimPar.fTFront->McMomentumZ(),
		      fSimPar.fTFront->Mass());
      }
      else {
//-----------------------------------------------------------------------------
// pathologicalcase when an upstream simulated MC particle is reconstructed 
// as the downstream one and there is no hit... efficiency doens't make sense
// just make sure the code doesnt' crash
//-----------------------------------------------------------------------------
	vdmom = fSimPar.fParticle->fStartMom;
      }

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

		    if ((103.5 < track->fP) && (track->fP < 105)) {
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

//_____________________________________________________________________________
void TTrackAnaModule::FillHistograms() {

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

  if ((cl_e > fMinETrig) && (cl0->Time() > 500)) {
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

    FillTrackHistograms(fHist.fTrack[0],trk);

    if (trk->fIDWord == 0) {
					// track passes selection "C" 

      FillTrackHistograms(fHist.fTrack[1],trk);

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
	FillTrackHistograms(fHist.fTrack[2],trk);
      }
      else {
					// events without a cluster
	FillTrackHistograms(fHist.fTrack[3],trk);
      }
//-----------------------------------------------------------------------------
// events with a good track, reconstructed clusters but without a match
//-----------------------------------------------------------------------------
      if ((fNClusters > 0) && (trk->NClusters() == 0)) {
	FillTrackHistograms(fHist.fTrack[4],trk);
      }
//-----------------------------------------------------------------------------
// TRK 5 : events with a good track, reconstructed clusters but without a match
//-----------------------------------------------------------------------------
      if ((trk->fVMaxEp) && (fabs(trk->fVMaxEp->fDt) < 2.5)) {
	FillTrackHistograms(fHist.fTrack[5],trk);
      }
//-----------------------------------------------------------------------------
// TRK 8: good track, |xslope| < 3
//-----------------------------------------------------------------------------
      if (fabs(trk->XSlope()) < 3.) {
	FillTrackHistograms(fHist.fTrack[8],trk);
      }
    }
//-----------------------------------------------------------------------------
// TRK 6 : events with a track and a cluster E/P > 0.4
// TRK 7 : events with a "Set C" track and a cluster E/P > 0.4
//-----------------------------------------------------------------------------
    if ((trk->Ep() > 0.4) && ( trk->Ep() < 1.2)) {
      FillTrackHistograms(fHist.fTrack[6],trk);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[7],trk);
      }
    }
//----------------------------------------------------------------------------
//TRK  9: events with track and with no EM cluster      with E < 60 MeV
//TRK 10: events with good track and with no EM cluster with E < 60 MeV
//TRK 10: events with good track and with no EM cluster with E < 60 MeV
//----------------------------------------------------------------------------
    if (cl_e > 60) {
      FillTrackHistograms(fHist.fTrack[9],trk);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[10],trk);
      }
    }
//-----------------------------------------------------------------------------
// TRK 11: tracks with P > 103.5 MeV
//-----------------------------------------------------------------------------
    if (trk->P() > 103.5) FillTrackHistograms(fHist.fTrack[11],trk);
//-----------------------------------------------------------------------------
// TRK_12: tracks with fcon < 1e-4
// TRK_13: "Set C" tracks with 100 <= P < 110 
// TRK_14: tracks with fcon < 1e-2
//-----------------------------------------------------------------------------
    if (trk->fFitCons < 1.e-4) FillTrackHistograms(fHist.fTrack[12],trk);

    if ((trk->fIDWord == 0) && (trk->P() >= 100.) && (trk->P() < 110.)) {
      FillTrackHistograms(fHist.fTrack[13],trk);
    }

    if (trk->fFitCons < 1.e-2) FillTrackHistograms(fHist.fTrack[14],trk);

    TStnTrack::InterData_t*    vt = trk->fVMinS;  // track-only
//-----------------------------------------------------------------------------
// TRK_15: all tracks which have intersection with the 1st disk
// TRK_16: all tracks which have intersection with the 2nd disk
// TRK_17: all tracks which do not have intersections with the calorimeter
//-----------------------------------------------------------------------------
    if (vt) {
      if      (vt->fID == 0) FillTrackHistograms(fHist.fTrack[15],trk);
      else if (vt->fID == 1) FillTrackHistograms(fHist.fTrack[16],trk);
    }
    else {
      FillTrackHistograms(fHist.fTrack[17],trk);
    }
//-----------------------------------------------------------------------------
// TRK_18: BEST_ID tracks with T0 > 700
//-----------------------------------------------------------------------------
    if ((trk->fIDWord == 0) && (trk->T0() > 700)) {
      FillTrackHistograms(fHist.fTrack[18],trk);
    }
//-----------------------------------------------------------------------------
// TRK_19: BEST_ID tracks with an associated cluster
//-----------------------------------------------------------------------------
    if ((trk->fIDWord == 0) && (trk->Ep() > 0)) {
      FillTrackHistograms(fHist.fTrack[19],trk);

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
// TRK_20: BEST_ID tracks intersecting DISK=0
// TRK_21: BEST_ID tracks intersecting DISK=1
//-----------------------------------------------------------------------------
    if (vt && (trk->fIDWord == 0)) {
      if      (vt->fID == 0) FillTrackHistograms(fHist.fTrack[20],trk);
      else if (vt->fID == 1) FillTrackHistograms(fHist.fTrack[21],trk);
    }
//-----------------------------------------------------------------------------
// TRK 22: BEST_ID tracks with an associated cluster and chi2(match) < 100
// TRK 23: BEST_ID tracks with an associated cluster and chi2(match) < 100 and LLHR(cal) > 0
//         this is interesting to see which muons are getting misidentified
//-----------------------------------------------------------------------------
    if ((trk->fIDWord == 0) && (tp->fEp > 0) && (tp->fChi2Tcm < 100.)) {

      FillTrackHistograms(fHist.fTrack[22],trk);
      if (trk->LogLHRCal() > 0) {
	FillTrackHistograms(fHist.fTrack[23],trk);
	
	if (GetDebugBit(36) && (trk->fP < 80.)) {
	  GetHeaderBlock()->Print(Form(" bit:036 trk p = %10.3f E/P = %10.3f", trk->fP,tp->fEp));
	}
      }
      else {
//-----------------------------------------------------------------------------
// TRK_24: BEST_ID  tracks with an associated cluster and chi2(match) < 100 and LLHR(cal) < 0
//         this set allows to see which electrons are getting misidentified
//-----------------------------------------------------------------------------
	FillTrackHistograms(fHist.fTrack[24],trk);
      }
//-----------------------------------------------------------------------------
// TRK_25: BEST_ID tracks, 100 < P < 110, 0 < E/p < 1.15,  |dt_corr| < 3, chi2(match) < 100
//-----------------------------------------------------------------------------
//      double dt_corr = trk->Dt(); // -1.;
      if ( (tp->fDt > fMinDtTcm) && (tp->fDt < fMaxDtTcm) && 
	   (tp->fEp <      1.15) && 
	   (trk->fP >      100.) && (trk->fP <       110)   ) {
	FillTrackHistograms(fHist.fTrack[25],trk);
//-----------------------------------------------------------------------------
// more details on the calorimeter-based likelihood 
// TRK_26 : TRK 25 events with LLHR_CAL > 0 ( interesting for muons)
// TRK_27 : TRK 25 events with LLHR_CAL < 0 ( interesting for electrons)
//-----------------------------------------------------------------------------
	double llhr_cal = trk->LogLHRCal();
	if (llhr_cal > 0) {
	  FillTrackHistograms(fHist.fTrack[26],trk);

	  if (llhr_cal > 5.) {
	    if (GetDebugBit(37)) {
	      GetHeaderBlock()->Print(Form("bit:037 trk p = %10.3f E/P = %10.3f", trk->fP,tp->fEp));
	    }
	  }
	}
	else {
	  FillTrackHistograms(fHist.fTrack[27],trk);
	}
      }

    }
//-----------------------------------------------------------------------------
// TRK_28 : events with a BEST_ID track and a cluster E/P > 1.1
//-----------------------------------------------------------------------------
    if (trk->fIDWord == 0) {
      if (tp->fEp > 1.1) {
	FillTrackHistograms(fHist.fTrack[28],trk);

	if (GetDebugBit(28)) {
	  GetHeaderBlock()->Print(Form(" bit:028 LLHR(CAL) = %10.3f ep = %10.3f dt = %10.3f",
				       trk->LogLHRCal(), trk->Ep(), trk->Dt()));
	}
      }
    }
//-----------------------------------------------------------------------------
// TRK_29: BEST_ID track 100 < P < 110, 0 < E/P < 1.15              : next to TRK_13
// TRK_32: BEST_ID track 100 < P < 110, 0 < E/P < 1.15, chi2tcm<100 : next to TRK_29
//-----------------------------------------------------------------------------
    if (trk->fIDWord == 0) {
      if ((tp->fEp > 0) && (trk->P() > 100.) && (trk->P() < 110.) && (tp->fEp < 1.15)) {
	FillTrackHistograms(fHist.fTrack[29],trk);

	if (tp->fChi2Tcm < 100.) {
	  FillTrackHistograms(fHist.fTrack[32],trk);
	}
      }
    }
//-----------------------------------------------------------------------------
// TRK_30: tracks with >= 25 hits
// TRK_31: tracks with >= 25 hits and Chi2/Ndof < 3
//-----------------------------------------------------------------------------
    if (trk->NActive() >= 25) {
      FillTrackHistograms(fHist.fTrack[30],trk);
      if (trk->Chi2Dof() < 3) {
	FillTrackHistograms(fHist.fTrack[31],trk);
      }
    }
//-----------------------------------------------------------------------------
// TRK_33: "DaveTrkQual" Q > 0.4 tracks
//-----------------------------------------------------------------------------
    if (tp->fIDWord[fBestID] == 0) {
      FillTrackHistograms(fHist.fTrack[33],trk);
    }
//-----------------------------------------------------------------------------
// TRK_34: "DaveTrkQual" Q > 0.1 tracks
//-----------------------------------------------------------------------------
    if (tp->fIDWord[1] == 0) {
      FillTrackHistograms(fHist.fTrack[34],trk);
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
      FillTrackHistograms(fHist.fTrack[40],trk);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[41],trk);
	// print run/event numbers :
	if (GetDebugBit(6)) {
	  double ep = trk->Ep();
	  if ((ep > 0.8) && (ep < 1.1)) {
	    GetHeaderBlock()->Print(Form(" bit:006 trk_41: track E/P = %8.3f",trk->Ep()));
	  }
	}
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[42],trk);
      }
    }
    else if (alg_mask == 2) {
//-----------------------------------------------------------------------------
// CalPatRec-only tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[50],trk);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[51],trk);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[52],trk);
      }
    }
    else if (alg_mask == 3) {
//-----------------------------------------------------------------------------
// TrkPatRec+CalPatRec tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[60],trk);
      if (trk->fIDWord == 0) {
	FillTrackHistograms(fHist.fTrack[61],trk);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[62],trk);
      }
    }
//-----------------------------------------------------------------------------
// TRK_63 : fIDWORD[bestID] == 0, TRkPatRec
// TRK_64 : fIDWORD[bestID] == 0, CalPatRec
//-----------------------------------------------------------------------------
    if (tp->fIDWord[fBestID] == 0) {
      if      (best_alg == 0) FillTrackHistograms(fHist.fTrack[63],trk);
      else if (best_alg == 1) FillTrackHistograms(fHist.fTrack[64],trk);
    }
//-----------------------------------------------------------------------------
// TRK_71: SetC tracks  103.5 < p < 105 : DIO studies
//-----------------------------------------------------------------------------
    if ((trk->fIDWord == 0) && (trk->fP > 103.5) && (trk->fP < 105.)) {
      FillTrackHistograms(fHist.fTrack[71],trk);
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
// assume less than 20 tracks 
//-----------------------------------------------------------------------------
int TTrackAnaModule::InitTrackPar(TStnTrackBlock*   TrackBlock  , 
				  TStnClusterBlock* ClusterBlock, 
				  TrackPar_t*       TrackPar    ) {
  TrackPar_t*           tp;
  TStnTrack*            track;
  int                   id_word, icorr (-1);
  double                xs;
  TEmuLogLH::PidData_t  dat;
//-----------------------------------------------------------------------------
// momentum corrections for TrkPatRec and CalPatRec
//-----------------------------------------------------------------------------
  const double kMomentumCorr[2] = { 0.049, 0.020 };
  const double kDtTcmCorr   [2] = { 0.22 , -0.30 }; // ns, sign: fit peak positions

  const char* block_name = TrackBlock->GetNode()->GetName();

  if (fApplyCorrections != 0) {
    if      (strcmp(block_name,"TrkPatRec"    ) == 0) icorr = 0;
    else if (strcmp(block_name,"CalPatRec"    ) == 0) icorr = 1;
    else if (strcmp(block_name,"TrackBlock"   ) == 0) icorr = 2;
    else if (strcmp(block_name,"TrackBlockDmm") == 0) icorr = -1;
    else if (strcmp(block_name,"TrackBlockUmp") == 0) icorr = -1;
    else {
      icorr = -999;
      Error("TTrackCompModule::InitTrackPar","IN TROUBLE");
      return -1;
    }
  }
//-----------------------------------------------------------------------------
// loop over tracks, assume 
//-----------------------------------------------------------------------------
  int ntrk = TrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    tp             = TrackPar+itrk;
    track          = TrackBlock->Track(itrk);

    for (int idd=0; idd<fNID; idd++) {
      tp->fIDWord[idd] = fTrackID[idd]->IDWord(track);
    }

    id_word        = tp->fIDWord[fBestID];
    track->fIDWord = id_word;

    if (id_word == 0) {
      fNGoodTracks += 1;
      if ((track->fVMaxEp != NULL) && (fabs(track->fVMaxEp->fDt) < 2.5)) {
	fNMatchedTracks += 1;
      }
    }
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
    if (icorr == 2) icorr = track->BestAlg();

    tp->fP = track->fP0;
    if (icorr >= 0) tp->fP  += kMomentumCorr[icorr];		// correcting

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

    tp->fDtZ0    = -1.e6;
    tp->fDtBack  = -1.e6;
    if (fSimPar.fTMid ) tp->fDtZ0   = track->T0()-fSimPar.fTMid->Time();
    if (fSimPar.fTBack) tp->fDtBack = track->TBack()-fSimPar.fTBack->Time();
//-----------------------------------------------------------------------------
// track residuals 
// 2016-10-26: replase VMaxEp with VMinS
//-----------------------------------------------------------------------------
    TStnTrack::InterData_t*  vr = track->fVMinS; 
    double    nx, ny;

    tp->fDiskID    = -1;
    tp->fEcl       = -1.e6;
    tp->fEp        = -1.e6;
    tp->fEDiskP    = -1.e6;

    double ecal    = fDiskCalorimeter->Energy();
    tp->fECalP     = ecal/track->fP2;

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
      tp->fEcl     = vr->fEnergy;
      tp->fEp      = tp->fEcl/track->fP2;
      tp->fDx      = vr->fDx;
      tp->fDy      = vr->fDy;
      tp->fDz      = vr->fDz;

      double edisk(-1.);
      if (fDiskCalorimeter->NDisks() > 0) edisk = fDiskCalorimeter->Disk(vr->fID)->Energy();
      tp->fEDiskP  = edisk/track->fP2;
//-----------------------------------------------------------------------------
// correct TrkpatRec and CalPatRec tracks differently
//-----------------------------------------------------------------------------
      tp->fDt  = vr->fDt;
      if (icorr >= 0) tp->fDt  -= kDtTcmCorr[icorr];

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
// it track hits are stored, calculate track mean time averaged over the hits 
//-----------------------------------------------------------------------------
    tp->fTMean = -1.e6;
    if (fTrackStrawHitBlock->NTracks() > 0) {
    
      int nhits = fTrackStrawHitBlock->NTrackHits(itrk);

					// calculate the mean time
      float tmean(0);
      for (int i=0; i<nhits; i++) {
	TTrackStrawHitData* hit = fTrackStrawHitBlock->Hit(itrk,i);
	tmean += hit->Time();
      }

      tp->fTMean  = tmean/nhits;
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
int TTrackAnaModule::Event(int ientry) {

  TDiskCalorimeter::GeomData_t disk_geom;

  fTrackBlock->GetEntry(ientry);
  fTrackStrawHitBlock->GetEntry(ientry);
  fClusterBlock->GetEntry(ientry);
  fStrawDataBlock->GetEntry(ientry);
  fCalDataBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fVdetBlock->GetEntry(ientry);
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
// use the first fit
//-----------------------------------------------------------------------------
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  int nvdhits = fVdetBlock->NHits();
  for (int i=0; i<nvdhits; i++) {
    TVdetHitData* vdhit = fVdetBlock->Hit(i);
    if (vdhit->PdgCode() == fSimp->fPdgCode) {
      if ((vdhit->Index() == 13) || (vdhit->Index() == 14)) {
	if (fDirection*vdhit->McMomentumZ() > 0) {
	  if (fSimPar.fTFront == 0) fSimPar.fTFront = vdhit;
	}
      }
      else if ((vdhit->Index() == 11) || (vdhit->Index() == 12)) {
	if (fDirection*vdhit->McMomentumZ() > 0) {
	  if (fSimPar.fTMid == 0) fSimPar.fTMid = vdhit;
	}
      }
      else if (vdhit->Index() == mu2e::VirtualDetectorId::TT_Back) {
	if (fDirection*vdhit->McMomentumZ() > 0) {
	  if (fSimPar.fTBack == NULL) fSimPar.fTBack = vdhit;
	}
      }
    }
  }

  if (fDiskCalorimeter->Initialized() == 0) {
    disk_geom.fNDisks = fCalDataBlock->NDisks();

    for (int i=0; i<disk_geom.fNDisks; i++) {
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

  fNGoodTracks    = 0;
  fNMatchedTracks = 0;
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

  InitTrackPar(fTrackBlock,fClusterBlock,fTrackPar);
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
void TTrackAnaModule::Debug() {

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
// bit 1: IDWord =0 0 tracks
//-----------------------------------------------------------------------------
    if (GetDebugBit(1) == 1) {
      if (trk->fIDWord == 0) {
	GetHeaderBlock()->Print(Form("bit_001: IDWord=0 p = %10.3f",
				     trk->Momentum()->P()));
      }
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
int TTrackAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TTrackAnaModule::Test001() {
}

