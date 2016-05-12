//////////////////////////////////////////////////////////////////////////////
// base class for track analysis modules - want unpacking to be done the same
// way
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
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------

#include "ana/TTrackAnaModuleBase.hh"

ClassImp(TTrackAnaModuleBase)
//-----------------------------------------------------------------------------
TTrackAnaModuleBase::TTrackAnaModuleBase(const char* name, const char* title):
  TStnModule(name,title)
{
//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
//-----------------------------------------------------------------------------
  fPdgCode       = 11;
  fGeneratorCode = 2;			// conversionGun, 28:StoppedParticleReactionGun
  fLogLH         = NULL;
}

//-----------------------------------------------------------------------------
TTrackAnaModuleBase::~TTrackAnaModuleBase() {
  for (int i=0; i<fNID; i++) delete fTrackID[i];

  if (fLogLH) delete fLogLH;
}



//_____________________________________________________________________________
int TTrackAnaModuleBase::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//-----------------------------------------------------------------------------
// fill efficiency histograms : need 10 histogram sets
// pitch = 1./tan(dip)
//-----------------------------------------------------------------------------
void TTrackAnaModuleBase::FillEfficiencyHistograms(TStnTrackBlock*  TrackBlock, 
						   TStnTrackID*     TrackID   , 
						   int              HistSet   ) {
  if (fSimp->NStrawHits() >= 20) {
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
		
		if (((id_word & TStnTrackID::kD1Bit) == 0) && 
		    ((id_word & TStnTrackID::kD1Bit) == 0)    ) {
		  
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

//-----------------------------------------------------------------------------
// assume less than 20 tracks 
//-----------------------------------------------------------------------------
int TTrackAnaModule::InitTrackPar(TStnTrackBlock*   TrackBlock  , 
				  TStnClusterBlock* ClusterBlock, 
				  TrackPar_t*       TrackPar    ) {
  TrackPar_t*           tp;
  TStnTrack*            track;
  int                   id_word, icorr;
  double                xs;
  TEmuLogLH::PidData_t  dat;
//-----------------------------------------------------------------------------
// momentum corrections for TrkPatRec and CalPatRec
//-----------------------------------------------------------------------------
  const double kMomentumCorr[2] = { 0.049, 0.020 };
  const double kDtTcmCorr   [2] = { 0.22 , -0.30 }; // ns, sign: fit peak positions

  const char* block_name = TrackBlock->GetNode()->GetName();

  if      (strcmp(block_name,"TrkPatRec"    ) == 0) icorr = 0;
  else if (strcmp(block_name,"CalPatRec"    ) == 0) icorr = 1;
  else if (strcmp(block_name,"TrackBlock"   ) == 0) icorr = 2;
  else if (strcmp(block_name,"TrackBlockDmm") == 0) icorr = -1;
  else {
    icorr = -999;
    Error("TTrackCompModule::InitTrackPar","IN TROUBLE");
    return -1;
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

    tp->fP = track->fP2;
    if (icorr >= 0) tp->fP  += kMomentumCorr[icorr];		// correcting

    tp->fDpF   = track->fP2    -track->fPFront;
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

