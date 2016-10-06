///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TROOT.h"
#include <cmath>
#include <cassert>

#include "murat/plot/smooth.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

#include "obj/TMu2eDatasets.hh"
#include "obj/TMu2eChannel.hh"

ClassImp(TMu2eChannel)

//-----------------------------------------------------------------------------
// static variables
//-----------------------------------------------------------------------------
TH1*    TMu2eChannel::g_h1(0);
TH1*    TMu2eChannel::g_h2(0);
TH1*    TMu2eChannel::g_h3(0);

				// cross section in picobarns
const double xs_wln = 2690.;
const double xs_zll =  255.;
//-----------------------------------------------------------------------------
// processes used in this analysis: first goes data, the last one - MC for the 
// signal process
//-----------------------------------------------------------------------------
//   {"tmet"   , "data"}, 
//   {"wenu"   , "mc"  },
//   {"wmunu"  , "mc"  },
//   {"zee"    , "mc"  },
//   {"zmumu"  , "mc"  },
//   {"ztautau", "mc"  },
//   {"qcd"    , "data"},
//   {"zzx" , "mc"  }                    // W->tau nu MC has to be the last one
// };

//-----------------------------------------------------------------------------
// versions I have: ""             : default
//-----------------------------------------------------------------------------
TMu2eChannel::TMu2eChannel(const char* Version, int LumiBin, const char* Signal): 
  analysis ("zzx",Version,Signal,LumiBin) {

  //  double         q_tot, q_loose, q_tight, xsec;
  int            mcflag, /*ntotal_mc, srr,*/ ipp;
  //  char           norm_hist_name[200];
  //  char           fn[200];
  const char    *pname /*, *ptype*/;
  //  TH1*           norm_hist;
  aprocess*      p;

  analysis::fDataset = "zzx";
//-----------------------------------------------------------------------------
// proc[0] : data, rest - backgrounds
// process name already implies whether it is MC or data...
// MC signal is the last process
//-----------------------------------------------------------------------------
  ipp = 0;
  fProc[ipp].name = "data"     ; fProc[ipp].mcflag = 0; ipp++;  // data
  fProc[ipp].name = "dio"      ; fProc[ipp].mcflag = 1; ipp++;  // DIO
  fProc[ipp].name = "cosmics"  ; fProc[ipp].mcflag = 1; ipp++;  // cosmics
  fProc[ipp].name = "conv"     ; fProc[ipp].mcflag = 3; ipp++;  // mu --> e conversion
  fProc[ipp].name = ""         ;                                // END
//-----------------------------------------------------------------------------
// define datasets
//-----------------------------------------------------------------------------
  analysis::fDsMetadata = new TMu2eDatasets(Version,Signal,LumiBin);
  analysis::fDsMetadata->Print("");
  analysis::init();
//-----------------------------------------------------------------------------
// initialize "data" process, the data come from proc[0]
//-----------------------------------------------------------------------------
  run_range_dat_t* rr;
  for (int irr=0; irr<fNRunRanges; irr++) { 
    printf(" ---- TMu2eChannel: RR = %2i\n",irr);
    if (fDsMetadata->UsedRR(irr) == 0)      goto NEXT_RR;

    rr           = fRR+irr;
    rr->fIntLumi = fDsMetadata->Luminosity(irr,LumiBin);
    rr->fDat     = new aprocess(fProc[0].name.Data(),0,irr,this);
//-----------------------------------------------------------------------------
// define backgrounds, the signal process (W->tau nu) is the last one
// use bgr[i].proc_name == 0 to terminate the loop
//------------------------------------------------------------------------------
    rr->fBgr  = new TClonesArray("aprocess");
    rr->fNBgr = 0;

    for (int i=1; fProc[i].name != ""; i++) {
      pname  = fProc[i].name.Data();
      mcflag = fProc[i].mcflag;
      printf("i,pname,mcflag: %i %s %i\n",i,pname,mcflag);
//-----------------------------------------------------------------------------
// create new background process
// MC signal is the last among the background processes
// MC also has Version=ze0s01 etc...
//-----------------------------------------------------------------------------
      p  = new ((*rr->fBgr)[i-1]) aprocess(pname,mcflag,irr,this);
      rr->fNBgr++;
      if (mcflag == 3) {
	rr->fSig = p;
      }
    }
  NEXT_RR:;
  }
}

//-----------------------------------------------------------------------------
TMu2eChannel::~TMu2eChannel() {
}


//-----------------------------------------------------------------------------
// produce 'path-qualified' name of a histogram in the input file
//-----------------------------------------------------------------------------
void TMu2eChannel::GetHistogramName(const char* HistSet     ,
				     int         RunRange    ,
				     int         Bin         ,
				     const char* FitHistName ,
				     char*       FullHistName) {

  sprintf(FullHistName,"%s_%i/%s",HistSet,Bin,FitHistName);
}

//-----------------------------------------------------------------------------
// produce 'path-qualified' name of a histogram in the input file
// histograms for "zee_fakes" and "zmm_fakes" are produced by TZAnaModule, 
// its standard name is "ZAna"
// all the rest histograms so far come from the module names "ZZAna"
// ... shame! this wasn't necessary...
//-----------------------------------------------------------------------------
int TMu2eChannel::GetHistogram(aprocess*       Process,
			       const char*     Module  ,
			       const char*     HistSet ,
			       int             RunRange,
			       int             Bin     ,
			       const char*     HistName,
			       int             Rebin   ,
			       TH1*&           Hist    ) {

  int rc(0);

  if (strcmp(Module,"Mu2eLimits") == 0) {
    rc = Process->get_h1(Module,HistSet,Bin,HistName,Rebin,Hist);
  }
  else {
    Error("GetHistogram",">>> not Mu2eLimits, return NULL");
    Hist = 0;
    rc   = -1;
  }
  
  return rc;
}

//-----------------------------------------------------------------------------
// standard naming convention for histograms used in fits...
//-----------------------------------------------------------------------------
void TMu2eChannel::get_hist_name(const char* Module       ,
				 const char* HistSet       ,
				 int         RunRange      ,
				 int         Bin           ,
				 const char* FitHistName   ,
				 int         Rebin         ,
				 char*       Name          ) {
  
  sprintf(Name,"%s/%s_%i/%s_rebin_%i",Module,HistSet,Bin,
	  FitHistName,
	  Rebin);
}


//-----------------------------------------------------------------------------
// we only need number of entries...
//-----------------------------------------------------------------------------
double TMu2eChannel::GetAcceptance(const aprocess*    Process, 
				    const char*        Module ,
				    const char*        HistSet,
				    int                RunRange,
				    int                Bin    ,
				    const char*        HistName) const {
  double  acc;
  TH1*    h1;
  char    name[100];
  int     bin;

  const char* fn = Process->GetHistFileName();

  bin = Bin;

  sprintf(name,"%s_%i/%s",HistSet,bin,HistName);
  h1 = ::gh1(fn,Module,name);

  acc = h1->GetEntries()/Process->GetNTotal()*Process->GetAccCorrFactor();

  return acc;
}

//-----------------------------------------------------------------------------
// set ID eff SF to 1
// moving towards having MC histograms filled with weights, such that scaling
// here would not be necessary...
//-----------------------------------------------------------------------------
double TMu2eChannel::GetIDEffSF(const aprocess*    Process,
				 const char*        Module ,
				 const char*        HistSet,
				 int                RunRange,
				 int                Bin    ,
				 const char*        HistName) const {
//-----------------------------------------------------------------------------
// use TCE-tagged Z's
//-----------------------------------------------------------------------------
//  double const tce_id_eff_sf [] = {
////     0.95668, 0.93807, 1.00080, 1.00252, 0.94399,
////     1.00190, 0.97465, 1.00928, 0.93759, 0.97651,
////     0.96009, 0.96709, 0.93094, 0.98974, 0.95648,
////     0.95917, 0.97634, 0.95651, 0.96108, 0.94285,
////     0.91004, 0.91891, 0.92991, 0.94990, 0.90399,
////     0.96241, 0.93492, -1
//
//    0.95003, 0.95630, 0.98153, 0.98441, 0.94983,
//    0.95589, 0.96353, 0.93320, 0.95009, 0.96872,
//    0.95951, 0.96990, 0.95133, 0.98277, 0.95768,
//    0.94282, 0.96033, 0.98182, 0.96765, 0.94764,
//    0.94279, 0.91435, 0.92258, 0.94499, 0.90771,
//    0.94421, 0.93988, -1
//  };
//
//  double const lce_id_eff_sf [] = {
////     1.02521, 0.96616, 0.99400, 1.02574, 0.96971,
////     1.03284, 1.00299, 1.02556, 0.96098, 1.00608,
////     0.99099, 0.97627, 0.95463, 1.00663, 0.94585,
////     1.00495, 0.98973, 0.97543, 0.98655, 0.97471,
////     0.96795, 0.95952, 0.96810, 0.98723, 0.95894,
////     0.99003, 0.96268, -1.
//
//    1.01591, 0.99298, 1.01401, 1.00342, 0.96718,
//    0.97666, 1.00017, 0.95334, 0.96519, 0.99902,
//    0.99116, 0.98526, 0.97738, 1.00927, 0.96164,
//    0.97271, 0.98037, 0.99113, 0.97901, 0.97595,
//    0.97725, 0.95297, 0.96772, 0.96827, 0.94892,
//    0.97400, 0.97805, -1
//  };
//
//  double const tpe_id_eff_sf [] = {
////     0.94111, 0.94235, 0.97723, 0.96251, 0.98077,
////     0.93700, 0.96179, 0.90483, 0.93254, 0.92605,
////     0.92560, 0.90626, 0.93250, 0.92177, 0.91125,
////     0.94640, 0.94197, 0.93215, 0.91651, 0.91872,
////     0.90157, 0.90522, 0.91689, 0.90166, 0.91225,
////     0.91662, 0.92496, -1
//
//    0.94498, 0.96032, 0.93960, 0.94549, 0.94610,
//    0.93112, 0.94418, 0.91310, 0.94314, 0.94356,
//    0.93918, 0.91467, 0.91241, 0.93347, 0.94947,
//    0.95464, 0.99323, 0.94817, 0.90907, 0.92039,
//    0.90081, 0.91988, 0.92822, 0.95227, 0.94105,
//    0.91409, 0.89860, -1
//  };

  double id_eff_sf = 1.;

  return id_eff_sf;
}


//-----------------------------------------------------------------------------
// plug tracking efficiency scale factor
//-----------------------------------------------------------------------------
double TMu2eChannel::GetTrEffSF(const aprocess* Process, 
				const char*     Module ,
				const char*     HistSet,
				int             RunRange,
				int             Bin     ,
				const char*     HistName) const {
  double trk_eff_sf;

  double const k_plug_trk_eff_sf[] = {
    0.97761, 0.98691, 0.98641, 0.99251, 0.97127,
    0.98635, 0.98073, 0.95097, 0.98956, 0.97041,
    0.97164, 0.97166, 0.98176, 0.96053, 0.97587,
    0.98664, 0.99632, 0.98849, 0.98850, 0.98636,
    0.98520, 0.99080, 0.99152, 0.98662, 0.99452,
    0.97999, 0.98060, -1
  };

  TString process_name  = Process->GetName();
  TString hist_set      = HistSet;

  trk_eff_sf = 1;

  if (process_name == "bhelzl") {
    if (hist_set == "zee_rr") {
      if (Bin == 3) {		        // TCE-TPE
	trk_eff_sf = k_plug_trk_eff_sf[RunRange];
      }
    }
    else if (hist_set == "zee") {
      if (Bin == 43) {		        // TCE-TPE
	trk_eff_sf = k_plug_trk_eff_sf[RunRange];
      }
    }
  }

  return trk_eff_sf;
}

//-----------------------------------------------------------------------------
// this routine returns total normalization for varios histograms 
// which are not normalized from the first principles...
//-----------------------------------------------------------------------------
double TMu2eChannel::GetFitQEvents(const aprocess*    Process, 
				   const char*        Module ,
				   const char*        HistSet,
				   int                RunRange,
				   int                Bin    ,
				   const char*        HistName) const {
  
  double fit_qevents(0);

  return fit_qevents;
}

