///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TMath.h"

#include "obj/analysis.hh"
#include "obj/aprocess.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

ClassImp(aprocess)

//-----------------------------------------------------------------------------
aprocess::aprocess(const char* Name):TNamed(Name,Name) {
  // a dataset may be a stripped one

  fHistFileName     = "";
  fRunRange         = -1;
  fDsMetadata       = NULL;
  fSrr              = -1;
  fMcFlag           = -1;
  fXSec             = 0;
  fColor            =  1;
  fStyle            = 1;
//-----------------------------------------------------------------------------
// number of events corresponding to the cross section
//-----------------------------------------------------------------------------
  fNTotal           = 0;
  fIntLumi          = 0;
  fAcc              = 0;
  //  fIdEff            = 1;
  fTrEff            = 1;
  fTightIdSF        = 1;
  fLooseIdSF        = 1;
  fTrEffSF          = 1;
  fNormHist         = 0;
  fQEvents          = 0;
  fAnalysis         = 0;
  //  fListOfHistograms = 0; 
  for (int i=0; i<kNBins; i++) { 
    fTight[i]      = 0; 
    fFitQEvents[i] = -1;
  }
}

//-----------------------------------------------------------------------------
aprocess::aprocess(const char* Name    , 
		   int         McFlag  , 
		   int         RunRange,
		   analysis*   Analysis) : TNamed(Name,Name) {
  //
  char           fn[200];
  TDsMetadata*   metadata;
  const char*    dsid;
  int            lumi_bin;
//-----------------------------------------------------------------------------
// fTitle = datasetID, same "zee" MC can use different MC datasets...
//-----------------------------------------------------------------------------
  fMcFlag           = McFlag;
  fRunRange         = RunRange;
  fAnalysis         = Analysis;
  
  if (fAnalysis) {
    lumi_bin          = fAnalysis->GetLumiBin();
    fDsMetadata       = Analysis->GetDsMetadata();

    metadata          = fDsMetadata->GetDsMetadata(Name);
    dsid              = metadata->GetTitle();
    SetTitle(dsid);

    fDsMetadata->GetHistFileName(Name,fMcFlag,RunRange,fn);
    fHistFileName     = fn;

    fSrr              = fDsMetadata->SubstituteRunRange(Name,RunRange);
    fIntLumi          = fDsMetadata->Luminosity(RunRange,lumi_bin);
    fXSec             = fDsMetadata->GetXSec(Name,McFlag);
//-----------------------------------------------------------------------------
// number of events corresponding to the cross section - important for MC
//-----------------------------------------------------------------------------
    fNTotal           = fDsMetadata->GetNTotalMc(Name,McFlag,fSrr);

    //    fListOfHistograms = new TObjArray(); 
    fColor            = fDsMetadata->GetHistColor(Name,McFlag);
    fStyle            = fDsMetadata->GetHistStyle(Name,McFlag);
  }

  for (int i=0; i<kNBins; i++) { 
    fTight[i]       = 0; 
    fFitQEvents[i]  = -1;
  }
//-----------------------------------------------------------------------------
// normalize MC histograms to the acc*lumi*cross_section from the 
// first principles, do not do it for the data
// however may want to scale some fake backgrounds... so fNormalize is not quite right...
//-----------------------------------------------------------------------------
  if (fMcFlag == 0) fNormalize = 0;
  else              fNormalize = 1;

  fAcc              = 0;
  fTrEff            = 1;
  fTightIdSF        = 1;
  fLooseIdSF        = 1;
  fTrEffSF          = 1;
  fNormHist         = 0;
  fQEvents          = 0;
}

//-----------------------------------------------------------------------------
aprocess::~aprocess() {
  //  delete fListOfHistograms;
}

//-----------------------------------------------------------------------------
// MC:   returned histogram is normalized according to the integrated luminosity
//       and cross-section
// data: returned histogram is normalized to a number of events returned by 
//       fAnalysis::get_nbgr-events
// returned histogram is cloned and is supposed to be owned by the caller
// may need to retrieve a histogram with a name, different from requested one
//-----------------------------------------------------------------------------
int aprocess::get_h1(const char* Module  , 
		     const char* HistSet , 
		     int         Bin     , 
		     const char* HistName, 
		     int         Rebin   ,
		     TH1*&       Hist    )  
{
  // note, that as overlaying 2D histograms doesn't make much sense, 
  // possibly rebinned
  // returned histogram is supposed to be owned by the caller
  
  TH1  *hist, *h1, *h2, *hnorm;

  char    name[500], nnorm[500], file_hist_name[200], job_hist_name[200];

  TString module_name;

  double  q, q1, e1, qent, sf(1.), qnorm, fit_qevents;
  double  acc, id_eff, id_eff_sf, tr_eff, tr_eff_sf;

  int     nbins, delete_h1(0);

  module_name = Module;
//-----------------------------------------------------------------------------
// process 'knows' about its run range...
// 'Module' is passed just for if's  .... file hist name does not include the Module name
//-----------------------------------------------------------------------------
  GetHistogramName(Module,HistSet,Bin,HistName,"file",file_hist_name);
//-----------------------------------------------------------------------------
// for data backgrounds check if a histogram with a given name has been assigned
// naming conventions are important to follow !
// example: name = 'TauAna/tau_1/ang_tr_rebin_5'
// however, name is the run-time name .... job hist name does include the module name
//-----------------------------------------------------------------------------
  GetHistogramName(Module,HistSet,Bin,HistName,"job",job_hist_name);

  sprintf(name,"%s_%s_rebin_%i",Module,job_hist_name,Rebin);

  h1 = 0; // (TH1*) fListOfHistograms->FindObject(name);

  if (h1) {
    Hist = (TH1*) h1->Clone();
  }
  else {
//-----------------------------------------------------------------------------
// fetch histogram from a file
//-----------------------------------------------------------------------------
    hist = ::gh1(fHistFileName,module_name,file_hist_name);
    if (hist == 0) {
      printf("aprocess::get_h1 ERROR 001: cant find %s/%s in %s, set Hist = 0 and EXIT\n",
	     module_name.Data(),file_hist_name,fHistFileName.Data());
      Hist = 0;
      return -1;
    }
    h1   = (TH1*) hist->Clone();
    delete_h1 = 1;
//-----------------------------------------------------------------------------
// define h1, ignore underflows and overflows, defined errors and normalize to
// right luminosity and efficiency
//-----------------------------------------------------------------------------
    if (fNormalize) {
//-----------------------------------------------------------------------------
// MC - (fXSec > 0) : do scaling
//-----------------------------------------------------------------------------
      fit_qevents = fAnalysis->GetFitQEvents(this,Module,HistSet,fRunRange,Bin,HistName);

      if (fit_qevents > 0) {
//-----------------------------------------------------------------------------
// if number of events has been determined by fit, use fit result
//-----------------------------------------------------------------------------
	fQEvents = fit_qevents;
	//qnorm    = hist->GetEntries();
	qnorm    = hist->Integral();// Aidan100714 histograms are now weighted

	if (fDebug[0] != 0) {
	  printf("[aprocess:get_h1] fQEvents=%10.3f qnorm=%10.3f\n",fQEvents,qnorm);
	}
      }
      else {
//-----------------------------------------------------------------------------
// otherwise determine normalization from the "first principles"
//-----------------------------------------------------------------------------
	acc       = fAnalysis->GetAcceptance(this,Module,HistSet,fRunRange,Bin,HistName);
	id_eff    = fAnalysis->GetIDEff     (this,Module,HistSet,fRunRange,Bin,HistName);
	id_eff_sf = fAnalysis->GetIDEffSF   (this,Module,HistSet,fRunRange,Bin,HistName);
	tr_eff    = fAnalysis->GetTrEff     (this,Module,HistSet,fRunRange,Bin,HistName);
	tr_eff_sf = fAnalysis->GetTrEffSF   (this,Module,HistSet,fRunRange,Bin,HistName);

	if (fDebug[0] != 0) {
	  printf("[aprocess:get_h1] acc, id_eff, id_eff_sf, tr_eff, tr_eff_sf = %10.3f %10.3f %10.3f %10.3f 510.3f\n", acc, id_eff, id_eff_sf, tr_eff, tr_eff_sf);
	}
	// fQEvents - number of events in the normalized histogram

	fQEvents  = fXSec*fIntLumi*acc*id_eff*tr_eff*id_eff_sf*tr_eff_sf; //Aidan: fine, effs=1 so use GetEntries
	if (fNormHist) {
	  qnorm     = fNormHist->GetEntries();
	}
	else {
	  qnorm     = hist->GetEntries();
	}

	if (fDebug[0] != 0) {
	  printf("[aprocess::get_h1] fQEvents, qnorm = %10.3f %10.3f\n",fQEvents,qnorm);
	}
      }

      if (qnorm > 0) sf = fQEvents/qnorm;
      else           { 
	printf("[aprocess::get_h1] WARNING: %s normalization histogram empty, set SF to 0!\n", 
	       GetName());
	sf = 0.;
      }
    }
    else {
//-----------------------------------------------------------------------------
// data histogram  - use AccScaleFactor
//-----------------------------------------------------------------------------
      sf = fDsMetadata->GetAccCorrFactor(GetName(),fMcFlag,fRunRange);
    }

    h1->Reset();
    nbins    = hist->GetNbinsX();

    if (fDebug[0] != 0) {
      printf("[aprocess::get_h1] sf = %10.3f\n",sf);
    }

    for (int i=1; i<=nbins; i++) {
      q  = hist->GetBinContent(i);
      q1 = q*sf;
      e1 = TMath::Sqrt(q)*sf;
      h1->SetBinContent(i,q1);
      h1->SetBinError  (i,e1);
    }
//-----------------------------------------------------------------------------
// return a clone of a histogram with rebinning factor appended to the name
// either of TH1::Clone and TH1::Rebin creates a new histogram
//-----------------------------------------------------------------------------
    if (Rebin == 1) Hist = (TH1*) h1->Clone(job_hist_name);
    else            Hist = h1->Rebin(Rebin,job_hist_name);
//-----------------------------------------------------------------------------
// finally make sure that we do not create multiple clones
//-----------------------------------------------------------------------------
    delete h1;
  }

  return 0;
}

//-----------------------------------------------------------------------------
// MC:   returned histogram is not normalized 
//       ANorm is the nrmalization factor
//       use to interface with MCLIMIT
// returned histogram is cloned and is supposed to be owned by the caller
// may need to retrieve a histogram with a name, different from requested one
//-----------------------------------------------------------------------------
int aprocess::get_h1_unorm(const char* Module  , 
			   const char* HistSet , 
			   int         Bin     , 
			   const char* HistName, 
			   int         Rebin   ,
			   TH1*&       Hist    ,
			   double*     Anorm)  
{
  // note, that as overlaying 2D histograms doesn't make much sense, 
  // possibly rebinned
  // returned histogram is supposed to be owned by the caller
  
  TH1  *hist, *h1, *h2, *hnorm;

  char    name[500], nnorm[500], file_hist_name[200], job_hist_name[200];

  TString module_name;

  double  q, q1, e1, qent, sf, qnorm, fit_qevents;
  double  acc, id_eff, id_eff_sf, tr_eff, tr_eff_sf;

  int     nbins, delete_h1(0);

  module_name = Module;
//-----------------------------------------------------------------------------
// process 'knows' about its run range...
// 'Module' is passed just for if's  .... file hist name does not include the Module name
//-----------------------------------------------------------------------------
  GetHistogramName(Module,HistSet,Bin,HistName,"file",file_hist_name);
//-----------------------------------------------------------------------------
// for data backgrounds check if a histogram with a given name has been assigned
// naming conventions are important to follow !
// example: name = 'TauAna/tau_1/ang_tr_rebin_5'
// however, name is the run-time name .... job hist name does include the module name
//-----------------------------------------------------------------------------
  GetHistogramName(Module,HistSet,Bin,HistName,"job",job_hist_name);
  sprintf(name,"%s_%s_rebin_%i",Module,job_hist_name,Rebin);
//-----------------------------------------------------------------------------
// fetch histogram from a file
//-----------------------------------------------------------------------------
  hist = ::gh1(fHistFileName,module_name,file_hist_name);
  if (hist == 0) {
    printf("aprocess::get_h1 ERROR 001: cant find %s/%s in %s, set Hist = 0 and EXIT\n",
	   module_name.Data(),file_hist_name,fHistFileName.Data());
    Hist = 0;
    return -1;
  }
  h1   = (TH1*) hist->Clone();
//-----------------------------------------------------------------------------
// define h1, ignore underflows and overflows, defined errors and normalize to
// right luminosity and efficiency
//-----------------------------------------------------------------------------
  fit_qevents = fAnalysis->GetFitQEvents(this,Module,HistSet,fRunRange,Bin,HistName);

  if (fit_qevents > 0) {
//-----------------------------------------------------------------------------
// if number of events has been determined by fit, use fit result
//-----------------------------------------------------------------------------
    fQEvents = fit_qevents;
    qnorm    = hist->Integral();// Aidan100714 histograms are now weighted
  }
  else {
//-----------------------------------------------------------------------------
// otherwise determine normalization from the "first principles"
//-----------------------------------------------------------------------------
    acc       = fAnalysis->GetAcceptance(this,Module,HistSet,fRunRange,Bin,HistName);
    id_eff    = fAnalysis->GetIDEff     (this,Module,HistSet,fRunRange,Bin,HistName);
    id_eff_sf = fAnalysis->GetIDEffSF   (this,Module,HistSet,fRunRange,Bin,HistName);
    tr_eff    = fAnalysis->GetTrEff     (this,Module,HistSet,fRunRange,Bin,HistName);
    tr_eff_sf = fAnalysis->GetTrEffSF   (this,Module,HistSet,fRunRange,Bin,HistName);

      // fQEvents - number of events in the normalized histogram

    fQEvents  = fXSec*fIntLumi*acc*id_eff*tr_eff*id_eff_sf*tr_eff_sf; //Aidan: fine, effs=1 so use GetEntries

    if (fNormHist) {
//-----------------------------------------------------------------------------
// fNormHist is defined for histograms corresponding to different jet 
// multiplicity bins etc
//-----------------------------------------------------------------------------
      qnorm     = fNormHist->GetEntries();
    }
    else {
      qnorm     = hist->GetEntries();
    }
  }
    
  if (qnorm > 0) *Anorm = fQEvents/qnorm;
  else           { 
    printf(" aprocess:get_h1 WARNING: %s normalization histogram empty, set SF to 0!\n", 
	   GetName());
    *Anorm = 0.;
  }
//-----------------------------------------------------------------------------
// return a clone of a histogram with rebinning factor appended to the name
// either of TH1::Clone and TH1::Rebin creates a new histogram
//-----------------------------------------------------------------------------
  if (Rebin == 1) Hist = (TH1*) h1->Clone(job_hist_name);
  else            Hist = h1->Rebin(Rebin,job_hist_name);
//-----------------------------------------------------------------------------
// finally delete h1 to ensure that we do not create multiple clones
//-----------------------------------------------------------------------------
  delete h1;

  return 0;
}

//-----------------------------------------------------------------------------
// produce 'path-qualified' name of a histogram in an input file
// default version, can be overloaded
// option = "job" or "file"
// for run range-dependent histograms use 'HistSet'="zee_rr_01" etc
//-----------------------------------------------------------------------------
void aprocess::GetHistogramName(const char* Module     ,
				const char* HistSet    ,
				int         Bin        ,
				const char* FitHistName,
				const char* Option,
				char*       FileHistName) const {
  if (strcmp(Option,"job") == 0) {
    sprintf(FileHistName,"%s_%i/%s",HistSet,Bin,FitHistName);
  }
  else if (strcmp(Option,"file") == 0) {
    sprintf(FileHistName,"%s_%i/%s",HistSet,Bin,FitHistName);
  }
}


//-----------------------------------------------------------------------------
// default: for backward compatibility
//-----------------------------------------------------------------------------
double aprocess::GetFitQEvents(const char*  Module ,
			       const char*  HistSet,
			       int          Bin    ,
			       const char*  HistName) const {

  double qev = fFitQEvents[Bin];
  return qev;
}


//-----------------------------------------------------------------------------
// default: for backward compatibility
//-----------------------------------------------------------------------------
double aprocess::GetIDEffSF(const char*  Module ,
			    const char*  HistSet,
			    int          Bin    ,
			    const char*  HistName) const {
  double id_eff_sf;

  if   (fTight[Bin] == 1) id_eff_sf = fTightIdSF;
  else                    id_eff_sf = fLooseIdSF;

  return id_eff_sf;
}


//-----------------------------------------------------------------------------
// default: for backward compatibility
//-----------------------------------------------------------------------------
double aprocess::GetTrEff  (const char*  Module ,
			    const char*  HistSet,
			    int          Bin    ,
			    const char*  HistName) const {
  return fTrEff;
}


//-----------------------------------------------------------------------------
// default: for backward compatibility
//-----------------------------------------------------------------------------
double aprocess::GetTrEffSF(const char*  Module ,
			    const char*  HistSet,
			    int          Bin    ,
			    const char*  HistName) const {
  return fTrEffSF;
}


//-----------------------------------------------------------------------------
// default: for backward compatibility
//-----------------------------------------------------------------------------
double aprocess::GetAcceptance(const char*  Module ,
			       const char*  HistSet,
			       int          Bin    ,
			       const char*  HistName) const {
  return fAcc;
}

//-----------------------------------------------------------------------------
// corrected acc = acc*corr_factor
// example: name='zee' title='ze0s01'
//-----------------------------------------------------------------------------
double aprocess::GetAccCorrFactor() const {
  double corr_factor;

  corr_factor = fDsMetadata->GetAccCorrFactor(GetName(),fMcFlag,fRunRange);

  return corr_factor;
}


//-----------------------------------------------------------------------------
// want to be able to redefine on the fly
//-----------------------------------------------------------------------------
int aprocess::GetColor() const {
  int color;

  color = fDsMetadata->GetHistColor(GetName(),fMcFlag);

  return color;
}

//-----------------------------------------------------------------------------
int aprocess::GetStyle() const {
  int style;

  style = fDsMetadata->GetHistStyle(GetName(),fMcFlag);

  return style;
}


//-----------------------------------------------------------------------------
const char* aprocess::GetLabel() const {

  return fDsMetadata->GetDsMetadata(GetName())->GetLabel();

}


//-----------------------------------------------------------------------------
// print information about a process
//-----------------------------------------------------------------------------
void aprocess::Print(const char*  Module  , 
		     const char*  HistSet ,
		     int          Bin     ,
		     const char*  HistName,
		     const char*  Opt     ) const {

  TString opt = Opt;

  if (opt.Index("banner") >= 0) {
    printf("------------------------------------------------------------------");
    printf("------------------------------------------------------------------\n");
    printf("process MC      xsec  intlumi      ACC    SF(l)  SF(t) ");
    printf("TR_EFF SF(TR)  N(total)  fQEvents   hist_file \n");
  } 

  printf("%-8s %1i", GetName(),fMcFlag);

  printf(" %9.2f %8.3f %10.3e",
	 fXSec,
	 fIntLumi,
	 fAnalysis->GetAcceptance(this,Module,HistSet,fRunRange,Bin,HistName));

				// meaning tight and loose....
  printf(" %6.4f %6.4f",
	 fAnalysis->GetIDEffSF(this,Module,HistSet,fRunRange,Bin,HistName),
	 fAnalysis->GetIDEffSF(this,Module,HistSet,fRunRange,Bin,HistName));

  printf(" %6.4f %6.4f %9.0f %9.2f",
	 fAnalysis->GetTrEff(this,Module,HistSet,fRunRange,Bin,HistName),
	 fAnalysis->GetTrEffSF(this,Module,HistSet,fRunRange,Bin,HistName),
	 fNTotal,
	 fAnalysis->GetFitQEvents(this,Module,HistSet,fRunRange,Bin,HistName));

  printf(" %s\n",fHistFileName.Data());
  
}


//-----------------------------------------------------------------------------
void aprocess::Print(const char*  Opt) const {
  printf("PROCESS : %s:%s\n",GetName(),GetTitle());
}
