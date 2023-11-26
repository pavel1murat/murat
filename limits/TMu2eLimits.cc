///////////////////////////////////////////////////////////////////////////////
// learn to use Tom's MCLIMIT
// needs libmclimit.so to be loaded, this is why it is a separate subpackage
// channel names: 
//
// 1. zz4l_mass
///////////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TObjString.h"
#include "TCanvas.h"

#include "math.h"
#include "string.h"
#include "limits/TMu2eLimits.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

ClassImp(TMu2eLimits)

//-----------------------------------------------------------------------------
TMu2eLimits::TMu2eLimits() {
  fChannelData = NULL;
}

//-----------------------------------------------------------------------------
int TMu2eLimits::InitChannelData(channel_data_t*& Data) {
//   TMu2eLimits::channel_data_t channel_data[] = {
// //------------------------------------------------------------------------------------------
// // define channels, channel with the name=0 - the last one
// // 
// //  ana init    bit      name           Module    HistSet HistBin  HistName  xmin xmax Rebin
// //------------------------------------------------------------------------------------------
//     { 0,  0,      0,   "test"    , "Mu2eLimits" ,  ""    ,  49,   "p"   ,    0.,    -1.,  1},
//     { 0,  0,   1<<0,   "mu2e"    , "Mu2eLimits" ,  "trk" ,  01,   "p"   ,   101.,   106., 1},
//     { 0,  0,      0,        0    , 0            ,  0     ,   0,        0,    0.,    -1.,  0}
//   };

  Data = new channel_data_t[3];
  // dont forget to initialize
  return 0;
}

//-----------------------------------------------------------------------------
// Mode =  0       : test
//      =  1       : mu2e
//
// Mass: mass of the expected resonance - chooses the MC histograms file
// as now we have several signals, need to be able to specify signal
//-----------------------------------------------------------------------------
TMu2eLimits::TMu2eLimits(int Mode, channel_data_t* ChannelData, const char* Signal) {
  // Mode defines which channels are used
  // Mass is the mass of a Mu2e resonance
  // ChannelData specifies which histograms to use for a given channel
  // Signal defines pt0 or pt100

					// 1. create a new mclimit_csm object
  fMcLimit   = new mclimit_csm();

  //  TH1::SetDefaultSumw2(1);

  fNullHyp   = 0;
  fNullHypPe = 0;
  fTestHyp   = 0;
  fTestHypPe = 0;
  fBgrLevel  = 0.01;
  fNBinsTest = 100; // 500;
					// 2013-12-21: use first lumi bin
  fMu2e      = new TMu2eChannel   ("mu2e_signal",0,Signal);

  fChannelName = new TObjArray();
  fModuleName  = new TObjArray();
  fHistName    = new TObjArray();
				        // count number of channels

  if (ChannelData == 0) {
    InitChannelData(fChannelData);
  }
  else                  fChannelData = ChannelData;

  fNChannels = 0;
  for (int i=0; fChannelData[i].name != 0; i++) {

    if      (strcmp(fChannelData[i].name,"test"   ) == 0) {
      fChannelData[i].ana           = fMu2e;
      fChannelData[i].init_function = TMu2eLimits::InitTestChannel;
    }
    else if (strcmp(fChannelData[i].name,"mu2e"   ) == 0) {
      fChannelData[i].ana           = fMu2e;
      fChannelData[i].init_function = TMu2eLimits::InitMu2eChannel;
    }
//-----------------------------------------------------------------------------
// selected channels are defined by Mode and their masks
//-----------------------------------------------------------------------------
    if (((Mode ==0) && (fChannelData[i].mask == 0)) || ((Mode & fChannelData[i].mask) != 0)) {
      printf(" >>>>> adding channel %s\n",fChannelData[i].name);
      fChannelName->Add(new TObjString(fChannelData[i].name));
      fModuleName ->Add(new TObjString(fChannelData[i].module_name) );
      fHistName   ->Add(new TObjString(fChannelData[i].hist_name)   );
      fNChannels++;
    }
  }

  fNPseudoExp   = 10000;
  Init();
}



//-----------------------------------------------------------------------------
TMu2eLimits::~TMu2eLimits() {

  if (fNullHyp  ) delete fNullHyp;
  if (fNullHypPe) delete fNullHypPe;

  if (fTestHyp  ) delete fTestHyp;
  if (fTestHypPe) delete fTestHypPe;

  delete fMcLimit;
}



//-----------------------------------------------------------------------------
// channels: 'zll_zll_mass' : mass distribution in 4l channel
//
// in each channel background templates come from Z+2jets,ttbar,WZ,ZZ,WW 
// 5 background templates per channel in total
// signal templates come from H(300) MC
//-----------------------------------------------------------------------------
int TMu2eLimits::Init() {

  int                rebin_factor;
  const char         *channel_name, *module_name, *hist_file;
  char               hist_name[100];
  TH1F*              hist(NULL);
  channel_data_t*    ch_data;
  aprocess           *dat; //, *bgr;

//-----------------------------------------------------------------------------
// create null and test hypotheses
//-----------------------------------------------------------------------------
  fNullHyp = new csm_model();
  fTestHyp = new csm_model();

  for (int i=0; i<fNChannels; i++) {
    channel_name = ((TObjString*) fChannelName->At(i))->String().Data();
    ch_data      = GetChannelData(channel_name);
//-----------------------------------------------------------------------------
// start from the data histogram, so far work only with 1D histograms
// how to make a 1D histogram for Mu2e ?
// assume only one run range
//-----------------------------------------------------------------------------
    if (strcmp(channel_name,"test") == 0) {
//-----------------------------------------------------------------------------
// test mode: data histogram
//-----------------------------------------------------------------------------
      hist = new TH1F("test_data","test data", fNBinsTest,101.,106.);
      for (int i=1; i<=fNBinsTest; i++) {
	hist->SetBinContent(i,fBgrLevel);
	hist->SetBinError  (i,fBgrLevel/10.);
      }
    }
    else if  (strcmp(channel_name,"mu2e") == 0) {
      module_name  = ch_data->module_name;
      sprintf(hist_name,"%s_%i/%s",
	      ch_data->hist_set,
	      ch_data->bin,
	      ch_data->hist_name);

      dat          = ch_data->ana->GetData(0);  // assume only one run range
      hist_file    = dat->GetHistFileName();

      if (hist_name[0] != 0) {
//-----------------------------------------------------------------------------
// so far only one Mu2e channel - momentum distribution
//-----------------------------------------------------------------------------
	hist = (TH1F*) gh1(hist_file,module_name,hist_name)->Clone();

	if (ch_data->xmin < ch_data->xmax) {
	  int nb = hist->GetNbinsX();
	  for (int i=1; i<=nb; i++) {
	    double x = hist->GetBinCenter(i);
	    if ((x < ch_data->xmin) || (x > ch_data->xmax)) {
	      hist->SetBinContent(i,0.);
	      hist->SetBinError  (i,0.);
	    }
	  }
	}

	rebin_factor = ch_data->rebin;
	if (rebin_factor != 1) {
	  hist->Rebin(rebin_factor);
	}
      }
    }

    fMcLimit->set_datahist(hist,channel_name);

    hist->Delete();
//-----------------------------------------------------------------------------
// for each channel establish null_hyp model
//-----------------------------------------------------------------------------
    fChannelModel[i] = new csm_channel_model(channel_name);

    if      (strcmp(channel_name,"test"   ) == 0)   { 
      ch_data->init_function(this, fChannelModel[i]);
      fNullHyp->add_chanmodel(fChannelModel[i],channel_name);
      fTestHyp->add_chanmodel(fChannelModel[i],channel_name);
      AddTestSignal(fTestHyp);
    }
    else {
//-----------------------------------------------------------------------------
// initialize fChanModel[i]
//-----------------------------------------------------------------------------
      ch_data->init_function(this, fChannelModel[i]);
      
      fNullHyp->add_chanmodel(fChannelModel[i],channel_name);
      fTestHyp->add_chanmodel(fChannelModel[i],channel_name);

      AddSignal(channel_name,fTestHyp);
    }
  }
//-----------------------------------------------------------------------------
// lock luminosity nuisance parameters for all models - later... if needed
//-----------------------------------------------------------------------------
//   char* pname[10];
//   pname[0] = "lumi";
//   int np

  fNullHypPe = (csm_model*) fNullHyp->Clone();
  fTestHypPe = (csm_model*) fTestHyp->Clone();

  fMcLimit->set_null_hypothesis(fNullHyp);
  fMcLimit->set_null_hypothesis_pe(fNullHypPe);
  fMcLimit->set_test_hypothesis(fTestHyp);
  fMcLimit->set_test_hypothesis_pe(fTestHypPe);

  return 0;
}


//-----------------------------------------------------------------------------
// add signal from a given 'Process' to a channel 'ChannelName' of a 'Model'
//-----------------------------------------------------------------------------
int   TMu2eLimits::RunPseudoExperiments(int NPseudoExp) {
//-----------------------------------------------------------------------------
// prepare to run pseudoexperiments
//-----------------------------------------------------------------------------
  if (NPseudoExp > 0) SetNPseudoExp(NPseudoExp);

  fMcLimit->setpxprintflag(true);
  fMcLimit->set_npe(fNPseudoExp);
  printf("# pseudoexperiments started\n");
  fMcLimit->run_pseudoexperiments();
  printf("# pseudoexperiments ended\n");

  // fTestHyp->print();

  printf("# Getting results\n");

  double tsobs  = fMcLimit->ts();
  double tsbm1  = fMcLimit->tsbm1();
  double tsbmed = fMcLimit->tsbmed();
  double tsbp1  = fMcLimit->tsbp1();
  double tssmed = fMcLimit->tssmed();
  double tssp1  = fMcLimit->tssp1();
  double tssm1  = fMcLimit->tssm1();

  printf("ts   : %15.8f\n",tsobs);
  printf("tsbm2: %15.8f\n",fMcLimit->tsbm2());
  printf("tsbm1: %15.8f\n",tsbm1);
  printf("tsbmed: %15.8f\n",tsbmed);
  printf("tsbp1: %15.8f\n",tsbp1);
  printf("tsbp2: %15.8f\n",fMcLimit->tsbp2());
  printf("tssm2: %15.8f\n",fMcLimit->tssm2());
  printf("tssm1:  %15.8f\n",tssm1);
  printf("tssmed: %15.8f\n",tssmed);
  printf("tssp1:  %15.8f\n",tssp1);
  printf("tssp2:  %15.8f\n",fMcLimit->tssp2());
  
  printf("CLs  :  %15.8f\n",fMcLimit->cls());
  printf("CLb  :  %15.8f\n",fMcLimit->clb());
  printf("1-CLb:  %15.8f\n",fMcLimit->omclb());
  printf("CLsb :  %15.8f\n",fMcLimit->clsb());
  
  printf("CLs -2sigma (bkg): %15.8f\n",fMcLimit->clsexpbm2());
  printf("CLs -1sigma (bkg): %15.8f\n",fMcLimit->clsexpbm1());
  printf("CLs median  (bkg): %15.8f\n",fMcLimit->clsexpbmed());
  printf("CLs +1sigma (bkg): %15.8f\n",fMcLimit->clsexpbp1());
  printf("CLs +2sigma (bkg): %15.8f\n",fMcLimit->clsexpbp2());
  
  printf("1-CLb -2sigma (sig): %15.8f\n",fMcLimit->omclbexpsm2());
  printf("1-CLb -1sigma (sig): %15.8f\n",fMcLimit->omclbexpsm1());

  double pval=fMcLimit->omclbexpsmed();
  printf("1-CLb median  (sig): %15.8f\n",pval);
  printf("CLb median    (sig): %15.8f\n",(1.-pval));
  printf("1-CLb +1sigma (sig): %15.8f\n",fMcLimit->omclbexpsp1());
  printf("1-CLb +2sigma (sig): %15.8f\n",fMcLimit->omclbexpsp2());
  
  double backSD = fabs(tsbp1-tsbm1)/2;
  double sigSD  = fabs(tssp1-tssm1)/2;
  
  double backMean = tsbmed;
  double sigMean  = tssmed;
  
  double backMeanErr = backSD/ sqrt(fNPseudoExp);
  double sigMeanErr  = sigSD / sqrt(fNPseudoExp);

  double backSDErr = backSD/sqrt(2*fNPseudoExp);
  double sigSDErr  = sigSD /sqrt(2*fNPseudoExp);

  double fom = (backMean - sigMean) / sqrt(backSD * backSD + sigSD * sigSD);

  double factorInFront = 1 / (backSD * backSD + sigSD * sigSD);
  double term1         = backMeanErr + 2 * backSD * backSDErr * fom;
  double term2         = sigMeanErr + 2 * sigSD * sigSDErr * fom;

  double fomErr        = factorInFront * sqrt(term1 * term1 + term2 * term2);

  printf("Figure of merit: %g +/- %g\n",fom,fomErr);
  
  /*
    printf("1-CLbw: " << fMcLimit->omclbw());
    printf("1-CLbw -2sigma (sig): " << fMcLimit->omclbexpsm2w());
    printf("1-CLbw -1sigma (sig): " << fMcLimit->omclbexpsm1w());
    printf("1-CLbw median  (sig): " << fMcLimit->omclbexpsmedw());
    printf("1-CLbw +1sigma (sig): " << fMcLimit->omclbexpsp1w());
    printf("1-CLbw +2sigma (sig): " << fMcLimit->omclbexpsp2w());
  */

  return 0;
}


//-----------------------------------------------------------------------------
// add signal from a given 'Process' to a channel 'ChannelName' of a 'Model'
//-----------------------------------------------------------------------------
int   TMu2eLimits::AddSignal(const char* ChannelName, csm_model*  Model) {

  TH1             /**h1(0), *h2(0),*/  *h11(0) ; //, *h12, *h15;
  int              np;
  const char      *pnames [10];
  double           par_sf_lo    [10];
  double           par_sf_hi    [10];
  TH1*             hist_shape_lo[10];
  double           nsig_shape_lo[10];
  TH1*             hist_shape_hi[10];
  double           nsig_shape_hi[10];
  double           sf;
  //  const char       *pname;
  //  char             signal_name[100];

  channel_data_t*  ch_data;
  aprocess         /* *process,*/ *signal;
  analysis         *ana;
//-----------------------------------------------------------------------------
// for Mu2e  
//-----------------------------------------------------------------------------
  ch_data = GetChannelData(ChannelName);
//-----------------------------------------------------------------------------
// all these datasets should correspond to the same luminosity and have 
// the same number of run ranges defined


//-----------------------------------------------------------------------------
  ana = ch_data->ana;
					// round-off down to 1 GeV

  //  sprintf(signal_name,"gzz%i",(int)fMass);

  signal = ana->GetSignalMC(0); // new aprocess(signal_name,1,0,ana); 
//-----------------------------------------------------------------------------
// this is needed for the histogram to be normalized
//-----------------------------------------------------------------------------
  signal->get_h1_unorm(ch_data->module_name,
		       ch_data->hist_set,
		       ch_data->bin,
		       ch_data->hist_name,
		       ch_data->rebin,
		       h11,
		       &sf);
//-----------------------------------------------------------------------------
// kludge : replace number of entries in the histogram by the fit result
//-----------------------------------------------------------------------------
  // double xxx [] = {
  //   250.,   220.,   
  //   300.,   201.,   
  //   325.,   236.,   
  //   350.,   238.,   
  //   400.,   241.,   
  //   500.,   225.,   
  //   600.,   222.,   
  //   700.,   210.,   
  //   800.,   233.,   
  //   900.,   223.,   
  //  1000.,   157.,
  //     -1.
  // };

  //  double coef [2] = { 246.9, -0.0542 };

  if (strcmp(ChannelName,"mu2e") == 0) {

    //    int imass(0);

    // for (int i=0; i<11; i++) {
    //   if ( fabs(xxx[2*i]-fMass) < 0.01) {
    // 	imass = i;
    //   }
    // }

    // // sf = sf * xxx[2*imass+1] / (coef[0] + coef[1]*fMass);

    // h11->Scale((coef[0] + coef[1]*fMass)/xxx[2*imass+1]);
    
  }
//-----------------------------------------------------------------------------
// optionally, assuming 1D histograms, define a used range 
//-----------------------------------------------------------------------------
  if (ch_data->xmin < ch_data->xmax) {
    int nb = h11->GetNbinsX();
    for (int i=1; i<=nb; i++) {
      double x = h11->GetBinCenter(i);
      if ((x < ch_data->xmin) || (x > ch_data->xmax)) {
	h11->SetBinContent(i,0.);
	h11->SetBinError  (i,0.);
      }
    }
  }
				// an option: for Mu2e truncate the signal distribution

  if (strcmp(ChannelName,"mu2e") == 0) {
    int nb = h11->GetNbinsX();
    for (int i=nb; i>=0; i--) {
      
      // double x = h11->GetBinCenter(i);
      //      if ((x < fMass-150.) || (x > fMass+150)) {
      // if ((x < fMass-100.) || (x > fMass+100)) {
      // 	h11->SetBinContent(i,0.);
      // 	h11->SetBinError  (i,0.);
      // }
    }
  }
//-----------------------------------------------------------------------------
// consider only one systematic uncertainty, corresponding to the luminosity
// histogram is already normalized to the yield ("X-section")
//-----------------------------------------------------------------------------
  np               = 1 ; // 1;
  pnames[0]        = "lumi";
  par_sf_lo[0]     = -0.1;
  par_sf_hi[0]     =  0.1;
  hist_shape_lo[0] =  NULL;
  hist_shape_hi[0] =  NULL;
  nsig_shape_lo[0] = 1.;
  nsig_shape_hi[0] = 1.;

  Model->add_template(h11,
		      sf,
		      np,
		      pnames,
		      par_sf_lo,
		      par_sf_hi,
		      hist_shape_lo,
		      nsig_shape_lo,
		      hist_shape_hi,
		      nsig_shape_hi,
		      // 1,  // Poisson errors
		      2,                 // use error bars stored in the histogram
		      1,                 // this is signal
		      ChannelName);
  return 0;
}



//-----------------------------------------------------------------------------
int TMu2eLimits::Poisson95CL(int NPseudoExp, int PrintPxFlag) {

  double chi2, s95;

  if (NPseudoExp > 0) SetNPseudoExp(NPseudoExp);

  fMcLimit->setpxprintflag(PrintPxFlag);

  chi2 = 1.; // fMcLimit->chisquared();
  s95  = fMcLimit->s95();

  printf("chi2 = %10.3f, s95=%10.5f\n",chi2,s95);
	
  return 0;
}

//-----------------------------------------------------------------------------
int TMu2eLimits::Bayes95CL(int NPseudoExp, int PrintPxFlag) {

  double chi2, s95, s95_err, sm2,sm1,smed,sp1,sp2;

  if (NPseudoExp > 0) SetNPseudoExp(NPseudoExp);

  fMcLimit->setpxprintflag(PrintPxFlag);
  fMcLimit->set_npe(fNPseudoExp);

  chi2 = 1.; // fMcLimit->chisquared();

  fMcLimit->bayes_heinrich_withexpect(0.95,&s95,&s95_err,
				      fNPseudoExp,
				      &sm2,&sm1,&smed,&sp1,&sp2);

  printf("chi2 = %10.3f, OBSERVED: s95=%10.5f +/- %10.5f ",chi2,s95,s95_err);
  printf("EXPECTED: sm2=%10.5f sm1=%10.5f smed=%10.5f  smp1=%10.5f smp2=%10.5f\n",
	 sm2,sm1,smed,sp1,sp2);
	
  return 0;
}

//-----------------------------------------------------------------------------
// CL = 0.95, 0.90 etc
//-----------------------------------------------------------------------------
int TMu2eLimits::BayesCL(double CL, int NPseudoExp, int PrintPxFlag) {

  double chi2, s, s_err, sm2,sm1,smed,sp1,sp2;

  if (NPseudoExp > 0) SetNPseudoExp(NPseudoExp);

  fMcLimit->setpxprintflag(PrintPxFlag);
  fMcLimit->set_npe(fNPseudoExp);

  chi2 = 1.; // fMcLimit->chisquared();

  fMcLimit->bayes_heinrich_withexpect(CL,&s,&s_err,
				      fNPseudoExp,
				      &sm2,&sm1,&smed,&sp1,&sp2);

  printf("chi2 = %10.3f, CL= %8.3f OBSERVED: s=%10.5f +/- %10.5f ",chi2,CL,s,s_err);
  printf("EXPECTED: sm2=%10.5f sm1=%10.5f smed=%10.5f  smp1=%10.5f smp2=%10.5f\n",
	 sm2,sm1,smed,sp1,sp2);
	
  return 0;
}

//-----------------------------------------------------------------------------
int TMu2eLimits::PlotHypothesis(const char* ChannelName    , 
				const char* Hypothesis     , 
				int         CreateNewCanvas) {

  char name[1000], message[1000];

  TCanvas* c;
  if (CreateNewCanvas) {
    sprintf(name,"%s:%s",ChannelName,Hypothesis);
    c = new TCanvas(name,name,500,500);
    c->cd();
  }

  TH1F* hist = (TH1F*) fMcLimit->get_datahist(ChannelName);

  if (hist != 0) {
    if (strcmp(Hypothesis,"null") == 0) {
      fNullHypPe->plotwithdata(ChannelName,hist);
    }
    else if (strcmp(Hypothesis,"test") == 0) {
      fTestHypPe->plotwithdata(ChannelName,hist);
    }
  }
  else {
    sprintf(message,"HISTOGRAM for ChannelName=%s Hypothesis=%s DOESN\'T EXIST",
	    ChannelName,Hypothesis);
    Error("PlotHypothesis",message);
  }

  return 0;
}


//-----------------------------------------------------------------------------
TMu2eLimits::channel_data_t*  TMu2eLimits::GetChannelData(const char* ChannelName) {
  
  channel_data_t*   dat(0);

  for (int i=0; fChannelData[i].name != 0; i++) {
    if      (strcmp(fChannelData[i].name,ChannelName) == 0) {
      dat = &fChannelData[i];
      break;
    }
  }
 
  if (dat == 0) {
    Error("GetChannelData",Form("%s NOT FOUND, RETURN NULL POINTER",ChannelName));
  }

  return dat;
}
