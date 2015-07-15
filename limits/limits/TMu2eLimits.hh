#ifndef murat_limits_TMu2eLimits
#define murat_limits_TMu2eLimits

#include "TNamed.h"
#include "TString.h"
#include "TH1.h"

#include "murat/mclimit/mclimit_csm.h"

#include "murat/obj/TMu2eChannel.hh"
#include "murat/obj/aprocess.hh"

class TMu2eLimits: public TNamed {

public:

//-----------------------------------------------------------------------------
  struct channel_data_t {
    analysis* ana;
    int          (*init_function) (TMu2eLimits* , csm_channel_model* ChannelModel);
    int          mask;
    char*        name;
    char*        module_name;
    char*        hist_set;
    int          bin;
    const char*  hist_name;
    float        xmin;
    float        xmax;
    int          rebin;
  };

  mclimit_csm*    fMcLimit;
  csm_model*      fNullHyp;
  csm_model*      fNullHypPe;
  csm_model*      fTestHyp;
  csm_model*      fTestHypPe;

  csm_channel_model*  fChannelModel[100];          // to begin with

  int             fNChannels;

  channel_data_t* fChannelData;

  TObjArray*      fChannelName;
  TObjArray*      fModuleName;
  TObjArray*      fHistName;

  TMu2eChannel*    fMu2e;

  int              fNPseudoExp;              // number of pseudo-experiments

//   TH1*            fDataHist;
//   TH1*            fBgrHist;
//   TH1*            fSigHist;

  double          fXMax;
  double          fXMin;
  int             fNBins;

  double          fBgrLevel;               // for test mode
  int             fNBinsTest;		   // for test mode
//-----------------------------------------------------------------------------
// constructors and destructor, allow to redefine channels from the script
//-----------------------------------------------------------------------------
  TMu2eLimits();
  TMu2eLimits(int Mode, channel_data_t* ChannelData = 0, const char* Signal = "gpt000");

  ~TMu2eLimits();

  int   Init();

  static int   InitTestChannel    (TMu2eLimits* Tzz, csm_channel_model* ChannelModel);
  static int   InitMu2eChannel    (TMu2eLimits* Tzz, csm_channel_model* ChannelModel);
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TMu2eChannel* GetMu2eChannel() { return fMu2e; }
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  int   InitChannelData(channel_data_t*& ChannelData);
  int   AddSignal    (const char* ChannelName, csm_model* Model);
  int   AddTestSignal(csm_model* Model);

  void  SetNPseudoExp(int NPseudoExp) {
    fNPseudoExp = NPseudoExp;
    fMcLimit->set_npe(NPseudoExp);
  }

  int   RunPseudoExperiments(int NPseudoExp = -1);
  int   Poisson95CL         (int NPseudoExp = -1, int PrintPxFlag = 0);
  int   Bayes95CL           (int NPseudoExp = -1, int PrintPxFlag = 0);
  int   BayesCL             (double CL, int NPseudoExp = -1, int PrintPxFlag = 0);

				// Hypothesis = "null" or "test"

  int   PlotHypothesis(const char* ChannelName        , 
		       const char* Hypothesis         , 
		       int         CreateNewCanvas = 0);

  channel_data_t*   GetChannelData(const char* ChannelName);

  ClassDef(TMu2eLimits,0) 

};

#endif
