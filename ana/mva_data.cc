///////////////////////////////////////////////////////////////////////////////
// framework
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

// Xerces XML Parser
#include <xercesc/dom/DOM.hpp>

#include "Mu2eUtilities/inc/MVATools.hh"

#include <stdlib.h>
#include "murat/ana/mva_data.hh"
#include "TString.h"

using std::string;
using std::vector;
//-----------------------------------------------------------------------------
// v573  trained MVA's - compared to Dave's TrkQual with 0.4 cut
// encoding of the last part of the name: (as integer)
// 1. logfcons used for training:  weight_mode
// 2. chi2d    used for training: 100+weight_mode
// weight_mode: see 
//-----------------------------------------------------------------------------

mva_data::data_t  mva_calpatrec_e11s5731_000 = {
  "calpatrec_e11s5731_logfcons_000",
  "CalPatRec/data/v5_7_7/MLP_weights_logfcons_0_uni.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_000.hist",
  0.60, 12
};

mva_data::data_t  mva_calpatrec_e11s5731_001 = {
  "calpatrec_e11s5731_logfcons_001",
  "CalPatRec/data/v5_7_7/MLP_weights_logfcons_1_lin.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_001.hist",
  0.65, 13
};

mva_data::data_t  mva_calpatrec_e11s5731_002 = {
  "calpatrec_e11s5731_logfcons_002",
  "CalPatRec/data/v5_7_7/MLP_weights_logfcons_2_exp.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_002.hist",
  0.85, 17
};

mva_data::data_t  mva_calpatrec_e11s5731_003 = {
  "calpatrec_e11s5731_logfcons_003",
  "CalPatRec/data/v5_7_7/MLP_weights_logfcons_3_pol.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_003.hist",
  0.85, 17  // not sure
};

mva_data::data_t  mva_calpatrec_e11s5731_004 = {
  "calpatrec_e11s5731_logfcons_004",
  "CalPatRec/data/v5_7_7/MLP_weights_logfcons_4_exp.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_004.hist",
  0.85, 17  // not sure
};

mva_data::data_t  mva_calpatrec_e11s5731_100 = {
  "calpatrec_e11s5731_chi2d_100",
  "CalPatRec/data/v5_7_7/MLP_weights_chi2d_0_uni.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_100.hist",
  0.65, 13  // not sure
};

mva_data::data_t  mva_calpatrec_e11s5731_101 = {
  "calpatrec_e11s5731_chi2d_101",
  "CalPatRec/data/v5_7_7/MLP_weights_chi2d_1_lin.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_101.hist",
  0.65, 13  // not sure
};

mva_data::data_t  mva_calpatrec_e11s5731_102 = {
  "calpatrec_e11s5731_chi2d_102",
  "CalPatRec/data/v5_7_7/MLP_weights_chi2d_2_exp.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_102.hist",
  0.65, 13  // not sure
};

mva_data::data_t  mva_calpatrec_e11s5731_103 = {
  "calpatrec_e11s5731_chi2d_103",
  "CalPatRec/data/v5_7_7/MLP_weights_chi2d_3_pol.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_103.hist",
  0.65, 13  // not sure
};

mva_data::data_t  mva_calpatrec_e11s5731_104 = {
  "calpatrec_e11s5731_chi2d_104",
  "CalPatRec/data/v5_7_7/MLP_weights_chi2d_4_exp.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_104.hist",
  0.65, 13  // not sure
};

mva_data::data_t  mva_calpatrec_e11s5731_202 = {
  "calpatrec_e11s5731_chi2d_202",
  "CalPatRec/data/v5_7_7/MLP_weights_chi2d_2_exp.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_202.hist",
  0.65, 13  // not sure
};

mva_data::data_t  mva_calpatrec_e11s5731_204 = {
  "calpatrec_e11s5731_chi2d_204",
  "../../alaha/dev/TrkQualPExp4Weights/TMVAClassification_MLP.weights.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_204.hist",
  0.65, 13  // not sure
};

mva_data::data_t  mva_trkpatrec_dave_002 = {
  "trkpatrec_dave_logfcons_002",
  "AnalysisConditions/weights/TrkQual.weights.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_002.hist",              // not correct
  0.40, 8
};

mva_data::data_t  mva_trkpatrec_e115731_001 = {
  "trkpatrec_e11s5731_logfcons_001",
  "CalPatRec/data/v5_7_7/MLP_weights_trkpatrec_logfcons_1_uni.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_001_000.hist",
  0.40, 8
};

//-----------------------------------------------------------------------------
// 2020-12-22 back to training, default settings:
//
// TCut signal_cuts("(tdip>0.5)&&(tdip<1.0)&&(t0err<5)&&(fabs(p-pmc)<0.25)");
// TCut bkg_cuts   ("(tdip>0.5)&&(tdip<1.0)&&(t0err<5)&&(p-pmc)>0.7"     );   or 0.6
//
// TString training_opt = "nTrain_Signal=50000:nTrain_Background=20000";
// training_opt        += ":nTest_Signal=50000:nTest_Background=20000";
//-----------------------------------------------------------------------------
mva_data::data_t  mva_fele2s51b1_0070 = {
  "fele2s51b1_0070",                                                         // name
  "su2020/data/trk_qual_mva/MLP_weights_0070.xml",                           // location of XML weights
  "$MU2E_HIST/su2020/su2020.fele2s51b1.track_comp_use_mva_0070.hist",        // histograms
  0.40, 8                                                                    // cut value, ID #
};


mva_data::data_t  mva_fele2s51b1_0060 = {
  "fele2s51b1_060",                                                          // name
  "su2020/data/trk_qual_mva/MLP_weights_0060.xml",                           // location of XML weights
  "$MU2E_HIST/su2020/su2020.fele2s51b1.track_comp_use_mva_0060.hist",        // histograms
  0.80, 16                                                                   // cut value, ID #
};

mva_data::data_t  mva_fele2s51b1_1055 = {                                    // 
  "fele2s51b1_1055",                                                         // name
  "su2020/data/trk_qual_mva/MLP_weights_1055.xml",                           // location of XML weights
  "$MU2E_HIST/su2020/su2020.fele2s51b1.track_comp_use_mva_1055.hist",        // histograms
  0.20, 4                                                                    // cut value, ID #
};

mva_data::data_t  mva_fele2s51b1_1060 = {                                    // 
  "fele2s51b1_1060",                                                         // name
  "su2020/data/trk_qual_mva/MLP_weights_1060.xml",                           // location of XML weights
  "$MU2E_HIST/su2020/su2020.fele2s51b1.track_comp_use_mva_1060.hist",        // histograms
  0.20, 4                                                                    // cut value, ID #
};

mva_data::data_t  mva_fele2s51b1_1070 = {
  "fele2s51b1_1070",                                                         // name
  "su2020/data/trk_qual_mva/MLP_weights_1070.xml",                           // location of XML weights
  "$MU2E_HIST/su2020/su2020.fele2s51b1.track_comp_use_mva_1070.hist",        // histograms
  0.20, 4                                                                    // cut value, ID #
};


//-----------------------------------------------------------------------------
mva_data::~mva_data() {
}

//-----------------------------------------------------------------------------
mva_data::mva_data() {
}

//-----------------------------------------------------------------------------
// CalPatRec ANN training sets
// ---------------------------
// Training Codes: 
//
// 0060 : PAR dPf > 0.60
// 0070 : PAR dPf > 0.70
// 1060 : DAR dPf > 0.60
// 1070 : DAR dPf > 0.70
//-----------------------------------------------------------------------------
mva_data::mva_data(const char* Dataset, int TrainingCode) {

  int error(0);

  fTrainingCode = TrainingCode;
  
  TString ds = Dataset;
  ds.ToUpper();

  if (ds == "E11S5731") {
    if      (TrainingCode ==   0) fData = mva_calpatrec_e11s5731_000;
    else if (TrainingCode ==   1) fData = mva_calpatrec_e11s5731_001;
    else if (TrainingCode ==   2) fData = mva_calpatrec_e11s5731_002;
    else if (TrainingCode ==   3) fData = mva_calpatrec_e11s5731_003;
    else if (TrainingCode ==   4) fData = mva_calpatrec_e11s5731_004;
    else if (TrainingCode == 100) fData = mva_calpatrec_e11s5731_100;
    else if (TrainingCode == 101) fData = mva_calpatrec_e11s5731_101;
    else if (TrainingCode == 102) fData = mva_calpatrec_e11s5731_102;
    else if (TrainingCode == 103) fData = mva_calpatrec_e11s5731_103;
    else if (TrainingCode == 104) fData = mva_calpatrec_e11s5731_104;
    else if (TrainingCode == 202) fData = mva_calpatrec_e11s5731_202;
    else if (TrainingCode == 204) fData = mva_calpatrec_e11s5731_204;
    else                  { error = 1; }
  }
  else if (ds = "FELE2S51B1") {
    if      (TrainingCode ==   60) fData = mva_fele2s51b1_0060;  // PAR, logfcons, uniform weight, dPf > 0.6
    else if (TrainingCode ==   70) fData = mva_fele2s51b1_0070;  // PAR, logfcons, uniform weight, dPf > 0.7
    else if (TrainingCode == 1055) fData = mva_fele2s51b1_1055;  // DAR, logfcons, w= 1/dPf      , dPf > 0.5
    else if (TrainingCode == 1060) fData = mva_fele2s51b1_1060;  // DAR, logfcons, uniform weight, dPf > 0.6
    else if (TrainingCode == 1070) fData = mva_fele2s51b1_1070;  // DAR, logfcons, uniform weight, dPf > 0.7
    else                  { error = 1; }
  }
  else                    { error = 1; }

  if (error != 0) {
    printf(" >>> ERROR in mva_data::mva_data(const char*,const char*,int): algorithm : %s training code: %5i. BAIL OUT\n",
	   Dataset,TrainingCode);
  }
}



//-----------------------------------------------------------------------------
int mva_data::Init() {
  string s1(this->XmlWeightsFile());

  fhicl::ParameterSet pset;

  printf(">>> [mva_data::Init] Init MVA from %s\n",s1.data());
  pset.put<string>("MVAWeights",s1);

  fMva = new mu2e::MVATools(pset);

  fMva->initMVA();
  fMva->showMVA();

  return 0;
}
