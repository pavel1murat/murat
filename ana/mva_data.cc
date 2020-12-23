//


#include <stdlib.h>
#include "murat/ana/mva_data.hh"
#include "TString.h"

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
// back to training
//-----------------------------------------------------------------------------
mva_data::data_t  mva_tpr_fele2s51b1_000 = {
  "tpr_fele2s51b1_000",                                                     // name
  "su2020/data/trk_qual_mva/MLP_weights_tpr_logfcons_000.xml",             // location of XML weights
  "$MU2E_HIST/su2020/su2020.fele2s51b1.track_comp_use_mva_0000.hist",        // histograms
  0.40, 8                                                                   // cut value, ID #
};


mva_data::data_t  mva_tpr_fele2s51b1_001 = {
  "tpr_fele2s51b1_001",                                                     // name
  "su2020/data/trk_qual_mva/MLP_weights_tpr_logfcons_001.xml",             // location of XML weights
  "$MU2E_HIST/su2020/su2020.fele2s51b1.track_comp_use_mva_0001.hist",        // histograms
  0.80, 16                                                                   // cut value, ID #
};

mva_data::data_t  mva_tpr_fele2s51b1_002 = {
  "tpr_fele2s51b1_002",                                                     // name
  "su2020/data/trk_qual_mva/MLP_weights_tpr_logfcons_002.xml",             // location of XML weights
  "$MU2E_HIST/su2020/su2020.fele2s51b1.track_comp_use_mva_0002.hist",        // histograms
  0.80, 16                                                                   // cut value, ID #
};

mva_data::data_t  mva_tpr_fele2s51b1_003 = {
  "tpr_fele2s51b1_001",                                                     // name
  "su2020/data/trk_qual_mva/MLP_weights_tpr_logfcons_003.xml",             // location of XML weights
  "$MU2E_HIST/su2020/su2020.fele2s51b1.track_comp_use_mva_0003.hist",        // histograms
  0.80, 16                                                                  // cut value, ID #
};

mva_data::data_t  mva_cpr_fele2s51b1_000 = {
  "cpr_fele2s51b1_000",                                                     // name
  "su2020/data/trk_qual_mva/MLP_weights_cpr_logfcons_000.xml",             // location of XML weights
  "$MU2E_HIST/su2020/su2020.fele2s51b1.track_comp_use_mva_1000.hist",        // histograms
  0.40, 8                                                                   // cut value, ID #
};


mva_data::data_t  mva_cpr_fele2s51b1_001 = {
  "cpr_fele2s51b1_001",                                                     // name
  "su2020/data/trk_qual_mva/MLP_weights_cpr_logfcons_001.xml",             // location of XML weights
  "$MU2E_HIST/su2020/su2020.fele2s51b1.track_comp_use_mva_1001.hist",        // histograms
  0.40, 8                                                                   // cut value, ID #
};

mva_data::data_t  mva_cpr_fele2s51b1_002 = {
  "cpr_fele2s51b1_002",                                                     // name
  "su2020/data/trk_qual_mva/MLP_weights_cpr_logfcons_002.xml",             // location of XML weights
  "$MU2E_HIST/su2020/su2020.fele2s51b1.track_comp_use_mva_1002.hist",        // histograms
  0.40, 8                                                                   // cut value, ID #
};

mva_data::data_t  mva_cpr_fele2s51b1_003 = {
  "cpr_fele2s51b1_003",                                                     // name
  "su2020/data/trk_qual_mva/MLP_weights_cpr_logfcons_003.xml",             // location of XML weights
  "$MU2E_HIST/su2020/su2020.fele2s51b1.track_comp_use_mva_1003.hist",        // histograms
  0.80, 16                                                                  // cut value, ID #
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
// 000-004 : log(fcons) with different weights 
// 101-104 : chi2d      with different weights
// 202     : chi2       training by Arpan
//-----------------------------------------------------------------------------
mva_data::mva_data(const char* TrkRecAlgorithm, const char* Dataset, int Type) {

  int error(0);
  
  TString ds = Dataset;
  ds.ToUpper();

  TString algo = TrkRecAlgorithm;
  algo.ToUpper();

  if (algo == "CPR") {
    if (ds == "E11S5731") {
      if      (Type ==   0) fData = mva_calpatrec_e11s5731_000;
      else if (Type ==   1) fData = mva_calpatrec_e11s5731_001;
      else if (Type ==   2) fData = mva_calpatrec_e11s5731_002;
      else if (Type ==   3) fData = mva_calpatrec_e11s5731_003;
      else if (Type ==   4) fData = mva_calpatrec_e11s5731_004;
      else if (Type == 100) fData = mva_calpatrec_e11s5731_100;
      else if (Type == 101) fData = mva_calpatrec_e11s5731_101;
      else if (Type == 102) fData = mva_calpatrec_e11s5731_102;
      else if (Type == 103) fData = mva_calpatrec_e11s5731_103;
      else if (Type == 104) fData = mva_calpatrec_e11s5731_104;
      else if (Type == 202) fData = mva_calpatrec_e11s5731_202;
      else if (Type == 204) fData = mva_calpatrec_e11s5731_204;
      else                  { error = 1; }
    }
    else if (ds = "FELE2S51B1") {
      if      (Type ==   0) fData = mva_cpr_fele2s51b1_000;  // logfcons, uniform weight
      else if (Type ==   1) fData = mva_cpr_fele2s51b1_001;  // logfcons, linear  weight
      else if (Type ==   2) fData = mva_cpr_fele2s51b1_002;  // logfcons, linear  weight
      else if (Type ==   3) fData = mva_cpr_fele2s51b1_003;  // logfcons, linear  weight
      else                  { error = 1; }
    }
    else                    { error = 1; }
  }
  else if (algo == "TPR") {
    printf(" ------ mva_data::mva_data tpr\n");
    if (ds == "DAVE") {
      if       (Type == 2)  fData = mva_trkpatrec_dave_002;
      else                  { error = 1; }
    }
    else if (ds == "E11S5731") {
      if       (Type == 1)  fData = mva_trkpatrec_e115731_001;
      else                  { error = 1; }
    }
    else if (ds == "FELE2S51B1") {
      printf(" ------ mva_data::mva_data fele2s51b1\n");
      if       (Type == 0)   fData = mva_tpr_fele2s51b1_000;  // logfcons, uniform weight
      else if  (Type == 1)   fData = mva_tpr_fele2s51b1_001;  // logfcons, linear  weight
      else if  (Type == 2)   fData = mva_tpr_fele2s51b1_002;  // logfcons, linear  weight
      else if  (Type == 3)   fData = mva_tpr_fele2s51b1_003;  // logfcons, linear  weight
      else                   { error = 1; }
    }
    else                     { error = 1; }
  }
  else                       { error = 1; }

  if (error != 0) {
    printf(" >>> ERROR in mva_data::mva_data(const char*,const char*,int): algorithm : %s dataset: %s type: %5i. BAIL OUT\n",
	   TrkRecAlgorithm,Dataset,Type);
  }
}

