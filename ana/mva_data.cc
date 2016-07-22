//


#include <stdlib.h>
#include "murat/ana/mva_data.hh"
#include "TString.h"

//-----------------------------------------------------------------------------
// v573  trained MVA's - compared to Dave's TrkQual with 0.4 cut
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
  "TrkDiag/test/TrkQual.weights.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_002.hist",
  0.40, 8
};

mva_data::data_t  mva_trkpatrec_e115731_001 = {
  "trkpatrec_e11s5731_logfcons_001",
  "CalPatRec/data/v5_7_7/MLP_weights_trkpatrec_logfcons_1_uni.xml",
  "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_001_000.hist",
  0.40, 8
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

  TString ds = Dataset;
  ds.ToUpper();

  TString algo = TrkRecAlgorithm;
  algo.ToUpper();

  if (algo == "CALPATREC") {
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
    }
  }
  else if (algo == "TRKPATREC") {
    if (ds == "DAVE") {
      if (Type == 2)  fData = mva_trkpatrec_dave_002;
    }
    else if (ds = "E11S5731") {
      if (Type == 1) fData = mva_trkpatrec_e115731_001;
    }
  }
  else {
    printf(" >>> ERROR in mva_data::mva_data(const char*,int): unknown algorithm : %s. BAIL OUT\n",
	   TrkRecAlgorithm);
  }
}

