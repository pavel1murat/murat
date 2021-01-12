///////////////////////////////////////////////////////////////////////////////
// note, that "andy" MVA's have less variables and have not been tested
///////////////////////////////////////////////////////////////////////////////
#include "murat/ana/mva_data.hh"
#include "TString.h"
//-----------------------------------------------------------------------------
// 2020-12-22 back to training, default settings:
//
// TCut signal_cuts("(tdip>0.5)&&(tdip<1.0)&&(t0err<5)&&(fabs(p-pmc)<0.25)");
// TCut bkg_cuts   ("(tdip>0.5)&&(tdip<1.0)&&(t0err<5)&&(p-pmc)>0.7"     );   or 0.6
//
// TString training_opt = "nTrain_Signal=50000:nTrain_Background=20000";
// training_opt        += ":nTest_Signal=50000:nTest_Background=20000";
//-----------------------------------------------------------------------------
mva_data::data_t  trq_mva_andy_0000 = {
  "andy_0000",                                                               // name
  "nactive:nafract:log10(fcons):momerr:t0err:nza_o_na:nma_o_nm",             // assume all floats
  "AnalysisConditions/weights/TrkQualNeg.weights.xml",                       // location of XML weights
  0.80                                                                       // default cut value
};

mva_data::data_t  trq_mva_andy_0100 = {
  "andy_0100",                                                               // name
  "nactive:nafract:log10(fcons):momerr:t0err:nza_o_na:nma_o_nm",             // assume all floats
  "AnalysisConditions/weights/TrkQualPos.weights.xml",                       // location of XML weights
  0.80                                                                       // default cut value
};

mva_data::data_t  trq_mva_fele2s51b1_0070 = {
  "fele2s51b1_0070",                                                         // name
  "nactive:nafract:log10(fcons):momerr:t0err:nda_o_na:nza_o_na:nma_o_nm",    // assume all floats
  "su2020/data/trq_mva/MLP_weights_0070.xml",                                // location of XML weights
  0.40                                                                       // default cut value
};

mva_data::data_t  trq_mva_fele2s51b1_0060 = {
  "fele2s51b1_060",                                                          // name
  "nactive:nafract:log10(fcons):momerr:t0err:nda_o_na:nza_o_na:nma_o_nm",    // assume all floats
  "su2020/data/trq_mva/MLP_weights_0060.xml",                                // location of XML weights
  0.80                                                                       // default cut value
};

mva_data::data_t  trq_mva_fele2s51b1_1055 = {                                    // 
  "fele2s51b1_1055",                                                         // name
  "nactive:nafract:log10(fcons):momerr:t0err:nda_o_na:nza_o_na:nma_o_nm",    // assume all floats
  "su2020/data/trq_mva/MLP_weights_1055.xml",                                // location of XML weights
  0.20                                                                       // default cut value
};

mva_data::data_t  trq_mva_fele2s51b1_1060 = {                                    // 
  "fele2s51b1_1060",                                                         // name
  "nactive:nafract:log10(fcons):momerr:t0err:nda_o_na:nza_o_na:nma_o_nm",    // assume all floats
  "su2020/data/trq_mva/MLP_weights_1060.xml",                                // location of XML weights
  0.20                                                                       // default cut value
};

mva_data::data_t  trq_mva_fele2s51b1_1070 = {
  "fele2s51b1_1070",                                                         // name
  "nactive:nafract:log10(fcons):momerr:t0err:nda_o_na:nza_o_na:nma_o_nm",    // assume all floats
  "su2020/data/trq_mva/MLP_weights_1070.xml",                                // location of XML weights
  0.20                                                                       // default default cut value
};

mva_data::data_t  trq_mva_fpos2s51b1_0150 = {
  "fpos2s51b1_0150",                                                         // positrons, PAR dpF>0.5 MeV/c
  "nactive:nafract:log10(fcons):momerr:t0err:nda_o_na:nza_o_na:nma_o_nm",    // assume all floats
  "su2020/data/trq_mva/MLP_weights_0150.xml",                                // location of XML weights
  0.20                                                                       // default default cut value
};

mva_data::data_t  trq_mva_fpos2s51b1_0170 = {
  "fpos2s51b1_0170",                                                         // positrons, PAR dpF>0.7 MeV/c
  "nactive:nafract:log10(fcons):momerr:t0err:nda_o_na:nza_o_na:nma_o_nm",    // assume all floats
  "su2020/data/trq_mva/MLP_weights_0170.xml",                                // location of XML weights
  0.20                                                                       // default default cut value
};

mva_data::data_t  trq_mva_fpos2s51b1_1150 = {
  "fpos2s51b1_1150",                                                         // positrons, DAR dpF>0.5 MeV/c
  "nactive:nafract:log10(fcons):momerr:t0err:nda_o_na:nza_o_na:nma_o_nm",    // assume all floats
  "su2020/data/trq_mva/MLP_weights_1150.xml",                                // location of XML weights
  0.20                                                                       // default default cut value
};

mva_data::data_t  trq_mva_fpos2s51b1_1170 = {
  "fpos2s51b1_1170",                                                         // positrons, DAR dpF>0.7 MeV/c
  "nactive:nafract:log10(fcons):momerr:t0err:nda_o_na:nza_o_na:nma_o_nm",    // assume all floats
  "su2020/data/trq_mva/MLP_weights_1170.xml",                                // location of XML weights
  0.20                                                                       // default default cut value
};

//------------------------------------------------------------------------------
// PID MVA's
//-----------------------------------------------------------------------------
mva_data::data_t  pid_mva_ele00s61b0_1000 = {
  "pid_ele00s61b0_1000",                                                     // name
  "ecl/p:ncr:seedfr:ele_dt:ele_dz:ele_dr:ele_path",                          // assume all floats
  "su2020/data/pid_mva/pid_MLP_weights_1000.xml",                            // location of XML weights
  0.50                                                                       // default cut value
};

//-----------------------------------------------------------------------------
// PAR ANN training sets
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

  if      (ds == "FELE2S51B1") {
    if      (TrainingCode ==   60) fData = trq_mva_fele2s51b1_0060;  // PAR, logfcons, uniform weight, dPf > 0.6
    else if (TrainingCode ==   70) fData = trq_mva_fele2s51b1_0070;  // PAR, logfcons, uniform weight, dPf > 0.7
    else if (TrainingCode == 1055) fData = trq_mva_fele2s51b1_1055;  // DAR, logfcons, w= 1/dPf      , dPf > 0.5
    else if (TrainingCode == 1060) fData = trq_mva_fele2s51b1_1060;  // DAR, logfcons, uniform weight, dPf > 0.6
    else if (TrainingCode == 1070) fData = trq_mva_fele2s51b1_1070;  // DAR, logfcons, uniform weight, dPf > 0.7
    else                           { error = 1; }
  }
  else if (ds == "FPOS2S51B1") {
    if      (TrainingCode ==  150) fData = trq_mva_fpos2s51b1_0150;  // PAR, logfcons, uniform weight, dPf > 0.5
    else if (TrainingCode ==  170) fData = trq_mva_fpos2s51b1_0170;  // PAR, logfcons, uniform weight, dPf > 0.7
    else if (TrainingCode == 1150) fData = trq_mva_fpos2s51b1_1150;  // DAR, logfcons, uniform weight, dPf > 0.5
    else if (TrainingCode == 1170) fData = trq_mva_fpos2s51b1_1170;  // DAR, logfcons, uniform weight, dPf > 0.7
    else                           { error = 1; }
  }
  else if (ds == "ELE00S61B0") {
    if      (TrainingCode == 1000) fData = pid_mva_ele00s61b0_1000;  // DAR PID MVA 
    else                           { error = 1; }
  }
  else if (ds == "ANDY") {
    if      (TrainingCode ==    0) fData = trq_mva_andy_0000;        // PAR TRQ MVA (negative)
    if      (TrainingCode ==  100) fData = trq_mva_andy_0100;        // PAR TRQ MVA (positive)
    else                           { error = 1; }
  }
  else                             { error = 1; }

  if (error == 0) {
    fMva = new TMVA::Reader("!Color:!Silent");

    TString s(VarNames());

    TObjArray* x = s.Tokenize(":");

    int n = x->GetEntries();

    for (int i=0; i<n; i++) {
      TObjString* name = (TObjString*) x->At(i);
      fMva->AddVariable(name->String().Data(),&fVar[i]);
    }

    fMva->BookMVA(Name(),XmlWeightsFile());
  }
  else {
    printf(" >>> ERROR in mva_data::mva_data(const char*,int): algorithm : %s training code: %5i. BAIL OUT\n",
	   Dataset,TrainingCode);
  }
}

//-----------------------------------------------------------------------------
mva_data::mva_data() {
}

//-----------------------------------------------------------------------------
mva_data::~mva_data() {
}

