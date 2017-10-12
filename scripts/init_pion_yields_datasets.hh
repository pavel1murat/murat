//
#ifndef __murat_scripts_init_pion_yields_datasets__
#define __murat_scripts_init_pion_yields_datasets__

#include "murat/scripts/dataset.hh"
//-----------------------------------------------------------------------------
// datasets
//-----------------------------------------------------------------------------
dataset_t  d_pion_yields_622_0003_ftfp_bert_atl;
dataset_t  d_pion_yields_622_0004_ftfp_bert_hp;
dataset_t  d_pion_yields_622_0010_qgsp_bert;
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void init_pion_yields_datasets() {
  const char* HistDir    = "/projects/hist/mu2e/v6_1_4";

  d_pion_yields_622_0003_ftfp_bert_atl.fName = "622_0003_ftfp_bert_atl";
  d_pion_yields_622_0003_ftfp_bert_atl.fFn    = Form("%s/pion_yields.622_0003_ftfp_bert_atl.g4s2_mubeam.bflash_ana.hist",HistDir);
  d_pion_yields_622_0003_ftfp_bert_atl.fLabel = "FTFP_BERT_ATL";
  d_pion_yields_622_0003_ftfp_bert_atl.fNPOT  = 5000.;

  d_pion_yields_622_0004_ftfp_bert_hp.fName = "622_0004_ftfp_bert_hp";
  d_pion_yields_622_0004_ftfp_bert_hp.fFn    = Form("%s/pion_yields.622_0004_ftfp_bert_hp.g4s2_mubeam.bflash_ana.hist",HistDir);
  d_pion_yields_622_0004_ftfp_bert_hp.fLabel = "FTFP_BERT_HP";
  d_pion_yields_622_0004_ftfp_bert_hp.fNPOT  = 5000.;

  d_pion_yields_622_0010_qgsp_bert.fName = "622_0010_qgsp_bert";
  d_pion_yields_622_0010_qgsp_bert.fFn    = Form("%s/pion_yields.622_0010_qgsp_bert.g4s2_mubeam.bflash_ana.hist",HistDir);
  d_pion_yields_622_0010_qgsp_bert.fLabel = "Shift20dd";
  d_pion_yields_622_0010_qgsp_bert.fNPOT  = 5000.;
}

#endif
