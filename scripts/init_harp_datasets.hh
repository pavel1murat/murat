//
#ifndef __murat_scripts_init_harp_datasets__
#define __murat_scripts_init_harp_datasets__

#include "murat/scripts/dataset.hh"
//-----------------------------------------------------------------------------
// "simulation" of the HARP "setup"
//-----------------------------------------------------------------------------
dataset_t  d_harp_622_0001_ftfp_bert;
dataset_t  d_harp_622_0002_missing;
dataset_t  d_harp_622_0003_ftfp_bert_atl;
dataset_t  d_harp_622_0004_ftfp_bert_hp;
dataset_t  d_harp_622_0005_ftfp_inclxx;
dataset_t  d_harp_622_0006_ftfp_inclxx_hp;
dataset_t  d_harp_622_0007_ftf_bic;
dataset_t  d_harp_622_0008_missing;
dataset_t  d_harp_622_0009_qbbc;
dataset_t  d_harp_622_0010_qgsp_bert;
dataset_t  d_harp_622_0011_qgsp_bert_hp;
dataset_t  d_harp_622_0012_qgsp_bic;
dataset_t  d_harp_622_0013_qgsp_bic_hp;
dataset_t  d_harp_622_0014_missing;
dataset_t  d_harp_622_0015_qgsp_ftfp_bert;
dataset_t  d_harp_622_0016_qgsp_inclxx;
dataset_t  d_harp_622_0017_qgsp_inclxx_hp;
dataset_t  d_harp_622_0018_qgs_bic;
dataset_t  d_harp_622_0019_shielding;
dataset_t  d_harp_622_0020_missing;
dataset_t  d_harp_622_0021_shieldingm;
dataset_t  d_harp_622_0022_nubeam;
//-----------------------------------------------------------------------------
void init_harp_datasets() {
  const char* HistDir    = "/projects/hist/mu2e/v6_1_4";
  
  d_harp_622_0001_ftfp_bert.fName = "harp_622_0001_ftfp_bert";
  d_harp_622_0001_ftfp_bert.fFn    = Form("%s/harp_622_0001.ftfp_bert_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0001_ftfp_bert.fLabel = "p+Ta, 8 GeV/c, FTFP_BERT";
  d_harp_622_0001_ftfp_bert.fNPOT  = 100000.;

  d_harp_622_0002_missing.fName = "harp_622_0001_ftfp_bert";
  d_harp_622_0002_missing.fFn    = Form("%s/harp_622_0001.ftfp_bert_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0002_missing.fLabel = "p+Ta, 8 GeV/c, FTFP_BERT";
  d_harp_622_0002_missing.fNPOT  = -1;

  d_harp_622_0003_ftfp_bert_atl.fName = "harp_622_0003_ftfp_bert_atl";
  d_harp_622_0003_ftfp_bert_atl.fFn    = Form("%s/harp_622_0003.ftfp_bert_atl_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0003_ftfp_bert_atl.fLabel = "p+Ta, 8 GeV/c, FTFP_BERT_ATL";
  d_harp_622_0003_ftfp_bert_atl.fNPOT  = 100000.;

  d_harp_622_0004_ftfp_bert_hp.fName = "harp_622_0004_ftfp_bert_hp";
  d_harp_622_0004_ftfp_bert_hp.fFn    = Form("%s/harp_622_0004.ftfp_bert_hp_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0004_ftfp_bert_hp.fLabel = "p+Ta, 8 GeV/c, FTFP_BERT_HP";
  d_harp_622_0004_ftfp_bert_hp.fNPOT  = 100000.;

  d_harp_622_0005_ftfp_inclxx.fName = "harp_622_0005_ftfp_inclxx";
  d_harp_622_0005_ftfp_inclxx.fFn    = Form("%s/harp_622_0005.ftfp_inclxx_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0005_ftfp_inclxx.fLabel = "p+Ta, 8 GeV/c, FTFP_INCLXX";
  d_harp_622_0005_ftfp_inclxx.fNPOT  = 100000.;

  d_harp_622_0006_ftfp_inclxx_hp.fName = "harp_622_0006_ftfp_inclxx_hp";
  d_harp_622_0006_ftfp_inclxx_hp.fFn    = Form("%s/harp_622_0006.ftfp_inclxx_hp_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0006_ftfp_inclxx_hp.fLabel = "p+Ta, 8 GeV/c, FTFP_INCLXX_HP";
  d_harp_622_0006_ftfp_inclxx_hp.fNPOT  = 100000.;

  d_harp_622_0007_ftf_bic.fName = "harp_622_0007_ftf_bic";
  d_harp_622_0007_ftf_bic.fFn    = Form("%s/harp_622_0007.ftf_bic_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0007_ftf_bic.fLabel = "p+Ta, 8 GeV/c, FTF_BIC";
  d_harp_622_0007_ftf_bic.fNPOT  = 100000.;

  d_harp_622_0008_missing.fName = "harp_622_0008_missing";
  d_harp_622_0008_missing.fFn    = Form("%s/harp_622_0008.ftfp_bert_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0008_missing.fLabel = "p+Ta, 8 GeV/c, MISSING";
  d_harp_622_0008_missing.fNPOT  = -1;

  d_harp_622_0009_qbbc.fName = "harp_622_0009_qbbc";
  d_harp_622_0009_qbbc.fFn    = Form("%s/harp_622_0009.qbbc_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0009_qbbc.fLabel = "p+Ta, 8 GeV/c, QBBC";
  d_harp_622_0009_qbbc.fNPOT  = 100000.;

  d_harp_622_0010_qgsp_bert.fName = "harp_622_0010_qgsp_bert";
  d_harp_622_0010_qgsp_bert.fFn    = Form("%s/harp_622_0010.qgsp_bert_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0010_qgsp_bert.fLabel = "p+Ta, 8 GeV/c, QGSP_BERT";
  d_harp_622_0010_qgsp_bert.fNPOT  = 100000.;

  d_harp_622_0011_qgsp_bert_hp.fName = "harp_622_0011_qgsp_bert";
  d_harp_622_0011_qgsp_bert_hp.fFn    = Form("%s/harp_622_0011.qgsp_bert_hp_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0011_qgsp_bert_hp.fLabel = "p+Ta, 8 GeV/c, QGSP_BERT_HP";
  d_harp_622_0011_qgsp_bert_hp.fNPOT  = 100000.;

  d_harp_622_0012_qgsp_bic.fName = "harp_622_0012_qgsp_bic";
  d_harp_622_0012_qgsp_bic.fFn    = Form("%s/harp_622_0012.qgsp_bic_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0012_qgsp_bic.fLabel = "p+Ta, 8 GeV/c, QGSP_BIC";
  d_harp_622_0012_qgsp_bic.fNPOT  = 100000.;

  d_harp_622_0013_qgsp_bic_hp.fName = "harp_622_0013_qgsp_bic";
  d_harp_622_0013_qgsp_bic_hp.fFn    = Form("%s/harp_622_0013.qgsp_bic_hp_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0013_qgsp_bic_hp.fLabel = "p+Ta, 8 GeV/c, QGSP_BIC_HP";
  d_harp_622_0013_qgsp_bic_hp.fNPOT  = 100000.;

  d_harp_622_0014_missing.fName = "harp_622_0014_missing";
  d_harp_622_0014_missing.fFn    = Form("%s/harp_622_0014.ftfp_bert_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0014_missing.fLabel = "p+Ta, 8 GeV/c, MISSING";
  d_harp_622_0014_missing.fNPOT  = -1;

  d_harp_622_0015_qgsp_ftfp_bert.fName = "harp_622_0015_qgsp_ftfp_bert";
  d_harp_622_0015_qgsp_ftfp_bert.fFn    = Form("%s/harp_622_0015.qgsp_ftfp_bert_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0015_qgsp_ftfp_bert.fLabel = "p+Ta, 8 GeV/c, QGSP_FTFP_BERT";
  d_harp_622_0015_qgsp_ftfp_bert.fNPOT  = 100000.;

  d_harp_622_0016_qgsp_inclxx.fName = "harp_622_0016_qgsp_inclxx";
  d_harp_622_0016_qgsp_inclxx.fFn    = Form("%s/harp_622_0016.qgsp_inclxx_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0016_qgsp_inclxx.fLabel = "p+Ta, 8 GeV/c, QGSP_INCLXX";
  d_harp_622_0016_qgsp_inclxx.fNPOT  = 100000.;

  d_harp_622_0017_qgsp_inclxx_hp.fName = "harp_622_0017_qgsp_inclxx_hp";
  d_harp_622_0017_qgsp_inclxx_hp.fFn    = Form("%s/harp_622_0017.qgsp_inclxx_hp_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0017_qgsp_inclxx_hp.fLabel = "p+Ta, 8 GeV/c, QGSP_INCLXX_HP";
  d_harp_622_0017_qgsp_inclxx_hp.fNPOT  = 100000.;

  d_harp_622_0018_qgs_bic.fName = "harp_622_0018_qgs_bic";
  d_harp_622_0018_qgs_bic.fFn    = Form("%s/harp_622_0018.qgs_bic_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0018_qgs_bic.fLabel = "p+Ta, 8 GeV/c, QGS_BIC";
  d_harp_622_0018_qgs_bic.fNPOT  = 100000.;

  d_harp_622_0019_shielding.fName  = "harp_622_0019_shielding";
  d_harp_622_0019_shielding.fFn    = Form("%s/harp_622_0019.shielding_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0019_shielding.fLabel = "p+Ta, 8 GeV/c, Shielding";
  d_harp_622_0019_shielding.fNPOT  = 100000.;

  d_harp_622_0020_missing.fName = "harp_622_0020_missing";
  d_harp_622_0020_missing.fFn    = Form("%s/harp_622_0020.ftfp_bert_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0020_missing.fLabel = "p+Ta, 8 GeV/c, MISSING";
  d_harp_622_0020_missing.fNPOT  = -1;

  d_harp_622_0021_shieldingm.fName = "harp_622_0021_shieldingm";
  d_harp_622_0021_shieldingm.fFn    = Form("%s/harp_622_0021.shieldingm_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0021_shieldingm.fLabel = "p+Ta, 8 GeV/c, ShieldingM";
  d_harp_622_0021_shieldingm.fNPOT  = 100000.;

  d_harp_622_0022_nubeam.fName  = "harp_622_0022_nubeam";
  d_harp_622_0022_nubeam.fFn    = Form("%s/harp_622_0022.nubeam_g4study2.g4val_ana.hist",HistDir);
  d_harp_622_0022_nubeam.fLabel = "p+Ta, 8 GeV/c, NuBeam";
  d_harp_622_0022_nubeam.fNPOT  = 100000.;
}

#endif
