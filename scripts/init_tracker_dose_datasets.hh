//
#ifndef __murat_scripts_init_tracker_dose_datasets__
#define __murat_scripts_init_tracker_dose_datasets__

#include "murat/scripts/dataset.hh"
//-----------------------------------------------------------------------------
// datasets
//-----------------------------------------------------------------------------
dataset_t  cd3_beam_cs1_mubeam_0506a_0000;
dataset_t  d_621_0001_cal_dose;                // former "tracker_dose_FLASH"

dataset_t  g4s2_622_0000_mubeam;
dataset_t  g4s2_622_0001_mubeam;

dataset_t  g4s3_622_0000_tgtstops; // output of stage 2 - step point MC's
dataset_t  g4s3_622_0000_ootstops; // output of stage 2 - step point MC's
dataset_t  g4s3_622_0001_tgtstops; // output of stage 2 - step point MC's
dataset_t  g4s3_622_0001_ootstops; // output of stage 2 - step point MC's

dataset_t  d_cd3_0001_g4s1_mubeam_stnmaker; // see catalogs for more information about the datasets
dataset_t  d_622_0009_g4s1_mubeam_stnmaker;
dataset_t  d_622_0010_g4s1_mubeam_stnmaker;

dataset_t  d_622_0005_cal_dose;
dataset_t  d_622_0006_cal_dose;
dataset_t  d_622_0007_cal_dose;
dataset_t  d_622_0008_cal_dose;
dataset_t  d_622_0011_cal_dose;
dataset_t  d_622_0012_cal_dose;
dataset_t  d_622_0013_cal_dose;
dataset_t  d_622_0014_cal_dose;

dataset_t  d_ts3_tooth_622_0001_g4s3_tgtstops;
//-----------------------------------------------------------------------------
void init_tracker_dose_datasets() {

  const char* HistDir    = "/projects/hist/mu2e/v6_1_4";
  
  cd3_beam_cs1_mubeam_0506a_0000.fName  = "g4s1_01_mubeam_0506a_0000";
  cd3_beam_cs1_mubeam_0506a_0000.fFn    = Form("%s/cd3_beam_cs1_mubeam_0506a_0000.stn_01.bflash_ana.hist",HistDir);
  cd3_beam_cs1_mubeam_0506a_0000.fLabel = "stage1";
  cd3_beam_cs1_mubeam_0506a_0000.fNPOT  = -1;
//-----------------------------------------------------------------------------
// v621 datasets (first stage of the rad dose studies)
//-----------------------------------------------------------------------------
  d_621_0001_cal_dose.fName  = "d_621_0001_cal_dose";
  d_621_0001_cal_dose.fFn    = Form("%s/tracker_dose_FLASH.hist",HistDir);
  d_621_0001_cal_dose.fLabel = "FLASH, default, Giani";
  d_621_0001_cal_dose.fNPOT  = 5.1e9;
//-----------------------------------------------------------------------------
// v622 datasets (second stage of the rad dose studies)
//-----------------------------------------------------------------------------
  g4s2_622_0000_mubeam.fName  = "g4s2_622_0000_mubeam";
  g4s2_622_0000_mubeam.fFn    = Form("%s/g4s2_622_0000_mubeam.spmc_ana.hist",HistDir);
  g4s2_622_0000_mubeam.fLabel = "80x40x40";
  g4s2_622_0000_mubeam.fNPOT  = 15600000;

  g4s2_622_0001_mubeam.fName  = "g4s2_622_0001_mubeam";
  g4s2_622_0001_mubeam.fFn    = Form("%s/g4s2_622_0001_mubeam.spmc_ana.hist",HistDir);
  g4s2_622_0001_mubeam.fLabel = "no tooth";
  g4s2_622_0001_mubeam.fNPOT  = 21200000;

  g4s3_622_0000_tgtstops.fName  = "g4s3_622_0000_tgtstops";
  g4s3_622_0000_tgtstops.fFn    = Form("%s/g4s3_622_0000_tgtstops.mustop_ana.hist",HistDir);
  g4s3_622_0000_tgtstops.fLabel = "80x40x40";

  g4s3_622_0000_ootstops.fName  = "g4s3_622_0000_ootstops";
  g4s3_622_0000_ootstops.fFn    = Form("%s/g4s3_622_0000_ootstops.mustop_ana.hist",HistDir);
  g4s3_622_0000_ootstops.fLabel = "80x40x40";

  g4s3_622_0001_tgtstops.fName  = "g4s3_622_0001_tgtstops";
  g4s3_622_0001_tgtstops.fFn    = Form("%s/g4s3_622_0001_tgtstops.mustop_ana.hist",HistDir);
  g4s3_622_0001_tgtstops.fLabel = "no tooth";

  g4s3_622_0001_ootstops.fName  = "g4s3_622_0001_ootstops";
  g4s3_622_0001_ootstops.fFn    = Form("%s/g4s3_622_0001_ootstops.mustop_ana.hist",HistDir);
  g4s3_622_0001_ootstops.fLabel = "no tooth";

  d_cd3_0001_g4s1_mubeam_stnmaker.fName  = "cd3_0001_g4s1_mubeam_stnmaker";
  d_cd3_0001_g4s1_mubeam_stnmaker.fFn    = Form("%s/cd3_0001.g4s1_mubeam_stnmaker.bflash_ana.hist",HistDir);
  d_cd3_0001_g4s1_mubeam_stnmaker.fLabel = "CD3";
  d_cd3_0001_g4s1_mubeam_stnmaker.fNPOT  = 1.e6;

  d_622_0009_g4s1_mubeam_stnmaker.fName  = "622_0009_g4s1_mubeam_stnmaker";
  d_622_0009_g4s1_mubeam_stnmaker.fFn    = Form("%s/622_0009.g4s1_mubeam_stnmaker.bflash_ana.hist",HistDir);
  d_622_0009_g4s1_mubeam_stnmaker.fLabel = "QGSP_BERT";
  d_622_0009_g4s1_mubeam_stnmaker.fNPOT  = 5.e3;

  d_622_0010_g4s1_mubeam_stnmaker.fName  = "622_0010_g4s1_mubeam_stnmaker";
  d_622_0010_g4s1_mubeam_stnmaker.fFn    = Form("%s/622_0010.g4s1_mubeam_stnmaker.bflash_ana.hist",HistDir);
  d_622_0010_g4s1_mubeam_stnmaker.fLabel = "ShieldingM";
  d_622_0010_g4s1_mubeam_stnmaker.fNPOT  = 5.e3;

  d_622_0010_g4s1_mubeam_stnmaker.fName  = "622_0010_g4s1_mubeam_stnmaker";
  d_622_0010_g4s1_mubeam_stnmaker.fFn    = Form("%s/622_0010.g4s1_mubeam_stnmaker.bflash_ana.hist",HistDir);
  d_622_0010_g4s1_mubeam_stnmaker.fLabel = "ShieldingM";
  d_622_0010_g4s1_mubeam_stnmaker.fNPOT  = 5.e3;

  d_622_0005_cal_dose.fName  = "622_0005_cal_dose";
  d_622_0005_cal_dose.fFn    = Form("%s/ts3_tooth.622_0005.cal_dose.hist",HistDir);
  d_622_0005_cal_dose.fLabel = "default_geometry";
  d_622_0005_cal_dose.fNPOT  = 61.2e6;

  d_622_0006_cal_dose.fName  = "622_0006_cal_dose";
  d_622_0006_cal_dose.fFn    = Form("%s/ts3_tooth.622_0006.cal_dose.hist",HistDir);
  d_622_0006_cal_dose.fLabel = "target hole R=21.6 mm";
  d_622_0006_cal_dose.fNPOT  = 99.6e6;

  d_622_0007_cal_dose.fName  = "622_0007_cal_dose";
  d_622_0007_cal_dose.fFn    = Form("%s/ts3_tooth.622_0007.cal_dose.hist",HistDir);
  d_622_0007_cal_dose.fLabel = "80x40x10 mm TS3 tooth";
  d_622_0007_cal_dose.fNPOT  = 100e6;

  d_622_0008_cal_dose.fName  = "622_0008_cal_dose";
  d_622_0008_cal_dose.fFn    = Form("%s/ts3_tooth.622_0008.cal_dose.hist",HistDir);
  d_622_0008_cal_dose.fLabel = "TS3 tooth + hole";
  d_622_0008_cal_dose.fNPOT  = 100e6;

  d_622_0011_cal_dose.fName  = "622_0011_cal_dose";
  d_622_0011_cal_dose.fFn    = Form("%s/ts3_tooth.622_0011.cal_dose.hist",HistDir);
  d_622_0011_cal_dose.fLabel = "Shift20uu";
  d_622_0011_cal_dose.fNPOT  = 6.0e6;

  d_622_0012_cal_dose.fName  = "622_0012_cal_dose";
  d_622_0012_cal_dose.fFn    = Form("%s/ts3_tooth.622_0012.cal_dose.hist",HistDir);
  d_622_0012_cal_dose.fLabel = "Shift20ud";
  d_622_0012_cal_dose.fNPOT  = 5.91e6;

  d_622_0013_cal_dose.fName  = "622_0013_cal_dose";
  d_622_0013_cal_dose.fFn    = Form("%s/ts3_tooth.622_0013.cal_dose.hist",HistDir);
  d_622_0013_cal_dose.fLabel = "Shift20du";
  d_622_0013_cal_dose.fNPOT  = 5.91e6;

  d_622_0014_cal_dose.fName  = "622_0014_cal_dose";
  d_622_0014_cal_dose.fFn    = Form("%s/ts3_tooth.622_0014.cal_dose.hist",HistDir);
  d_622_0014_cal_dose.fLabel = "Shift20dd";
  d_622_0014_cal_dose.fNPOT  = 5.97e6;

  d_ts3_tooth_622_0001_g4s3_tgtstops.fName  = "ts3_tooth_622_0001_g4s3_tgtstops";
  d_ts3_tooth_622_0001_g4s3_tgtstops.fFn    = Form("%s/ts3_tooth.622_0001.g4s3_tgtstops.bflash_ana_spmc.hist",HistDir);
  d_ts3_tooth_622_0001_g4s3_tgtstops.fLabel = "622_0001.g4s3_tgtstops.bflash_ana_spmc";
  d_ts3_tooth_622_0001_g4s3_tgtstops.fNPOT  = -1;
}


#endif
