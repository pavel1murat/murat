//
#ifndef __murat_scripts_init_ts_warm_bore_datasets__
#define __murat_scripts_init_ts_warm_bore_datasets__

#include "murat/scripts/dataset.hh"

namespace ts_warm_bore {
//-----------------------------------------------------------------------------
// "simulation" of the TS_WARM_BORE "setup"
//-----------------------------------------------------------------------------
  dataset_t  d_711_1000_s1_mubeam;   // 4.96M events, nominal geometry
  dataset_t  d_711_1000_s2_mubeam;   // 4.96M events, nominal geometry
  dataset_t  d_711_1000_s3_tgtstops; // 4.96M events, nominal geometry
  dataset_t  d_711_1000_s3_ootstops; // 4.96M events, nominal geometry

  dataset_t  d_711_1001_s1_mubeam;   // 5M events, worst case geometry
  dataset_t  d_711_1001_s2_mubeam;   // 5M events, worst case geometry
  dataset_t  d_711_1001_s3_tgtstops; // 5M events, worst case geometry
  dataset_t  d_711_1001_s3_ootstops; // 5M events, worst case geometry

  dataset_t  d_711_1005_s1_mubeam;   // 4.98M events, mitigated geometry
  dataset_t  d_711_1005_s2_mubeam;   // 4.98M events, mitigated geometry
  dataset_t  d_711_1005_s3_tgtstops; // 4.98M events, mitigated geometry
  dataset_t  d_711_1005_s3_ootstops; // 4.98M events, mitigated geometry

  dataset_t  d_711_1006_s1_mubeam;   // 4.98M events, nominal geometry, new target with fins
  dataset_t  d_711_1006_s2_mubeam;   // 4.98M events, nominal geometry, new target with fins
  dataset_t  d_711_1006_s3_tgtstops; // 4.98M events, nominal geometry, new target with fins
  dataset_t  d_711_1006_s3_ootstops; // 4.98M events, nominal geometry, new target with fins

  dataset_t  d_711_1007_s1_mubeam;   // 4.98M events, nominal geometry, target with R=5 mm
  dataset_t  d_711_1007_s2_mubeam;   // 4.98M events, nominal geometry, target with R=5 mm
  dataset_t  d_711_1007_s3_tgtstops; // 4.98M events, nominal geometry, target with R=5 mm
  dataset_t  d_711_1007_s3_ootstops; // 4.98M events, nominal geometry, target with R=5 mm
//-----------------------------------------------------------------------------
  void init_datasets() {
    const char* HistDir    = "/projects/hist/mu2e/ts_warm_bore";
  
    d_711_1000_s1_mubeam.fName  = "ts_warm_bore.711_1000_s1_mubeam";
    d_711_1000_s1_mubeam.fFn    = Form("%s/ts_warm_bore.711_1000_s1_mubeam_stn.spmc_ana.hist",HistDir);
    d_711_1000_s1_mubeam.fLabel = "default TS geometry";
    d_711_1000_s1_mubeam.fNPOT  = 4960000.;

    d_711_1000_s2_mubeam.fName  = "ts_warm_bore.711_1000_s2_mubeam";
    d_711_1000_s2_mubeam.fFn    = Form("%s/ts_warm_bore.711_1000_s2_mubeam_stn.spmc_ana.hist",HistDir);
    d_711_1000_s2_mubeam.fLabel = "default TS geometry";
    d_711_1000_s2_mubeam.fNPOT  = 4960000.;

    d_711_1000_s3_tgtstops.fName  = "ts_warm_bore.711_1000_s3_tgtstops";
    d_711_1000_s3_tgtstops.fFn    = Form("%s/ts_warm_bore.711_1000_s3_tgtstops_stn.mustop_ana.hist",HistDir);
    d_711_1000_s3_tgtstops.fLabel = "default TS geometry";
    d_711_1000_s3_tgtstops.fNPOT  = 4960000.;

    d_711_1000_s3_ootstops.fName  = "ts_warm_bore.711_1000_s3_ootstops";
    d_711_1000_s3_ootstops.fFn    = Form("%s/ts_warm_bore.711_1000_s3_ootstops_stn.mustop_ana.hist",HistDir);
    d_711_1000_s3_ootstops.fLabel = "default TS geometry";
    d_711_1000_s3_ootstops.fNPOT  = 4960000.;

    d_711_1001_s1_mubeam.fName  = "ts_warm_bore.711_1001_s1_mubeam";
    d_711_1001_s1_mubeam.fFn    = Form("%s/ts_warm_bore.711_1001_s1_mubeam_stn.spmc_ana.hist",HistDir);
    d_711_1001_s1_mubeam.fLabel = "misaligned TS warm bores";
    d_711_1001_s1_mubeam.fNPOT  = 5000000.;

    d_711_1001_s2_mubeam.fName  = "ts_warm_bore.711_1001_s2_mubeam";
    d_711_1001_s2_mubeam.fFn    = Form("%s/ts_warm_bore.711_1001_s2_mubeam_stn.spmc_ana.hist",HistDir);
    d_711_1001_s2_mubeam.fLabel = "misaligned TS warm bores";
    d_711_1001_s2_mubeam.fNPOT  = 5000000.;

    d_711_1001_s3_tgtstops.fName  = "ts_warm_bore.711_1001_s3_tgtstops";
    d_711_1001_s3_tgtstops.fFn    = Form("%s/ts_warm_bore.711_1001_s3_tgtstops_stn.mustop_ana.hist",HistDir);
    d_711_1001_s3_tgtstops.fLabel = "misaligned TS warm bores";
    d_711_1001_s3_tgtstops.fNPOT  = 5000000.;

    d_711_1001_s3_ootstops.fName  = "ts_warm_bore.711_1001_s3_ootstops";
    d_711_1001_s3_ootstops.fFn    = Form("%s/ts_warm_bore.711_1001_s3_ootstops_stn.mustop_ana.hist",HistDir);
    d_711_1001_s3_ootstops.fLabel = "misaligned TS warm bores";
    d_711_1001_s3_ootstops.fNPOT  = 5000000.;

    d_711_1005_s1_mubeam.fName    = "ts_warm_bore.711_1005_s1_mubeam";
    d_711_1005_s1_mubeam.fFn      = Form("%s/ts_warm_bore.711_1005_s1_mubeam_stn.spmc_ana.hist",HistDir);
    d_711_1005_s1_mubeam.fLabel   = "adjusted TS warm bores";
    d_711_1005_s1_mubeam.fNPOT    = 4980000.;

    d_711_1005_s2_mubeam.fName    = "ts_warm_bore.711_1005_s2_mubeam";
    d_711_1005_s2_mubeam.fFn      = Form("%s/ts_warm_bore.711_1001_s2_mubeam_stn.spmc_ana.hist",HistDir);
    d_711_1005_s2_mubeam.fLabel   = "adjusted TS warm bores";
    d_711_1005_s2_mubeam.fNPOT    = 4980000.;

    d_711_1005_s3_tgtstops.fName  = "ts_warm_bore.711_1005_s3_tgtstops";
    d_711_1005_s3_tgtstops.fFn    = Form("%s/ts_warm_bore.711_1005_s3_tgtstops_stn.mustop_ana.hist",HistDir);
    d_711_1005_s3_tgtstops.fLabel = "adjusted TS warm bores";
    d_711_1005_s3_tgtstops.fNPOT  = 4980000.;

    d_711_1005_s3_ootstops.fName  = "ts_warm_bore.711_1005_s3_ootstops";
    d_711_1005_s3_ootstops.fFn    = Form("%s/ts_warm_bore.711_1005_s3_ootstops_stn.mustop_ana.hist",HistDir);
    d_711_1005_s3_ootstops.fLabel = "adjusted TS warm bores";
    d_711_1005_s3_ootstops.fNPOT  = 4980000.;

    d_711_1006_s1_mubeam.fName    = "ts_warm_bore.711_1006_s1_mubeam";
    d_711_1006_s1_mubeam.fFn      = Form("%s/ts_warm_bore.711_1006_s1_mubeam_stn.spmc_ana.hist",HistDir);
    d_711_1006_s1_mubeam.fLabel   = "target with fins";
    d_711_1006_s1_mubeam.fNPOT    = 4960000.;

    d_711_1006_s2_mubeam.fName    = "ts_warm_bore.711_1006_s2_mubeam";
    d_711_1006_s2_mubeam.fFn      = Form("%s/ts_warm_bore.711_1001_s2_mubeam_stn.spmc_ana.hist",HistDir);
    d_711_1006_s2_mubeam.fLabel   = "target with fins";
    d_711_1006_s2_mubeam.fNPOT    = 4960000.;

    d_711_1006_s3_tgtstops.fName  = "ts_warm_bore.711_1006_s3_tgtstops";
    d_711_1006_s3_tgtstops.fFn    = Form("%s/ts_warm_bore.711_1006_s3_tgtstops_stn.mustop_ana.hist",HistDir);
    d_711_1006_s3_tgtstops.fLabel = "target with fins";
    d_711_1006_s3_tgtstops.fNPOT  = 4960000.;

    d_711_1006_s3_ootstops.fName  = "ts_warm_bore.711_1006_s3_ootstops";
    d_711_1006_s3_ootstops.fFn    = Form("%s/ts_warm_bore.711_1006_s3_ootstops_stn.mustop_ana.hist",HistDir);
    d_711_1006_s3_ootstops.fLabel = "target with fins";
    d_711_1006_s3_ootstops.fNPOT  = 4960000.;

    d_711_1007_s1_mubeam.fName    = "ts_warm_bore.711_1007_s1_mubeam";
    d_711_1007_s1_mubeam.fFn      = Form("%s/ts_warm_bore.711_1007_s1_mubeam_stn.spmc_ana.hist",HistDir);
    d_711_1007_s1_mubeam.fLabel   = "target R=5 mm";
    d_711_1007_s1_mubeam.fNPOT    = 4960000.;

    d_711_1007_s2_mubeam.fName    = "ts_warm_bore.711_1007_s2_mubeam";
    d_711_1007_s2_mubeam.fFn      = Form("%s/ts_warm_bore.711_1001_s2_mubeam_stn.spmc_ana.hist",HistDir);
    d_711_1007_s2_mubeam.fLabel   = "target R=5 mm";
    d_711_1007_s2_mubeam.fNPOT    = 4960000.;

    d_711_1007_s3_tgtstops.fName  = "ts_warm_bore.711_1007_s3_tgtstops";
    d_711_1007_s3_tgtstops.fFn    = Form("%s/ts_warm_bore.711_1007_s3_tgtstops_stn.mustop_ana.hist",HistDir);
    d_711_1007_s3_tgtstops.fLabel = "target R=5 mm";
    d_711_1007_s3_tgtstops.fNPOT  = 4960000.;

    d_711_1007_s3_ootstops.fName  = "ts_warm_bore.711_1007_s3_ootstops";
    d_711_1007_s3_ootstops.fFn    = Form("%s/ts_warm_bore.711_1007_s3_ootstops_stn.mustop_ana.hist",HistDir);
    d_711_1007_s3_ootstops.fLabel = "target R=5 mm";
    d_711_1007_s3_ootstops.fNPOT  = 4960000.;
  }
}
#endif
