//
#ifndef __murat_scripts_datasets_hh__
#define __murat_scripts_datasets_hh__

  // const char* e41s5721_track_ana  = "~/hist/mu2e/v5_7_0/e41s5721.track_ana.hist" ; // HEE+BGR , matcorr in CalPatRec
  // const char* e41s5721_track_comp = "~/hist/mu2e/v5_7_0/e41s5721.track_comp.hist"; // HEE+BGR , matcorr in CalPatRec

  // const char* e42s5721_track_ana  = "~/hist/mu2e/v5_7_0/e42s5721.track_ana.hist" ; // CE+BGR  , matcorr in CalPatRec
  // const char* e42s5721_track_comp = "~/hist/mu2e/v5_7_0/e42s5721.track_comp.hist"; // CE+BGR  , matcorr in CalPatRec

  struct dataset_t {
    const char *name;
    const char *fn_track_comp; // hist file name for "track_comp" job
    const char *fn_track_ana;  // hist file name for "track_ana"  job
    const char *label;         // dataset label to be printed, so far, CE" or "HEE"
    int          type;         //  0:  data, 1:MC
    int          ngen;	       //  number of generated events
  };

//-----------------------------------------------------------------------------
// v600 v0 datasets 
//-----------------------------------------------------------------------------
dataset_t g50s6000 = { "g50s6000",                       // DIO LL B=0.5 T
		       "~/hist/mu2e/v6_0_0/g50s6000.track_comp_use_mva.hist",
		       "~/hist/mu2e/v6_0_0/g50s6000.track_ana.hist",   // missing 
		       "DIO-LL-B=0.5T",
		       1,
		       10000000};

dataset_t d50s6000 = { "d50s6000",                       // DIO LO B=0.5 T
		       "~/hist/mu2e/v6_0_0/d50s6000.track_comp_use_mva.hist",
		       "~/hist/mu2e/v6_0_0/d50s6000.track_ana.hist",   // missing 
		       "DIO-LO-B=0.5T",
		       1,
		       10000000};

dataset_t e31s6000 = { "e31s6000",                       // CE LL + MIX-CD3 x1
		       "~/hist/mu2e/v6_0_0/e31s6000.track_comp_use_mva.hist",
		       "~/hist/mu2e/v6_0_0/e12s5740.track_ana.hist",   // missing 
		       "CE-photos+mixcd3.x1",
		       1,
		       475000};

dataset_t e30s6000 = { "e30s6000",                       // CE LL
		       "~/hist/mu2e/v6_0_0/e30s6000.track_comp_use_mva.hist",
		       "~/hist/mu2e/v6_0_0/e30s6000.track_ana.hist",   // missing 
		       "CE-photos",
		       1,
		       900000};

dataset_t e21s6000 = { "e21s6000",
		       "~/hist/mu2e/v6_0_0/e21s6000.track_comp_use_mva.hist",
		       "~/hist/mu2e/v6_0_0/e21s6000.track_ana.hist",
		       "CE+cd3-detmix-cut.v566b",
		       1,
		       900000};

dataset_t e00s6000 = { "e00s6000",                       // E105 .. TCM
		       "~/hist/mu2e/v6_0_0/e00s6000.track_comp_use_mva.hist",
		       "~/hist/mu2e/v6_0_0/e00s6000.track_ana.hist",   // missing 
		       "E105-v600",
		       1,
		       1000000};
//-----------------------------------------------------------------------------
// 574 v0 datasets 
//-----------------------------------------------------------------------------
dataset_t e12s5740 = { "e12s5740",
		       "~/hist/mu2e/v5_7_0/e12s5740.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/e12s5740.track_ana.hist",
		       "HEE+cd3-detmix-cut.v566b.x2",
		       1,
		       1310000};
//-----------------------------------------------------------------------------
// 573 v1 datasets (include Giani's pattern recognition cleanup)
//-----------------------------------------------------------------------------
dataset_t m02s5731 = { "m02s5731",
		       "~/hist/mu2e/v5_7_0/mo2s5731.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/m02s5731.track_ana.hist",
		       "M105+cd3-detmix-cut.v566b.x2",
		       1,
		       460000};

dataset_t e22s5731 = { "e22s5731",
		       "~/hist/mu2e/v5_7_0/e22s5731.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/e22s5731.track_ana.hist",
		       "CE+cd3-detmix-cut.v566b.x2",
		       1,
		       860000};

dataset_t e21s5731 = { "e21s5731",
		       "~/hist/mu2e/v5_7_0/e21s5731.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/e21s5731.track_ana.hist",
		       "CE+cd3-detmix-cut.v566b",
		       1,
		       890000};

dataset_t e11s5731 = { "e11s5731",
		       "~/hist/mu2e/v5_7_0/e11s5731.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/e11s5731.track_ana.hist",
		       "HEE+cd3-detmix-cut.v566b",
		       1,
		       2290000};
//-----------------------------------------------------------------------------
// 573 v0 datasets
//-----------------------------------------------------------------------------
dataset_t e00s5730 = { "e00s5730",
		       "~/hist/mu2e/v5_7_0/e00s5730.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/e00s5730.track_ana.hist",
		       "E105",
		       1,
		       5000000};

dataset_t e40s5720 = { "e40s5720",
		       "~/hist/mu2e/v5_7_0/e40s5720.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/e40s5720.track_ana.hist",
		       "E105",
			 1,
		       2500000};

dataset_t e41s5721 = { "e41s5721",
		       "~/hist/mu2e/v5_7_0/e41s5721.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/e41s5721.track_ana.hist",
		       "HEE+cd3-detmix-cut.v566b",
		       1,
		       5110000};

dataset_t e42s5721 = { "e42s5721",
		       "~/hist/mu2e/v5_7_0/e42s5721.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/e42s5721.track_ana.hist",
		       "CE+cd3-detmix-cut.v566b",
		       1,
		       880000};

dataset_t m40s5720 = { "m40s5720",
		       "~/hist/mu2e/v5_7_0/m40s5720.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/m40s5720.track_ana.hist",
		       "M105",
		       1,
		       2500000};

dataset_t m40s5721 = { "m40s5721",
		       "~/hist/mu2e/v5_7_0/m40s5721.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/m40s5721.track_ana.hist",
		       "M105+cd3-detmix-cut.v566b",
		       1,
		       1990000};

#endif
