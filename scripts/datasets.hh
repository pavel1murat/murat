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
// 573 v1 datasets (include Giani's pattern recognition cleanup)
//-----------------------------------------------------------------------------
dataset_t e21s5731 = { "e21s5731",
		       "~/hist/mu2e/v5_7_0/e21s5731.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/e21s5731.track_ana.hist",
		       "CE+MIXCD3-cut-v2",
		       1,
		       960000};

dataset_t e21s5731 = { "e21s5731",
		       "~/hist/mu2e/v5_7_0/e11s5731.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/e11s5731.track_ana.hist",
		       "HEE+MIXCD3-cut-v2",
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
		       "HEE",
		       1,
		       5110000};

dataset_t e42s5721 = { "e42s5721",
		       "~/hist/mu2e/v5_7_0/e42s5721.track_comp.hist",
		       "~/hist/mu2e/v5_7_0/e42s5721.track_ana.hist",
		       "CE",
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
		       "M105+MIXCD3-cut-v2",
		       1,
		       1990000};

#endif
