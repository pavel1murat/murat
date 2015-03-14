//
#include "limits/TMu2eLimits.hh"
// #include "murat/ana/TAnaUtils.hh"

//-----------------------------------------------------------------------------
int TMu2eLimits::InitMu2eChannel(TMu2eLimits* Tzz, csm_channel_model* Channel) {
// in the 4L mass channel add 3 histograms: M(ZeeZee), M(ZeeZmm) and M(ZmmZmm)xs

  TH1    *h11(0);
  int    np;
//-----------------------------------------------------------------------------
// 
// everything else should be much lower 
//-----------------------------------------------------------------------------
  const char        *name;
  channel_data_t*   ch_data;
  const char        *zz_pnames      [10];
  double            zz_par_sf_lo    [10];
  double            zz_par_sf_hi    [10];
  TH1*              zz_hist_shape_lo[10];
  double            zz_nsig_shape_lo[10];
  TH1*              zz_hist_shape_hi[10];
  double            zz_nsig_shape_hi[10];
  double            sf;
  aprocess*         proc;
  analysis*         ana;

  np                  = 1;           // 1;
  zz_pnames[0]        = "lumi";
  zz_par_sf_lo[0]     = -0.1;
  zz_par_sf_hi[0]     =  0.1;
  zz_hist_shape_lo[0] =  NULL;
  zz_hist_shape_hi[0] =  NULL;
  zz_nsig_shape_lo[0] = 1.;
  zz_nsig_shape_hi[0] = 1.;
//-----------------------------------------------------------------------------
// background histogram for mi ---> e conversion
//-----------------------------------------------------------------------------
  name     = Channel->GetName();
  ch_data  = Tzz->GetChannelData(name);
  ana      = ch_data->ana;

  int nproc = ana->GetNProcesses(0);

  for (int i=0; i<nproc; i++) {
    proc = ana->GetProcess(i,0);
					// skip MC signal processes
    if (proc->GetMcFlag() != 3) {
//-----------------------------------------------------------------------------
// this call returns 'Poisson-like' histogram in 'h11' and a scale factor 
// in 'sf'
//-----------------------------------------------------------------------------
      proc->get_h1_unorm(ch_data->module_name,
			 ch_data->hist_set,
			 ch_data->bin,
			 ch_data->hist_name,
			 ch_data->rebin,
			 h11,
			 &sf);
//-----------------------------------------------------------------------------
// optionally, assuming 1D histograms, define a used range 
//-----------------------------------------------------------------------------
      if (ch_data->xmin < ch_data->xmax) {
	int nb = h11->GetNbinsX();
	for (int i=1; i<=nb; i++) {
	  double x = h11->GetBinCenter(i);
	  if ((x < ch_data->xmin) || (x > ch_data->xmax)) {
	    h11->SetBinContent(i,0.);
	    h11->SetBinError  (i,0.);
	  }
	}
      }
//-----------------------------------------------------------------------------
// a smeared and weighted histogram is used as input, also use errors stored 
// in the histogram
//-----------------------------------------------------------------------------
      Channel->add_template(h11,
			    sf,
			    np,
			    zz_pnames,
			    zz_par_sf_lo,
			    zz_par_sf_hi,
			    zz_hist_shape_lo,
			    zz_nsig_shape_lo,
			    zz_hist_shape_hi,
			    zz_nsig_shape_hi,
			    // 1,                  // use poisson errors 
			    2,                  // use errors stored in the histogram
			    0);                 // this is background
    }
  }

  return 0;
}

