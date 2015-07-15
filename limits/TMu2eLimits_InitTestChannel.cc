//
#include "limits/TMu2eLimits.hh"
// #include "murat/ana/TAnaUtils.hh"

//-----------------------------------------------------------------------------
int TMu2eLimits::InitTestChannel(TMu2eLimits* Tzz, csm_channel_model* Channel) {
// test channel

  TH1    *h11; //, *h12, *h15;  
  int    np;
//-----------------------------------------------------------------------------
// in 4L channel consider only SM ZZ process 
// everything else should be much lower 
//-----------------------------------------------------------------------------
  const char  *zz_pnames [10];
  double zz_par_sf_lo    [10];
  double zz_par_sf_hi    [10];
  TH1*   zz_hist_shape_lo[10];
  double zz_nsig_shape_lo[10];
  TH1*   zz_hist_shape_hi[10];
  double zz_nsig_shape_hi[10];
  double zz_hist_sf;

  zz_hist_sf          = 1.;
  np                  = 0; // 1;
  zz_pnames[0]        = "lumi";
  zz_par_sf_lo[0]     = -0.01;
  zz_par_sf_hi[0]     =  0.01;
  zz_hist_shape_lo[0] =  NULL;
  zz_hist_shape_hi[0] =  NULL;
  zz_nsig_shape_lo[0] = 1.;
  zz_nsig_shape_hi[0] = 1.;
//-----------------------------------------------------------------------------
// test background histogram: add it to null and test channles
//-----------------------------------------------------------------------------
  h11 = new TH1F("test_bgr","Test Background",Tzz->fNBinsTest,101.,106.);

  for (int i=1; i<=Tzz->fNBinsTest; i++) {
    h11->SetBinContent(i,Tzz->fBgrLevel);
    h11->SetBinError  (i,Tzz->fBgrLevel/10.);
  }

  Channel->add_template(h11,
			zz_hist_sf,
			np,
			zz_pnames,
			zz_par_sf_lo,
			zz_par_sf_hi,
			zz_hist_shape_lo,
			zz_nsig_shape_lo,
			zz_hist_shape_hi,
			zz_nsig_shape_hi,
			2,                  // use errors stored in the histogram
			0);                 // this is background

  return 0;
}



//-----------------------------------------------------------------------------
// add signal from a given 'Process' to a channel 'ChannelName' of a 'Model'
//-----------------------------------------------------------------------------
int   TMu2eLimits::AddTestSignal(csm_model*  Model) {

  TH1             /* *h1(0), *h2(0),*/ *h11(0); //, *h12, *h15;
  int              np;
  const char      *pnames [10];
  double           par_sf_lo    [10];
  double           par_sf_hi    [10];
  TH1*             hist_shape_lo[10];
  double           nsig_shape_lo[10];
  TH1*             hist_shape_hi[10];
  double           nsig_shape_hi[10];
  double           hist_sf;
  //  const char*      pname;

  //  channel_data_t*  ch_data;
  //  aprocess         *process, *signal;
  //  analysis         *ana;
//-----------------------------------------------------------------------------
// test signal histogram
//-----------------------------------------------------------------------------
  h11 = new TH1F("test_sig","Test Signal",fNBinsTest,101,106);
  for (int i=1; i<=fNBinsTest; i++) {
    if (i == 1) {
      h11->SetBinContent(i,10.);
      h11->SetBinError  (i,0.001);
    }
    else {
      h11->SetBinContent(i,0.);
      h11->SetBinError  (i,0.);
    }
  }

  hist_sf          = 1.; // already normalized to the luminosity
  np               = 0 ; // 1;
  pnames[0]        = "lumi";
  par_sf_lo[0]     = -0.1;
  par_sf_hi[0]     =  0.1;
  hist_shape_lo[0] =  NULL;
  hist_shape_hi[0] =  NULL;
  nsig_shape_lo[0] = 1.;
  nsig_shape_hi[0] = 1.;

  Model->add_template(h11,
		      hist_sf,
		      np,
		      pnames,
		      par_sf_lo,
		      par_sf_hi,
		      hist_shape_lo,
		      nsig_shape_lo,
		      hist_shape_hi,
		      nsig_shape_hi,
		      2,                 // use error bars stored in the histogram
		      1,                 // this is signal
		      "test");
  return 0;
}


