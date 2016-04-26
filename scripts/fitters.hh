//
#ifndef __murat_scripts_fitters_hh__
#define __murat_scripts_fitters_hh__

  struct fitter_t {
    const char *name;		// fitter name
    const char *script;         // name of teh script to load
  };

fitter_t fitter_asymm_gauss = { "asymm_gauss", "murat/scripts/fit_asymm_gauss.C" };

#endif
