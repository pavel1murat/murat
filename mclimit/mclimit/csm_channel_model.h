#ifndef csm_channel_model_h
#define csm_channel_model_h

#include "TNamed.h"
#include "TH1.h"

#include <vector>
using std::vector;

// use to steer the handling of bin-to-bin uncertainties in templates

//maximum number of iterations in the quadratic system solver
#define CSM_MAXITER 100
//if all rates change fractionally by less than this, then we declare
//the system to be solved.
#define PREC1 1.0e-8

#define CSM_NOBINERR 0
#define CSM_POISSON_BINERR 1
#define CSM_GAUSSIAN_BINERR 2

// horizontal is the interpolation style which calls csm_pvmorph (or csm_pvmorph_2d)
// vertical interpolation just does a linear interpolation bin-by-bin, but not
// letting any bin go below zero

typedef enum {
  CSM_INTERP_HORIZONTAL,
  CSM_INTERP_VERTICAL,
  CSM_INTERP_HORIZONTAL_EXTRAP,
  CSM_INTERP_VERTICAL_EXTRAP
} INTERPSTYLE;

struct svstruct_s {
  Int_t itemplate;   // which template histo this syst. variation applies to
  char *sysname;     // name of nuisance parameter this corresponds to
  Double_t sysfracl;
  Double_t sysfrach;
                     // if there is no shape uncertainty associated with this dependence of a particular
                     // template on a nuisance parameter, then these pointers should be set to zero.
  TH1 *lowshape;     // for shape uncertainty -- low histogram shape histo id
  TH1 *highshape;    // for shape uncertainty -- high histogram shape
  Double_t xsiglow;  // how many sigma low lowshape corresponds to (should be a negative number)
  Double_t xsighigh; // how many sigma high highshape corresponds to. (should be a positive number)
};

typedef struct svstruct_s svstruct;

class csm_model;
class csm;

// interpolate histogram with a pvmorph-style procedure, or with vertical interpolation.
// csm_interpolate_histogram interpolates the errors too, while csm_interpolate_histogram_noerr
// is a speedup which interpolates only the bin contents and not the errors

void csm_interpolate_histogram(TH1*,Double_t,TH1*,Double_t,TH1*,Double_t,INTERPSTYLE);

void csm_interpolate_histogram_noerr(TH1*,Double_t,TH1*,Double_t,TH1*,Double_t,INTERPSTYLE);

// version to be used with cascading shape errors -- needs a central shape, a varied shape,
// and a shape to apply the variations to (which may not be either of the above, but the result
// of a different shape variation) startshape.  The output is outshape.

// csm_interpolate histogram2 calls csm_interpolate_histogram3 twice, once for the
// bin contents, once for the errors.

void csm_interpolate_histogram2(TH1* central, Double_t paramcentral,
                                TH1* varied, Double_t paramvaried,
                                TH1* startshape, 
                                TH1* outshape,
                                Double_t param,
                                INTERPSTYLE istyle);

// here's a version that just interpolates the bin contents but not the errors
//  (twice as fast)

void csm_interpolate_histogram2_noerr(TH1* central, Double_t paramcentral,
                                TH1* varied, Double_t paramvaried,
                                TH1* startshape, 
                                TH1* outshape,
                                Double_t param,
                                INTERPSTYLE istyle);

// this routine just interpolates the histogram contents

void csm_interpolate_histogram3(TH1* central, Double_t paramcentral,
                                TH1* varied, Double_t paramvaried,
                                TH1* startshape, 
                                TH1* outshape,
                                Double_t param,
                                INTERPSTYLE istyle);

void csm_pvmc(Int_t nb, Double_t *dist1, Double_t *dist2, Double_t *dist3, Double_t *distn,
	      Double_t par1, Double_t par2, Double_t parn);

void csm_pvmc2d(Int_t nx, Int_t ny, Double_t *xydist1, 
                Double_t *xydist2, Double_t *xydist3, Double_t *xydistn,
                Double_t par1, Double_t par2, Double_t parn);

void csm_yproj(Int_t nx, Int_t ny, Double_t *xydist, Double_t *ydist);

void csm_ycont(Int_t ny, Double_t *ydist1, Double_t *ydist2,
               Double_t *ydist3, Double_t *ydistn,
               Double_t *alpha1, Double_t *alpha2, Double_t *alpha3);

void csm_ycontaux(Int_t ny, Double_t *y, Double_t *yn,
                  Double_t *alpha);

void csm_acnvec2(Double_t *vec, Int_t n);



class csm_channel_model: public TNamed {
 public:
  csm_channel_model(const char* Name="");
  ~csm_channel_model();
  void add_template( TH1 *,          // template histogram
		     Double_t,       // scale factor to multiply template by to compare w/ data
		     Int_t,          // number of nuisance parameters (Gaussian of unit width)
		     const char *[], // nuisance parameter names 
		     Double_t *,     // fractional uncertainty on sf due to each nuisance parameter -- low side
		     Double_t *, // fractional uncertainty on sf due to each nuisance parameter -- high side
		     TH1 *[],    // array of low histogram shapes, one for each nuisance param
		     Double_t *, // number of sigma low for each nuisance parameter shape variation
		     TH1 *[],    // array of high histogram shapes, one for each nuisance param
		     Double_t *, // number of sigma high for each shape variation
		     Int_t,      // MC statistical error steering (inupt).    Values and meaning:
                                  //    2:  Pay attention to error bars in each bin
                                  //    1:  Poisson errors with unit-weight entries (see below)
                                  //    0:  No bin-by-bin errors.
                                  //  Notes on Poisson (value=1):  There is a split between the
                                  // Poisson handling when the model is an ensemble model (for pseudoexperiment
                                  // generation) and when it is a test model (for fitting to data and pseudodata).
                                  //  In an ensemble model, the bin contents are treated
                                  // as Poisson means and the values are fluctuated according to Poisson statistics.
                                  // These random components are then used with their scale factors as components
                                  // of another Poisson mean for generating pseudodata. 
                                  // In the model used to fit to the data, these values are treated as integer
                                  // measurements from subsidiary experiments.
		     Int_t);     // Scale flag -- set to 1 if you want this template to be scaled in the s95
                                  // calculation, 0 if you don't.  It's intended that applications should set this to 1 for
                                  // signal histograms, and 0 for background histograms.


  csm_channel_model* clone() const;    // make an exact copy of this channel model -- all the internal
                                 // cloned histograms are cloned again.  Better than just
                                 // assigning a new model to this one because the
                                 // destructors won't delete the same memory twice.
  void nuisance_response(Int_t, char *[], Double_t []); // update the internal copies of the varied histograms
                                                        // and scale factors according to the nuisance parameters supplied
  void undo_nuisance_response();  // resets all the varied copies of the histogram templates and scale factors
                                   // to their original values

  void   print();                            // print out some details of the model
  void   plotwithdata(const TH1* Hist);      // compare data with a model.
  void   candcheck(TH1*);                    // print out high s/b candidates
  double kstest(TH1*);                       // get ROOT's raw KS prob
  double kstest_px(TH1*);                    // get ROOT's raw KS px prob
//-----------------------------------------------------------------------------
// adding two models together, and multiplying a model by a scalar
// note that all of these operations on models create a new model 
// (which must be cleaned up later)
//-----------------------------------------------------------------------------
  csm_channel_model* add(csm_channel_model&);// adds the components of two models together
  csm_channel_model* scale(Double_t);        // scales all parts of the model
  csm_channel_model* scalesignal(Double_t);  // scales only those parts of the model called "signal"

  csm_channel_model* scale_err(Double_t);    // scales rates up with scale factor, and scales
                                             // systematic errors down with scale factor.  n.b. --
                                             // MC statistical errors on model histos cannot be scaled

  Double_t chisquared1(const TH1* Hist);     // Inputs a pointer to data histogram -- 
					     // computes the chisquared of this model
                                             // compared to the supplied data histogram, 
					     // without minimizing over nuisance
					     // parameters, a la Tom Devlin's note

  Int_t checkneg();                          // check for negative bins
  
  void set_interpolation_style(INTERPSTYLE); // CSM_INTERP_HORIZONTAL or CSM_INTERP_VERTICAL
                                             // horizontal (csm_pvmorph) is the default 
					     // if this is not called.

  friend class     csm_model;                // so we can access systematic error names and limits
  friend class     mclimit_csm;

 private:
  vector<TH1*>     histotemplate;
  vector<TH1*>     histotemplate_varied;
  vector<Double_t> sft;
  vector<Double_t> sft_varied;
  vector<Int_t>    poissflag;
  vector<Int_t>    scaleflag;
  vector<svstruct> syserr;
  INTERPSTYLE      chan_istyle;

  ClassDef(csm_channel_model,0)
};

#endif
