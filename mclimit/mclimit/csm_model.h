#ifndef csm_model_h
#define csm_model_h

#include "TObject.h"
#include "TH1.h"
#include "TRandom3.h"

#include <vector>
using std::vector;

#include "csm_channel_model.h"


// Constraint equations between nuisance parameters.  Sometimes there are fewer degrees of freedom
// than there are nuisance parameters.  Example:  MET vs. Iso evaluation of the non-W
// background in a signal region, using the 4-sector method.  A*C/B=D relates four
// nuisance parameters down to three: D is computed from A, B, and C.  Note:  multiple constraints between nuisance
// parameters will not be solved for, they will just be evaluated in order.  In pseudexperiments,
// nuisance parameters that are constrained are not randomly chosen, but rather are calculated
// based on the other nuisance parameters.  The arbitrary function here means that no checking
// can be done to make sure that the computed nuisance parameters are physical.

struct npcstruct_s
{
  Int_t ninput;       // number of nuisance paramters input to the calculation of a constrained one
  char **pnameinput;  // the names of the input nuisance parameters
  char *pnameoutput;  // the name of the output parameter
  Double_t (*f)(Double_t*);  // a function taking an array of all the input parameters and computing the output parameter
};

typedef struct npcstruct_s npcstruct;


/* a full model for a particular channel.  csm_model is just a collection
   of channel models and associated names */


class csm_model: public TObject {
 public:
   csm_model();
   ~csm_model();
   void add_template( TH1 *,      //template histogram
	   	      Double_t,   //scale factor to multiply template by to compare w/ data
                      Int_t,      // number of nuisance parameters (Gaussian of unit width)
                      const char *[],   // nuisance parameter names 
                      Double_t *, // fractional uncertainty on sf due to each nuisance parameter -- low side
                      Double_t *, // fractional uncertainty on sf due to each nuisance parameter -- high side
                      TH1 *[],    // array of low histogram shapes, one for each nuisance param
                      Double_t *, // number of sigma low for each nuisance parameter shape variation (should be negative!)
                      TH1 *[],    // array of high histogram shapes, one for each nuisance param
		      Double_t *, // number of sigma high for each shape variation (should be positive!)
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
                      Int_t,      // Scale flag -- set to 1 if you want this template to be scaled in the s95
                                  // calculation, 0 if you don't.  It's intended that applications should set this to 1 for
                                  // signal histograms, and 0 for background histograms.
                      const char *);    // Channel name
//-----------------------------------------------------------------------------
// sets the interpolation style for a particlar channel 
// first arg: channel name
// second arg: interpolation style: CSM_INTERP_HORIZONTAL or CSM_INTERP_VERTICAL
//-----------------------------------------------------------------------------
   void set_interpolation_style(const char* ChannelName, INTERPSTYLE Style);  

   void add_chanmodel(csm_channel_model*,
                      const char*); // instead of adding a template at a time, let's add a whole channel's model

   void add_npcons(Int_t, char**, const char*, Double_t (*f)(Double_t*)); // for a nuisance parameter which
	    // can be computed given the values of other nuisance parameters, this allows
	    // the user to specify such a funtion.   Arguments;  number of n.p.'s the one to calculate
            // depends on, their names, the name of the n.p. to calculate, and a pointer to the
            // function which does the job.

   void add_npbounds(char *npname, Double_t lowbound, Double_t highbound);

   void plotwithdata(const char* Channel, const TH1* Hist); // plot a named channel's model with the supplied data histogram
   void candcheck(const char*,TH1*); // check candidates
   double kstest(const char*,TH1*); // get raw ROOT ks prob
   double kstest_px(const char*,TH1*); // get ROOT's ks test prob with px simulation (stat only)

   virtual TObject* Clone(const char* Opt = "") const; // make an exact copy of this model -- all the internal
                                                       // cloned histograms are cloned again.  Better than just
                                                       // assigning a new model to this one because the
                                                       // destructors won't delete the same memory twice.

   void print();                 // print out some details of the model
   void print(char *channame);   // just print out the piece corresponding to channame
   // adding two models together, and multiplying a model by a scalar

   // note that all of these operations on models create a new model (which must be cleaned up later)

   csm_model* add(csm_model&); // adds the components of two models together
   csm_model* scale(Double_t);        // scales all parts of the model
   csm_model* scalesignal(Double_t);  // scales only those parts of the model called "signal"
   csm_model* scale_err(Double_t);  // scales rates up with scale factor, and scales
                                            // systematic errors down with scale factor.  n.b. --
                                            // MC statistical errors on model histos cannot be scaled
   void nuisance_response(Int_t, char *[], Double_t []); // updates the fluctuated version of the histogram templates
                                 // and scale factors inside the channel model
                                 // according to the nuisance parameters provided.  Inputs: number of nuisance
                                 // parameters and their names
   void undo_nuisance_response();  // resets all the varied copies of the histogram templates and scale factors
                                   // to their original values
   void varysyst();   // randomly choose nuisance parameter values and call nuisance_response
   void print_nuisance_params();  // after calling varysyst, can print out the randomly chosen parameters and names
   void single_pseudoexperiment(TH1 *[]); // generate pseudodata for all the channels in this model.
   void list_nparams(vector<char *> *, vector<Double_t> *, vector<Double_t> *); // get a list of
      // unique nuisance parameter names for all the channels in this model and their most
      // restrictive lower and upper bounds.

   Double_t chisquared1(const TH1* []); // calls the chisquared routine for each channel model.

   friend class mclimit_csm;
   friend class csm;

 private:
  vector<char*> channame;
  vector<csm_channel_model*> chanmodel;
  Int_t lookup_add_channame(const char*);  // look up the channel name in the channame vector.  If it's
                                       // not there, add it.
  vector<npcstruct> npcm;  // constraint equations between nuisance parameters

  vector<char*> npbname;   // upper and lower bounds
  vector<Double_t> npbhigh;
  vector<Double_t> npblow;
  vector<Double_t> npvalp; // handles to stored lists of nuisance parameter names and values for
  vector<char*> npnp;

  TRandom3*     fRandom;

  ClassDef(csm_model,0)
};



#endif
