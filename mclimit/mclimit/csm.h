#ifndef csm_h
#define csm_h

#include "TObject.h"
#include "TH1.h"

#include <vector>
using std::vector;

#include "csm_model.h"

class csm : public TObject {
  public:
    	csm();
	~csm();
        void     set_htofit(const TH1* Hist, const char* Channel);
        void     set_modeltofit(csm_model* Model, const TH1** Hist); // a set of template histograms to fit to the data histograms
	Double_t chisquared();           // calculates chisquared
        Int_t ndof();                    // calculates an (approximate!) NDOF
	Int_t    getnparams();
        Double_t getparam(Int_t);
        Double_t getperror(Int_t);
        Double_t getcov(Int_t, Int_t);  // get an entry out of the covariance matrix.
        const char* getpname(Int_t);
        csm_model* getbestmodel(); // a model with the template histograms all normalized to the best fit values.
                                   // suitable for plotting.  This varies the input model (given in set_modeltofit)
                                   //  with the best fit nuisance parameters and returns a pointer to the same model
                                   // given in set_modeltofit.
        void plotcompare(char *); // make a stacked plot of data in the named channel against the central
                                              // value model provided.
        void setminuitmaxcalls(Int_t);  // Maximum number of function calls MINUIT is allowed to do per minimization
                                        // default: 500
        Int_t getminuitmaxcalls(); 
        void setminosmaxcalls(Int_t);  // Maximum number of function calls MINOS is allowed to do per parameter
                                        // default: 500
        Int_t getminosmaxcalls(); 
        void setminuitstepsize(Double_t);  // Initial step size for MINUIT fit parameters default: 0.1
        Double_t getminuitstepsize(); 
        void setminosflag(bool);   // true: call MINOS. Best to have printing set too if 
                                   // you're going to run MINOS.  False:  Do not call MINOS (default) 
        bool getminosflag(); 
	void setprintflag(bool);  // true:  let MINUIT print stuff out;  FALSE -- turn off MINUIT printing (default)
	bool getprintflag();

 private:
        vector<Double_t> fitparam;   // parameters of the fit
	vector<Double_t> fiterror;   // errors on fit parameters
        vector<char*> fitparamname; //  names of fit parameters
	Double_t *fitcov;           //  pointer to covariance matrix
	Int_t nfitcov;              //  size of covariance matrix
        Int_t minuitmaxcalls;       //  how many calls to the function Minuit is allowed to make
        Int_t minosmaxcalls;       //  how many calls to the function Minos is allowed to make
	Double_t minuitstepsize;    // initial step size for MINUIT parameters -- default 0.1
        bool minuitprintflag;       // true if we let MINUIT print stuff out.
        bool minosflag;             // true if we are to call MINOS (most useful when the printflag is on too)

	ClassDef(csm,0)
};

#endif
