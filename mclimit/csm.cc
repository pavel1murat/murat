///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cstdlib>
using namespace std;

#include "TMath.h"
#include "TMinuit.h"
#include "mclimit/csm.h"
#include "mclimit/csm_model.h"


ClassImp(csm)

static csm_model*     modeltofit;
static vector<char*>  npfitname;
static vector<const TH1*>   datatofit;
static vector<char*>  datatofitname;
static vector<Int_t>  constrainedfitparam;


void csm_minuit_fcn(Int_t &npar, double *gin, double &f,
                        double *par, Int_t iflag);


/*----------------------------------------------------------------------------*/
// Use TMinuit to minimize the chisquared in T. Devlin's note CDF 3126 wrt the
// nuisance parameters.
// updated here -- do a joint minimization over shared nuisance parameters
// for several histograms (channels).
// global (in this file) declarations are at the top of the source file
// constructor
//-----------------------------------------------------------------------------
csm::csm() {
  //  cout << "Chisquared Minimizer Constructor called\n";
  datatofit.clear();
  datatofitname.clear();
  constrainedfitparam.clear();
  npfitname.clear();
  fitcov=0;  // we have not yet allocated memory for the fit error matrix.
  nfitcov=0;
  minuitmaxcalls = 500;
  minosmaxcalls = 500;
  minuitstepsize = 0.1;
  minuitprintflag = 0;
  minosflag = 0;
}


//-----------------------------------------------------------------------------
//destuctor
//-----------------------------------------------------------------------------
csm::~csm() {
  Int_t i;
  //cout << "Chisquared Minimizer Destructor called\n";

  // clear out static global variables

  for (i=0;i<(Int_t) datatofit.size(); i++)
    {
      delete datatofit[i];
      delete[] datatofitname[i];
    }
  datatofit.clear();
  datatofitname.clear();

  // clear out allocated memory pointed to by our private members

  for (i=0;i<(Int_t) fitparamname.size();i++)
    {
      delete[] fitparamname[i];
    }
  if (fitcov) delete[] fitcov;
}



//-----------------------------------------------------------------------------
Double_t csm::chisquared() {
  vector<char*> npn;
  vector<Double_t> nplb;
  vector<Double_t> nphb;
  Double_t arglist[20];
  Int_t ierflag = 0;
  Int_t i,j,icons;
  Double_t cresult;
  Double_t param,paramerror;

  modeltofit->list_nparams(&npn, &nplb, &nphb);
  Int_t npns = npn.size();

  if (npns > 0)
    {
      
      TMinuit *mnp = new TMinuit(npns+1);

      mnp->SetFCN(csm_minuit_fcn);
      if (!minuitprintflag)
	{
          arglist[0] = -1;
          mnp->mnexcm("SET PRINT", arglist, 1, ierflag);
          mnp->mnexcm("SET NOW",arglist,1,ierflag);
	}

      arglist[0] = 2; 
      mnp->mnexcm("SET STRATEGY", arglist, 1, ierflag);	
      mnp->mnexcm("SET NOG",arglist,1,ierflag);  // no gradiants required of FCN


      mnp->SetMaxIterations(minuitmaxcalls); // doesn't seem to do anything in TMinuit
      char npname[10];

      for (i=0;i < (Int_t) npns;i++)
        {
          sprintf(npname,"np%d",i);
	  //cout << "setting minuit parameter: " << npname << " " << nplb[i] << " " << nphb[i] << endl;
	  TString npname2 = npname;
          mnp->mnparm(i,npname2,0.0,minuitstepsize,nplb[i],nphb[i],ierflag);
          icons = 1;
	  for (j=0;j<(Int_t) modeltofit->npcm.size();j++)
	    {
	      if (strcmp(modeltofit->npcm[j].pnameoutput,npn[i])==0)
		{
		  ierflag = mnp->FixParameter(i);
		  icons = 0;
		}
	    }
          char *s = new char[strlen(npn[i])+1];
          strcpy(s,npn[i]);
          npfitname.push_back(s);  // this copy is in static global storage so the minuit function knows about it
          if (strstr(npn[i],"UNCONSTRAINED") != 0)
  	    {
	      icons = 0;
	    }
          constrainedfitparam.push_back(icons);
        }

      arglist[0] = 1;
      mnp->mnexcm("SET ERR", arglist ,1,ierflag); 
      ierflag = 0;
      arglist[0] = minuitmaxcalls;  // here's where maxcalls makes a difference
      arglist[1] = 1.;

      //      mnp->mnexcm("SIMPLEX", arglist ,2,ierflag);
      //      mnp->mnexcm("MIGRAD", arglist ,2,ierflag);

      mnp->mnexcm("MINI", arglist ,2,ierflag);
      mnp->mnexcm("IMPROVE", arglist ,2,ierflag);

      if (minosflag) 
	{
          arglist[0] = minosmaxcalls;
          mnp->mnexcm("MINOS",arglist,1,ierflag);
	}

      //cout << "Number of function calls in Minuit: " << mnp->fNfcn << endl;

      // copy best fit parameters for outside use

      cresult = mnp->fAmin;
      fitparam.clear();
      fiterror.clear();
      for (i=0;i<(Int_t) fitparamname.size();i++)
        {
          delete[] fitparamname[i];
        }
      fitparamname.clear();

      // allocate memory for the covariance matrix only if we have to.
      // (re-use the old memory if it has the right size)

      if (nfitcov != npns)
	{
	  if (fitcov != 0)
	    { 
	      delete[] fitcov;
	    }
	  nfitcov = 0;
	}
      if (nfitcov == 0)
        { 
	  fitcov = new Double_t[npns*npns];
	  nfitcov = npns;
	}

      for (i=0;i < (Int_t) npns;i++)
        {
          mnp->GetParameter(i,param,paramerror);
          fitparam.push_back(param);
          fiterror.push_back(paramerror);
          char *s = new char[strlen(npn[i])+1];
          strcpy(s,npn[i]);
          fitparamname.push_back(s); // this copy's part of the class private members
	  mnp->mnemat(fitcov,npns);
        }

      delete mnp;
    }
  else
    {
      i = 0;
      csm_minuit_fcn(i,0,cresult,0,0);
    }
  if (cresult < 0)
    { 
      //cout << "chisquared less than zero: " << cresult << " setting it to zero" << endl;
      cresult = 0;
    }

  for (i=0;i<(Int_t) npfitname.size();i++)
    {
      delete[] npfitname[i];
    }
  npfitname.clear();
  constrainedfitparam.clear();

  return(cresult);
}


void csm::setminuitmaxcalls(Int_t maxcalls)
{
  minuitmaxcalls = maxcalls;
}
Int_t csm::getminuitmaxcalls()
{
  return(minuitmaxcalls);
}

void csm::setminosmaxcalls(Int_t maxcalls)
{
  minosmaxcalls = maxcalls;
}
Int_t csm::getminosmaxcalls()
{
  return(minosmaxcalls);
}

void csm::setminuitstepsize(Double_t stepsize)
{
  minuitstepsize = stepsize;
}
Double_t csm::getminuitstepsize()
{
  return(minuitstepsize);
}

void csm::setprintflag(bool pf)
{
  minuitprintflag = pf;
}
bool csm::getprintflag()
{
  return(minuitprintflag);
}

void csm::setminosflag(bool mf)
{
  minosflag = mf;
}
bool csm::getminosflag()
{
  return(minosflag);
}

//-----------------------------------------------------------------------------
// put in the data histograms in the same order we built up the model histograms
//-----------------------------------------------------------------------------
void csm::set_htofit(const TH1 *Hist, const char* Channel) {
  Int_t i,ifound,j,jfound;
  vector<char*>::iterator dni;
  vector<const TH1*>::iterator dfi;
  char *s;

  ifound = -1;
  jfound = -1;
  for (i=0; i < (Int_t) datatofitname.size(); i++) {
    j = (Int_t) strcmp(Channel,datatofitname[i]);
    if (j == 0) {
      ifound = i;
    }
    if (j>0 && jfound == -1) {
      jfound = i;
    }
  }
/* if the name isn't already in the list, add it to the vector of names and
   make a blank model for it too.  Put the new name in it sorted place, sorted
   by increasing sort order of the name strings.  If the name is on the 
   list, replace the existing data histogram with a clone of the one supplied. */
  
  if (ifound == -1) {
    s = new char[strlen(Channel)+1];
    strcpy(s,Channel);
    if (jfound == -1) {
      datatofitname.push_back(s);
      datatofit.push_back((TH1*) Hist->Clone());
    }
    else {
      dni = datatofitname.begin() + jfound;
      datatofitname.insert(dni,s);
      dfi = datatofit.begin() + jfound;
      datatofit.insert(dfi,(TH1*) Hist->Clone());
    }
  }
  else {
    delete datatofit[ifound];
    datatofit[ifound] = (TH1*) Hist->Clone();
  }
}

//-----------------------------------------------------------------------------
void csm::set_modeltofit(csm_model* Model, const TH1** Hist) {
  modeltofit = Model;

  // assume (check?) that the model channel names match up with the data channel names
  // define histograms to fit

  int nc = Model->channame.size();

  for (int i=0; i<nc; i++) {
    set_htofit(Hist[i],Model->channame[i]);
  }
}


//-----------------------------------------------------------------------------
csm_model* csm::getbestmodel() {
  Int_t i;

  // make a local array of pointers to nuisance parameter names
  
  char **fpnameloc = new char *[fitparamname.size()];
  for (i=0;i<(Int_t) fitparamname.size();i++)
    {
      fpnameloc[i] = fitparamname[i];
      //cout << "in getbestmodel, paramname: " << fpnameloc[i] << endl;
    }
  Double_t *parloc = new Double_t[fitparam.size()];
  for (i=0;i<(Int_t) fitparam.size();i++)
    {
      parloc[i] = fitparam[i];
      //cout << "in getbestmodel, param: " << parloc[i] << endl;
    }

  modeltofit->nuisance_response(fitparam.size(),fpnameloc,parloc);
  delete[] fpnameloc;
  delete[] parloc;
  return(modeltofit);
}

void csm::plotcompare(char *cname)
{
  Int_t i;
  for (i=0;i<(Int_t)datatofitname.size();i++)
    {
      if (strcmp(datatofitname[i],cname)==0)
	{
	  modeltofit->plotwithdata(cname,datatofit[i]);
	}
    }
}

// Number of degrees of freedom -- this is approximately true for
// large statistics (in fact, the whole chisquared idea is only approximately
// true in cases of large statistics where distributions are approximately Gaussian)
// Degrees of freedom "freeze out" as the expected number of events gets small
// (<5 or so).  A bin with no data and no expectation shouldn't contribute either
// to the chisquared or the number of degrees of freedom, and neither really should
// a bin with 1E-6 expected and no observed events.  This routine won't draw the
// line (and even interpolated histograms can have variable numbers of bins with
// zero expectation).  This routine's very naive and just counts bins, filled or not.

Int_t csm::ndof()
{
  Int_t ndofl,i;
  vector<char*> npn;
  vector<Double_t> nplb;
  vector<Double_t> nphb;

  ndofl = 0;
  for (i=0;i<(Int_t) datatofit.size();i++)
    {
      ndofl += datatofit[i]->GetNbinsX()*datatofit[i]->GetNbinsY();
    }
  modeltofit->list_nparams(&npn, &nplb, &nphb);
  ndofl -= npn.size();
  cout << "nDOF isn't very clearly defined here... todo" << endl;
  return(ndofl);
}


Int_t csm::getnparams()
{
  return(fitparam.size());
}

Double_t csm::getparam(Int_t iparam)
{
  return(fitparam[iparam]);
}

Double_t csm::getcov(Int_t iparam, Int_t jparam)
{
  Int_t nparams = fitparam.size();
  return(fitcov[iparam+nparams*jparam]);
}

Double_t csm::getperror(Int_t iparam)
{
  return(fiterror[iparam]);
}

const char* csm::getpname(Int_t iparam)
{
  return(fitparamname[iparam]);
}

//-----------------------------------------------------------------------------
// Call the individual channel chisquared calculators inside here.
// the parameters par are labeled by their names npfitname in the 
// static global vector.
//-----------------------------------------------------------------------------
void csm_minuit_fcn(Int_t &npar, Double_t */*gin*/, Double_t &f, Double_t *par, Int_t /*iflag*/)
{
  Int_t i;

  // adjust the model according to the nuisance paramters
  
  char **fpnameloc = new char *[npfitname.size()];
  for (i=0;i<(Int_t) npfitname.size();i++) {
    fpnameloc[i] = npfitname[i];
    // cout << "in minuit fit fcn: " << i << " " << npfitname[i] << endl;
  }

  modeltofit->nuisance_response(npfitname.size(),fpnameloc,par);

  const TH1** dfloc = new const TH1*[datatofit.size()];
  for (i=0;i<(Int_t) datatofit.size();i++) {
    dfloc[i] = datatofit[i];
  }

  // cout << "In minuit function: printing out the model" << endl;
  // modeltofit->print();

  f = modeltofit->chisquared1(dfloc);

  // Gaussian constraints for variables which are constrained.

  for (i=0;i<npar;i++) {
    if (constrainedfitparam[i] != 0) {
      //cout << "In fcn: " << i << " " << par[i] << endl;
      f += par[i]*par[i];
    }
  }

  // cout << "end of computation of f in minuit_fit_fcn: " << f << endl;
  // printf("%20.14f\n",f);

  delete[] fpnameloc;
  delete[] dfloc;
}

