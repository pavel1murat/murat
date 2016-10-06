///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <limits>
#include <stdlib.h>
#include <math.h>

using namespace std;

#include "TRandom.h"
#include "mclimit/csm_model.h"

/*----------------------------------------------------------------------------*/
// update the internal representation of the nuisance response of the model.
// Use the method for doing this with each
// channel separately.

void csm_model::nuisance_response(Int_t nparams,
                                  char *paramname[],
                                  Double_t paramvalue[])
{
  Int_t i,nchans;
  Int_t ipar,icons,j,k,ifound,jfound;
  Double_t *cinput;

  /*
  cout << "in model nuisance response: " << endl;
  for (i=0;i < (Int_t) nparams; i++)
  {
    cout << "param: " << i << " name: " << paramname[i] << endl;
  }
  */

  // compute the nuisance parameters which are functions of the others
  // loop over constraints and replace the parameter values with the constrained ones
  // do not overwrite the input parameters but make our own copy

  Double_t *parloc = new Double_t[nparams];
  for (i=0;i<nparams;i++)
    {
      parloc[i] = paramvalue[i];
    }
  for (icons=0;icons<(Int_t)npcm.size();icons++)
    {
      jfound = 0;
      for (ipar=0;ipar<nparams;ipar++)
	{
	  //cout << "Comparing: " << npcm[icons].pnameoutput << " with " << paramname[ipar] << endl;
          if (strcmp(npcm[icons].pnameoutput,paramname[ipar])==0)
	    {
	      jfound = 1;
	      cinput = new Double_t[npcm[icons].ninput];
	      for (j=0;j<(Int_t)npcm[icons].ninput;j++)
	        {
	          ifound = 0;
	          for (k=0;k<nparams;k++)
		    {
		      if (strcmp(npcm[icons].pnameinput[j],paramname[k])==0)
		        {
		          cinput[j] = parloc[k];
		          ifound = 1;
		        }
		    }
	          if (ifound == 0)
		    {
		      cout << "Didn't find parameter name: " << 
                          npcm[icons].pnameinput[j] << 
                          " in the list of nuisance parameters" << endl;
		      exit(0);
		    }
	        }
	      parloc[ipar] = npcm[icons].f(cinput);
	      delete[] cinput;
	    }
	}
      // Not a disaster if we didn't find a parameter on the list.  Sometimes we
      // only fit for a subset of parameters and they aren't on the list.
      if (jfound == 0) {
	cout << "Constraint equation found for nuisance parameter not on our list: " 
             << npcm[icons].pnameoutput << endl;
      //  exit(0);
      }
    }

  nchans = (Int_t) channame.size();
  for (i=0;i< nchans; i++)
    {
      chanmodel[i]->nuisance_response(nparams,paramname,parloc);
    }
  delete[] parloc;
}

/*----------------------------------------------------------------------------*/

void csm_model::undo_nuisance_response()
{
  Int_t i,nchans;
  nchans = (Int_t) channame.size();
  for (i=0;i< nchans; i++)
    {
      chanmodel[i]->undo_nuisance_response();
    }
}

//----------------------------------------------------------------------------
// collect all nuisance parameter names and upper and lower bounds for this model
// Where the bounds come from -- do not allow extrapolation on histogram shapes.  
// (all histogram extrapolation should be done and verified by the user)
// Also do not allow any contribution to signal or background to go negative.
//-----------------------------------------------------------------------------
void csm_model::list_nparams(vector<char*> *npn, vector<Double_t> *nplb, vector<Double_t> *nphb) {

  Int_t              ifound;
  csm_channel_model* cm;
  Double_t           nplb_tmp, nphb_tmp, a, b, c, disc, xp, xm, xht, xlt;

  npn->clear();
  nplb->clear();
  nphb->clear();

  for (uint i=0; i<channame.size(); i++) {
    cm = chanmodel[i];
    for (uint j=0; j<cm->syserr.size(); j++) {
      //cout << "sys error item channel: " << i << 
      //" error index: " << j << " " << cm->syserr[j].sysname << " " <<
      //(cm->syserr[j].highshape != 0) << " " <<  
      //(cm->syserr[j].lowshape != 0) << " " <<
      //cm->syserr[j].xsiglow << " " << cm->syserr[j].xsighigh << endl;  
      
      // question -- do we need to consider nuisance parameter variations beyond 20 sigma?
      // probably not if we only need 5-sigma discovery significance.
      
      nplb_tmp = -20;
      nphb_tmp = 20;

      // Require the user to supply shape variations out to the number of sigma
      // we will investigate here.  This program won't do shape extrapolations internally,
      // but the csm_pvmorph subroutine supplied will in fact extrapolate.  Users should
      // look at what they get when extrapolating histograms, though -- check and validate.
      
      // 26 Feb 2008 -- allow shape extrapolations of templates beyond the provided ranges.

      if (cm->chan_istyle != CSM_INTERP_HORIZONTAL_EXTRAP && cm->chan_istyle != CSM_INTERP_VERTICAL_EXTRAP) {
	if (cm->syserr[j].lowshape != 0) { 
	  nplb_tmp = max(nplb_tmp,cm->syserr[j].xsiglow); 
	}
	if (cm->syserr[j].highshape != 0){ 
	  nphb_tmp = min(nphb_tmp,cm->syserr[j].xsighigh); 
	}
      }

      // limit the nuisance paramters also so that individual scale factors do not go negative.
      // There's protection in the fit function, but we need the pseudoexperiments also to
      // be sensible -- this is the equivalent (using the asymmetric errors supplied) of the
      // truncated Gaussian
      
      a = (cm->syserr[j].sysfrach + cm->syserr[j].sysfracl)/2.0;
      b = (cm->syserr[j].sysfrach - cm->syserr[j].sysfracl)/2.0;
      c = 1;
      if (a == 0) {
//-----------------------------------------------------------------------------
// symmetric error
//-----------------------------------------------------------------------------
	if (b > 0) { 
	  nplb_tmp = max(nplb_tmp,-1.0/b); 
	}
	if (b < 0) { 
	  nphb_tmp = min(nphb_tmp,-1.0/b); 
	}
      }
      else {
//-----------------------------------------------------------------------------
// asymmetric error
//-----------------------------------------------------------------------------
	disc = b*b - 4.0*a*c;
	if (disc > 0) { 
	  xp = (-b + sqrt(disc))/(2.0*a);
	  xm = (-b - sqrt(disc))/(2.0*a);
	  xht = max(xp,xm);
	  xlt = min(xp,xm); 
	  // we know that a nuisance parameter value of 0 has a non-negative prediction,
	  // but the choice of which of these two solutions to a quadratic to take depends
	  // on which side of zero they are on.
	  if (xht < 0) {
	    nplb_tmp = max(nplb_tmp,xht);
	  }
	  else if (xlt > 0) {
	    nphb_tmp = min(nphb_tmp,xlt);
	  }
	  else {
	    nphb_tmp = min(nphb_tmp,xht);
	    nplb_tmp = max(nplb_tmp,xlt);
	  }
	}
      }
      
      ifound = -1;
      for (uint k=0;k<npn->size();k++) {
	if (strcmp(cm->syserr[j].sysname,(*npn)[k]) == 0) { 
	  ifound = k; 
	}
      }
      if (ifound == -1) {
	npn->push_back(cm->syserr[j].sysname);
	nplb->push_back(nplb_tmp);
	nphb->push_back(nphb_tmp);
	//cout << "sysname: " << cm->syserr[j].sysname << " assigned ranges: " << nplb_tmp << " " << nphb_tmp << endl;
      }
      else {
	(*nplb)[ifound] = max((*nplb)[ifound],nplb_tmp);
	(*nphb)[ifound] = min((*nphb)[ifound],nphb_tmp); 
	//cout << "sysname: " << cm->syserr[j].sysname << " reassigned ranges: " << nplb_tmp << " " << nphb_tmp <<  " " <<
	// (*nplb)[ifound] << " " << (*nphb)[ifound] << endl;
      }
    }
  }
  // add in user-specified bounds (27 Feb 2008)

  for (uint k=0; k<npbname.size(); k++) {
    for (uint j=0; j< npn->size(); j++) {
      if (strcmp(npbname[k],(*npn)[j])==0) {
	(*nplb)[j] = max((*nplb)[j],npblow[k]);
	(*nphb)[j] = min((*nphb)[j],npbhigh[k]); 
      }
    }
  }
}

//----------------------------------------------------------------------------
// a splitoff from single_pseudoexperiment -- just vary the templates but do  
// not generate pseudodata.   Useful for interfacing with Joel's program      
//----------------------------------------------------------------------------
void csm_model::varysyst() {
  vector<Double_t> nplb;
  vector<Double_t> nphb;
  Double_t xval;

 // systematically fluctuate our model

  list_nparams(&npnp, &nplb, &nphb);   // this clears and fills the vectors of name pointers and bounds
  //cout << " in pe: npn.size " << npn.size() << endl;
  npvalp.clear();

  for (uint i=0; i<npnp.size(); i++) {
    do {
      xval = fRandom->Gaus(0,1);
      //cout << i << " " << nplb[i] << " " << xval << " " << nphb[i] << endl;
    } while (xval < nplb[i] || xval > nphb[i]);
    npvalp.push_back(xval);
  }
  nuisance_response(npnp.size(),&(npnp[0]),&(npvalp[0]));
}


//-----------------------------------------------------------------------------
// prints out names and values of nuisance parameters generated in the latest 
// call to varysyst.
//-----------------------------------------------------------------------------
void csm_model::print_nuisance_params()
{
  cout << "Nuisance parameter listing: " << endl;
  for (int i=0;i<(int) npnp.size();i++)
    {
      cout << i << " " << npnp[i] << " " << npvalp[i] << endl; 
    }
}

//-----------------------------------------------------------------------------
// Generate a single pseudoexperiment from a model -- fluctuate all nuisance parameters
// with their uncertainties -- the pseudodata histograms are in the same order as 
// the channels in the model description, with the same binning assumed.  The psuedodata
// histograms should be allocated in the calling routine.  That way the histograms don't
// have to be continually created and destroyed for each pseudoexperiment but can be
// re-used.
//-----------------------------------------------------------------------------
void csm_model::single_pseudoexperiment(TH1 *pseudodata[]) {
  Int_t ichan,itpl,ibinx,ibiny,nbinsx,nbinsy,nchans,ntemplates;
  csm_channel_model* cm;
  Double_t bintot;
  TH1* ht;
  Double_t r;

  // call nuisance_response with random nuisance parameters

  varysyst();
 
  // generate random pseudodata.  Randomly fluctuate the Poisson subsidiary
  // experiments (a "systematic effect") to figure out what the proper mean
  // is for the main experiment.

  nchans = channame.size();
  for (ichan=0;ichan<nchans;ichan++) {
    cm = chanmodel[ichan];
    nbinsx = cm->histotemplate[0]->GetNbinsX();
    nbinsy = cm->histotemplate[0]->GetNbinsY(); 
    for (ibinx=0;ibinx<nbinsx;ibinx++) {
      for (ibiny=0;ibiny<nbinsy;ibiny++) {
	bintot = 0;
	ntemplates = (Int_t) cm->histotemplate.size();
	for (itpl=0; itpl<ntemplates; itpl++) {
	  ht = cm->histotemplate_varied[itpl];

	  if (nbinsy == 1) r = ht->GetBinContent(ibinx+1); 
	  else             r = ht->GetBinContent(ibinx+1,ibiny+1); 

	  if (cm->poissflag[itpl] == CSM_POISSON_BINERR){ 
	    r = fRandom->Poisson(r); 
	  }
	  else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR) { 
	    double histerr,edraw;
	    if (nbinsy==1)  histerr = ht->GetBinError(ibinx+1);
	    else            histerr = ht->GetBinError(ibinx+1,ibiny+1);

	    do { 
	      edraw = fRandom->Gaus(0,histerr); 
	    } while (edraw+r<r*1E-6);           // don't let it hit zero or go negative.
	    r += edraw;
	  }
	  r *= cm->sft_varied[itpl];
	  bintot += r;
	}
	r = fRandom->Poisson(bintot);

	if (nbinsy == 1) pseudodata[ichan]->SetBinContent(ibinx+1,r); 
	else             pseudodata[ichan]->SetBinContent(ibinx+1,ibiny+1,r); 
      } // end loop over binsy
    } // end loop over binsx
  } // end loop over channels
}


/*----------------------------------------------------------------------------*/

// A model is a collection of channel models and names

csm_model::csm_model() {
  // use own random number generator
  fRandom = new TRandom3();
}

/*----------------------------------------------------------------------------*/

csm_model::~csm_model()
{
  Int_t i,j;
  for (i=0; i < (Int_t) channame.size(); i++)
    {
      delete[] channame[i];
      delete chanmodel[i];
    }
  for (i=0;i<(Int_t) npcm.size();i++)
    {
      for (j=0;j<npcm[i].ninput;j++)
	{
	  delete[] npcm[i].pnameinput[j];
	}
      delete[] npcm[i].pnameinput;
      delete[] npcm[i].pnameoutput;
    }
  for (i=0;i<(Int_t) npbname.size();i++)
    {
      delete[] npbname[i];
    }

  /* The vectors themselves are deleted when the class instance is deleted */

  delete fRandom;
}

/*----------------------------------------------------------------------------*/

void csm_model::add_template(TH1 *template_hist, //Poisson or non-Poisson histogram
                                Double_t sf,        //scale factor to multiply template by to compare w/ data 
                                                    //(e.g., (data_lum/MC_lum) for a MC Poisson histogram
                                Int_t nnp,          // number of nuisance parameters (Gaussian of unit width)
                                const char* npname[],     // nuisance parameter names 
                                Double_t *nps_low,  // fractional uncertainty on sf due to each nuisance parameter -- low side
                                Double_t *nps_high, // fractional uncertainty on sf due to each nuisance parameter -- high side
		                                    // typically nps_low and nps_high are input with opposite signs -- if opposite
                                                    // variations of the nuisance parameter create opposite changes in sf.  The sign
                                                    // is retained in the calculation in case both variations of a nuisance parameter
                                                    // shift the normalization in the same way (either both + or both -)
                                TH1 *lowshape[],    // array of low hisogram shapes, one for each nuisance param (null if no shape error)
                                Double_t *lowsigma, // number of sigma low for each nuisance parameter shape variation
                                TH1 *highshape[],   // array of high histogram shapes, one for each nuisance param (null if no shape error)
	                        Double_t *highsigma, // number of sigma high for each shape variation
                                Int_t pflag,         // Poisson flag -- 1 if Poisson, 0 of not.  2 if Gaussian error from the histo contents
                                Int_t sflag,         // scale flag -- 1 if signal, 0 if background (for use with s95 calculator)
                                const char *cname)
{
  Int_t i;
  i = lookup_add_channame(cname);
  chanmodel[i]->add_template(template_hist,sf,nnp,npname,nps_low,
                             nps_high,lowshape,lowsigma,highshape,
                             highsigma,pflag,sflag);
}

/*----------------------------------------------------------------------------*/
// add a whole channel's model to the total set of models.

void csm_model::add_chanmodel(csm_channel_model *cm, const char *cname)
{
  Int_t ichan;

  ichan = lookup_add_channame(cname);
  chanmodel[ichan] = cm->clone();
}

void csm_model::add_npbounds(char *pname, Double_t lowbound, Double_t highbound)
{
  char *s = new char[strlen(pname)+1];
  strcpy(s,pname);
  npbname.push_back(s);
  npblow.push_back(lowbound);
  npbhigh.push_back(highbound);
}

/*----------------------------------------------------------------------------*/
// add a constraint function between nuisance parameters.  Make our own copies of all
// the names and the function pointer.

void csm_model::add_npcons(Int_t nparin,char **parin, const char *parout, Double_t (*f)(Double_t*))
{
  Int_t i,j;
  npcstruct npc;
  char *s;

  for (i=0;i<(Int_t) npcm.size();i++)
    {
      if (strcmp(parout,npcm[i].pnameoutput)==0)
	{
	  cout << "Warning: Two constraint functions for the same nuisance parameter: " << parout << " defined" << endl;
          exit(0); // bad enough to crash
	}
      for (j=0;j<npcm[i].ninput;j++)
	{
	  if (strcmp(parout,npcm[i].pnameinput[j]) == 0)
	    {
	      cout << "Warning: nuisance parameter: " << npcm[i].pnameinput[j] << " depends on " << parout << endl;
              cout << "but " << parout << "is computed itsef by a constraint after " << npcm[i].pnameinput[j] << "is computed." << endl;
              exit(0); // bad enough to crash
	    }
	}
    }

  npc.ninput = nparin;
  npc.pnameinput = new char*[nparin];
  for (i=0;i<nparin;i++)
    {
      s = new char[strlen(parin[i])+1];
      strcpy(s,parin[i]);
      if (strcmp(s,parout)==0)
	{
	  cout << "Constraint function for nuisance parameter: " << s 
               << " depends on nuisance parameter " << s << endl;
	  exit(0);
	}
      npc.pnameinput[i] = s;
    }
  s = new char[strlen(parout)+1];
  strcpy(s,parout);
  npc.pnameoutput = s;
  npc.f = f;
  npcm.push_back(npc);
}

/*----------------------------------------------------------------------------*/
TObject* csm_model::Clone(const char* Opt) const {
  Int_t i;
  csm_model* mclone = new csm_model;

  for (i=0;i < (Int_t) channame.size(); i++)
    {
      mclone->add_chanmodel(chanmodel[i],channame[i]);
    } 
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      mclone->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  for (i=0;i < (Int_t) npbname.size(); i++)
    {
      mclone->add_npbounds(npbname[i],npblow[i],npbhigh[i]);
    }
  return(mclone);
}

/*----------------------------------------------------------------------------*/
// includes all the new contributing histograms, and also collects together all constraint
// relationships between nuisance parameters.

csm_model* csm_model::add(csm_model &a)
{
  Int_t i;
  csm_model* mclone = (csm_model*) a.Clone();
  for (i=0; i < (Int_t) channame.size(); i++)
    {
      mclone->add_chanmodel(chanmodel[i],channame[i]);
    }
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      mclone->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  for (i=0;i < (Int_t) npbname.size(); i++)
    {
      mclone->add_npbounds(npbname[i],npblow[i],npbhigh[i]);
    }
  return(mclone);
}

/*----------------------------------------------------------------------------*/
csm_model* csm_model::scale(Double_t coefficient)
{
  Int_t i;
  csm_channel_model* scmodel;
  csm_model* smodel = new csm_model;

  for (i=0; i< (Int_t) channame.size(); i++)
    {
      scmodel = chanmodel[i]->scale(coefficient);
      smodel->add_chanmodel(scmodel,channame[i]);
      delete scmodel;
    }
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      smodel->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  for (i=0;i < (Int_t) npbname.size(); i++)
    {
      smodel->add_npbounds(npbname[i],npblow[i],npbhigh[i]);
    }
  return(smodel);
}

/*----------------------------------------------------------------------------*/
csm_model* csm_model::scalesignal(Double_t coefficient) {
  Int_t i;
  csm_channel_model* scmodel;
  csm_model* smodel = new csm_model;

  for (i=0; i< (Int_t) channame.size(); i++) {
    scmodel = chanmodel[i]->scalesignal(coefficient);
    smodel->add_chanmodel(scmodel,channame[i]);
    delete scmodel;
  }
  for (i=0;i < (Int_t) npcm.size(); i++) {
    smodel->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
  }
  for (i=0;i < (Int_t) npbname.size(); i++) {
    smodel->add_npbounds(npbname[i],npblow[i],npbhigh[i]);
  }
  return(smodel);
}

/*----------------------------------------------------------------------------*/
csm_model* csm_model::scale_err(Double_t coefficient) {
  Int_t i;
  csm_channel_model* scmodel;
  csm_model* smodel = new csm_model;

  for (i=0; i< (Int_t) channame.size(); i++)
    {
      scmodel = chanmodel[i]->scale_err(coefficient);
      smodel->add_chanmodel(scmodel,channame[i]);
      delete scmodel;
    }
  for (i=0;i < (Int_t) npcm.size(); i++)
    {
      smodel->add_npcons(npcm[i].ninput,npcm[i].pnameinput,npcm[i].pnameoutput,npcm[i].f);
    }
  for (i=0;i < (Int_t) npbname.size(); i++)
    {
      smodel->add_npbounds(npbname[i],npblow[i],npbhigh[i]);
    }
  return(smodel);
}

/*----------------------------------------------------------------------------*/
/* Build the list of channel models and channel names that is sorted by channel */
/* name during the building process.  Return the vector index to use to refer */
/* to this particular channel. */

Int_t csm_model::lookup_add_channame(const char* cname) {
  Int_t i,ifound,j,jfound;
  char *s;
  csm_channel_model *cm; 
  vector<char*>::iterator cni;
  vector<csm_channel_model*>::iterator cmi;

  ifound = -1;
  jfound = -1;
  for (i=0; i < (Int_t) channame.size(); i++)
    {
      j = (Int_t) strcmp(cname,channame[i]);
      if (j == 0)
	{
	  ifound = i;
	}
      if (j>0 && jfound == -1)
	{
	  jfound = i;
	}
    }
  /* if the name isn't already in the list, add it to the vector of names and
     make a blank model for it too.  Put the new name in it sorted place, sorted
     by increasing sort order of the name strings */

  if (ifound == -1)
    {
      s = new char[strlen(cname)+1];
      cm = new csm_channel_model;
      strcpy(s,cname);
      if (jfound == -1)
	{
          ifound = channame.size();
          channame.push_back(s);
          chanmodel.push_back(cm);
	}
      else
	{
	  ifound = jfound;
	  cni = channame.begin() + jfound;
	  channame.insert(cni,s);
	  cmi = chanmodel.begin() + jfound;
	  chanmodel.insert(cmi,cm);
	}
    }

  return(ifound);
}

//-----------------------------------------------------------------------------
void csm_model::print() {
  Int_t i,j;

  cout << "csm_model::print  -- printing out model information" << endl;
  for (i=0;i<(Int_t) channame.size();i++)
    {
      cout << "Channel: " << i << " Name: " << channame[i] << endl;
      chanmodel[i]->print();
    }
  for (i=0;i<(Int_t) npcm.size();i++)
    {
      cout << "-------------------" << endl;
      cout << "Constraint equation:  " << npcm[i].pnameoutput << " is computed from " << endl;
      for (j=0;j<npcm[i].ninput;j++)
	{
	  cout << npcm[i].pnameinput[j] << endl;
	}
    }
  cout << "-------------------" << endl;
  for (i=0;i<(Int_t) npbname.size(); i++)
    {
      cout << "NP Bounds.  Name, lowbound highbound: " << npbname[i] << " " << npblow[i] << " " << npbhigh[i] << endl; 
    }
}


//-----------------------------------------------------------------------------
// print out just a piece of a model, indexed by channel name.
void csm_model::print(char *channame) {
  Int_t i = lookup_add_channame(channame);
  cout << "Printing One Channel: " << channame << endl;
  chanmodel[i]->print();
}


/*------------------------------------------------------------------------*/
// and a method to allow an object of type csm_model to plot up one of its
// channels with some data compared.

void csm_model::plotwithdata(const char* cname, const TH1* dh)
{
  Int_t i;
  for (i=0;i<(Int_t)channame.size();i++)
    {
      if (strcmp(cname,channame[i])==0)
	{
	  chanmodel[i]->plotwithdata(dh);
	}
    }
}
 
/*------------------------------------------------------------------------*/

// check candidates

void csm_model::candcheck(const char* cname, TH1* dh)
{
  Int_t i;
  for (i=0;i<(Int_t)channame.size();i++)
    {
      if (strcmp(cname,channame[i])==0)
	{
	  chanmodel[i]->candcheck(dh);
	}
    }
}
 
/*------------------------------------------------------------------------*/

double csm_model::kstest(const char* cname, TH1* dh)
{
  Int_t i;
  double ksresult=0; 
  for (i=0;i<(Int_t)channame.size();i++)
    {
      if (strcmp(cname,channame[i])==0)
	{
	  ksresult = chanmodel[i]->kstest(dh);
	}
    }
  return(ksresult);
}
 
/*------------------------------------------------------------------------*/
double csm_model::kstest_px(const char* cname, TH1* dh) {
  Int_t i;
  double ksresult=0;
  for (i=0;i<(Int_t)channame.size();i++)
    {
      if (strcmp(cname,channame[i])==0)
	{
	  ksresult = chanmodel[i]->kstest_px(dh);
	}
    }
  return(ksresult);
}
 
//-----------------------------------------------------------------------------
Double_t csm_model::chisquared1(const TH1 **dh) {
  Int_t i;
  Double_t cs;
  cs = 0;
  for (i=0;i<(Int_t)chanmodel.size();i++) {
    cs += chanmodel[i]->chisquared1(dh[i]);
  }
  return(cs);
}

/*------------------------------------------------------------------------*/
//  Set the interpolation style for a particular channel.  Two methods    
//  one for channel models, and one if you just have a pointer to a csm_model
/*------------------------------------------------------------------------*/
void csm_model::set_interpolation_style(const char *cname, INTERPSTYLE istyle) {
  Int_t i;
  for (i=0; i<(Int_t) channame.size(); i++) {
    if (strcmp(channame[i],cname)==0) {
      chanmodel[i]->set_interpolation_style(istyle);
    }
  }
}

