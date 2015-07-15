///////////////////////////////////////////////////////////////////////////////
#include "mclimit/csm_channel_model.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <limits>
#include <iostream>
#include <stddef.h>
#include <algorithm>
#include <assert.h>

using namespace std;

#include "THStack.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "TMatrixT.h"

ClassImp(csm_channel_model)

/*----------------------------------------------------------------------------*/
/* A channel model is a sum of template histograms along with systematic errors */
// constructor
//-----------------------------------------------------------------------------
csm_channel_model::csm_channel_model(const char* Name): TNamed(Name,Name) {
  chan_istyle = CSM_INTERP_HORIZONTAL;  //  defaults to csm_pvmorph interpolation
}

/*----------------------------------------------------------------------------*/
// destructor
//-----------------------------------------------------------------------------
csm_channel_model::~csm_channel_model() {
  Int_t i;

  //cout << "Called csm_channel_model destructor" << endl;

  // deallocate memory used to save sytematic error names

  for (i=0;i < (Int_t) syserr.size();i++)
    {
      delete[] syserr[i].sysname;
    }
  // deallocate cloned input histograms
  for (i=0;i < (Int_t) histotemplate.size();i++)
    {
      delete histotemplate[i];
      delete histotemplate_varied[i];
    }
  for (i=0;i < (Int_t) syserr.size();i++)
    {
      if (syserr[i].lowshape !=0)
	{
	  delete syserr[i].lowshape;
	}
      if (syserr[i].highshape !=0)
	{
	  delete syserr[i].highshape;
	}
    }
}


void csm_channel_model::print()
{
  Int_t i,j;
  Double_t ssum = 0;
  Double_t bsum = 0;
  Double_t ssumv = 0;
  Double_t bsumv = 0;
  Double_t central_integral=0;

  cout << "Begin-----------------csm_channel_model::print()------------" << endl;

  for(i=0;i < (Int_t) histotemplate.size();i++)
    {
      cout << endl;
      cout << "Template " << i << endl;
      cout << "  Histogram name: " << histotemplate[i]->GetName();
      cout << "  Histogram title: " << histotemplate[i]->GetTitle();
      cout << "  sft: " << sft[i] << endl;
      cout << "  sft_varied: " << sft_varied[i] << endl;
      cout << "  poissflag: " << poissflag[i] << endl;
      cout << "  signalflag: " << scaleflag[i] << endl;
      central_integral = histotemplate[i]->Integral();
      cout << "  Integral: " << central_integral << endl;
      cout << "  Scaled Integral: " << histotemplate[i]->Integral()*sft[i] << endl;
      cout << "  Scaled Integral with all syst: " << histotemplate_varied[i]->Integral()*sft_varied[i] << endl;
      cout << "  Template bbeta: " << (histotemplate_varied[i]->Integral()*sft_varied[i])/(histotemplate[i]->Integral()*sft[i]) << endl;
      if (scaleflag[i]) 
	{
	  ssum += histotemplate[i]->Integral()*sft[i];
	  ssumv += histotemplate_varied[i]->Integral()*sft_varied[i];
	}
      else
	{
	  bsum += histotemplate[i]->Integral()*sft[i];
	  bsumv += histotemplate_varied[i]->Integral()*sft_varied[i];
	}

      //histotemplate[i]->Print("all");
      Double_t errtotup = 0;
      Double_t errtotdown = 0;
      Double_t uperrloc,downerrloc,hsi,hsl,fracerr;

      for (j=0;j< (Int_t) syserr.size();j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      cout << "Syst: " << syserr[j].sysname << endl;
	      uperrloc = syserr[j].sysfrach;
	      downerrloc = syserr[j].sysfracl;
	      cout << "  Up rate error: " << uperrloc << endl;
	      cout << "  Down rate error: " << downerrloc << endl;

	      if (syserr[j].highshape != 0)
		{
		  hsi = syserr[j].highshape->Integral();
		  cout << "  Up shape error provided sigma: " << syserr[j].xsighigh << 
                          " integral: " << hsi << endl;
		  fracerr = (hsi-central_integral)/central_integral;
		  cout << "    Fractional error due to the shape: " << fracerr << endl;
		  uperrloc += fracerr;
		  cout << "    Total up rate error incl. shape: " << uperrloc << endl;
		  //syserr[j].highshape->Print("all");
		}
	      if (syserr[j].lowshape != 0)
		{
		  hsl = syserr[j].lowshape->Integral();
		  cout << "  Down shape error provided sigma: " << syserr[j].xsiglow << 
                          " integral: " << hsl << endl;

		  fracerr = (hsl-central_integral)/central_integral;
		  cout << "    Fractional error due to the shape: " << fracerr << endl;
		  downerrloc += fracerr;
		  cout << "    Total down rate error incl. shape: " << downerrloc << endl;

		  //syserr[j].lowshape->Print("all");
		}

	      errtotup += uperrloc*uperrloc;
	      errtotdown += downerrloc*downerrloc;

	    }
	}
      errtotup = sqrt(errtotup);
      errtotdown = sqrt(errtotdown);
      cout << "Total relative error on this template (up): " << errtotup << endl;
      cout << "Total relative error on this template (down): " << errtotdown << endl;
    }

  cout << "Total signal this channel in this model: " << ssum << endl;
  cout << "Total background this channel in this model: " << bsum << endl;
  cout << "Syst. Varied Total signal this channel in this model: " << ssumv << endl;
  cout << "Syst. Varied Total background this channel in this model: " << bsumv << endl;
  cout << "End-------------------csm_channel_model::print()------------" << endl;
}
/*----------------------------------------------------------------------------*/
void csm_channel_model::add_template(TH1 *template_hist, 
				     Double_t sf, 
				     Int_t nnp, 
				     const char* npname[],
				     Double_t *nps_low,  
				     Double_t *nps_high, 
				     TH1 *lowshape[],   
				     Double_t *lowsigma,  
				     TH1 *highshape[],   
				     Double_t *highsigma, 
				     Int_t pflag,         
				     Int_t sflag)
{
//-----------------------------------------------------------------------------
// template_hist: Poisson or non-Poisson histogram
// sf           : scale factor to multiply by to compare with the data for a MC Poisson histogram
// nnp          : number of nuisance parameters (gaussian with unit width)
// npname[]     : nuisance parameter names
// nps_low[]    : fractional uncertainty on SF due to each nuisance parameter - low side
// nps_high[]   : fractional uncertainty on SF due to each nuisance parameter - high side
//
//                typically nps_low and nps_high are input with opposite signs -- if opposite
//                variations of the nuisance parameter create opposite changes in sf.  The sign
//                is retained in the calculation in case both variations of a nuisance parameter
//                shift the normalization in the same way (either both + or both -)
// lowshape[]   : array of low hisogram shapes, one for each nuisance param (null if no shape error)
// lowsigma     : number of sigma low for each nuisance parameter shape variation
// highshape[]  : array of high histogram shapes, one for each nuisance param (null if no shape error)
// highsigma    : number of sigma high for each shape variation
// pflag        : Poisson flag -- 1 if Poisson, 0 of not.  2 if Gaussian error from the histo contents
// sflag        : scale flag -- 1 if signal, 0 if background (for use with s95 calculator)
//-----------------------------------------------------------------------------
  int i;
  svstruct ses;
  char *s;

  sft.push_back(sf);
  sft_varied.push_back(sf);
  poissflag.push_back(pflag);
  scaleflag.push_back(sflag);
  for (i=0;i<nnp;i++) {
    s = new char[strlen(npname[i])+1];
    strcpy(s,npname[i]);
    ses.sysname = s;
    ses.itemplate = histotemplate.size();
    ses.sysfracl = nps_low[i];
    ses.sysfrach = nps_high[i];
    if (lowshape[i] !=0) {
      ses.lowshape = (TH1*) lowshape[i]->Clone();
    }
    else {
      ses.lowshape = 0;
    }
    if (highshape[i] !=0) {
      ses.highshape = (TH1*) highshape[i]->Clone();
    }
    else {
      ses.highshape = 0;
    }
    ses.xsiglow = lowsigma[i];
    ses.xsighigh = highsigma[i];
    syserr.push_back(ses);

    if (ses.highshape != 0) {
      if (ses.highshape->GetNbinsX() != template_hist->GetNbinsX()) {
	cout << "Chisquared minmization:  histo template and high shape have different bin counts." << endl;
	cout << template_hist->GetNbinsX() << " != " << ses.highshape->GetNbinsX() << endl;
	exit(0);
      }
    }
    if (ses.lowshape != 0) {
      if (ses.lowshape->GetNbinsX() != template_hist->GetNbinsX()) {
	cout << "Chisquared minmization:  histo template and low shape have different bin counts." << endl;
	cout <<  template_hist->GetNbinsX() << " != " << ses.lowshape->GetNbinsX() << endl;
	exit(0);
      }
    }
  }
  histotemplate.push_back((TH1*) template_hist->Clone());
  histotemplate_varied.push_back((TH1*) template_hist->Clone());
  //cout << "model::add_template: " << histotemplate.size() << endl;
  //gDirectory->ls();
}

//-----------------------------------------------------------------------------
// make a copy of this model by adding the templates over again.
// that way the clone can be deleted by itself, and the destructor
// won't try to delete allocated memory twice
//-----------------------------------------------------------------------------
csm_channel_model* csm_channel_model::Clone() {
  Int_t     i,j,nnp,ntemplates,nsys;
  Double_t *nps_low   = new Double_t[syserr.size()];
  Double_t *nps_high  = new Double_t[syserr.size()];
  Double_t *lowsigma  = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];

  TH1 **lowshape      = new TH1 *[syserr.size()];
  TH1 **highshape     = new TH1 *[syserr.size()];
  const char **ename  = new const char *[syserr.size()];

  csm_channel_model* mclone = new csm_channel_model(this->GetName());

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i < ntemplates;i++) {
    nnp = 0;
    for (j=0;j < nsys;j++) {
      if (syserr[j].itemplate == i) {
	nps_low[nnp] = syserr[j].sysfracl;
	nps_high[nnp] = syserr[j].sysfrach;
	lowshape[nnp] = syserr[j].lowshape;
	highshape[nnp] = syserr[j].highshape;
	lowsigma[nnp] = syserr[j].xsiglow;
	highsigma[nnp] = syserr[j].xsighigh;
	ename[nnp] = syserr[j].sysname;
	nnp++;
      }
    }
      //cout << "Model clone adding template " << i << endl;
    mclone->add_template(histotemplate[i],sft[i],nnp,ename,
			 nps_low,nps_high,lowshape,lowsigma,
			 highshape,highsigma,poissflag[i],scaleflag[i]);
  }

  mclone->chan_istyle = chan_istyle;

  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;

  return mclone;
}

/*----------------------------------------------------------------------------*/

// addition of two models makes a new model with the sum of the templates

csm_channel_model* csm_channel_model::add(csm_channel_model &a)
{
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1 **lowshape = new TH1*[syserr.size()];
  TH1 **highshape = new TH1*[syserr.size()];
  const char **ename = new const char*[syserr.size()];

  csm_channel_model* mclone = a.Clone();

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i < ntemplates;i++)
    {
      nnp = 0;
      for (j=0;j < nsys;j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      nps_low[nnp] = syserr[j].sysfracl;
	      nps_high[nnp] = syserr[j].sysfrach;
	      lowshape[nnp] = syserr[j].lowshape;
	      highshape[nnp] = syserr[j].highshape;
	      lowsigma[nnp] = syserr[j].xsiglow;
	      highsigma[nnp] = syserr[j].xsighigh;
	      ename[nnp] = syserr[j].sysname;
	      nnp++;
	    }
	}
      mclone->add_template(histotemplate[i],sft[i],nnp,ename,
                           nps_low,nps_high,lowshape,lowsigma,
                           highshape,highsigma,poissflag[i],scaleflag[i]);
    }
  mclone->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return(mclone);
}

/*----------------------------------------------------------------------------*/
// multiplication of a model and a scalar
// all templates are scaled up 
//-----------------------------------------------------------------------------
csm_channel_model* csm_channel_model::scale(Double_t coefficient) {
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1 **lowshape = new TH1*[syserr.size()];
  TH1 **highshape = new TH1*[syserr.size()];
  const char **ename = new const char*[syserr.size()];

  csm_channel_model* smodel = new csm_channel_model;

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i < ntemplates;i++) {
    nnp = 0;
    for (j=0;j < nsys;j++) {
      if (syserr[j].itemplate == i) {
	nps_low[nnp] = syserr[j].sysfracl;
	nps_high[nnp] = syserr[j].sysfrach;
	lowshape[nnp] = syserr[j].lowshape;
	highshape[nnp] = syserr[j].highshape;
	lowsigma[nnp] = syserr[j].xsiglow;
	highsigma[nnp] = syserr[j].xsighigh;
	ename[nnp] = syserr[j].sysname;
	nnp++;
      }
    }
    smodel->add_template(histotemplate[i],coefficient*sft[i],nnp,ename,
			 nps_low,nps_high,lowshape,lowsigma,
			 highshape,highsigma,poissflag[i],scaleflag[i]);
  }
  smodel->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return smodel;
}

/*----------------------------------------------------------------------------*/
// multiplication of a model and a scalar -- scale the systematic
// errors down with 1/sqrt(coefficient)
//-----------------------------------------------------------------------------
csm_channel_model* csm_channel_model::scale_err(Double_t coefficient) {
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t escale;
  Double_t *nps_low = new Double_t[syserr.size()];
  Double_t *nps_high = new Double_t[syserr.size()];
  Double_t *lowsigma = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1 **lowshape = new TH1*[syserr.size()];
  TH1 **highshape = new TH1*[syserr.size()];
  const char **ename = new const char*[syserr.size()];

  csm_channel_model* smodel = new csm_channel_model;

  escale = 1.0/sqrt(coefficient);

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0;i< ntemplates;i++)
    {
      nnp = 0;
      for (j=0;j< nsys;j++)
	{
	  if (syserr[j].itemplate == i)
	    {
	      nps_low[nnp] = syserr[j].sysfracl*escale;
	      nps_high[nnp] = syserr[j].sysfrach*escale;
	      lowshape[nnp] = syserr[j].lowshape;
	      highshape[nnp] = syserr[j].highshape;
	      lowsigma[nnp] = syserr[j].xsiglow/escale;
	      highsigma[nnp] = syserr[j].xsighigh/escale;
	      ename[nnp] = syserr[j].sysname;
	      nnp++;
	    }
	}
      smodel->add_template(histotemplate[i],coefficient*sft[i],nnp,ename,
                           nps_low,nps_high,lowshape,lowsigma,
                           highshape,highsigma,poissflag[i],scaleflag[i]);
    }
  smodel->chan_istyle = chan_istyle;
  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;
  return(smodel);
}

//-----------------------------------------------------------------------------
// multiplication of only parts of a model and a scalar
// only signal is scaled up
//-----------------------------------------------------------------------------
csm_channel_model* csm_channel_model::scalesignal(Double_t coefficient) {
  Int_t i,j,nnp,ntemplates,nsys;
  Double_t *nps_low   = new Double_t[syserr.size()];
  Double_t *nps_high  = new Double_t[syserr.size()];
  Double_t *lowsigma  = new Double_t[syserr.size()];
  Double_t *highsigma = new Double_t[syserr.size()];
  TH1 **lowshape      = new TH1*    [syserr.size()];
  TH1 **highshape     = new TH1*    [syserr.size()];
  const char **ename  = new const char*[syserr.size()];
  Double_t sc1;

  csm_channel_model* smodel = new csm_channel_model;

  ntemplates = (Int_t) histotemplate.size();
  nsys = (Int_t) syserr.size();
  for (i=0; i<ntemplates; i++) {
    nnp = 0;
    for (j=0;j < nsys;j++) {
      if (syserr[j].itemplate == i) {
	nps_low[nnp] = syserr[j].sysfracl;
	nps_high[nnp] = syserr[j].sysfrach;
	lowshape[nnp] = syserr[j].lowshape;
	highshape[nnp] = syserr[j].highshape;
	lowsigma[nnp] = syserr[j].xsiglow;
	highsigma[nnp] = syserr[j].xsighigh;
	ename[nnp] = syserr[j].sysname;
	nnp++;
      }
    }
    sc1 = sft[i];
    if (scaleflag[i] != 0) {
      sc1 = coefficient*sft[i];
    }
    smodel->add_template(histotemplate[i],sc1,nnp,ename,
			 nps_low,nps_high,lowshape,lowsigma,
			 highshape,highsigma,poissflag[i],scaleflag[i]);
  }
  smodel->chan_istyle = chan_istyle;

  delete[] ename;
  delete[] lowshape;
  delete[] highshape;
  delete[] lowsigma;
  delete[] highsigma;
  delete[] nps_low;
  delete[] nps_high;

  return smodel;
}

/*----------------------------------------------------------------------------*/

// Create a fluctuated channel model -- input a list of nuisance parameter names
// and values, and return a pointer to a new channel model which has 
// responded to those nuisance parameters.  Be sure to delete it when done.
// todo -- make sure that the shape errors accumulate.  Suggestion of John
// Zhou: average all shape error interpolations.  --  Better compounded interpolations
// introduced in Spring 2007

void csm_channel_model::nuisance_response(Int_t nparams,
                                          char *paramname[],
                                          Double_t paramvalue[])
{
  Int_t i,j,itpl,nsys; //,ntemplates;

  /*
  cout << "in channel nuisance response: " << endl;
  for (i=0;i < (Int_t) syserr.size(); i++)
  {
    cout << "error source: " << i << " name: " << syserr[i].sysname << endl;
  }
  for (i=0;i < (Int_t) nparams; i++)
  {
    cout << "param: " << i << " name: " << paramname[i] << endl;
  }
  */

  undo_nuisance_response();
  //  ntemplates = (Int_t) histotemplate.size();

  TH1* hcl = 0;

  nsys = (Int_t) syserr.size();
  for (i=0;i < nsys; i++)
    {
      for (j=0; j<nparams; j++)
	{
	  if (strcmp(syserr[i].sysname,paramname[j]) == 0)
	    {
	      itpl = syserr[i].itemplate;
	      sft_varied[itpl] *= max(1E-8,( 
                  (syserr[i].sysfrach+syserr[i].sysfracl)*paramvalue[j]*paramvalue[j]/2.0 +
                  (syserr[i].sysfrach-syserr[i].sysfracl)*paramvalue[j]/2.0 + 1.0));
	      if (paramvalue[j]>0)
		{
		  if (syserr[i].highshape != 0)
		    {
		      if (hcl == 0)
			{
                           hcl = (TH1*) histotemplate[itpl]->Clone();
			}
		      else
			{
			  hcl->Reset();
			  hcl->Add(histotemplate_varied[itpl],1.0);
			}
		      if (poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{
		          csm_interpolate_histogram2(histotemplate[itpl],0.0,
                                                 syserr[i].highshape,syserr[i].xsighigh,
                                                 hcl,histotemplate_varied[itpl],paramvalue[j],chan_istyle);
			}
		      else
			{
		          csm_interpolate_histogram2_noerr(histotemplate[itpl],0.0,
                                                 syserr[i].highshape,syserr[i].xsighigh,
                                                 hcl,histotemplate_varied[itpl],paramvalue[j],chan_istyle);
			}

		      //cout << "did a +interpolation " << i << " " << j << " param: " << paramvalue[j] << endl;
		    }
		}
	      else
		{
		  if (syserr[i].lowshape != 0)
		    {
		      if (hcl == 0)
			{
                           hcl = (TH1*) histotemplate[itpl]->Clone();
			}
		      else
			{
			  hcl->Reset();
			  hcl->Add(histotemplate_varied[itpl],1.0);
			}
		      if (poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{
		          csm_interpolate_histogram2(histotemplate[itpl],0.0,
                                                syserr[i].lowshape, -syserr[i].xsiglow,
                                                hcl,histotemplate_varied[itpl],-paramvalue[j],chan_istyle);
			}
		      else
			{
		          csm_interpolate_histogram2_noerr(histotemplate[itpl],0.0,
                                                syserr[i].lowshape, -syserr[i].xsiglow,
                                                hcl,histotemplate_varied[itpl],-paramvalue[j],chan_istyle);
			}
		      //cout << "did a -interpolation " << i << " " << j << " param: " << paramvalue[j] <<  endl;
		    }
		}
	    }
	}
    }
  if (hcl != 0)
    {
      delete hcl;
    }
}

/*----------------------------------------------------------------------------*/
// resets all the varied histograms and scales to their unvaried states.

void csm_channel_model::undo_nuisance_response()
{
  Int_t i,ntemplates,nbinsx,nbinsy,ix,iy;

  ntemplates = (Int_t) histotemplate.size();
  for (i=0;i<ntemplates;i++)
    {
      sft_varied[i] = sft[i];
      nbinsx = histotemplate[i]->GetNbinsX();
      nbinsy = histotemplate[i]->GetNbinsY();
      if (nbinsy==1)
	{
          for (ix=1;ix<=nbinsx;ix++)
	    {
	      histotemplate_varied[i]->SetBinContent(ix,histotemplate[i]->GetBinContent(ix));
	      histotemplate_varied[i]->SetBinError(ix,histotemplate[i]->GetBinError(ix));
	    }
	}
      else
	{
          for (ix=1;ix<=nbinsx;ix++)
	    {
	      for (iy=1;iy<=nbinsy;iy++)
		{
	          histotemplate_varied[i]->SetBinContent(ix,iy,histotemplate[i]->GetBinContent(ix,iy));
	          histotemplate_varied[i]->SetBinError(ix,iy,histotemplate[i]->GetBinError(ix,iy));
		}
	    }
	}
    }
}

/*----------------------------------------------------------------------------*/
// check to see if any bin has a total negative prediction in this channel

Int_t csm_channel_model::checkneg()
{
  cout << "csm_channel_model::checkneg() to be written" << endl;
  return(0);
}


/*------------------------------------------------------------------------*/
// make a plot of the results, along with some data
void csm_channel_model::plotwithdata(const TH1* Hist)
{
  Int_t i,ntemplates,nbinsy;
  Double_t stackmax,datamax,plotmax;

  TH1* dh = (TH1*) Hist->Clone();

  THStack *hs = new THStack("hs",dh->GetTitle());
  ntemplates = (Int_t) histotemplate.size();
  TLegend *slegend = (TLegend*) new TLegend(0.7,0.6,0.89,0.89);

  for (i=0;i<ntemplates;i++)
    {
      TH1* htl = (TH1*) histotemplate_varied[i]->Clone();
      htl->Scale(sft_varied[i]);
      htl->SetFillColor(i+40);
      htl->SetFillStyle(1001);
      hs->Add(htl);
    }

  TList *hlist = hs->GetHists();
  TObjLink *lnk = hlist->LastLink();          
  while (lnk)
    {  slegend->AddEntry(lnk->GetObject(),lnk->GetObject()->GetName(),"F");
       lnk = lnk->Prev();                       
    }     
  // make sure the plot is big enough to fit the data, the model stack,
  // and the data error bars with a little room to spare
  stackmax = hs->GetMaximum();
  datamax = dh->GetMaximum();
  nbinsy = dh->GetNbinsY();
  datamax += sqrt(datamax);
  plotmax = max(datamax,stackmax);
  hs->SetMaximum(plotmax);
  dh->SetMaximum(plotmax);
  hs->SetMinimum(0);
  dh->SetMinimum(0);
  if (nbinsy==1)
    {
      hs->Draw("HIST");
      dh->SetMarkerStyle(20);
      dh->SetMarkerColor(kBlack);
      dh->DrawCopy("E0SAME");
    }
  else
    {
      hs->Draw();
      dh->SetMarkerStyle(20);
      dh->SetMarkerColor(kBlack);
      dh->DrawCopy("LEGO,SAME");
    }
  slegend->AddEntry(dh,dh->GetName(),"P");
  slegend->SetHeader(dh->GetTitle());
  slegend->Draw();
					// memory cleanup
  delete dh;
}

void csm_channel_model::candcheck(TH1 *dh)
{
  cout << dh->GetTitle() << " Candidate Check " << endl;

  Double_t sumsb = 0;
  Double_t ssum = 0;
  Double_t bsum = 0;
  Double_t dsum = 0;
  Int_t nbinsx = histotemplate[0]->GetNbinsX();
  Int_t nbinsy = histotemplate[0]->GetNbinsY();
  if (nbinsx != dh->GetNbinsX())
    {
      cout << "data histogram and model template have different numbers of x bins: " <<
	nbinsx << " != " << dh->GetNbinsX() << endl;
      return;
    }
  if (nbinsy != dh->GetNbinsY())
    {
      cout << "data histogram and model template have different numbers of y bins: " <<
	nbinsy << " != " << dh->GetNbinsY() << endl;
      return;
    }

  Int_t ibinx,ibiny;

  // Double_t hcb[nbinsx][nbinsy];
  // Double_t hcs[nbinsx][nbinsy];

  TMatrixT<double> hcb(nbinsy,nbinsx);
  TMatrixT<double> hcs(nbinsy,nbinsx);

  for (ibinx=0;ibinx<nbinsx;ibinx++)
    {
      for (ibiny=0;ibiny<nbinsy;ibiny++)
	{
	  hcb[ibinx][ibiny] = 0;
	  hcs[ibinx][ibiny] = 0;
	}
    }

  for(Int_t ic=0;ic < (Int_t) histotemplate.size();ic++)
    {
      for (ibinx=0;ibinx<nbinsx;ibinx++)
	{
	  for (ibiny=0;ibiny<nbinsy;ibiny++)
	    {
	      if (scaleflag[ic])
		{
		  hcs[ibinx][ibiny] += sft_varied[ic]*histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1);
		  if (hcs[ibinx][ibiny] < 0) 
		    {
		      cout << "Negative signal expectation (" << ibinx+1 << "," << ibiny+1 << "):" << hcs[ibinx][ibiny] << endl;
		    }
		}
	      else
		{
		  hcb[ibinx][ibiny] += sft_varied[ic]*histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1);
		  if (hcb[ibinx][ibiny] < 0) 
		    {
		      cout << "Negative background expectation (" << ibinx+1 << "," << ibiny+1 << "):" << hcb[ibinx][ibiny] << endl;
		    }
		}
	    }
	}
    }

  for (ibinx=0;ibinx<nbinsx;ibinx++)
    {
      for (ibiny=0;ibiny<nbinsy;ibiny++)
	{
	  ssum += hcs[ibinx][ibiny];
	  bsum += hcb[ibinx][ibiny];
          Double_t dc = dh->GetBinContent(ibinx+1,ibiny+1);
	  //cout << "*@*" << hcs[ibinx][ibiny] << " " << hcb[ibinx][ibiny] << " " << dc << endl;
	  dsum += dc;
	  if (hcs[ibinx][ibiny]>0 && hcb[ibinx][ibiny]<=0)
	    {
	      cout << "Null background with expected signal: (" << ibinx+1 << "," << ibiny+1 
                   << ") Cands: " << dc << " Signal: " << hcs[ibinx][ibiny] << endl;
	    }
	  else
	    {
	      if (dc > 0)
		{
		  Double_t sbratio = hcs[ibinx][ibiny]/hcb[ibinx][ibiny];
		  sumsb += dc*sbratio;
		  if (sbratio>0.3)
		    {
		      cout << "High s/b candidate(s): (" << ibinx+1 << "," << ibiny+1 << ") cands: " << dc << " s/b: " << sbratio << 
			" s: " << hcs[ibinx][ibiny] << " b: " << hcb[ibinx][ibiny] << endl;
		    }
		}
	    }
	}
    }
  cout << "S/B sum over all candidates: " << sumsb << endl;
  cout << "S sum over all bins: " << ssum << endl;
  cout << "B sum over all bins: " << bsum << endl;
  cout << "D sum over all bins: " << dsum << endl;
}

double csm_channel_model::kstest(TH1* dh)
{
  Int_t i,ntemplates;
  ntemplates = (Int_t) histotemplate.size();
  double tout;

  TH1* hsum = (TH1*) histotemplate_varied[0]->Clone();
  hsum->Sumw2();
  hsum->Reset();


  for (i=0;i<ntemplates;i++)
    {
      TH1* htl = (TH1*) histotemplate_varied[i];
      hsum->Add(htl,sft_varied[i]);
    }

  tout = hsum->KolmogorovTest(dh);
  delete hsum;
  return(tout);
}

//-----------------------------------------------------------------------------
double csm_channel_model::kstest_px(TH1* dh)
{
  Int_t i,ntemplates;
  ntemplates = (Int_t) histotemplate.size();
  double tout;

  TH1* hsum = (TH1*) histotemplate_varied[0]->Clone();
  hsum->Sumw2();
  hsum->Reset();


  for (i=0;i<ntemplates;i++)
    {
      TH1* htl = (TH1*) histotemplate_varied[i];
      hsum->Add(htl,sft_varied[i]);
    }

  tout = hsum->KolmogorovTest(dh,"X");
  delete hsum;
  return(tout);
}

/*------------------------------------------------------------------------*/

/* chisquared1 evaluates a chisquared function in the style of T. Devlin's CDF 3126, eq's 9 and 10
   The signal is a sum of signal contributions and the background is a sum of
   background contributions.  This chisquared function is meant to be minimized 
   with respect to the free nuisance parameters 

   This version does one 1D or 2D histogram at a time.

   This version allows for multiple sources of signal and multiple sources of background,
   some of each of which are estimated using finite MC or data statistics in each bin.
   This routine does not distinguish between a signal source and a background source --
   finding the chisquared of a data distribution to a sum of models does not need a distinction
   at this level.  Instead, one may compare the chisquared of the same data against collections
   of models that include signals and those that do not include signals, calling this routine
   twice (or more times).

   CDF 3126 describes how to minimize the chisquared function over each bin's uncertain
   Poisson-constrained rates.  When multiple sources are allowed to be estimated from Poisson
   distributed subsidiary measurements, the quadratic polynomial to be solved for turns
   into a system of quadratic equations which is solved here iteratively.

   This function is meant to be part of a MINUIT minimization over the nuisance
   parameters.


            input: TH1 *dh -- data histogram to compare the channel's model against
            
            output:  chi squared, the function value.

   Update 5 July 2006 -- Reading Barlow and Beeston about bins with zero MC prediction in one
   or more source.  Take the one with the strongest contribution (here taken from the normalization
   scale factors), and set the others to zero when solving the n coupled quadratic equations.

   Update 8 Dec, 2007 -- put in the error bars on the template histograms when the flag is
   CSM_GAUSSIAN_BINERR (new feature) as if they were Poisson (i.e., no different terms in
   the likelihood function -- they're probably Poisson underneath anyhow, but more often, they
   are a complicated mixture of MC or data events with different weights, and all treatments
   of them are approximations).  

*/

Double_t csm_channel_model::chisquared1(const TH1 *dh)
{
  Double_t chi2;
  Int_t ip1,ic,ibinx,ibiny,iter,iprec;
  Double_t A,B,C,D;
  Double_t csum,cpsum,gbc,gbe;
  Int_t nbinsx,nbinsy;
  Double_t nsubs;
  Int_t dtb;  // data observed in a single bin
  Int_t nc;

  // number of template histograms
  nc = (Int_t) histotemplate.size();

  Int_t lpoissflag[nc];  // local Poisson flag -- if error is zero for a template in a bin
                         // reclassify it as no bin error

  // allocate rho1 and rho2 for all templates, even though we're only going to need
  // them for the Poisson-distributed ones
  // push them on the stack -- 8 dec 2007

  Double_t rho1[nc];
  Double_t rho2[nc];
  Int_t zlist[nc];
  Double_t sfgp[nc];  // scale factor needed to approximate Gaussian errors as Poisson

  nbinsx = dh->GetNbinsX();
  nbinsy = dh->GetNbinsY();

  chi2 = 0;

  for (ibinx=0;ibinx<nbinsx;ibinx++) {
    for (ibiny=0;ibiny<nbinsy;ibiny++) {
//-----------------------------------------------------------------------------
// bin content of the data histogram
//-----------------------------------------------------------------------------
      if (nbinsy==1) dtb = (Int_t) dh->GetBinContent(ibinx+1);
      else           dtb = (Int_t) dh->GetBinContent(ibinx+1,ibiny+1);

      //cout << "In chi2calc: " << ibinx << " " << ibiny << " " << dtb << endl;

      // if we are told to pay attention to the error but it's zero, reclassify it locally
      // as a zero-error bin in this template

      for (ic=0; ic<nc; ic++) {
	sfgp[ic] = 1.0;
	lpoissflag[ic] = poissflag[ic];
	if (nbinsy == 1) { 
	  gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1); 
	  gbe = histotemplate_varied[ic]->GetBinError(ibinx+1); 
	}
	else {
	  gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); 
	  gbe = histotemplate_varied[ic]->GetBinError(ibinx+1,ibiny+1); 
	}
	if (poissflag[ic] == CSM_GAUSSIAN_BINERR && gbe == 0) { 
	  lpoissflag[ic] = CSM_NOBINERR; 
	}
	if (lpoissflag[ic] == CSM_GAUSSIAN_BINERR) {
	  sfgp[ic] = gbc/(gbe*gbe);
	  lpoissflag[ic] = CSM_POISSON_BINERR;
	}
      }
//-----------------------------------------------------------------------------
// the sum of zero-bin-error contributions, varied by the nuisance parameters 
//-----------------------------------------------------------------------------
      csum = 0;
      for (ic=0;ic<nc;ic++) {
	if (lpoissflag[ic] == CSM_NOBINERR) {
	  if (nbinsy == 1) gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1);
	  else             gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1);
	  csum += gbc*sft_varied[ic];
	}
      }
	 
      if (csum < 0) { 
	chi2 += 1E10; 
      }

      /* solve for the rho's in each bin for each source of Poisson-estimated model rate
	 rho1 is the current estimate used to compute the rho2's.  On each iteration,
	 copy the previous iteration's rho2 into the rho1 array and re-solve for rho2.
	 start with nominal central values from the subsidiary measurements */

      int    haszero =  0;
      int    im1     = -1;
      double xm1     =  0;

      for (ic=0;ic<nc;ic++) { 
	rho1[ic] = 0;
	rho2[ic] = 0;
	if (lpoissflag[ic] == CSM_POISSON_BINERR) {

	  if (nbinsy == 1) gbc = max(0.,histotemplate_varied[ic]->GetBinContent(ibinx+1)); 
	  else             gbc = max(0.,histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1)); 

	  rho2[ic] = gbc*sft_varied[ic];
	  rho1[ic] = rho2[ic];
	  if (gbc == 0 || sft_varied[ic] == 0) { 
	    haszero   = 1;
	    zlist[ic] = 1;
	    if (sft_varied[ic] > xm1) { 
	      xm1 = sft_varied[ic]; 
	      im1 = ic;
	    }
	  }
	  else {
	    zlist[ic] = 0;
	  }
	}
      }

      if (haszero != 0 && im1 > -1) {
	zlist[im1] = 0;
      }

      for (iter=0;iter<CSM_MAXITER;iter++) {
	for (ic=0;ic<nc;ic++) { 
	  if (lpoissflag[ic] == CSM_POISSON_BINERR) {
	    rho1[ic] = rho2[ic];
	    if (zlist[ic] == 1) { 
	      rho1[ic] = 0;
	    }
	  }
	}
	  
	for (ic=0; ic<nc; ic++) {
	  if (lpoissflag[ic] == CSM_POISSON_BINERR) {
	    if (zlist[ic] == 0) {
	      D = csum;
	      for (ip1=0;ip1<nc;ip1++) {
		if ( (lpoissflag[ip1]==CSM_POISSON_BINERR) && ip1 != ic) D += rho1[ip1]; 
	      }

	      A = 1.0 + sfgp[ic]/sft_varied[ic];

	      if (nbinsy == 1) gbc = sfgp[ic]*max(0.,histotemplate_varied[ic]->GetBinContent(ibinx+1)); 
	      else             gbc = sfgp[ic]*max(0.,histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1)); 

	      B        = A*D - dtb - gbc;
	      C        = -gbc*D;
	      rho2[ic] = (-B + sqrt(B*B - 4.0*A*C))/(2.0*A);
	      //cout << "ABC: " << A << " " << B << " " << C << endl;
	    }
	    else { 
//-----------------------------------------------------------------------------
// a la Barlow and Beeston, set only one prediction to nonzero if we have 
// zero MC -- the "strongest" one among all the contributions with zero MC prediction
//-----------------------------------------------------------------------------
	      rho2[ic] = 0;
	    }
	  }
	}

	iprec = 0;

	for (ic=0; ic<nc; ic++) {
	  if (lpoissflag[ic] == CSM_POISSON_BINERR) {
	    if (fabs(rho1[ic]) < PREC1) { 
	      if (fabs(rho2[ic]-rho1[ic]) > PREC1) { 
		iprec = 1;
		break;
	      }
	    }
	    else {
	      if ( fabs((rho2[ic]-rho1[ic])/rho1[ic])>PREC1 ) {
		iprec = 1;
		break;
	      }
	    }
	  }
	}

	if (iprec == 0) break;
      }  /* end loop over iterations to compute the rho's.  rho2 is the computed array */

      /*
	if (CSM_DEBUGPRINT >0 && iprec ==1)
	{
	// cout << "csm_chisquared1: iterations failed to converge " << endl;
	cout << "In chi2calc: " << ibinx << " " << ibiny << " " << dtb << endl;
	for (ic=0;ic<nc;ic++)
	{
	if (nbinsy == 1)
	{ gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1); }
	else
	{ gbc = histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); }
	if (lpoissflag[ic] == CSM_POISSON_BINERR)
	{
	cout << "Poisson contrib " << ic << " " << gbc << " " << sft_varied[ic] << endl;
	}
	else
	{
	cout << "Non-poisson contrib " << ic << " " << gbc << " " << sft_varied[ic] << endl;
	}
	}
	}
      */
//-----------------------------------------------------------------------------
// When the iterations fail to converge, it is usually an oscillatory
// solution.  Pick the rho1 or the rho2 array which minimizes chisquare
// first compute the chisquare using the rho2 array, and if we need to,
// redo it with the rho1 array, and pick the smaller of the two.
//-----------------------------------------------------------------------------
      Double_t chi2a = chi2;

      cpsum = csum;
					// in this loop cpsum will get incremented only 
					// in case of the poisson errors. what is csum ?
      for (ic=0;ic<nc;ic++) {
	if (lpoissflag[ic] == CSM_POISSON_BINERR) {
	  cpsum += rho2[ic];
	}
      }

      // 2010-11-16 P.Murat	  if (cpsum > 0) {
      if (cpsum > 0) {
	chi2a += cpsum;
	chi2a -= dtb;
	if (dtb > 0) { 
	  chi2a -= dtb*log(cpsum/((Double_t) dtb));
	}
      }
      else if (dtb > 0) { 
//-----------------------------------------------------------------------------
// dtb > 0 but cpsum <= 0 - what does this mean? should be some inconsistency
//-----------------------------------------------------------------------------
	chi2a += 1E10; 
      }

      for (ic=0;ic<nc;ic++) {
	if ( (lpoissflag[ic] == CSM_POISSON_BINERR) && sft_varied[ic] > 0) {
	  if (nbinsy == 1) {
	    nsubs = sfgp[ic]*histotemplate_varied[ic]->GetBinContent(ibinx+1); 
	  }
	  else {
	    nsubs = sfgp[ic]*histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); 
	  }
	  
#ifdef DEBUGPRINTC2

	  Double_t c2cont = ( rho2[ic]*sfgp[ic]/sft_varied[ic] - nsubs);
	  if (nsubs > 0 && rho2[ic] > 0) {
	    c2cont -= ((Double_t) nsubs)*log(rho2[ic]*sfgp[ic]/(sft_varied[ic]*((Double_t) nsubs)));
	  }
	  if (c2cont > 0.0) {
	    //  if (ibinx == 9 && ibiny == 4)
	    cout << "in chi2calc: " << ibinx << " " << ibiny << " " << ic << " " << 
	      rho2[ic] << " " << sft_varied[ic] << " " << nsubs << " " << c2cont;
	    if (c2cont>0.01) {
	      cout << "*" << endl;
	    }
	    else {
	      cout << endl;
	    }
	  }
#endif
	  chi2a += (rho2[ic]*sfgp[ic]/sft_varied[ic] - nsubs);
	  if (nsubs > 0 && rho2[ic] > 0) {
	    chi2a -= ((Double_t) nsubs)*log(rho2[ic]*sfgp[ic]/(sft_varied[ic]*((Double_t) nsubs)));
	  }
	}
      }

      Double_t chi2b = chi2;

      if (iprec == 1) {
	cpsum = csum;
	for (ic=0;ic<nc;ic++) {
	  if (lpoissflag[ic] == CSM_POISSON_BINERR) {
	    cpsum += rho1[ic];
	  }
	}

	if (cpsum > 0) {
	  chi2b += cpsum;
	  chi2b -= dtb;
	  if (dtb>0) {
	    chi2b -= dtb*log(cpsum/((Double_t) dtb));
	  }
	}
	else if (dtb > 0) { 
	  // again, cpsum <= 0 for non-zero data bin - what does it mean?
	  chi2b += 1E10; 
	}

	for (ic=0;ic<nc;ic++) {
	  if ( (lpoissflag[ic] == CSM_POISSON_BINERR) && sft_varied[ic] > 0) {
	    if (nbinsy == 1) { 
	      nsubs = sfgp[ic]*histotemplate_varied[ic]->GetBinContent(ibinx+1); 
	    }
	    else { 
	      nsubs = sfgp[ic]*histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); 
	    }

	    /*
	      cout << "in chi2calc: " << ibinx << " " << ibiny << " " << ic << " " << 
	      rho1[ic] << " " << sft_varied[ic] << " " << nsubs << " " sfgp[ic] << endl;
	    */
	    chi2b += (rho1[ic]*sfgp[ic]/sft_varied[ic] - nsubs);
	    if (nsubs > 0 && rho1[ic] > 0) {
	      chi2b -= ((Double_t) nsubs)*log(rho1[ic]*sfgp[ic]/(sft_varied[ic]*((Double_t) nsubs)));
	    }
	  }
	}
      }
      else {
	chi2b = chi2a;
      }

      chi2 = min(chi2a,chi2b);
    } /* end loop over binsy */
  }   /* end loop over binsx */

  chi2 *= 2.0;

  //cout << "chisquared calc: " << chi2 << endl;
  //if ( !(chi2<0 || chi2>=0))
  //  {
  //    cout << "Bad chi2: " << endl;
  //    exit(0);
  //  }

  return chi2;

}

//-----------------------------------------------------------------------------
void csm_channel_model::set_interpolation_style(INTERPSTYLE istyle) {
  chan_istyle = istyle;
}



/*------------------------------------------------------------------------*/

// interpolate 1D histograms and 2D histograms
// histo a corresponds to parameter xa, histo b corresponds to xb.
// xc is input, and histogram c is the interpolated output

// new version -- rely on the more general interpolator with three inputs, but reduce the argument
// count for backward compatibility

// csm_interpolate_histogram interpolates the bin contents and errors

void csm_interpolate_histogram(TH1* a, Double_t xa, 
                               TH1* b, Double_t xb,
                               TH1* c, Double_t xc,
                               INTERPSTYLE istyle)
{
  csm_interpolate_histogram2(a,xa,b,xb,a,c,xc,istyle);
}

// csm_interpolate_histogram_noerr interpolates just the bin contents but not the errors

void csm_interpolate_histogram_noerr(TH1* a, Double_t xa, 
                               TH1* b, Double_t xb,
                               TH1* c, Double_t xc,
                               INTERPSTYLE istyle)
{
  csm_interpolate_histogram2_noerr(a,xa,b,xb,a,c,xc,istyle);
}

// interpolate 1D histograms and 2D histograms
// histo a corresponds to parameter xa, histo b corresponds to xb.
// xc is input, and histogram c is the interpolated output
// d is the histogram to apply the shift given by a and b to, for compounded interpolations.

// approximate attempt to interpolate the uncertainties too.  Problem is, an interpolated
// histogram is a long sum of pieces interpolated from the same central value histogram,
// and thus the errors are correlated in interesting ways.
// A subterfuge -- jut linearly interpolate the errors in the same way that the
// bin contents are linearly interpolated.  It's not a full error propagation.  Halfway interpolations
// really are averages of statistically uncertain histograms, and thus the error on the average should
// be a bit better than the error on either one.  But itnterpolate again, and correlations have to be
// taken into account to do it right.
// we've also lost at this point whether the errors need to be interpolated, but let's
// do them for all histograms.
// speedup 9 Dec 2007 -- avoid cloning TH1's as this is slow

void csm_interpolate_histogram2(TH1* a, Double_t xa, 
                                TH1* b, Double_t xb,
				TH1* d,
                                TH1* c, Double_t xc,
                                INTERPSTYLE istyle)
{
  Int_t i,j;
  Double_t xtmp;
  Int_t nbinsa = a->GetNbinsX();
  Int_t nbinsb = b->GetNbinsX();
  Int_t nbinsc = c->GetNbinsX();
  Int_t nbinsd = d->GetNbinsX();
  Int_t nbinsya = a->GetNbinsY();
  Int_t nbinsyb = b->GetNbinsY();
  Int_t nbinsyc = c->GetNbinsY();
  Int_t nbinsyd = d->GetNbinsY();

  if (nbinsa != nbinsb)
    {
      cout << "nbins mismatch1 in csm_interpolate_histogram2: " << nbinsa << " " << nbinsb << endl;
    }
  if (nbinsb != nbinsc)
    {
      cout << "nbins mismatch2 in csm_interpolate_histogram2: " << nbinsb << " " << nbinsc << endl;
    }
  if (nbinsc != nbinsd)
    {
      cout << "nbins mismatch3 in csm_interpolate_histogram2: " << nbinsc << " " << nbinsd << endl;
    }
  if (nbinsya != nbinsyb)
    {
      cout << "nbinsy mismatch1 in csm_interpolate_histogram2: " << nbinsya << " " << nbinsyb << endl;
    }
  if (nbinsyb != nbinsyc)
    {
      cout << "nbinsy mismatch2 in csm_interpolate_histogram2 " << nbinsyb << " " << nbinsyc << endl;
    }
  if (nbinsyc != nbinsyd)
    {
      cout << "nbinsy mismatch3 in csm_interpolate_histogram2: " << nbinsyc << " " << nbinsyd << endl;
    }

  if (xb == xa)
    {
      cout << "xb == xa in csm_interpolate_histogram2 " << xa << endl;
      cout << "fatal error -- exiting." << endl;
      exit(0);
    }

  // interpolate contents

  csm_interpolate_histogram3(a,xa,b,xb,d,c,xc,istyle);

  // swap errors and contents and interpolate again  (approximate method for evaluating
  // errors on interpolated histograms)
  // be careful to swap only once, even if some pointers are repeated.

  if (nbinsya == 1)
    {
      for (i=1;i<=nbinsa;i++)
	{
	  xtmp = a->GetBinContent(i);
	  a->SetBinContent(i,a->GetBinError(i));
	  a->SetBinError(i,xtmp);
	  if (a != b)
	    {
	      xtmp = b->GetBinContent(i);
	      b->SetBinContent(i,b->GetBinError(i));
	      b->SetBinError(i,xtmp);
	    }
	  // c is the output histogram -- hopefully it is not the same as one of the input histograms
          xtmp = c->GetBinContent(i);
          // c->SetBinContent(i,c->GetBinError(i));
	  c->SetBinError(i,xtmp);

	  if (a != d && b != d)
	    {
	      xtmp = d->GetBinContent(i);
	      d->SetBinContent(i,d->GetBinError(i));
	      d->SetBinError(i,xtmp);
	    }
	}
    }
  else
    {
      for (i=1;i<=nbinsa;i++)
	{
	  for (j=1;j<=nbinsya;j++)
	    {
	       xtmp = a->GetBinContent(i,j);
	       a->SetBinContent(i,j,a->GetBinError(i,j));
	       a->SetBinError(i,j,xtmp);
	       if (a != b)
		 {
	           xtmp = b->GetBinContent(i,j);
	           b->SetBinContent(i,j,b->GetBinError(i,j));
	           b->SetBinError(i,j,xtmp);
		 }
	       xtmp = c->GetBinContent(i,j);
	       //c->SetBinContent(i,j,c->GetBinError(i,j));
	       c->SetBinError(i,j,xtmp);
	       if (a != d && b != d)
		 {
	           xtmp = d->GetBinContent(i,j);
	           d->SetBinContent(i,j,d->GetBinError(i,j));
	           d->SetBinError(i,j,xtmp);
		 }
	    }
	}
    }

  // interpolate the errors now and swap them back -- put the
  // original histograms back together again too

  csm_interpolate_histogram3(a,xa,b,xb,d,c,xc,istyle);

  if (nbinsya == 1)
    {
      for (i=1;i<=nbinsa;i++)
	{
	  xtmp = a->GetBinContent(i);
	  a->SetBinContent(i,a->GetBinError(i));
	  a->SetBinError(i,xtmp);
	  if (a != b)
	    {
	      xtmp = b->GetBinContent(i);
	      b->SetBinContent(i,b->GetBinError(i));
	      b->SetBinError(i,xtmp);
	    }
	  xtmp = c->GetBinContent(i);
	  c->SetBinContent(i,c->GetBinError(i));
	  c->SetBinError(i,xtmp);
	  if (a != d && b != d)
	    {
	      xtmp = d->GetBinContent(i);
	      d->SetBinContent(i,d->GetBinError(i));
	      d->SetBinError(i,xtmp);
	    }
	}
    }
  else
    {
      for (i=1;i<=nbinsa;i++)
	{
	  for (j=1;j<=nbinsya;j++)
	    {
	       xtmp = a->GetBinContent(i,j);
	       a->SetBinContent(i,j,a->GetBinError(i,j));
	       a->SetBinError(i,j,xtmp);
	       if (a != b)
		 {
	           xtmp = b->GetBinContent(i,j);
	           b->SetBinContent(i,j,b->GetBinError(i,j));
	           b->SetBinError(i,j,xtmp);
		 }
	       xtmp = c->GetBinContent(i,j);
	       c->SetBinContent(i,j,c->GetBinError(i,j));
	       c->SetBinError(i,j,xtmp);
	       if (a != d && b != d)
		 {
	           xtmp = d->GetBinContent(i,j);
	           d->SetBinContent(i,j,d->GetBinError(i,j));
	           d->SetBinError(i,j,xtmp);
		 }
	    }
	}
    }
}

void csm_interpolate_histogram2_noerr(TH1* a, Double_t xa, 
                                      TH1* b, Double_t xb,
				      TH1* d,
                                      TH1* c, Double_t xc,
                                      INTERPSTYLE istyle)
{
  Int_t nbinsa = a->GetNbinsX();
  Int_t nbinsb = b->GetNbinsX();
  Int_t nbinsc = c->GetNbinsX();
  Int_t nbinsd = d->GetNbinsX();
  Int_t nbinsya = a->GetNbinsY();
  Int_t nbinsyb = b->GetNbinsY();
  Int_t nbinsyc = c->GetNbinsY();
  Int_t nbinsyd = d->GetNbinsY();

  if (nbinsa != nbinsb)
    {
      cout << "nbins mismatch1 in csm_interpolate_histogram2_noerr: " << nbinsa << " " << nbinsb << endl;
    }
  if (nbinsb != nbinsc)
    {
      cout << "nbins mismatch2 in csm_interpolate_histogram2_noerr: " << nbinsb << " " << nbinsc << endl;
    }
  if (nbinsc != nbinsd)
    {
      cout << "nbins mismatch3 in csm_interpolate_histogram2_noerr: " << nbinsc << " " << nbinsd << endl;
    }
  if (nbinsya != nbinsyb)
    {
      cout << "nbinsy mismatch1 in csm_interpolate_histogram2_noerr: " << nbinsya << " " << nbinsyb << endl;
    }
  if (nbinsyb != nbinsyc)
    {
      cout << "nbinsy mismatch2 in csm_interpolate_histogram2_noerr " << nbinsyb << " " << nbinsyc << endl;
    }
  if (nbinsyc != nbinsyd)
    {
      cout << "nbinsy mismatch3 in csm_interpolate_histogram2_noerr: " << nbinsyc << " " << nbinsyd << endl;
    }

  if (xb == xa)
    {
      cout << "xb == xa in csm_interpolate_histogram2_noerr " << xa << endl;
      cout << "fatal error -- exiting." << endl;
      exit(0);
    }

  // interpolate just the bin contents

  csm_interpolate_histogram3(a,xa,b,xb,d,c,xc,istyle);

}


void csm_interpolate_histogram3(TH1* a, Double_t xa, 
                                TH1* b, Double_t xb,
				TH1* d,
                                TH1* c, Double_t xc,
                                INTERPSTYLE istyle)
{
  Double_t hnorma,hnormb,hnormc,hnormd,hnormci;
  Int_t i,j;
  Double_t gbc;

  Int_t nbinsa = a->GetNbinsX();
  Int_t nbinsb = b->GetNbinsX();
  Int_t nbinsc = c->GetNbinsX();
  Int_t nbinsd = d->GetNbinsX();
  Int_t nbinsya = a->GetNbinsY();
  Int_t nbinsyb = b->GetNbinsY();
  Int_t nbinsyc = c->GetNbinsY();
  Int_t nbinsyd = d->GetNbinsY();

  if (a->Integral()<=0 || b->Integral()<=0)
    { 
      for (i=1;i<=nbinsc;i++)
	{
	   for (j=1;j<=nbinsyc;j++)
	     { c->SetBinContent(i,j,0);
	     }
	} 
      //c->Reset();
      return;
    }
    
  if (nbinsya == 1)
    {
      Double_t *dista = new Double_t[nbinsa];
      Double_t *distb = new Double_t[nbinsb];
      Double_t *distc = new Double_t[nbinsc];
      Double_t *distd = new Double_t[nbinsd];

      hnorma = 0;
      hnormb = 0;
      hnormd = 0;
      for (i=0;i<nbinsa;i++)
        { dista[i] = a->GetBinContent(i+1); hnorma += dista[i]; }
      for (i=0;i<nbinsb;i++)
        { distb[i] = b->GetBinContent(i+1); hnormb += distb[i]; }
      for (i=0;i<nbinsd;i++)
        { distd[i] = d->GetBinContent(i+1); hnormd += distd[i]; }

      hnormc = hnorma + (xc-xa)*(hnormb-hnorma)/(xb-xa);
      // linearly interpolate the normalization between the central value and
      // the varied template.
      hnormc = hnormd*(hnormc/hnorma); // scale the normalization with the new template

      if (istyle == CSM_INTERP_HORIZONTAL || istyle == CSM_INTERP_HORIZONTAL_EXTRAP)
	{
           csm_pvmc(nbinsa,dista,distb,distd,distc,xa,xb,xc);
           hnormci = 0;
           for (i=0;i<nbinsc;i++) { hnormci += distc[i]; }

           for (i=0;i<nbinsc;i++)
           {
             c->SetBinContent(i+1,distc[i]*hnormc/hnormci);
           }
	}
      else if (istyle == CSM_INTERP_VERTICAL || istyle == CSM_INTERP_VERTICAL_EXTRAP)
	{
	  for (i=0;i<nbinsa;i++)
	    {
	      gbc = distd[i] + ((xc-xa)/(xb-xa))*(distb[i]-dista[i]);
	      if (gbc < 0) 
		{
		  gbc = 0;
		}

	      c->SetBinContent(i+1,gbc);
	    }
	}
      else
	{
	  cout << "csm_interpolate_histogram: unknown interpolation style " << istyle << endl;
	  exit(0);
	}

      //cout << xa << " " << xb << " " << xc << endl;

      delete[] dista;
      delete[] distb;
      delete[] distc;
      delete[] distd;
    }
  else         // 2d case
    {
      Double_t *distxya = new Double_t[nbinsa*nbinsya];
      Double_t *distxyb = new Double_t[nbinsb*nbinsyb];
      Double_t *distxyc = new Double_t[nbinsc*nbinsyc];
      Double_t *distxyd = new Double_t[nbinsd*nbinsyd];

      hnorma = 0;
      for (j=0;j<nbinsya;j++)
	{
	  for (i=0;i<nbinsa;i++)
	    {
	      gbc = a->GetBinContent(i+1,j+1);
	      distxya[i+nbinsa*j] = gbc;
	      hnorma += gbc;
	    }
	}

      hnormb = 0;
      for (j=0;j<nbinsyb;j++)
	{
	  for (i=0;i<nbinsb;i++)
	    {
	      gbc = b->GetBinContent(i+1,j+1);
	      distxyb[i+nbinsb*j] = gbc;
	      hnormb += gbc;
	    }
	}

      hnormd = 0;
      for (j=0;j<nbinsyd;j++)
	{
	  for (i=0;i<nbinsd;i++)
	    {
	      gbc = d->GetBinContent(i+1,j+1);
	      distxyd[i+nbinsb*j] = gbc;
	      hnormd += gbc;
	    }
	}

      hnormc = hnorma + (xc-xa)*(hnormb-hnorma)/(xb-xa);
      // linearly interpolate the normalization between the central value and
      // the varied template.
      hnormc = hnormd*(hnormc/hnorma); // scale the normalization with the new template

      if (istyle == CSM_INTERP_HORIZONTAL || istyle == CSM_INTERP_HORIZONTAL_EXTRAP)
	{
          csm_pvmc2d(nbinsa,nbinsya,
                     distxya,
                     distxyb,
                     distxyd,
		     distxyc,
                     xa, xb, xc);

          hnormci = 0;
          for (j=0;j<nbinsyc;j++)
	    {
              for (i=0;i<nbinsc;i++)
	        {
	          hnormci += distxyc[i+nbinsc*j];
	        }
	    }
          for (j=0;j<nbinsyc;j++)
	    {
              for (i=0;i<nbinsc;i++)
	        {
	          c->SetBinContent(i+1,j+1,distxyc[i+nbinsc*j]*hnormc/hnormci);
	        }
	    }
	}
      else if (istyle == CSM_INTERP_VERTICAL || istyle == CSM_INTERP_VERTICAL_EXTRAP)
	{
          for (j=0;j<nbinsyc;j++)
	    {
              for (i=0;i<nbinsc;i++)
	        {
		  gbc = distxyd[i+nbinsc*j] + ((xc-xa)/(xb-xa))*(distxyb[i+nbinsc*j]-distxya[i+nbinsc*j]);
		  if (gbc < 0)
		    {
		      gbc = 0;
		    }
	          c->SetBinContent(i+1,j+1,gbc);
	        }
	    }
	}
      else
	{
	  cout << "csm_interpolate_histogram: unknown interpolation style " << istyle << endl;
	  exit(0);
	}

      delete[] distxya;
      delete[] distxyb;
      delete[] distxyc;
      delete[] distxyd;
    }

}

/*------------------------------------------------------------------------*/

/* compounded interpolation -- dist1 = central value shape, dist2 = syst. varied shape,
  dist3 = shape to distort (may be the result of previous distortions for compounded shape
  variations), distn = resultant shape.  par1 = value of parameter for dist1.  par2 = value of
  parameter (like # of sigma) for dist2.  parn = value of parameter for the output histogram
  Built on the idea of d_pvmorph, but generalized a bit.  Returns a null histogram if any of
  the three input histograms has zero or negative total sum.
*/

//#define DEBUGPVMC

void csm_pvmc(Int_t nb, Double_t *dist1, Double_t *dist2, Double_t *dist3, Double_t *distn,
	      Double_t par1, Double_t par2, Double_t parn)
{
  Int_t nb3 = nb*3 + 3;
  Int_t i,j,k,k1,k2;
  Double_t total;
  Double_t wta,wtb;
  Double_t yd[nb3];
  Double_t id[nb3];
  Double_t xd[nb3];
  Double_t xdis[nb3];
  Double_t ydis[nb3];
  Double_t ydi[nb + 1];
  Int_t idx[nb3];
  Int_t ifirst;
  Double_t x1l,x2l,x3l,y1l,y2l,y3l;
  Double_t x1,x2,x3,y1,y2,y3;
  Double_t xloc,yloc,x1i,x2i,x3i;

  // default output -- empty distribution

  for (i=0;i<nb;i++)
    {
      distn[i] = 0.0;
    }

  // default index list

  for (i=0;i<nb3;i++)
    { idx[i] = i;}

  // parameter weights

  if (par2 != par1) 
    {
      wta = 1. - (parn-par1)/(par2-par1);
      wtb = 1. + (parn-par2)/(par2-par1);
    }
  else
    {
      wta = 0.5;
      wtb = 0.5;
    }

  // suppress warning messages in case of extrapolations

  //  if ( (parn>par1 && parn>par2) || (parn<par1 && parn<par2) )
  //  {
  //    cout << "CSM_PVMC: Histogram Extrapolation: " << parn << 
  //            " is not between " << par1 << " and " << par2 << endl;
  //  }

  // Fill cumulative distribution arrays -- squeeze out repeated entries
  // due to empty bins

  // The first point in the cumulative distributions has zero integral and
  // starts at the left-hand edge of the first bin with any value in it.
  // The id array says which distribution it came from, and the
  // xd array gives the x value at which the cumulative distribution is evaluated
  // (at the right-hand edge of the bin)

  j = 0;
  total = 0;
  for (i=0;i<nb;i++)
    {
      if (dist1[i] < 0)
	{ 
	  cout << "Negative bin entry found in dist1 in csm_pvmc" << endl;
	  cout << i << " " << dist1[i] << endl;
	  exit(0);
	}
      total += dist1[i];
    }
  if (total <= 0) return;

  yd[j] = 0;
  id[j] = 1;
  j++;
  ifirst = 1;
  for (i=0;i<nb;i++)
    {
      if (dist1[i] > 0)
	{ 
	  if (ifirst==1)
	    {
	      ifirst = 0;
	      xd[j-1] = (Double_t) i;
	    }
	  yd[j] = yd[j-1] + dist1[i]/total;
	  id[j] = 1;
	  xd[j] = (Double_t) i+1;
	  j++;
	}
    }

  total = 0;
  for (i=0;i<nb;i++)
    {
      if (dist2[i] < 0)
	{ 
	  cout << "Negative bin entry found in dist2 in csm_pvmc" << endl;
	  cout << i << " " << dist2[i] << endl;
	  exit(0);
	}
      total += dist2[i];
    }
  if (total <= 0) return;
  yd[j] = 0;
  id[j] = 2;
  j++;
  ifirst = 1;
  for (i=0;i<nb;i++)
    {
      if (dist2[i]>0)
	{
	  if (ifirst==1)
	    {
	      ifirst = 0;
	      xd[j-1] = (Double_t) i;
	    }
          yd[j] = yd[j-1] + dist2[i]/total;
	  id[j] = 2;
	  xd[j] = (Double_t) i+1;
	  j++;
	}
    }

  total = 0;
  for (i=0;i<nb;i++)
    {
      if (dist3[i] < 0)
	{ 
	  cout << "Negative bin entry found in dist3 in csm_pvmc" << endl;
	  cout << i << " " << dist3[i] << endl;
	  exit(0);
	}
      total += dist3[i];
    }
  if (total <= 0) return;
  yd[j] = 0;
  id[j] = 3;
  j++;
  ifirst = 1;
  for (i=0;i<nb;i++)
    {
      if (dist3[i]>0)
	{
	  if (ifirst==1)
	    {
	      ifirst = 0;
	      xd[j-1] = (Double_t) i;
	    }
          yd[j] = yd[j-1] + dist3[i]/total;
	  id[j] = 3;
	  xd[j] = (Double_t) i+1;
	  j++;
	}
    }

  // Sort all of the edges of the cumulative distribution functions
  // j is the number of entries in the yd, xd and id arrays

  TMath::Sort(j,yd,idx,0);

#ifdef DEBUGPVMC
  for (i=0;i<j;i++)
    {
      cout << i << " " << xd[i] << " " << " " << yd[i] << " " << id[i] << endl;
    }
  cout << "Sort index" << endl;
  for (i=0;i<j;i++)
    {
      cout << i << " " << idx[i] << endl;
    }
#endif

  x1l = 0;
  x2l = 0;
  x3l = 0;

  y1l = 0;
  y2l = 0;
  y3l = 0;

  x1 = 0;
  x2 = 0;
  x3 = 0;

  y1 = 0;
  y2 = 0;
  y3 = 0;

  // the three lowest points in the sort should all have zero integral --
  // interpolate the x's of these

  for (i=0;i<3;i++)
    {
       if ( id[idx[i]] == 1 )
         {
           x1 = xd[idx[i]];
           y1 = yd[idx[i]]; // should be zero
         }
      else if ( id[idx[i]] == 2 )
         {
           x2 = xd[idx[i]];
           y2 = yd[idx[i]];  // should be zero
         }
      else if ( id[idx[i]] == 3 )
         {
           x3 = xd[idx[i]];
           y3 = yd[idx[i]];  // should be zero
         }
    }
  // don't have the other ends of the line segments yet -- find them as we go along.

#ifdef DEBUGPVMC
  cout << "first bins: " << x1 << " " << x2 << " " << x3 << endl;
#endif
  y1l = y1;
  y2l = y2;
  y3l = y3;
  x1l = x1;
  x2l = x2;
  x3l = x3;

  // first point on interpolated curve -- zero integral.

  k = 0;
  xdis[k] = wta*x1l + wtb*x2l - x1l + x3l;
  xdis[k] = min((Double_t) (nb+1),max(0.0,xdis[k]));
  ydis[k] = 0;

#ifdef DEBUGPVMC
  cout << "first point: " << xdis[0] << " " << ydis[0] << endl;
#endif

  for (i=3;i<j;i++)
    {
      xloc = xd[idx[i]];
      yloc = yd[idx[i]];

      if (id[idx[i]] == 1 )
	{
	  x1l = x1;
	  y1l = y1;
	  x1 = xloc;
	  y1 = yloc;

	  if (yloc>y2)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==2)
		    {
		      y2l = y2;
		      x2l = x2;
		      y2 = yd[idx[k1]];
		      x2 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	  if (yloc>y3)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==3)
		    {
		      y3l = y3;
		      x3l = x3;
		      y3 = yd[idx[k1]];
		      x3 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	}
      else if (id[idx[i]] == 2 )
	{
	  x2l = x2;
	  y2l = y2;
	  x2 = xloc;
	  y2 = yloc;

	  if (yloc>y1)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==1)
		    {
		      y1l = y1;
		      x1l = x1;
		      y1 = yd[idx[k1]];
		      x1 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	  if (yloc>y3)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==3)
		    {
		      y3l = y3;
		      x3l = x3;
		      y3 = yd[idx[k1]];
		      x3 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	}
      else if (id[idx[i]] == 3 )
	{
	  x3l = x3;
	  y3l = y3;
	  x3 = xloc;
	  y3 = yloc;

	  if (yloc>y2)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==2)
		    {
		      y2l = y2;
		      x2l = x2;
		      y2 = yd[idx[k1]];
		      x2 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	  if (yloc>y1)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==1)
		    {
		      y1l = y1;
		      x1l = x1;
		      y1 = yd[idx[k1]];
		      x1 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	}

      if (yloc>ydis[k] && ydis[k] < 0.999999999)
        {

#ifdef DEBUGPVMC
	  cout << "Interpolating: " << x1 << " " << x2 << " " << x3 << endl;
	  cout << "Interpolating: " << x1l << " " << x2l << " " << x3l << endl;
	  cout << "Interpolating: " << y1 << " " << y2 << " " << y3 << endl;
	  cout << "Interpolating: " << y1l << " " << y2l << " " << y3l << endl;
#endif

	  k++;
	  ydis[k] = yloc;
	  if (yloc == y1l)
	    {
	      x1i = x1l;
	    }
	  else
	    {
	      x1i = x1l + (yloc-y1l)*(x1-x1l)/(y1-y1l);
	    }
	  if (yloc == y2l)
	    {
	      x2i = x2l;
	    }
	  else
	    {
	      x2i = x2l + (yloc-y2l)*(x2-x2l)/(y2-y2l);
	    }
	  if (yloc == y3l)
	    {
	      x3i = x3l;
	    }
	  else
	    {
	      x3i = x3l + (yloc-y3l)*(x3-x3l)/(y3-y3l);
	    }
	  xdis[k] = x3i + wta*x1i + wtb*x2i - x1i;
          xdis[k] = min((Double_t) (nb+1),max(0.0,xdis[k]));
	  if (xdis[k]<xdis[k-1])
	    {
	      k--;
	      ydis[k] = yloc;
	    }

#ifdef DEBUGPVMC
	  cout << "point: " << k << endl;
	  cout << "x1i, x2i, x3i: " << x1i << " " << x2i << " " << x3i << endl;
	  cout << "Interpolated: " << xdis[k] << " " << yloc << endl;
#endif
	}
    }

#ifdef DEBUGPVMC
  for (i=0;i<=k;i++)
    {
      cout << "IC before bin: " << i << " " << xdis[i] << " " << ydis[i] << endl;
    }
#endif


  // k is the index of the last entry in the xdis, ydis interpolated array.
  // find the places where the piecewise linear interpolated cumulative distribution
  // crosses the bin edges  the index on ydi is the low bin edge.

  ydi[0] = 0.0;
  for (i=1;i<(1 + nb);i++)
    {
      ydi[i] = 1.0;
    }

  Int_t k2last;
  k1 = 0;
  for (i=0;i<=k;i++)
    {
      k2last = k1;
      for (k2=k1+1;k2<(nb+1);k2++)
	{
	  if ( (Double_t) k2 < xdis[i] )
	    {
	      if (i==0)
		{
		  ydi[k2] = 0;
		}
	      else
		{
	           ydi[k2] = ydis[i-1] + ( (Double_t) k2  - xdis[i-1] )*
                                         (ydis[i]-ydis[i-1])/(xdis[i]-xdis[i-1]);
#ifdef DEBUGPVMC
		   cout << "filling bins: " << i << " " << k1 << " " << k2 << " " << ydi[k2] << endl; 
#endif

		}
	       k2last = k2;
	    } 
	  if ( (Double_t) k2 > xdis[i] ) break;
	}
      k1 = k2last;
    }

#ifdef DEBUGPVMC
  for (i=0;i<(nb+1);i++)
    {
      cout << "interp. cumulative: " << i << " " << ydi[i] << endl;
    }
#endif

  // differentiate to get the output distn

  for (i=0;i<(nb);i++)
    {
      distn[i] = ydi[i+1] - ydi[i]; 
    }
}

/*------------------------------------------------------------------------*/


/*  
 Re-coded version of d_pvmorph_2d from Alex Read.  C version from Tom Junk
 Added feature of compounding shape variations as systematic errors.
 February 2007

......Do a linear interpolation between three two-dimensional
      probability distributions (scatterplots) as a function
      of the characteristic parameter of the distribution, for use
      in both interpolation and in application of systematic uncertainties.
      xydist1 is the "central value" histogram
      xydist2 is the "systematically varied" histogram
      xydist3 is the histogram to apply the variation to
      xydistn is the output histogram.  See csm_pvmc
        for compounded systematic variation applicaiton in 1D
        (used repeatedly in here).

      This is a generalization of csm_pvmc. The 2d distribution
      can be move around and be stretched or squeezed in two
      dimenions but finite rotations (changes in the correlation)
      are poorly approximated.

 nx        : Number of x-bins in the input and output distributions.
 ny        : Number of y-bins in the input and output distributions.
 xydist1,xydist2,xydist3,xydistn
           : Bin contents of scatterplots. The arrays should be
             packed with the index running fastest over the x
             dimension.
             Contents are in xydist[ix+nx*iy]
 par1,par2,parn     : Values of the linear parameters that characterise the
                      histograms (e.g. the Higgs mass).

 Output: xydistn.  Same binning as xydist1,xydist2,xydist3.
                   Its memory must be allocated outside of the routine

*/

//#define DEBUGPVMC2D

void csm_pvmc2d(Int_t nx, Int_t ny, Double_t *xydist1, 
                Double_t *xydist2, Double_t *xydist3, Double_t *xydistn,
                Double_t par1, Double_t par2, Double_t parn)
{
  Double_t ydist1[ny],ydist2[ny],ydist3[ny],ydistn[ny];
  Double_t xtemp1[nx],xtemp2[nx],xtemp3[nx],xtempn[nx];
  Double_t alpha1[ny*ny],alpha2[ny*ny],alpha3[ny*ny];
  Int_t i,j,k;

  // Project xydist1,2,3 onto the y axis and normalize

  csm_yproj(nx,ny,xydist1,ydist1);
  csm_yproj(nx,ny,xydist2,ydist2);
  csm_yproj(nx,ny,xydist3,ydist3);

  // Interpolate the y-projections 

  csm_pvmc(ny,ydist1,ydist2,ydist3,ydistn,par1,par2,parn);

#ifdef DEBUGPVMC2D
  for (i=0;i<ny;i++)
    {
      cout << "iy: " << i << " " << ydist1[i] << " " << ydist2[i] << " " << ydistn[i] << endl;
    }
#endif

  // Find out which y bins of histograms 1,2,3 contribute
  // to each bin of the interpolated distribution ydistn

  csm_ycont(ny,ydist1,ydist2,ydist3,ydistn,alpha1,alpha2,alpha3);

  // Extract the x-distributions in the y-slice determined above
  // and interpolate them

  for (i=0;i<ny;i++)  // loop over resulting bins
    {
      for (k=0;k<nx;k++) xtemp1[k] = 0;
      for (k=0;k<nx;k++) xtemp2[k] = 0;
      for (k=0;k<nx;k++) xtemp3[k] = 0;
      for (j=0;j<ny;j++) // loop over contributing bins
	{
	  for (k=0;k<nx;k++)
	    {
	      xtemp1[k] += alpha1[j+ny*i]*xydist1[k+nx*j];
	      xtemp2[k] += alpha2[j+ny*i]*xydist2[k+nx*j];
	      xtemp3[k] += alpha3[j+ny*i]*xydist3[k+nx*j];
	    }
	}
      // Interpolate the x distributions
      csm_pvmc(nx,xtemp1,xtemp2,xtemp3,xtempn,par1,par2,parn);

      // Insert the interpolated x distribution into the final output dist
      for (k=0;k<nx;k++) xydistn[k+nx*i] = xtempn[k]*ydistn[i];
    }
}

/*
 Re-coded version of d_ypvscat from Alex Read.  C version from Tom Junk
 Project a scatterplot onto the y-axis.The
 projection is normalized so that the sum of the bin contents is 1.0.

 nx,ny    : Number of bins in the scatterplot for the x and y coordindates.
            The projection is done onto <ny> bins.
 xydist   : The 2-dimensional array of the probabilities
 ydist    : The 1-dimensional array of the 2d probabilities projected onto
            the y-axis.

 Inputs : nx,ny,xydist
 Outputs: ydist (ny is unchanged from input to output)
*/

void csm_yproj(Int_t nx, Int_t ny, Double_t *xydist, Double_t *ydist)
{
  Int_t i,j;
  Double_t total;

  for (i=0;i<ny;i++) ydist[i] = 0;
  total = 0;

  for (i=0;i<ny;i++)
    {
      for (j=0;j<nx;j++)
	{
          ydist[i] += xydist[j+nx*i];
	}
      total += ydist[i]; 
    }

  if (total>0)
    {
      for (i=0;i<ny;i++) ydist[i] /= total;
    }
}

/*
Recoded d_getycont -- original by Alex Read, recoded by Tom Junk
February 2007

<ydist1> and <ydist2> and <ydist3>
are the projections on the y-axis of three
scatterplots which are going to be interpolated. <ydistn> is
the interpolated 1d distribution which represent the projection
of the interpolated scatterplot on the y-axis. This routine determines
which bins of <ydist1> and <ydist2> and <ydist3>
contribute and by what amount to
each bin of <ydistn>. This information is used in csm_pvmc2d to 
determine the input distributions in the x-direction of each
y-bin: these are then interpolated and accumulated in the interpolated
2d distribution.

Inputs : ny,ydist1,ydist2,ydist3,ydistn
Outputs: alpha1,alpha2,alpha3

alpha1[iyc+ny*iy] encodes the contribution of bin iyc in ydist1
to to bin iy in ydistn

*/

void csm_ycont(Int_t ny, Double_t *ydist1, Double_t *ydist2,
               Double_t *ydist3, Double_t *ydistn,
               Double_t *alpha1, Double_t *alpha2, Double_t *alpha3)
{
  Double_t y[ny+1];
  Double_t yn[ny+1];
  Int_t i;

  // Make arrays to describe the straight-line approximations
  // to the four cumulative distributions, y1,y2,y3,yn 
  // Make sure to start out with a point at 0

  yn[0] = 0;
  for (i=0;i<ny;i++) yn[i+1] = ydistn[i];
  csm_acnvec2(yn,ny+1);

  y[0] = 0;
  for (i=0;i<ny;i++) y[i+1] = ydist1[i];
  csm_acnvec2(y,ny+1);
  csm_ycontaux(ny,y,yn,alpha1);
#ifdef DEBUGPVMC2D
  for (i=0;i<ny+1;i++)
    {
      cout << "getting alpha1: " << i << " " << y[i] << " " << yn[i] << endl;
    }
  Int_t j; 
  for (i=0;i<ny;i++)
    {
      for (j=0;j<ny;j++)
	{
	   cout << i << " " << j << " " << alpha1[i+ny*j] << endl;
	}
    }
  
#endif

  y[0] = 0;
  for (i=0;i<ny;i++) y[i+1] = ydist2[i];
  csm_acnvec2(y,ny+1);
  csm_ycontaux(ny,y,yn,alpha2);

  y[0] = 0;
  for (i=0;i<ny;i++) y[i+1] = ydist3[i];
  csm_acnvec2(y,ny+1);
  csm_ycontaux(ny,y,yn,alpha3);
}

void csm_ycontaux(Int_t ny, Double_t *y, Double_t *yn,
                  Double_t *alpha)
{
  Int_t ny2;
  Int_t i,j;

  ny2 = ny*ny;

  // clear out the alpha array

  for (i=0;i<ny2;i++) alpha[i] = 0;

  // loop over bins and see what fraction each contributes

  for (i=0;i<ny;i++) // interpolated histogram
    {
      for (j=0;j<ny;j++) // contributing histogram bin
	{
          if (y[j+1]-y[j]>0)
	    {
	      // first case -- contributing bin entirely contained
	      // within the interpolated output bin.
	      if (y[j]>=yn[i] && y[j]<yn[i+1] && 
		  y[j+1]>=yn[i] && y[j+1]<yn[i+1])
		{
		  alpha[j+ny*i] = 1;
		}
	      // second case -- interpolated output bin is entirely
	      // contained within the contributing bin
	      else if (y[j]<yn[i] && y[j+1] >= yn[i+1])
		{
		  alpha[j+ny*i] = (yn[i+1]-yn[i])/(y[j+1]-y[j]);
		}
	      // third case -- contributing bin straddles the
	      // left edge of the interpolated output bin but ends inside
	      // the bin
	      else if (y[j]<yn[i] && y[j+1]>=yn[i] && y[j+1]<yn[i+1])
		{
		  alpha[j+ny*i] = (y[j+1]-yn[i])/(y[j+1]-y[j]);
		}
	      // fourth case -- contributing bin straddles the
	      // right edge of the interpolated output bin but starts inside
	      // the output bin
	      else if (y[j]>=yn[i] && y[j]<yn[i+1] && y[j+1]>=yn[i+1])
		{
		  alpha[j+ny*i] = (yn[i+1]-y[j])/(y[j+1]-y[j]);
		}
	      // non-overlapping case -- do nothing.
	      // save some time if we're beyond the edge
	      if (y[j]>yn[i+1]) break; 
	    }
	}
    } 
}

// Integrate and normalize a vector -- pretty much the same as
// csm_acnvec

void csm_acnvec2(Double_t *vec, Int_t n)
{
  Int_t i;
  Double_t tot;
  for (i=1;i<n;i++)
    {
      vec[i] += vec[i-1];
    }
  tot = vec[n-1];
  if (tot > 0)
    {
      for (i=0;i<n;i++) vec[i] /= tot;
    }
}


