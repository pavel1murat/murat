///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TROOT.h"

#include "math.h"
#include "string.h"
#include "limits/TLHChannel.hh"

ClassImp(TLHChannel)

//-----------------------------------------------------------------------------
TLHChannel::TLHChannel(): TNamed() {
  fNExpected       = 0;
  fProbHist        = 0;
  fListOfProcesses = new TObjArray();
}

//-----------------------------------------------------------------------------
TLHChannel::TLHChannel(const char* Name, const char* Title): TNamed(Name,Title) {
  fNExpected       = 0;
  fProbHist        = 0;
  fListOfProcesses = new TObjArray();
}



//-----------------------------------------------------------------------------
TLHChannel::~TLHChannel() {
  delete fListOfProcesses;
}

//-----------------------------------------------------------------------------
int TLHChannel::Init() {

  TLHProcess   *p;
  char         name[200];
  double       nexp, qint;
  TH1*         h;
  TObject*     o;

  int nb = fListOfProcesses->GetEntriesFast();

  fNExpected = 0;

  for (int i=0; i<nb; i++) {
    p           = (TLHProcess*) fListOfProcesses->At(i);
    h           = p->GetProbHist ();
    nexp        = p->GetNExpected();
    fNExpected += nexp;
    if (fProbHist == 0) {
      sprintf(name,"%s_channel_bgr_prob_hist",GetName());

      while (o = gROOT->FindObject(name)) delete o;
      fProbHist = (TH1*) h->Clone(name);
      fProbHist->Scale(nexp);
    }
    else {
      fProbHist->Add(h,nexp);
    }
  }
//-----------------------------------------------------------------------------
// finally, normalize both histograms to the unit integral
// we need combined probability distribution to plot likelihood...
//-----------------------------------------------------------------------------
  qint = fProbHist->Integral();
  fProbHist->Scale(1./qint);

  return 0;
}

