///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TObjString.h"
#include "TCanvas.h"

#include "math.h"
#include "string.h"
#include "limits/TLHProcess.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

ClassImp(TLHProcess)

//-----------------------------------------------------------------------------
TLHProcess::TLHProcess(): TNamed() {
}

//-----------------------------------------------------------------------------
TLHProcess::TLHProcess(const char* Name, const char* Title, double NExpected): 
  TNamed(Name,Title) {
  fNExpected = NExpected;
}

//-----------------------------------------------------------------------------
TLHProcess::~TLHProcess() {
}


//-----------------------------------------------------------------------------
// a histogram could be either 1- or 2-dimensional
//-----------------------------------------------------------------------------
void TLHProcess::SetProbHist(const char* Filename, 
			     const char* Module  , 
			     const char* Hist    ,
			     int         NDim    ) {
  double qent;
  TFile* f;

  if (Module) {
    if (NDim == 1) {
      Error("SetProbHist","1D Not implemented yet");
    }
    else {
      fProbHist = gh2(Filename,Module,Hist);
    }
  }
  else {
					// Module == 0
    f = TFile::Open(Filename);
    if (NDim == 1) {
      Error("SetProbHist","1D Not implemented yet");
    }
    else {
      fProbHist = (TH2F*) f->Get(Hist);
    }
  }
//-----------------------------------------------------------------------------
// normalize fProdHist to intergal of 1
//-----------------------------------------------------------------------------
  qent  = fProbHist->Integral();
  //    fProbHist->Rebin2D(2,2);
  fProbHist->Scale(1./qent);
}


