///////////////////////////////////////////////////////////////////////////////



#include "Stntuple/val/stntuple_val_functions.hh"
#include "murat/ana/prob_dist.hh"


int prob_dist::fgIndex(0);

//-----------------------------------------------------------------------------
prob_dist::prob_dist() {
  h     = NULL;
  hprob = NULL;
}


//-----------------------------------------------------------------------------
prob_dist::prob_dist(TH1F* Hist) {

  h     = (TH1*) Hist->Clone(Form("hprob_dist_h_%i"  ,fgIndex));
  hprob = (TH1F*) h->Clone(Form("hprob_dist_hprob_%i",fgIndex));
  fgIndex += 1;

  hprob->Reset();

  int nb = h->GetNbinsX();
  
  double anorm = h->Integral(1,nb);

  for (int i=1; i<=nb; i++) {
    double prob = h->Integral(1,i)/anorm;
    hprob->SetBinContent(i,prob);
  }
}


//-----------------------------------------------------------------------------
prob_dist::prob_dist(const char* Fn, const char* Folder, const char* Hist) {

  h     = (TH1*) gh1(Fn,Folder,Hist)->Clone(Form("hprob_dist_h_%i"  ,fgIndex));
  hprob = (TH1F*) h->Clone(Form("hprob_dist_h_%i"  ,fgIndex));
  fgIndex += 1;

  hprob->Reset();

  int nb = h->GetNbinsX();
  
  double anorm = h->Integral(1,nb);

  for (int i=1; i<=nb; i++) {
    double prob = h->Integral(1,i)/anorm;
    hprob->SetBinContent(i,prob);
  }
}


//-----------------------------------------------------------------------------
prob_dist::~prob_dist() {
}


//-----------------------------------------------------------------------------
double prob_dist::prob(double X) {
  double f(0);

  int nb = h->GetNbinsX();
					// assume all bins are the same
  double binw = h->GetBinWidth(1);

  if      (X <  hprob->GetBinCenter( 1)-binw/2) return 0.;
  else if (X >= hprob->GetBinCenter(nb)+binw/2) return 1.;
  
  for (int i=1; i<=nb; i++) {
    double x = h->GetBinCenter(i);
    if (x+binw/2 > X) {
      f = hprob->GetBinContent(i);
      break;
    }
  }

  return f;
}
