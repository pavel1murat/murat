///////////////////////////////////////////////////////////////////////////////
// plot chi2/Ndf curve for Ndf degrees of freedom
///////////////////////////////////////////////////////////////////////////////
#ifndef __plot_chi2__
#define __plot_chi2__

void plot_chi2(double XMin, double XMax, int Ndf, const char* Opt, double Norm) {

  double xmin = XMin;
  double xmax = XMax;
  int nbins   = 1000;
  
  TH1F* h = new TH1F("h1","h1",nbins,xmin,xmax);

  double bin = (xmax-xmin)/nbins;

  for (int i=1; i<=nbins; i++) {
    double x = (xmin+bin*(i-0.5))*Ndf;
    double y = ROOT::Math::chisquared_pdf(x,Ndf);

    printf("x,y: %12.5e %12.5e\n",x,y);

    h->SetBinContent(i,y);
    h->SetBinError  (i,0);
  }

  h->SetStats(0);
  h->DrawNormalized(Opt,Norm);
}

#endif
