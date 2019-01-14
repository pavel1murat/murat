///////////////////////////////////////////////////////////////////////////////
// cb2: Crystall Ball function with polynomial tails on both sides
///////////////////////////////////////////////////////////////////////////////

namespace {
  TH1F* _Hist;
  TF1*  _Func;
};

//-----------------------------------------------------------------------------
double cb2_gauss(double* X, double* P) {
  double x0, f(0), sig, xdx;

  x0         = P[1];
  sig        = P[2];
  xdx        = (X[0]-x0)/sig;

  f = P[0]*TMath::Exp(-xdx*xdx/2.);
  
  return f;
}

//-----------------------------------------------------------------------------
double cb2_f(double* X, double* P) {
  double x0, x1, f, sig, alpha1, abs_alpha1, alpha2, abs_alpha2, n1, n2, a, b, dx0, dx1, xdx;

  x0         = P[1];
  sig        = P[2];
  alpha1     = P[3];
  n1         = P[4];
  alpha2     = P[5];
  n2         = P[6];

  dx0        = X[0]-x0;
  xdx        = dx0/sig;

  abs_alpha1 = fabs(alpha1);
  abs_alpha2 = fabs(alpha2);

  if (xdx < -abs_alpha1) {
    a = TMath::Power((n1/abs_alpha1),n1)*TMath::Exp(-(alpha1*alpha1)/2);
    b = n1/abs_alpha1-abs_alpha1;
    f = P[0]*a*TMath::Power((b-xdx),-n1);
  }
  else if (xdx < abs_alpha2) { 
					// gauss
    f = P[0]*TMath::Exp(-xdx*xdx/2.);
  }
  else {
    a = TMath::Power((n2/abs_alpha2),n2)*TMath::Exp(-(alpha2*alpha2)/2);
    b = -n2/abs_alpha2+abs_alpha2;
    f = P[0]*a*TMath::Power((xdx-b),-n2);
  }
  
  return f;
}

//-----------------------------------------------------------------------------
void cb2_init_parameters(TF1* F, 
			 double P0, double P1, double P2, double P3, 
			 double P4, double P5, double P6) {
  F->SetParameter(0,P0);
  F->SetParameter(1,P1);
  F->SetParameter(2,P2);
  F->SetParameter(3,P3);
  F->SetParameter(4,P4);
  F->SetParameter(5,P5);
  F->SetParameter(6,P6);
}

//-----------------------------------------------------------------------------
void cb2_fit(TH1F* Hist, TF1* F, double XMin, double XMax) {
  Hist->GetXaxis()->SetRangeUser(XMin,XMax);
  Hist->Fit(F,"L","",XMin,XMax);
}

//-----------------------------------------------------------------------------
void cb2_create_fit_function(TF1*& F, double X0, double XMin, double XMax) {

  _Func = new TF1("cb2_f",cb2_f,XMin,XMax,7);

  F->SetParName(0,"anorm");
  F->SetParName(1,"mean");
  F->SetParName(2,"sigma");
  F->SetParName(3,"alpha1");
  F->SetParName(4,"N1");
  F->SetParName(5,"alpha2");
  F->SetParName(6,"N2");

  F->SetParLimits(2,0.,1.);
  //  F->SetParLimits(5,0.,5.);
  //  F->SetParLimits(6,0.,5.);

  F->SetLineWidth(1);
}

//-----------------------------------------------------------------------------
void cb2_fit_crystal_ball(TH1* Hist, double X0, double XMin, double XMax, double Sigma = -1.) {

  TH1F* h      = (TH1F*) Hist;
  double anorm = h->GetEntries()/10;

  printf("anorm = %12.5e\n",anorm);

  cb2_create_fit_function(_Func,X0,XMin,XMax);
  cb2_init_parameters    (_Func,anorm,X0,0.180,1,.4,1.,4);

//   h->Draw();
//   _Func->Draw("same");

  if (Sigma > 0) { // effectively, fix Sigma
    _Func->SetParameter(2,Sigma);
    _Func->SetParLimits(2,Sigma*(1+1.e-5),Sigma*(1+1.e-5));
  }

  cb2_fit                (h,_Func,XMin,XMax);

  TF1* f2 = new TF1("cb2_gauss",cb2_gauss,XMin,XMax,3);
  f2->SetParameter(0,_Func->GetParameter(0));
  f2->SetParameter(1,_Func->GetParameter(1));
  f2->SetParameter(2,_Func->GetParameter(2));

  TH1F* h_cb2_gaus = (TH1F*) h->Clone("h_cb2_gaus");
  h_cb2_gaus->Reset();

  for (int i=1; i<h->GetNbinsX(); i++) {
    float x = h->GetBinCenter(i);
    float y = f2->Eval(x);
    h_cb2_gaus->SetBinContent(i,y);
    h_cb2_gaus->SetBinError  (i,0);
  }
  
  h_cb2_gaus->SetLineColor(kBlue+3);
  h_cb2_gaus->SetFillColor(kBlue+3);
  h_cb2_gaus->SetFillStyle(3003);
  
  h_cb2_gaus->Draw("same");
//-----------------------------------------------------------------------------
// to calculate the background fraction, use region from alpha2 to XMax
//-----------------------------------------------------------------------------
  double total = Hist->Integral(1,Hist->GetNbinsX());

  double x0     = _Func->GetParameter("mean"  );
  double sigma  = _Func->GetParameter("sigma" );
  double alpha2 = _Func->GetParameter("alpha2");
  double x1     = x0+fabs(sigma*alpha2);

  printf("x0 = %10.4f, sigma = %10.4f alpha2 = %10.4f\n",x0,sigma, alpha2);
 
  double htail = (_Func->Integral(x1,XMax)-f2->Integral(x1,XMax))/Hist->GetBinWidth(1);

  printf("total = %12.5f, HTail = %12.5f, htail/total : %12.5e\n",total,htail,htail/total);
}

//-----------------------------------------------------------------------------
void cb2_fit_crystal_ball(const char* File, const char* Module, const char* Hist, 
			  double X0, double XMin, double XMax, double Sigma = -1.) {

  _Hist = (TH1F*) gh1(File,Module,Hist)->Clone("h_cb2_fit_crystal_ball");

  cb2_fit_crystal_ball(_Hist,X0,XMin,XMax,Sigma);
}

