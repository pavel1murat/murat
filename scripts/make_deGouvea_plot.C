///////////////////////////////////////////////////////////////////////////////
// ROOT macro corresponding to Bob's doc-db 4314 11/2016 matematica notebook
// 2016-11-24 P.Murat
// call :
// -------
//              make_deGouvea_plot(Print)
//
// Print = 0: no printout (default)
//       = 1: print .pdf file
///////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include "TF1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"

#include "TLatex.h"
//-----------------------------------------------------------------------------
void draw_label_abs(const char* Text, double X1, double Y1, double FontSize, double Font = 52) {

  TLatex* label = new TLatex(X1,Y1,Text);
  label->SetNDC(kFALSE);
  label->SetTextSize(FontSize);
  label->SetTextFont(Font);
  label->Draw();
}

double mmu   = 0.10566;        // muon mass
double alpha = 1./137.03599;   // em alpha
double Gf    = 1.16637e-5;     // Fermi constant
double pi    = M_PI;

double MEG_LIMIT = 4.2e-13;    // MEG'2016 limit on BR(mu --> e gamma)

double Dr;
double Dl;

//------------------------------------------------------------------------------
double MuEGEnergyScale(double b) {

  double x = 3*alpha/(32*pi)*(Dr*Dr+Dl*Dl)/b;
  double f = (1./1000)*pow(x,1/4.);

  return f;
}

//------------------------------------------------------------------------------
double Mu2eEnergyScale(double b, double N) {
  
  double mmu4 = pow(mmu,4);
  double x    = 3*alpha/(32*pi)*(Dr*Dr+Dl*Dl)/b;
  double f    = (1./1000)*pow(Gf*Gf*mmu4/(96*pi*pi*pi*alpha)*3*1.e12*N*x,1/4.);

  //  printf("mmu, mmu4, x, f : %14.7f %14.7f  %14.7f %14.7f \n",mmu, mmu4,x,f);

  return f;
}

//-----------------------------------------------------------------------------
double Mu2e_II(double Kappa) {

  double powerBetterForMu2e2 = 10;
  double sf = pow(powerBetterForMu2e2,1./4);
  double f  = Mu2eEnergyScale(6e-18,1.1)*1./(1+Kappa) + 8500.*sf*Kappa/(1+Kappa);
  return f;
}


//-----------------------------------------------------------------------------
double Mu2e(double Kappa) {
  double f = Mu2eEnergyScale(6e-17,1.1)*1./(1+Kappa) + 8500.*Kappa/(1+Kappa);
  return f;
}

//-----------------------------------------------------------------------------
// 1.25 vs 1.1 - is itAu vs Al ?
//-----------------------------------------------------------------------------
double Sindrum(double Kappa) {
  double f = Mu2eEnergyScale(7e-13,1.25)*1./(1+Kappa) + 1100.*Kappa/(1+Kappa);
  return f;
}


//-----------------------------------------------------------------------------
double Meg(double Kappa, double Limit) {
  double f = MuEGEnergyScale(Limit)*(1./(1+Kappa));
  return f;
}


//-----------------------------------------------------------------------------
void make_deGouvea_plot(int Print = 0) {

  Dr = 16*sqrt(2)*pi*pi/Gf;
  Dl = Dr;

  printf (" Dr, Df = %14.7g %14.7g\n",Dr,Dl);
//-----------------------------------------------------------------------------
// set X-coordinates
//-----------------------------------------------------------------------------
  double x[1000], y[1000];

  int imin = -3;
  int imax =  3;
  
  x[0] = 1.*pow(10,imin);

  int loc = 1;
  for (int i=imin; i<imax; i++) {
    for (int j=1; j<10; j++) {
      x[loc] = j*pow(10,i);
      loc++;
    }
  }

  x[loc] = x[loc-1];
  loc++;

  TCanvas* c = new TCanvas("c","c",1200,800);
  //  c->SetPad(0,1000,99,50000);

  c->Draw();
  y[0]     = 1.; // Mu2e(x[0])*1.e-1;
//-----------------------------------------------------------------------------
// Mu2e-II
//-----------------------------------------------------------------------------
  y[loc-1] = y[0];
  for (int i=1; i<loc-1; i++) y[i] = Mu2e_II(x[i]);

  TGraph* g_mu2e_II = new TGraph(loc,x,y);

  g_mu2e_II->SetLineColor(kRed+2);
  g_mu2e_II->SetTitle("");
  g_mu2e_II->GetYaxis()->SetTitle("#Lambda, TeV");
  g_mu2e_II->GetYaxis()->SetTitleOffset(0.5);
  g_mu2e_II->GetXaxis()->SetTitle("k");
  g_mu2e_II->GetXaxis()->SetTitleOffset(0.5);
  //  g_mu2e_II->GetXaxis()->SetRangeUser(1.2e-3,999);
  g_mu2e_II->SetMinimum(100.);
  g_mu2e_II->SetMaximum(50000.);
  g_mu2e_II->SetLineStyle(2);
  g_mu2e_II->Draw("al");
//-----------------------------------------------------------------------------
// Mu2e
//-----------------------------------------------------------------------------
//  printf("Mu2E(0) = %14.7g %14.7g\n",Mu2e(0),Mu2e(1.e4)); 

  for (int i=1; i<loc-1; i++) y[i] = Mu2e(x[i]);
  
  TGraph* g_mu2e = new TGraph(loc,x,y);
  g_mu2e->SetTitle("");
  g_mu2e->GetYaxis()->SetTitle("#Lambda, TeV");

  g_mu2e->GetXaxis()->SetTitle("k");
  g_mu2e->SetLineColor(kBlue+3);

  g_mu2e->Draw("l,same");
//-----------------------------------------------------------------------------
// SINDRUM-II
//-----------------------------------------------------------------------------
  y[loc-1] = y[0];
  for (int i=1; i<loc-1; i++) y[i] = Sindrum(x[i]);

  TGraph* g_sindrum = new TGraph(loc,x,y);

  g_sindrum->SetLineColor(kRed+2);
  g_sindrum->Draw("l,same");
  g_sindrum->SetFillStyle(3005);
  g_sindrum->SetFillColor(kRed-9);
  g_sindrum->Draw("f,same");
//-----------------------------------------------------------------------------
// MEG
//-----------------------------------------------------------------------------
  y[loc-1] = y[0];
  for (int i=1; i<loc-1; i++) y[i] = Meg(x[i], MEG_LIMIT);

  TGraph* g_meg = new TGraph(loc,x,y);

  g_meg->Draw("l,same");
  g_meg->SetFillStyle(3004);
  g_meg->SetFillColor(kCyan-9);
  g_meg->Draw("f,same");
//-----------------------------------------------------------------------------
// MEG upgrade - x10 sensitivity
//-----------------------------------------------------------------------------
  y[loc-1] = y[0];
  for (int i=1; i<loc-1; i++) y[i] = Meg(x[i], MEG_LIMIT/10.);

  TGraph* g_meg_upgrade = new TGraph(loc,x,y);

  g_meg_upgrade->SetLineStyle(2);
  g_meg_upgrade->Draw("l,same");
//-----------------------------------------------------------------------------
// labels
//-----------------------------------------------------------------------------
  draw_label_abs("R_{#mue}^{Au}(#muN #rightarrow eN) < 7 #times 10^{-13}", 7,800,0.025,52);
  draw_label_abs("Excluded #muN #rightarrow eN"                          ,10,400,0.025,52);
  draw_label_abs("SINDRUM-II"                                            ,15,300,0.025,52);

  draw_label_abs("BR(#mu #rightarrow e#gamma) < 4.2 #times 10^{-13}",1.2e-3,550,0.025,52);
  draw_label_abs("Excluded #mu #rightarrow e#gamma"                 ,1.5e-3,400,0.025,52);
  draw_label_abs("MEG"                                              ,4.5e-3,300,0.025,52);

  draw_label_abs("R_{#mue}^{Al}(#muN #rightarrow eN) < 6 #times 10^{-17}", 7, 5500,0.025,52);
  draw_label_abs("R_{#mue}^{Al}(#muN #rightarrow eN) < 6 #times 10^{-18}", 7,11000,0.025,52);

  draw_label_abs("BR(#mu #rightarrow e#gamma) < 4.2 #times 10^{-14}",1.2e-3,1050,0.025,52);

  draw_label_abs("all limits at 90\% CL",1.2e-3,4.e4,0.025,52);

  c->SetLogx(1);
  c->SetLogy(1);
  c->Modified();
  c->Update();

  if (Print != 0) c->Print("deGouvea_plot.pdf");
}
