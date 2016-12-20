//
#include "TGraph.h"
#include "TLegend.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TArrow.h"
#include "TPave.h"
#include "TText.h"

TGraph  *gr_pion_arrivals(0), *gr_pion_arrivals_2(0);
TGraph  *gr_muon_arrivals(0), *gr_muon_arrivals_2(0);
TGraph  *gr_muon_decays  (0), *gr_muon_decays_2  (0);

//-----------------------------------------------------------------------------
int make_pion_arrival_graph(TGraph*& Graph, double TimeOffset) {
  int npt(31);

  double data [] = {
    67.    ,0.         ,
    78.4419,0.00268220 ,
    87.2960,0.00591928 ,
    98.3637,0.00996563 ,
    106.115,0.0136075  ,
    112.773,0.0192731  ,
    124.955,0.0250732  ,
    136.039,0.0322224  ,
    144.905,0.0377529  ,
    153.765,0.0420692  ,
    164.837,0.0469250  ,
    173.683,0.0486781  ,
    183.638,0.0512406  ,
    194.682,0.0507001  ,
    202.408,0.0496202  ,
    216.748,0.0458417  ,
    225.569,0.0427381  ,
    233.269,0.0368015  ,
    245.388,0.0308646  ,
    254.193,0.0247931  ,
    264.102,0.0187215  ,
    274.017,0.0137291  ,
    282.833,0.00981608 ,
    294.968,0.00671224 ,
    306.002,0.00455284 ,
    313.727,0.00306825 ,
    325.871,0.00185312 ,
    334.703,0.00090806 ,
    344.644,0.00077236 ,
    353.479,0.00050184 ,
    363.   , 0
  };

  double x[1000], y[1000];

  for (int i=0; i<npt; i++) {
    x[i] = data[2*i  ]+TimeOffset;
    y[i] = data[2*i+1];
  }

  if (Graph) delete Graph;
  
  Graph = new TGraph(npt,x,y);
  Graph->SetLineColor(kBlue+3);
  Graph->SetLineWidth(2);
  Graph->SetFillColor(0);

  //  gr_pion_arrivals->Draw("al");

  return 0;
}

//-----------------------------------------------------------------------------
int make_muon_arrival_graph(TGraph*& Graph, double TimeOffset) {
  int npt(41);

  double data [] = {
    72.    ,0.           ,
    89.4794,0.00106243   ,
    104.955,0.00308482   ,
    114.908,0.00524255   ,
    124.861,0.00753519   ,
    134.818,0.0103675    ,
    143.668,0.0127951    ,
    154.733,0.0164367    ,
    164.696,0.0204832    ,
    172.445,0.0235854    ,
    183.509,0.0268223    ,
    194.564,0.0285752    ,
    204.519,0.0311377    ,
    214.468,0.0324860    ,
    224.417,0.0338343    ,
    233.258,0.0346430    ,
    246.510,0.0339674    ,
    255.338,0.0323478    ,
    263.059,0.0301887    ,
    275.203,0.0289736    ,
    287.337,0.0258697    ,
    293.960,0.0249248    ,
    304.991,0.0220909    ,
    314.922,0.0202014    ,
    325.955,0.0176373    ,
    334.782,0.0157479    ,
    344.710,0.0131839    ,
    354.644,0.0118340    ,
    364.578,0.0103492    ,
    373.407,0.00899945   ,
    384.445,0.00751459   ,
    394.381,0.00643454   ,
    403.214,0.00562439   ,
    414.256,0.00481407   ,
    426.401,0.00373384   ,
    436.337,0.00278870   ,
    458.426,0.00211241   ,
    500.394,0.000760005  ,
    541.263,0.000486947  ,
    583.237,0.000213801  ,
    645.095,0.
  };

  double x[1000], y[1000];

  for (int i=0; i<npt; i++) {
    x[i] = data[2*i  ]+TimeOffset;
    y[i] = data[2*i+1];
  }

  if (Graph) delete Graph;
  
  Graph = new TGraph(npt,x,y);
  Graph->SetLineColor(kBlue-3);
  Graph->SetFillColor(kBlue-3);
  Graph->SetFillStyle(3004);
  Graph->SetLineWidth(2);

  //  gr_muon_arrivals->Draw("al");

  return 0;
}

//-----------------------------------------------------------------------------
int make_muon_decay_graph(TGraph*& Graph, double TimeOffset) {
  int npt(34);

  double data [] = {
    70.6974,0.0         ,
    87.2666,0.000388067 ,
    108.257,0.000926032 ,
    134.771,0.00146356  ,
    164.598,0.00200082  ,
    197.742,0.00294254  ,
    222.046,0.00348025  ,
    262.924,0.00482608  ,
    289.438,0.00549851  ,
    317.059,0.00657558  ,
    351.302,0.00657286  ,
    389.966,0.00710943  ,
    439.672,0.00670076  ,
    489.378,0.00642700  ,
    541.293,0.00615306  ,
    602.043,0.00547371  ,
    661.692,0.00533406  ,
    722.443,0.00478961  ,
    796.449,0.00424411  ,
    891.444,0.00396676  ,
    963.241,0.00342143  ,
    1059.34,0.00300908  ,
    1146.60,0.00259743  ,
    1225.03,0.00259121  ,
    1290.20,0.00218131  ,
    1365.31,0.00217535  ,
    1434.90,0.00203492  ,
    1502.28,0.00189467  ,
    1563.04,0.00162003  ,
    1653.61,0.00134303  ,
    1711.05,0.00133847  ,
    1759.66,0.00119970  ,
    1818.20,0.000925243 ,
    2000.  ,0.0 
  };

  double x[1000], y[1000];

  for (int i=0; i<npt; i++) {
    x[i] = data[2*i  ]+TimeOffset;
    y[i] = data[2*i+1];
  }

  if (Graph) delete Graph;
  
  Graph = new TGraph(npt,x,y);
  Graph->SetLineColor(kRed+3);
  Graph->SetFillColor(kRed+3);
  Graph->SetFillStyle(3005);
  Graph->SetLineWidth(2);

  //  gr_muon_arrivals->Draw("al");

  return 0;
}

//-----------------------------------------------------------------------------
int plot_arrival_times(int Print = 0) {

  TCanvas* c = new TCanvas("c_arrivals","arrivals",1400,500);

  TH1F* h = new TH1F("h_arrivals","",2500,-200,2300);
  h->SetMaximum(0.07);
  h->SetStats(0);
  TAxis* ax = h->GetXaxis();
  ax->SetTitle("Time, ns");
  ax->SetLabelSize(0.05);
  ax->SetTitleSize(0.05);
  //  ax->SetTitleOffset(0.05);

  TAxis* ay = h->GetYaxis();
  ay->SetLabelSize(0.05);
  ay->SetTitleSize(0.05);
  ay->SetTitleOffset(0.6);
  ay->SetTitle("dN/dt, arbitrary normalization");
  
  h->Draw();

  make_pion_arrival_graph(gr_pion_arrivals,0);
  gr_pion_arrivals->Draw("same");

  make_pion_arrival_graph(gr_pion_arrivals_2,1695.);
  gr_pion_arrivals_2->Draw("same");

  make_muon_arrival_graph(gr_muon_arrivals,0);
  gr_muon_arrivals->Draw("l,same");
  gr_muon_arrivals->Draw("f,same");

  make_muon_arrival_graph(gr_muon_arrivals_2,1695);
  gr_muon_arrivals_2->Draw("l,same");
  gr_muon_arrivals_2->Draw("f,same");

  make_muon_decay_graph(gr_muon_decays,0);
  gr_muon_decays->Draw("l,same");
  gr_muon_decays->Draw("f,same");

  make_muon_decay_graph(gr_muon_decays_2,1695);
  gr_muon_decays_2->Draw("l,same");
  gr_muon_decays_2->Draw("f,same");

  TLegend* leg = new TLegend(0.4,0.7,0.7,0.85);
  leg->SetBorderSize(0);

  leg->AddEntry(gr_pion_arrivals,"pion arrival to stopping target","f");
  leg->AddEntry(gr_muon_arrivals,"muon arrival to stopping target","f");
  leg->AddEntry(gr_muon_decays  ,"muon decays"  ,"f");

  leg->Draw();

  // a1 and a2: beam arrival times, a3: start of the measurement
  
  TArrow* a1 = new TArrow(0,0.02,0,0.0,0.015);
  a1->SetLineWidth(2);
  a1->Draw();

  TArrow* a2 = new TArrow(1695.,0.02,1695,0.0,0.015);
  a2->SetLineWidth(2);
  a2->Draw();
  
  TPave* pave = new TPave(700.,0,1695,0.01,1);
  pave->SetFillColor(2);
  pave->SetFillStyle(3004);
  pave->Draw();

  TText* txt1 = new TText(900.,0.011,"measurement time window");
  txt1->SetTextFont(42);
  txt1->Draw();

  TText* txt2 = new TText(-50.,0.022,"beam");
  txt2->SetTextFont(42);
  txt2->Draw();
			  
  TText* txt3 = new TText(-50+1695.,0.022,"beam");
  //  txt3->SetNDC(kFALSE);
  txt3->SetTextFont(42);
  txt3->Draw();

  if (Print != 0) {
    c->Print("arrival_time_stopping_target.eps");
  }
			  
  return 0;
}
