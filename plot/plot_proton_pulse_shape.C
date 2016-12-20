///////////////////////////////////////////////////////////////////////////////
// draw N pulses
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include "TObject.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TInterpreter.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"

#include "murat/scripts/plot_utilities.hh"
#include "murat/scripts/create_ntuple.C"


TCanvas* c;
TNtuple* nt_pulse (NULL);
TGraph*  g_pulse[100]; // assume NPulses < 100

//-----------------------------------------------------------------------------
// 
int make_proton_pulse_graph(TGraph*& Graph, int PulseNumber) {
  double x[10000], y[10000];
  
  if (nt_pulse == NULL) create_ntuple("ConditionsService/data/potTimingDistribution_avg.txt",nt_pulse);

  float var[100];

  TBranch* b1 = nt_pulse->GetBranch("var1");
  TBranch* b2 = nt_pulse->GetBranch("var2");

  nt_pulse->SetBranchAddress("var0",&var[0]);
  nt_pulse->SetBranchAddress("var1",&var[1]);

  //  b2->SetBranchAddress(&var[1]);

  // nt_pulse->Draw("var1:var0");
  // nt_pulse->SetMarkerStyle(1);

  int nent = nt_pulse->GetEntries();

  for (int ientry=0; ientry<nent; ientry++) {
    nt_pulse->GetEntry(ientry);

    x[ientry] = var[0]+PulseNumber*1695;
    y[ientry] = var[1];
  }

  if (Graph != NULL) delete Graph;
  Graph = new TGraph(nent,x,y);
  return 0;
};


//-----------------------------------------------------------------------------
// drawing arrivals makes sense for  two proton pulses - then arrivals show up
// in between
//-----------------------------------------------------------------------------
void plot_proton_pulse_shape(int NPulses = 1, int DrawArrivals = 0, int Print = 0) {

  c = new TCanvas("c_pulse","pulse",1500,400);
  //  c->SetLogy(1);

  TH1F* h = new TH1F("h1","",600+1695*(NPulses-1),-300,300+1695*(NPulses-1));

  h->SetStats(0);
  h->SetMaximum(0.03);
  h->Draw();
  
  for (int i=0; i<NPulses; i++) g_pulse[i] = NULL;

  for (int i=0; i<NPulses; i++) {
    make_proton_pulse_graph(g_pulse[i],i);
    
    g_pulse[i]->SetFillStyle(1001);
    g_pulse[i]->SetFillColor(602);
    
    if (i == 0) {
      g_pulse[0]->SetTitle("");
      g_pulse[0]->GetXaxis()->SetTitle("Time, ns");
      g_pulse[0]->GetXaxis()->SetRangeUser(-300,300+1695*(NPulses-1));
      g_pulse[0]->GetXaxis()->SetLabelSize(0.06);
      g_pulse[0]->GetXaxis()->SetTitleSize(0.06);
      g_pulse[0]->GetXaxis()->SetTitleOffset(0.4);
    }
    g_pulse[i]->Draw("l,fill");
  }

  draw_label_ndc("Mu2e proton beam time structure",0.2,0.85,0.06,52);

  if (DrawArrivvals != 0) {
    g_pion_
  }

  if (Print != 0) {
    char name[100];
    sprintf(name,"proton_pulse_%i.eps",NPulses);
    c->Print(name);
  }
}


//-----------------------------------------------------------------------------
void plot_proton_pulse_with_ac_dipole(float Scale = 50, int Print = 0) {

  float x[10000], y[10000];

  TCanvas* c = new TCanvas("c_pulse","pulse",1300,700);
  c->SetLogy(1);
  
  TNtuple* nt_pulse (NULL);
  TNtuple* ac_dipole(NULL);

  create_ntuple("ConditionsService/data/potTimingDistribution_avg.txt",nt_pulse);

  create_ntuple("ConditionsService/data/acDipoleTransmissionFunction.txt",ac_dipole);

  float var[100];
//-----------------------------------------------------------------------------
// read proton pulse shape
//-----------------------------------------------------------------------------
  nt_pulse->SetBranchAddress("var0",&var[0]);
  nt_pulse->SetBranchAddress("var1",&var[1]);

  int nent = nt_pulse->GetEntries();

  for (int ientry=0; ientry<nent; ientry++) {
    nt_pulse->GetEntry(ientry);

    x[ientry] = var[0];
    y[ientry] = var[1]*Scale;
  }

  TGraph* gr_pulse = new TGraph(nent,x,y);
  gr_pulse->SetFillStyle(1001);
  gr_pulse->SetFillColor(624);
  gr_pulse->Draw("al,fill");

  gr_pulse->SetTitle("");
  gr_pulse->GetXaxis()->SetTitle("Time, ns");
  gr_pulse->GetXaxis()->SetLabelSize(0.05);
  gr_pulse->GetXaxis()->SetTitleSize(0.05);
  gr_pulse->GetXaxis()->SetTitleOffset(1);
  gr_pulse->GetYaxis()->SetRangeUser(1.e-8,10);
//-----------------------------------------------------------------------------
// read transmission function
//-----------------------------------------------------------------------------
  ac_dipole->SetBranchAddress("var0",&var[0]);
  ac_dipole->SetBranchAddress("var1",&var[1]);

  nent = ac_dipole->GetEntries();

  for (int ientry=0; ientry<nent; ientry++) {
    ac_dipole->GetEntry(ientry);

    x[ientry] = var[0];
    y[ientry] = var[1];
  }

  TGraph* gr_ac_dipole = new TGraph(nent,x,y);
  gr_ac_dipole->SetLineWidth(2.);
  gr_ac_dipole->Draw("same");


  draw_label_ndc("proton pulse",0.45,0.3,0.05,52);
  draw_label_ndc("AC dipole transmission",0.65,0.65,0.05,52);

  if (Print != 0) {
    char name[100];
    sprintf(name,"proton_pulse_with_ac_dipole.eps");
    c->Print(name);
  }
}

