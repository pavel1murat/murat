///////////////////////////////////////////////////////////////////////////////
// plots for Mu2e-3722: Mu2e sensitivity studies
// -------------------------
// directory with the histogram files: by default: $WORK_DIR/results
// can be redefined with "mu2e.HistDir" in .rootrc
// 
// figures:
// ---------
// Fig.  1: expected momentum spectrum at 105 MeV - smearing
// Fig.  2: momentum spectrum of the
// fig.  3: track reco efficiency vs track momentum
//
// Fig. 11: plot lower limit
// Fig. 12: plot upper limit
//
// Fig. 21: MCLIMIT plot for the null hypothesis
// Fig. 22: MCLIMIT plot for the test hypothesis
//
// Fig.1080: reconstructed momentum spectrum of  80 MeV electrons
// Fig.1087: reconstructed momentum spectrum of  87 MeV electrons
// Fig.1090: reconstructed momentum spectrum of  90 MeV electrons
// Fig.1095: reconstructed momentum spectrum of  95 MeV electrons
// Fig.1100: reconstructed momentum spectrum of 100 MeV electrons
// Fig.1105: reconstructed momentum spectrum of 105 MeV electrons
// 
////////////////////////////////////////////////////////////////////////////////
#include "plot/TMu2e3722.hh"

#include "TPaveLabel.h"
#include "TArrow.h"
#include "TObjArray.h"
#include "TGraphErrors.h"
#include "TObjString.h"
#include <iostream>

ClassImp(TMu2e3722)

//_____________________________________________________________________________
TMu2e3722::TMu2e3722(int PlotMode, int BlessingMode): TPlotNote(PlotMode,BlessingMode) {

  const char* hist_dir;
  char         res_dir[500];
  
  fFiguresDir   = Form("%s/mu2e_3722",gEnv->GetValue("mu2e.Figures","."));
  fWorkDir      = gSystem->Getenv("WORK_DIR");

  sprintf(res_dir,"%s/results",fWorkDir.Data());
  hist_dir      = gEnv->GetValue("mu2e.HistDir",res_dir);

  sprintf(f_conv01,"%s/limits/conv01.mu2e_limits.hist", hist_dir); 
  sprintf(f_cosm01,"%s/limits/cosm01.mu2e_limits.hist", hist_dir); 
  sprintf(f_data01,"%s/limits/data01.mu2e_limits.hist", hist_dir); 
  sprintf(f_dio01 ,"%s/limits/dio01.mu2e_limits.hist" , hist_dir); 
  sprintf(f_dio02 ,"%s/limits/dio02.mu2e_limits.hist" , hist_dir); 

  sprintf(f_080,"%s/dev2/egun_080_080_001_disk_670_330_700_tcalm002.hist",hist_dir);
  sprintf(f_087,"%s/dev2/egun_087_087_001_disk_670_330_700_tcalm002.hist",hist_dir);
  sprintf(f_090,"%s/dev2/egun_090_090_001_disk_670_330_700_tcalm002.hist",hist_dir);
  sprintf(f_095,"%s/dev2/egun_095_095_001_disk_670_330_700_tcalm002.hist",hist_dir);
  sprintf(f_100,"%s/dev2/egun_100_100_001_disk_670_330_700_tcalm002.hist",hist_dir);
  sprintf(f_105,"%s/dev2/egun_105_105_001_disk_670_330_700_tcalm002.hist",hist_dir);

  sprintf(res_dir,"%s",fWorkDir.Data());
  hist_dir      = gEnv->GetValue("mu2e.HistDir",res_dir);

					// integral of the crystal ball function
  fNConvEle = 1.0;
//-----------------------------------------------------------------------------
// different settings (colors etc) for talk, note and paper modes
//-----------------------------------------------------------------------------
  // if (PlotMode == kNoteMode) {
  // }
  // else if (PlotMode == kTalkMode) {
  // }
  // else if (PlotMode == kPaperMode) {
  // }

};

//_____________________________________________________________________________
TMu2e3722::~TMu2e3722() {
};

//-----------------------------------------------------------------------------
const char* TMu2e3722::GetFiguresDir() {
  return fFiguresDir.Data();
}

void TMu2e3722::get_filename(int Figure, char* Filename) {

  if      (Figure ==   1) sprintf(Filename,"fig_%03i_%s",Figure,"momentum");
  else if (Figure == -12) sprintf(Filename,"fig_%03i_%s",Figure,"prt_energy");
  else {
//-----------------------------------------------------------------------------
// undefined, just figure
//-----------------------------------------------------------------------------
    sprintf(Filename,"fig_%i",Figure);
  }
}
//-----------------------------------------------------------------------------
// determine name of the .eps and .gif files in the ./figures directory
// today's default set is 'frr_03', all the plots should be coming from it
//-----------------------------------------------------------------------------
void TMu2e3722::remake_plots() {

  int fig [] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
		 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
		 21, 22, 23,
		 103, 113, 123,
		 -1
  };

  int figi;

  for (int i=0; fig[i]>0; i++) {
    figi = fig[i];
    plot (figi); 
    print(figi); 
  }
}


//-----------------------------------------------------------------------------
double TMu2e3722::f_crystal_ball(double* X, double* P) {
  double f, alpha, n, a, b;
  double dx = (X[0]-P[1])/P[2];

  alpha = P[3];
  n     = P[4];

  if (dx > -alpha) {
    f = P[0]*TMath::Exp(-dx*dx/2.);
  }
  else {
    a = pow((n/fabs(alpha)),n)*TMath::Exp(-alpha*alpha/2);
    b = n/fabs(alpha)-fabs(alpha);
    f = P[0]*a*pow((b-dx),-n);
  }

  return f;
}

//-----------------------------------------------------------------------------
double TMu2e3722::f_eff_mom(double* X, double* P) {
  double f, dx;

  dx = (X[0]-P[1])/P[2];

  if (dx < 0) {
    f = 0;
  }
  else {
    f = P[0]*(1-TMath::Exp(-dx));
  }

  return f;
}

//-----------------------------------------------------------------------------
// return the track reconstruction efficiency vs the track momentum
// the numbers are obtained with 'dev2' release (-D 2013-09-04) and corrected 
// for the finite stopping target size by the ratio of CE/singleE(105 MeV/c) 
// efficiencies 
//-----------------------------------------------------------------------------
void TMu2e3722::eff_vs_mom(TGraphErrors** G) {
  int const np = 6;

  double dio_sf = 4889./5591.;
					// momenta
  double p[np]  = {  80., 87.,     90.,   95.,  100.,  105.};
					// numbers of Set C tracks 
  double qn[np] = {   7., 1368., 2784., 4459., 5473., 5591.};

  double ex[np], y[np], ey[np];
					// total number of generated events  - 10000,
					// -0.6 < cos_the < 0.8 gives a factor of 1.4/2.0 = 0.7 
  for (int i=0; i<np; i++) {
    y [i] = qn[i]/(10000*2./1.4)*dio_sf;
    ex[i] = 0;
    ey[i] = 0.02; // y[i]/sqrt(qn[i]);
  }

  (*G) = new TGraphErrors(np,p,y,ex,ey);
}

//-----------------------------------------------------------------------------
// make asymmetric graphs for plotting the sensitivity limits
//-----------------------------------------------------------------------------
void TMu2e3722::make_graphs(double* Data, TGraphAsymmErrors **G1, TGraphAsymmErrors **G2) {

  int n = 0;

  double x[100], y[100], exl[100], exh[100], eyl1[100], eyh1[100], eyl2[100], eyh2[100];

  int loc;

  for (int i=0; Data[6*i] > 0;  i++) {
    loc = 6*i;

    x[i] = Data[loc];
    y[i] = Data[loc+3]*fNConvEle;

    if (Data[loc+6] > 0) {
      exh[i] = (Data[loc+6]-Data[loc])/2;
      if (loc == 0) {
	exl[i] = 0;
      }
      else {
	exl[i] = exh[i-1];
      }
    }
    else {
					// last point
      exh[i] = 0;
      exl[i] = exh[i-1];
    }

    exl[i] = 0;
    exh[i] = 0;

    eyl1[i] = 0; // y[i]-Data[loc+2];
    eyh1[i] = 0; // Data[loc+4]-y[i];

    eyl2[i] = 0; // y[i]-Data[loc+1];
    eyh2[i] = 0; // Data[loc+5]-y[i];

    n += 1;

    printf("%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",
	   x[i],y[i],exl[i],exh[i],eyl1[i],eyh1[i],eyl2[i],eyl2[i]);
  }

  *G1 = new TGraphAsymmErrors(n,x,y,exl,exh,eyl1,eyh1);
  *G2 = new TGraphAsymmErrors(n,x,y,exl,exh,eyl2,eyh2);

}

//-----------------------------------------------------------------------------
void TMu2e3722::plot(Int_t Figure, const char* CanvasName) {
//-----------------------------------------------------------------------------
//  make sure all the needed scripts are loaded
//-----------------------------------------------------------------------------
  char        macro[200];

  const char* script[] = { 
//    "plot/make_fake_rate_hist.C",
//    "plot/make_ratio_hist.C",
//    "plot/plot_ntrk10.C",
    0 
  };

  const char* work_dir = gSystem->Getenv("PWD");

  for (int i=0; script[i] != 0; i++) {
    sprintf(macro,"%s/murat/%s",work_dir,script[i]);
    if (! gInterpreter->IsLoaded(macro)) gInterpreter->LoadMacro(macro);
  }
//-----------------------------------------------------------------------------
  char        name[200], title[200];

  //  int old_opt_stat = gStyle->GetOptStat();
  
  if ((! CanvasName) || (CanvasName[0] == '0')) {
    sprintf(name,"fig_%i",Figure);
    sprintf(title,"figure_%i",Figure);
  }
  else {
    sprintf(name,"%s",CanvasName);
    sprintf(title,"figure %i",Figure);
  }
//-----------------------------------------------------------------------------
// Fig.  1: DIO energy (ibgr = 1)
//-----------------------------------------------------------------------------
  if (Figure == 1) {
    TH1D   *h1, *h2;
    fCanvas  = new_slide(name,title,1,1,1100,800);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    h1 = (TH1D*) gh1(f_dio01,"Mu2eLimits","trk_1/p");
    h2 = (TH1D*) gh1(f_dio02,"Mu2eLimits","trk_1/p");

    h2->GetYaxis()->SetRangeUser(1.e-8,0.2);
    h2->GetXaxis()->SetRangeUser(101,106);
    h2->Draw();
    h1->Draw("same");
  }
//-----------------------------------------------------------------------------
// Fig.  2: DIO vs converison electrons
//-----------------------------------------------------------------------------
  if (Figure == 2) {
    TH1D   *h1, *h2;
    fCanvas  = new_slide(name,title,1,1,1100,800);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    h1 = (TH1D*) gh1(f_dio02 ,"Mu2eLimits","trk_1/p");
    h2 = (TH1D*) gh1(f_conv01,"Mu2eLimits","trk_1/p");

    h2->GetYaxis()->SetRangeUser(1.e-8,0.5);
    h2->GetXaxis()->SetRangeUser(101,106);
    h1->Draw();
    h2->Draw("same");
  }
//-----------------------------------------------------------------------------
// Fig.  3: track reco efficiency vs momentum
//-----------------------------------------------------------------------------
  if (Figure == 3) {
    //    TH1D          *h1, *h2;
    TGraphErrors  *g;

    fCanvas  = new TCanvas(name,title,1000,700);

    fCanvas->cd();

    eff_vs_mom(&g);

    TF1* fit_eff = new TF1("fit_eff",f_eff_mom,94,106,3);
    
    fit_eff->SetParameter(0,0.40);
    fit_eff->SetParameter(1,80.);
    fit_eff->SetParameter(2,10.);

    g->Fit("fit_eff");
    g->SetTitle("Track reconstruction efficiency vs momentum");
    g->GetXaxis()->SetTitle("P (MeV/c)");
    g->Draw("ALP");
  }
//-----------------------------------------------------------------------------
// Fig.  11: sensitivity vs lower limit
//-----------------------------------------------------------------------------
  if (Figure == 11) {
    //    TH1D   *h1, *h2;
    fCanvas  = new TCanvas(name,title,1000,600);

    fCanvas->cd();

    double data[] = { // 90% CL lower limit
      102.00,   0.40753,   0.43430,   0.49896,   0.69210,   0.89419,
      102.20,   0.41419,   0.43821,   0.49108,   0.69656,   0.98080,
      102.40,   0.42184,   0.43731,   0.47825,   0.65203,   0.90417,
      102.60,   0.42726,   0.43842,   0.47994,   0.69928,   0.94021,
      102.80,   0.43035,   0.43718,   0.48746,   0.68653,   0.89646,
      103.00,   0.43606,   0.43606,   0.47875,   0.68411,   0.95873,
      103.10,   0.44287,   0.44287,   0.48535,   0.70436,   0.93870,
      103.20,   0.45240,   0.45240,   0.45611,   0.71296,   0.95328,
      103.25,   0.45839,   0.45839,   0.45839,   0.70574,   0.99245,
      103.25,   0.45830,   0.45830,   0.45830,   0.70180,   0.97913,
      103.30,   0.46507,   0.46507,   0.46507,   0.68547,   0.93740,
      103.50,   0.50366,   0.50366,   0.50366,   0.74841,   1.05128,
      103.60,   0.53342,   0.53342,   0.53342,   0.72942,   0.88497,
      103.80,   0.63784,   0.63784,   0.63784,   0.63784,   1.06116,
      104.00,   0.92844,   0.92844,   0.92844,   0.92844,   1.55243,
      -1
    };

    // double data[] = { // 95% CL lower limit
    //   102.00, 0.47058,   0.50766 ,  0.57184  ,  0.77512 ,  1.01867 ,
    //   102.50, 0.49714,   0.51434 ,  0.57499  ,  0.78270 ,  1.05134 ,
    //   102.60, 0.50084,   0.51467 ,  0.56311  ,  0.76529 ,  1.01099 ,
    //   102.80, 0.50560,   0.51429 ,  0.56871  ,  0.77358 ,  0.98065 ,
    //   102.90, 0.50667,   0.51352 ,  0.56259  ,  0.78019 ,  1.02713 ,
    //   103.00, 0.51121,   0.51121 ,  0.56219  ,  0.76254 ,  1.00466 ,
    //   103.01, 0.51184,   0.51184 ,  0.56894  ,  0.77200 ,  1.01177 ,
    //   103.20, 0.52989,   0.52989 ,  0.55239  ,  0.79894 ,  1.04101 ,
    //   103.21, 0.53118,   0.53118 ,  0.55290  ,  0.79898 ,  1.05356 ,
    //   103.25, 0.53674,   0.53674 ,  0.53674  ,  0.79886 ,  1.04385 ,
    //   103.30, 0.54461,   0.54461 ,  0.54461  ,  0.77570 ,  1.01970 ,
    //   103.40, 0.56400,   0.56400 ,  0.56400  ,  0.81076 ,  1.12310 ,
    //   103.50, 0.58968,   0.58968 ,  0.58968  ,  0.83669 ,  1.14230 ,
    //   103.60, 0.62450,   0.62450 ,  0.62450  ,  0.76985 ,  0.97574 ,
    //   -1,
    // };

    TGraphAsymmErrors  *g1, *g2;

    make_graphs(data,&g1,&g2);

    g1->SetTitle("Sensitivity vs P_{min}");
    g2->SetTitle("Sensitivity vs P_{min}");
    
    g1->GetXaxis()->SetTitle("P_{min} (MeV/c)");
    g2->GetXaxis()->SetTitle("P_{min} (MeV/c)");
    
    g1->GetYaxis()->SetTitle("Sensitivity, 90% CL (#times 10^{-16}),");
    g2->GetYaxis()->SetTitle("Sensitivity, 90% CL (#times 10^{-16}),");
    
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(1);
    
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(1);
    
    g2->SetFillColor(kBlue-7);
    g2->SetFillStyle(3001);
    //  g2->Draw("ALE3");
    g2->Draw("ALP");

    // g1->SetFillColor(kRed);
    // g1->SetFillStyle(3001);
    // g1->SetLineWidth(3);
    // g1->Draw("EL3,same");
  }
//-----------------------------------------------------------------------------
// Fig.  12: sensitivity vs upper limit (90% CL)
//-----------------------------------------------------------------------------
  if (Figure == 12) {
    TGraphAsymmErrors  *g1, *g2;
    //    TH1D               *h1, *h2;

    fCanvas  = new TCanvas(name,title,1000,600);

    fCanvas->cd();
    double data[] = {
      104.0,   0.87626,   0.87626,   0.87626,  1.32187,  1.81316,
      104.2,   0.57786,   0.57786,   0.57786,  0.84861,  1.17150,
      104.4,   0.47879,   0.47879,   0.47879,  0.69401,  0.94140,
      104.6,   0.45558,   0.45558,   0.45558,  0.69723,  0.89079,
      104.8,   0.45249,   0.45249,   0.45249,  0.70543,  0.94330,
      105.0,   0.45240,   0.45240,   0.45611,  0.71296,  0.95328,
      105.2,   0.45252,   0.45252,   0.45252,  0.70743,  0.98064,
      105.4,   0.44951,   0.45266,   0.45266,  0.70778,  0.90961,
      105.6,   0.44782,   0.45279,   0.45279,  0.69625,  0.94890,
      105.8,   0.44776,   0.45293,   0.45293,  0.69665,  0.89411,
      106.0,   0.44789,   0.45306,   0.45306,  0.70001,  0.91767,
      -1
    };

    // double data[] = {  // 95% CL
    //   104.0   , 0.95050,  0.95050,   1.09512,   1.47119,   1.98562,  
    //   104.2   , 0.64498,  0.64498,   0.72418,   0.97650,   1.20231,  
    //   104.4   , 0.53959,  0.53959,   0.58985,   0.81145,   1.07142,  
    //   104.5   , 0.52197,  0.52197,   0.57741,   0.79872,   1.01729,  
    //   104.7   , 0.51202,  0.51202,   0.56680,   0.77574,   1.08880,  
    //   104.8   , 0.51128,  0.51128,   0.56229,   0.77108,   1.01240,  
    //   105.0   , 0.51121,  0.51121,   0.56219,   0.76254,   1.00466,  
    //   105.5   , 0.51157,  0.51157,   0.57825,   0.79305,   1.01255,  
    //   106.0   , 0.51196,  0.51196,   0.56816,   0.77893,   1.05164,  
    //   106.5   , 0.51196,  0.51196,   0.56816,   0.77893,   1.05164,  
    //   -1
    // };

    make_graphs(data,&g1,&g2);

    g1->SetTitle("Sensitivity vs P_{max}");
    g2->SetTitle("Sensitivity vs P_{max}");
    
    g1->GetXaxis()->SetTitle("P_{max} (MeV/c)");
    g2->GetXaxis()->SetTitle("P_{max} (MeV/c)");
    
    g1->GetYaxis()->SetTitle("Sensitivity, 90% CL (#times 10^{-16}),");
    g2->GetYaxis()->SetTitle("Sensitivity, 90% CL (#times 10^{-16}),");
    
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(1);
    
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(1);
    
    g2->SetFillColor(kBlue-7);
    g2->SetFillStyle(3001);
    //  g2->Draw("ALE3");
    g2->Draw("ALP");

    // g1->SetFillColor(kRed);
    // g1->SetFillStyle(3001);
    // g1->SetLineWidth(3);
    // g1->Draw("EL3,same");
  }
//-----------------------------------------------------------------------------
// Fig.  21: 
// mu2e_limits(1,0.9,103.25,104.8,1000,0)
// gLMu2e->PlotHypothesis("mu2e","test")
//-----------------------------------------------------------------------------
  if (Figure == 21) {
    fCanvas  = new_slide(name,title,1,1,1100,700);
  }
//-----------------------------------------------------------------------------
// Fig.  22: 
// mu2e_limits(1,0.9,103.25,104.8,1000,0)
// gLMu2e->PlotHypothesis("mu2e","null")
//-----------------------------------------------------------------------------
  if (Figure == 22) {
    fCanvas  = new_slide(name,title,1,1,1100,700);
  }
//-----------------------------------------------------------------------------
// Fig.  1080: fit 80 MeV electrons
//-----------------------------------------------------------------------------
  if (Figure == 1080) {
    fCanvas  = new_slide(name,title,1,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    fH1[0] = (TH1F*) gh1(f_080,"TCalm002","trk_1/p2")->Clone("h1_0_fig_1080");
    
    double XMin, XMax;

    XMin =  65.;
    XMax =  87.;

    TF1* f = new TF1("f_crystal_ball",TMu2e3722::f_crystal_ball,90,100,5);
    
    f->SetParameter(0,1000);
    f->SetParameter(1,79.);
    f->SetParameter(2,0.3);
    f->SetParameter(3,2);
    f->SetParameter(4,1);

    f->SetParLimits(2,0,10);
    f->SetParLimits(4,0,10);

    fH1[0]->GetXaxis()->SetRangeUser(70,100);
    fH1[0]->Fit(f,"","",XMin,XMax);
  }
//-----------------------------------------------------------------------------
// Fig.  1087: fit 87 MeV electrons
//-----------------------------------------------------------------------------
  if (Figure == 1087) {
    fCanvas  = new_slide(name,title,1,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    fH1[0] = (TH1F*) gh1(f_087,"TCalm002","trk_1/p2")->Clone("h1_0_fig_1080");
    
    double XMin, XMax;

    XMin =  65.;
    XMax =  87.;

    TF1* f = new TF1("f_crystal_ball",TMu2e3722::f_crystal_ball,90,100,5);
    
    f->SetParameter(0,1000);
    f->SetParameter(1,86.);
    f->SetParameter(2,0.3);
    f->SetParameter(3,2);
    f->SetParameter(4,1);

    f->SetParLimits(2,0,10);
    f->SetParLimits(4,0,10);

    fH1[0]->GetXaxis()->SetRangeUser(77,107);
    fH1[0]->Fit(f,"","",XMin,XMax);
  }
//-----------------------------------------------------------------------------
// Fig.  1090: fit 90 MeV electrons
//-----------------------------------------------------------------------------
  if (Figure == 1090) {
    fCanvas  = new_slide(name,title,1,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    fH1[0] = (TH1F*) gh1(f_090,"TCalm002","trk_1/p2")->Clone("h1_0_fig_1090");
    
    double XMin, XMax;

    XMin =  80.;
    XMax =  97.;

    TF1* f = new TF1("f_crystal_ball",TMu2e3722::f_crystal_ball,80,110,5);
    
    f->SetParameter(0,1000);
    f->SetParameter(1,89.);
    f->SetParameter(2,0.3);
    f->SetParameter(3,2);
    f->SetParameter(4,1);

    f->SetParLimits(4,0,10);

    fH1[0]->GetXaxis()->SetRangeUser(80,110);
    fH1[0]->Fit(f,"","",XMin,XMax);
  }
//-----------------------------------------------------------------------------
// Fig.1095: fit 95 MeV electrons
//-----------------------------------------------------------------------------
  if (Figure == 1095) {
    fCanvas  = new_slide(name,title,1,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    fH1[0] = (TH1F*) gh1(f_095,"TCalm002","trk_1/p2")->Clone("h1_0_fig_1095");
    
    double XMin, XMax;

    XMin =  85.;
    XMax =  97.;

    TF1* f = new TF1("f_crystal_ball",TMu2e3722::f_crystal_ball,80,110,5);
    
    f->SetParameter(0,1000);
    f->SetParameter(1,94.);
    f->SetParameter(2,0.3);
    f->SetParameter(3,2);
    f->SetParameter(4,1);

    fH1[0]->GetXaxis()->SetRangeUser(80,110);
    fH1[0]->Fit(f,"","",XMin,XMax);
  }
//-----------------------------------------------------------------------------
// Fig.1100: fit 100 MeV electrons
//-----------------------------------------------------------------------------
  if (Figure == 1100) {
    fCanvas  = new_slide(name,title,1,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    fH1[0] = (TH1F*) gh1(f_100,"TCalm002","trk_1/p2")->Clone("h1_0_fig_1100");
    
    double XMin, XMax;

    XMin =  90.;
    XMax = 102.;

    TF1* f = new TF1("f_crystal_ball",TMu2e3722::f_crystal_ball,80,110,5);
    
    f->SetParameter(0,1000);
    f->SetParameter(1,99.);
    f->SetParameter(2,0.3);
    f->SetParameter(3,2);
    f->SetParameter(4,1);

    fH1[0]->GetXaxis()->SetRangeUser(85,115);
    fH1[0]->Fit(f,"","",XMin,XMax);
  }
//-----------------------------------------------------------------------------
// Fig.1105: fit 105 MeV electrons
//-----------------------------------------------------------------------------
  if (Figure == 1105) {
    fCanvas  = new_slide(name,title,1,1,1100,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fP1->cd(1);

    fH1[0] = (TH1F*) gh1(f_105,"TCalm002","trk_1/p2")->Clone("h1_0_fig_1105");
    
    double XMin, XMax;

    XMin =  95.;
    XMax = 107.;

    TF1* f = new TF1("f_crystal_ball",TMu2e3722::f_crystal_ball,80,110,5);
    
    f->SetParameter(0,1000);
    f->SetParameter(1,104.);
    f->SetParameter(2,0.3);
    f->SetParameter(3,2);
    f->SetParameter(4,1);

    fH1[0]->GetXaxis()->SetRangeUser(90,120);
    fH1[0]->Fit(f,"","",XMin,XMax);
  }
}

