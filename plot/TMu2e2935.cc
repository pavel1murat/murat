///////////////////////////////////////////////////////////////////////////////
// plots for mu2e-2935
// -------------------------
// directory with the histogram files: by default: $WORK_DIR/results
// can be redefined with "zzx.HistDir" in .rootrc
// 
// figures:
// ---------
// Fig.   1:  delta(T) for Set C tracks with clusters
// Fig.   2:  E/P      for Set C tracks with clusters 
//
// Fig.   3:  eff vs Z2, Mau7b  R(in) = 36 cm - optimization
// Fig.   4:  eff vs Z2, Mau7b  R(in) = 33 cm - optimization
//
// Fig.   5:  D0 for conversion electrons and cosmic muons
// Fig.   6:  Z0 for conversion electrons and cosmic muons
// Fig.   7:  occupancy for mixed events (vane-based MC) vs the calorimeter row 
// 
// Fig.  23: eff vs Z2, Mau8-uniform-in-DS, R(in) = 36 cm - optimization
// Fig.  24: eff vs Z2, Mau8-uniform-in-DS, R(in) = 33 cm
//
// Fig.  33: eff vs Z2, Mau8-full, R(in) = 36 cm - optimization
// Fig.  34: eff vs Z2, Mau8-full, R(in) = 33 cm - optimization
////////////////////////////////////////////////////////////////////////////////
#include "plot/TMu2e2935.hh"

#include "TPaveLabel.h"
#include "TArrow.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TProfile.h"
#include <iostream>

ClassImp(TMu2e2935)

//_____________________________________________________________________________
TMu2e2935::TMu2e2935(int PlotMode, int BlessingMode): TPlotNote(PlotMode,BlessingMode) {

  const char* hist_dir;
  char         res_dir[500];
  
  fFiguresDir   = Form("%s/mu2e_2935",gEnv->GetValue("mu2e.Figures","."));
  fWorkDir      = gSystem->Getenv("WORK_DIR");

  sprintf(res_dir,"%s/results",fWorkDir.Data());
  hist_dir      = gEnv->GetValue("mu2e.HistDir",res_dir);

  me000001_tcalm003 = Form("%s/v3_0_0/me000001.tcalm003.hist"     , hist_dir); 
  me000201_tcalm002 = Form("%s/v3_0_0/me000201_tcalm002.hist"     , hist_dir); 
  me000301_tcalm003 = Form("%s/v3_0_0/me000301.tcalm003_Mau8.hist", hist_dir); 
  me000401_tcalm003 = Form("%s/v3_0_0/me000401.tcalm003_Mau8.hist", hist_dir); 

  sprintf(res_dir,"%s",fWorkDir.Data());

  fIPlane1 = 10;

  fNPoints = 60;
  fStep    = 5.;   // 5cm

  fZMax    = fStep*fNPoints;
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
TMu2e2935::~TMu2e2935() {
};

//-----------------------------------------------------------------------------
const char* TMu2e2935::GetFiguresDir() {
  return fFiguresDir.Data();
}

void TMu2e2935::get_filename(int Figure, char* Filename) {

  if      (Figure == -1) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_dt");
  else if (Figure == -2) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_ep");
  else if (Figure == -3) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_dy");
  else if (Figure == -4) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_dz");
  else if (Figure == -5) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_d0");
  else if (Figure == -6) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_z0");
  else if (Figure == -100) sprintf(Filename,"fig_%03i_%s",Figure,"EoverP");
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
void TMu2e2935::remake_plots() {

  int fig [] = { 610, 

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
void TMu2e2935::plot(Int_t Figure, const char* CanvasName) {
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
//  Fig.1:  R=360
//-----------------------------------------------------------------------------
  if (Figure == 1) {

    int first_bin = 37;
    
    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    TH1F   *h_qtot, *h1; 
    TH2F   *h2f;  
    
    const int np = 60;
    float   eff [np], err[np], x[np], ex[np], eff2[np], err2[np];

    h_qtot  = gh1(me000001_tcalm003,"TCalm003","evt_0/ce_costh");
    double qtot = h_qtot->GetEntries();
    
    double step = 5.; // in cm
    double zmax = np*step;
//-----------------------------------------------------------------------------
//  efficiency of the first disk as function of its Z position
//-----------------------------------------------------------------------------
//    h_eff = new TH1F("h_eff","eff",1,0,150);

    for (int i=0; i<np; i++) {
      h1     = gh1(me000001_tcalm003,"TCalm003",Form("evt_0/rt_%02i",i));
      x  [i] = 5*i;
      ex[i]  = 0;
      eff[i] = h1->Integral(first_bin,100)/qtot;
      err[i] = sqrt(h1->Integral(first_bin,100))/qtot;
    }

    fP1->cd(1);
					// normalize axes

    h2f = new TH2F("h2f","h2f",1,0.,zmax,1,0.,1.);
    h2f->Draw();

    TGraphErrors* gre = new TGraphErrors(np, x,eff,ex,err);
    gre->Draw("APe,same");
    
    TF1 *f1 = new TF1("f1","[0]+[1]*sin([2]*x+[3])",0,zmax);
    f1->SetParameters(0.3,0.05,0.,0.);
    gre->Fit("f1","","ALP",0,zmax-10);
    
    double p2 = f1->GetParameter(2);
    double p3 = f1->GetParameter(3);
    
    double offset = -(p3-TMath::ACos(0))/p2;
    
    printf(" disk#1: offset = %10.3f\n",offset);
//-----------------------------------------------------------------------------
// position the second disk: evt_1 - events with R < 360 miss 
//-----------------------------------------------------------------------------
    fP1->cd(2);

    for (int i=0; i<np; i++) {
      h1     = gh1(me000001_tcalm003,"TCalm003",Form("evt_1/rt_%02i",i));
      x   [i] = 5*i;
      ex  [i]  = 0;
      eff2[i] = h1->Integral(first_bin,100)/qtot;
      err2[i] = sqrt(h1->Integral(first_bin,100))/qtot;
    }
  
    TGraphErrors* gre2 = new TGraphErrors(np, x,eff2,ex,err2);
    gre2->Draw("APe,same");
    
    f1->SetParameters(0.1,0.1,0.,0.);
    gre2->Fit("f1","","ALP",0,zmax-10.);
    
    p2     = f1->GetParameter(2);
    p3     = f1->GetParameter(3);
    
    offset = -(p3-TMath::ACos(0))/p2;
    
    printf(" disk#2: offset = %10.3f\n",offset);
  }
//-----------------------------------------------------------------------------
//  Fig.2: R=330
//-----------------------------------------------------------------------------
  if (Figure == 2) {

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    TH1F   *h_qtot, *h1; 
    TH2F   *h2f;  
    
    int     first_bin = 34; // R= 330

    const   int np    = 60;

    float   eff [np], err[np], x[np], ex[np], eff2[np], err2[np];

    h_qtot      = gh1(me000001_tcalm003,"TCalm003","evt_0/ce_costh");
    double qtot = h_qtot->GetEntries();
    
    double step = 5.; // in cm

    double zmax = np*step;
//-----------------------------------------------------------------------------
//  efficiency of the first disk
//-----------------------------------------------------------------------------
    for (int i=0; i<np; i++) {
      h1     = gh1(me000001_tcalm003,"TCalm003",Form("evt_0/rt_%02i",i));
      x  [i] = 5*i;
      ex[i]  = 0;
      eff[i] = h1->Integral(first_bin,100)/qtot;
      err[i] = sqrt(h1->Integral(first_bin,100))/qtot;
    }

    fP1->cd(1);
					// normalize axes

    h2f = new TH2F("h2f","h2f",1,0.,zmax,1,0.,1.);
    h2f->Draw();

    TGraphErrors* gre = new TGraphErrors(np, x,eff,ex,err);
    gre->Draw("APe,same");
    
    TF1 *f1 = new TF1("f1","[0]+[1]*sin([2]*x+[3])",0,zmax);
    f1->SetParameters(0.3,0.05,0.,0.);
    gre->Fit("f1","","ALP",0,zmax-10);
    
    double p2 = f1->GetParameter(2);
    double p3 = f1->GetParameter(3);
    
    double offset = -(p3-TMath::ACos(0))/p2;
    
    printf(" disk#1: offset = %10.3f\n",offset);
//-----------------------------------------------------------------------------
// position the second disk: evt_1 - events with R < 330 miss 
//-----------------------------------------------------------------------------
    fP1->cd(2);

    for (int i=0; i<np; i++) {
      h1     = gh1(me000001_tcalm003,"TCalm003",Form("evt_2/rt_%02i",i));
      x   [i] = 5*i;
      ex  [i]  = 0;
      eff2[i] = h1->Integral(first_bin,100)/qtot;
      err2[i] = sqrt(h1->Integral(first_bin,100))/qtot;
    }
  
    TGraphErrors* gre2 = new TGraphErrors(np, x,eff2,ex,err2);
    gre2->Draw("APe,same");
    
    f1->SetParameters(0.1,0.1,0.,0.);
    gre2->Fit("f1","","ALP",0,zmax-10.);
    
    p2 = f1->GetParameter(2);
    p3 = f1->GetParameter(3);
    
    offset = -(p3-TMath::ACos(0))/p2;
    
    printf(" disk#2: offset = %10.3f\n",offset);
  }
//-----------------------------------------------------------------------------
//  Fig.3: for a given first plane position , defined by fIPlane1 (0 to 59), 
//         plot efficiency of the second disk (left) 
//         and the total efficiency (right) as a function of the second disk position
//         R(in) = 36 cm
//-----------------------------------------------------------------------------
  if (Figure == 3) {

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    const int np = 60;

    //    int const kFirstBin=37;

    float   eff [np], err[np], x[np], ex[np], eff2[np], err2[np];
    //    TH1F    *h_qtot, *h0, *h1;
//-----------------------------------------------------------------------------
// normalize to the number of tracks - 1 track/event at max...
//-----------------------------------------------------------------------------
    calculate_eff_vs_z2(me000001_tcalm003,37,np,x,ex,eff,err,eff2,err2);

    fP1->cd(1);

    TH2F* h2f = new TH2F("h2f","h2f",1,0.,fZMax,1,0.,1.);
    h2f->Draw();
    
    TGraphErrors* gre = new TGraphErrors(np, x,eff,ex,err);
    gre->SetTitle("Disk 2 acceptance vs Z_2, Z_1 = 50 cm");
    gre->Draw("APe,same");
    
    gre->Fit("gaus","","ALP,same",fIPlane1*fStep+20,fIPlane1*fStep+110);

    fP1->cd(2);
    TGraphErrors* gre2 = new TGraphErrors(np, x,eff2,ex,err2);
    gre2->SetTitle("Total calorimeter acceptance vs Z_2, Z_1 = 50 cm");
    h2f->Draw();
    gre2->Draw("APe,same");
    gre2->Fit("gaus","","APe",fIPlane1*fStep+20,fIPlane1*fStep+110);
  }
//-----------------------------------------------------------------------------
//  Fig4: for a given first plane position , defined by fIPlane1 (0 to 59), R(in) = 33cm
//         plot efficiency of the second disk (left) 
//         and the total efficiency (right) as a function of the second disk position
//-----------------------------------------------------------------------------
  if (Figure == 4) {
    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    const int np = 60;

    float   eff [np], err[np], x[np], ex[np], eff2[np], err2[np];
//-----------------------------------------------------------------------------
// normalize to the number of tracks - 1 track/event at max...
//-----------------------------------------------------------------------------
    calculate_eff_vs_z2(me000001_tcalm003,34,np,x,ex,eff,err,eff2,err2);

    fP1->cd(1);

    TH2F* h2f = new TH2F("fig4_h2f","h2f",1,0.,fZMax,1,0.,1.);
    h2f->Draw();
    
    TGraphErrors* gre = new TGraphErrors(np, x,eff,ex,err);
    gre->SetTitle("Disk 2 acceptance vs Z_2, Z_1 = 50 cm");
    gre->Draw("APe,same");
    
    gre->Fit("gaus","","ALP,same",fIPlane1*fStep+20,fIPlane1*fStep+110);

    fP1->cd(2);
    TGraphErrors* gre2 = new TGraphErrors(np, x,eff2,ex,err2);
    gre2->SetTitle("Total calorimeter acceptance vs Z_2, Z_1 = 50 cm");
    h2f->Draw();
    gre2->Draw("APe,same");
    gre2->Fit("gaus","","APe",fIPlane1*fStep+20,fIPlane1*fStep+110);
  }
//-----------------------------------------------------------------------------
//  Fig.5: for a given Z1 - position of the first disk, find Z position of the
//         second disk maximizing the overall efficiency, plot result vs Z1
//         right plot - approximate distance between the disks (on 5cm grid)
//-----------------------------------------------------------------------------
  if (Figure == 5) {

    int const kFirstBin = 37;

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    const int np = 60;
    float   eff [np], err[np], x[np], ex[np], d2[np];
    float   eff_33[np], err_33[np] ; // , d2_33[np];
    double  eff1, eff2_max, eff_2, d2_max(-1.);
    TH1     *h0, *h_qtot, *h1;
//-----------------------------------------------------------------------------
// normalize to the number of tracks - 1 track/event at max...
//-----------------------------------------------------------------------------
    h_qtot  = gh1(me000001_tcalm003,"TCalm003","trk_0/tdip");

    double qtot = h_qtot->GetEntries();
    double step = 5.; // in cm

    //    double zmax = np*step;
    double zmax = 200.;

    for (int i1=0; i1<np; i1++) {
      h0       = gh1(me000001_tcalm003,"TCalm003",Form("evt_0/rt_%02i",i1));
      eff1     = h0->Integral(kFirstBin,100)/qtot;
      eff2_max = 0;

      for (int i2=i1+1; i2<np; i2++) {
	h1     = gh1(me000001_tcalm003,"TCalm003",Form("evt_0/rt0_%02i_%02i",i1,i2));
	eff_2  = h1->Integral(kFirstBin,100)/qtot;

	if (eff_2 > eff2_max) {
	  eff2_max  = eff_2;
	  d2_max    = (i2-i1)*5;
	}
      }
      eff[i1] = eff1+eff2_max;
      err[i1] = 0.01;
      x  [i1] = step*i1;
      ex [i1] = 0;
      d2 [i1] = d2_max;
    }
//-----------------------------------------------------------------------------
// now - efficiency for the 33 cm cut
//-----------------------------------------------------------------------------
    for (int i1=0; i1<np; i1++) {
      h0       = gh1(me000001_tcalm003,"TCalm003",Form("evt_0/rt_%02i",i1));
      eff1     = h0->Integral(34,100)/qtot;
      eff2_max = 0;

      for (int i2=i1+1; i2<np; i2++) {
	h1     = gh1(me000001_tcalm003,"TCalm003",Form("evt_0/rt33_0_%02i_%02i",i1,i2));
	eff_2  = h1->Integral(34,100)/qtot;

	if (eff_2 > eff2_max) {
	  eff2_max  = eff_2;
	  d2_max    = (i2-i1)*5;
	}
      }

      eff_33[i1] = eff1+eff2_max;
      err_33[i1] = 0.01;
      //      d2_33 [i1] = d2_max;
    }
    fP1->cd(1);

    TH2F* h2f = new TH2F("fig5_h2f","h2f",100,0.,zmax,1,0.8,1.2);
    h2f->GetXaxis()->SetRangeUser(0,zmax-1.e3);
    h2f->Draw();

    TGraphErrors* gre = new TGraphErrors(np, x,eff,ex,err);
    gre->SetTitle("Efficiency vs Z_{1}");
    gre->GetYaxis()->SetRangeUser(0.8,1.1);
    gre->SetMarkerStyle(20);
    gre->SetMarkerSize(1);
    gre->Draw("pe,same");
    
    TGraphErrors* gre_33 = new TGraphErrors(np, x,eff_33,ex,err_33);
    gre_33->SetTitle("Efficiency vs Z_{1}");
    gre_33->GetYaxis()->SetRangeUser(0.8,1.1);
    gre_33->SetMarkerStyle(24);
    gre_33->SetMarkerSize(1);
    gre_33->Draw("pe,same");
    
    fP1->cd(2);

    h2f = new TH2F("fig5_h2f2","h2f2",1,0.,zmax,1,0.,100.);
    h2f->GetXaxis()->SetRangeUser(0,zmax);
    h2f->Draw();

    TGraphErrors* gre2 = new TGraphErrors(np, x,d2,ex,err);
    gre2->SetTitle("Distance between the calorimeter disks");
    gre2->GetYaxis()->SetRangeUser(0.,100.);
    gre2->SetMarkerStyle(20);
    gre2->SetMarkerSize(1);
    gre2->Draw("APe,same");
  }
//-----------------------------------------------------------------------------
//  Fig.6:  R=33 cm
//-----------------------------------------------------------------------------
  if (Figure == 6) {
    int const kFirstBin = 34;

    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    const int np = 60;
    float   eff [np], err[np], x[np], ex[np], /* eff2[np], err2[np], */ d2[np];
    double  eff1, eff2_max, eff_2, d2_max(-1.);
    TH1     *h0, *h_qtot, *h1; // , *h2, *h3;
//-----------------------------------------------------------------------------
// normalize to the number of tracks - 1 track/event at max...
//-----------------------------------------------------------------------------
    h_qtot  = gh1(me000001_tcalm003,"TCalm003","trk_0/tdip");

    double qtot = h_qtot->GetEntries();
    double step = 5.; // in cm
					// above Z1 = 200cm everything becomes unreliable
    //    double zmax = np*step;
    double zmax = 200.;

    for (int i1=0; i1<np; i1++) {
//-----------------------------------------------------------------------------
// efficiency of the first disk
//-----------------------------------------------------------------------------
      h0       = gh1(me000001_tcalm003,"TCalm003",Form("evt_0/rt_%02i",i1));
      eff1     = h0->Integral(kFirstBin,100)/qtot;
      eff2_max = 0;

      for (int i2=i1+1; i2<np; i2++) {
	h1     = gh1(me000001_tcalm003,"TCalm003",Form("evt_0/rt33_0_%02i_%02i",i1,i2));
	eff_2  = h1->Integral(kFirstBin,100)/qtot;

	if (eff_2 > eff2_max) {
	  eff2_max  = eff_2;
	  d2_max    = (i2-i1)*5;
	}
      }

      eff[i1] = eff1+eff2_max;
      err[i1] = 0.01;
      x  [i1] = step*i1;
      ex [i1] = 0;
      d2 [i1] = d2_max;
    }

    fP1->cd(1);

    TH2F* h2f = new TH2F("fig5_h2f","h2f",1,0.,zmax,1,0.,1.);
    h2f->Draw();

    TGraphErrors* gre = new TGraphErrors(np, x,eff,ex,err);
    gre->SetTitle("Efficiency vs Z_{1}");
    gre->GetYaxis()->SetRangeUser(0.5,1.1);
    gre->SetMarkerStyle(20);
    gre->SetMarkerSize(1);
    gre->Draw("APe,same");
    
    fP1->cd(2);

    h2f = new TH2F("fig5_h2f2","h2f2",1,0.,zmax,1,0.,100.);
    h2f->GetXaxis()->SetRangeUser(0,zmax);
    h2f->Draw();

    TGraphErrors* gre2 = new TGraphErrors(np, x,d2,ex,err);
    gre2->SetTitle("Distance between the calorimeter disks");
    gre2->GetYaxis()->SetRangeUser(0.,100.);
    gre2->SetMarkerStyle(20);
    gre2->SetMarkerSize(1);
    gre2->Draw("APe,same");
  }
//-----------------------------------------------------------------------------
//  Fig.7 occupancy vs vane row for mixed events (DSID=me000201)
//-----------------------------------------------------------------------------
  if (Figure == 7) {
    TH2F* h2;
    TProfile* hpx;

    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    h2       = gh2(me000201_tcalm002,"TCalm002","evt_0/ncch_vs_row_1");

    fP1->cd(1);

    hpx = h2->ProfileX();
    hpx->SetStats(0);
    hpx->SetMinimum(0);
    hpx->SetMarkerStyle(20);
    hpx->SetTitle("N_{hits} vs Vane Row");
    hpx->GetXaxis()->SetRangeUser(0,10.99);
    hpx->GetXaxis()->SetTitle("crystal row, vane geometry");
    
    hpx->Draw();
  }
//-----------------------------------------------------------------------------
//  Fig.23: eff vs Z2, Mau8-uniform-in-DS, R(in) = 36 cm
//         for a given first plane position , defined by fIPlane1 (0 to 59), 
//         plot efficiency of the second disk (left) 
//         and the total efficiency (right) as a function of the second disk position
//-----------------------------------------------------------------------------
  if (Figure == 23) {

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    const int np = 60;

    float   eff [np], err[np], x[np], ex[np], eff2[np], err2[np];
//-----------------------------------------------------------------------------
// normalize to the number of tracks - 1 track/event at max...
//-----------------------------------------------------------------------------
    calculate_eff_vs_z2(me000301_tcalm003,37,np,x,ex,eff,err,eff2,err2);

    fP1->cd(1);

    double zmax = fZMax;
    TH2F* h2f = new TH2F("h2f","h2f",1,0.,zmax,1,0.,1.);
    h2f->Draw();
    
    TGraphErrors* gre = new TGraphErrors(np, x,eff,ex,err);
    gre->SetTitle("Disk 2 acceptance vs Z_2, Z_1 = 50 cm");
    gre->Draw("APe,same");
    
    gre->Fit("gaus","","ALP,same",fIPlane1*fStep+20,fIPlane1*fStep+110);

    fP1->cd(2);
    TGraphErrors* gre2 = new TGraphErrors(np, x,eff2,ex,err2);
    gre2->SetTitle("Total calorimeter acceptance vs Z_2, Z_1 = 50 cm");
    h2f->Draw();
    gre2->Draw("APe,same");
    gre2->Fit("gaus","","APe",fIPlane1*fStep+20,fIPlane1*fStep+110);
  }
//-----------------------------------------------------------------------------
//  Fig.24: eff vs Z2, Mau8-uniform-in-DS, R(in) = 33 cm
//         for a given first plane position , defined by fIPlane1 (0 to 59), 
//         plot efficiency of the second disk (left) 
//         and the total efficiency (right) as a function of the second disk position
//-----------------------------------------------------------------------------
  if (Figure == 24) {

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");

    const int np = 60;

    float   eff [np], err[np], x[np], ex[np], eff2[np], err2[np];
//-----------------------------------------------------------------------------
// normalize to the number of tracks - 1 track/event at max...
//-----------------------------------------------------------------------------
    int const kFirstBin=34;
    calculate_eff_vs_z2(me000301_tcalm003,kFirstBin,np,x,ex,eff,err,eff2,err2);

    fP1->cd(1);

    double zmax = fZMax;
    
    TH2F* h2f = new TH2F("h2f","h2f",1,0.,zmax,1,0.,1.);
    h2f->Draw();
    
    TGraphErrors* gre = new TGraphErrors(np, x,eff,ex,err);
    gre->SetTitle("Disk 2 acceptance vs Z_2, Z_1 = 50 cm");
    gre->Draw("APe,same");
    
    gre->Fit("gaus","","ALP,same",fIPlane1*fStep+20,fIPlane1*fStep+110);

    fP1->cd(2);
    TGraphErrors* gre2 = new TGraphErrors(np, x,eff2,ex,err2);
    gre2->SetTitle("Total calorimeter acceptance vs Z_2, Z_1 = 50 cm");
    h2f->Draw();
    gre2->Draw("APe");
    gre2->Fit("gaus","","APe",fIPlane1*fStep+20,fIPlane1*fStep+110);
  }
//-----------------------------------------------------------------------------
//  Fig 33: eff vs Z2, Mau8-full R(in) = 36 cm
//          for a given first plane position , defined by fIPlane1 (0 to 59), 
//          plot efficiency of the second disk (left) 
//          and the total efficiency (right) as a function of the second disk position
//-----------------------------------------------------------------------------
  if (Figure == 33) {
    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    const int np = 60;
    float     eff [np], err[np], x[np], ex[np], eff2[np], err2[np];
//-----------------------------------------------------------------------------
// normalize to the number of tracks - 1 track/event at max...
//-----------------------------------------------------------------------------
    calculate_eff_vs_z2(me000401_tcalm003,37,np,x,ex,eff,err,eff2,err2);
   
    fP1->cd(1);

    double zmax = fZMax;
    TH2F* h2f = new TH2F("fig4_h2f","h2f",1,0.,zmax,1,0.,1.);
    h2f->Draw();
    
    TGraphErrors* gre = new TGraphErrors(np,x,eff,ex,err);
    gre->SetTitle("Disk 2 acceptance vs Z_2, Z_1 = 50 cm");
    gre->Draw("APe,same");
    
    gre->Fit("gaus","","ALP,same",fIPlane1*fStep+20,fIPlane1*fStep+110);
    gre->SetTitle("Disk 2 acceptance vs Z_2, Z_1 = 50 cm");
    gre->Draw("APe,same");
    gre->Fit("gaus","","ALP,same",fIPlane1*fStep+20,fIPlane1*fStep+110);

    fP1->cd(2);
    TGraphErrors* gre2 = new TGraphErrors(np, x,eff2,ex,err2);
    gre2->SetTitle("Total calorimeter acceptance vs Z_2, Z_1 = 50 cm");
    h2f->Draw();
    gre2->Draw("APe");
    gre2->Fit("gaus","","APe",fIPlane1*fStep+20,fIPlane1*fStep+110);
  }
//-----------------------------------------------------------------------------
//  Fig 34: eff vs Z2, Mau8-uniform-in-DS, R(in) = 33 cm - first bin 34
//          for a given first plane position , defined by fIPlane1 (0 to 59), 
//          plot efficiency of the second disk (left) 
//          and the total efficiency (right) as a function of the second disk position
//-----------------------------------------------------------------------------
  if (Figure == 34) {

    fCanvas  = new_slide(name,title,2,1,1200,600);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    const int np = 60;

    float   eff [np], err[np], x[np], ex[np], eff2[np], err2[np];
//-----------------------------------------------------------------------------
// normalize to the number of tracks - 1 track/event at max...
//-----------------------------------------------------------------------------
    calculate_eff_vs_z2(me000401_tcalm003,34,np,x,ex,eff,err,eff2,err2);
   
    fP1->cd(1);

    double zmax = fZMax;
    TH2F* h2f = new TH2F("fig4_h2f","h2f",1,0.,zmax,1,0.,1.);
    h2f->Draw();
    
    TGraphErrors* gre = new TGraphErrors(np,x,eff,ex,err);
    gre->SetTitle("Disk 2 acceptance vs Z_2, Z_1 = 50 cm");
    gre->Draw("APe,same");
    gre->Fit("gaus","","ALP,same",fIPlane1*fStep+20,fIPlane1*fStep+110);

    fP1->cd(2);
    TGraphErrors* gre2 = new TGraphErrors(np, x,eff2,ex,err2);
    gre2->SetTitle("Total calorimeter acceptance vs Z_2, Z_1 = 50 cm");
    h2f->Draw();
    gre2->Draw("APe");
    gre2->Fit("gaus","","APe",fIPlane1*fStep+20,fIPlane1*fStep+110);
  }

//-----------------------------------------------------------------------------
//  Fig.100: distribution of the E/P ratio for downstream tracks
//           which are matched with an EMC with energy>10MeV and
//           which have a distance from the EMC (for each coordinate)
//           less than 15 mm (half crystal size)
//-----------------------------------------------------------------------------
  if (Figure == 100) {
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
  }

}



//-----------------------------------------------------------------------------
void TMu2e2935::calculate_eff_vs_z2(const char*   Filename, 
				    int FirstBin, int NPoints, 
				    float*   X , float* Ex  , 
				    float* Eff , float* Err ,
				    float* Eff2, float* Err2) {
  TH1F* h1;

  TH1F* h_qtot  = gh1(Filename,"TCalm003","trk_0/tdip");
  
  double qtot = h_qtot->GetEntries();

  TH1F* h0    = gh1(Filename,"TCalm003",Form("evt_0/rt_%02i",fIPlane1));
  
  double eff1 = h0->Integral(FirstBin,100)/qtot;
  
  for (int i2=0; i2<fIPlane1+1; i2++) {
    X [i2] = 5*i2;
    Ex[i2]  = 0;
    Eff[i2] = 0.;
    Err[i2] = 0.;
    
    Eff2[i2] = 0.;
    Err2[i2] = 0.;
  }

  for (int i2=fIPlane1+1; i2<NPoints; i2++) {
    h1     = gh1(Filename,"TCalm003",Form("evt_0/rt33_0_%02i_%02i",fIPlane1,i2));
    X [i2] = 5*i2;
    Ex[i2]  = 0;
    Eff[i2] = h1->Integral(FirstBin,100)/qtot;
    Err[i2] = sqrt(h1->Integral(FirstBin,100))/qtot;
    
    Eff2[i2] = h1->Integral(FirstBin,100)/qtot  + eff1;
    Err2[i2] = sqrt(h1->Integral(FirstBin,100)+h0->Integral(FirstBin,100))/qtot;
  }

  printf("eff(plane=%2i) = %8.3f\n",fIPlane1,eff1);  
}
