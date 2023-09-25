///////////////////////////////////////////////////////////////////////////////
// plots for cdf-10008
// -------------------------
// directory with the histogram files: by default: $WORK_DIR/results
// can be redefined with "zzx.HistDir" in .rootrc
// 
// figures:
// ---------
// Fig.   1:  not yet defined
//  Fig.11: Cluster time 
////////////////////////////////////////////////////////////////////////////////
#include "plot/TMu2e3189.hh"

#include "TPaveLabel.h"
#include "TArrow.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <iostream>

ClassImp(TMu2e3189)

//_____________________________________________________________________________
TMu2e3189::TMu2e3189(int PlotMode, int BlessingMode): TPlotNote(PlotMode,BlessingMode) {

  const char* hist_dir;
  char         res_dir[500];
  
  fFiguresDir   = Form("%s/mu2e_3189",gEnv->GetValue("mu2e.Figures","."));
  fWorkDir      = gSystem->Getenv("WORK_DIR");

  sprintf(res_dir,"%s/results",fWorkDir.Data());
  hist_dir      = gEnv->GetValue("mu2e.HistDir",res_dir);

  dio_50_105_50b = Form("%s/mu2e_3189/decay_in_orbit_50_105_50b.hist", hist_dir); 
  dio_00_050_80c = Form("%s/mu2e_3189/decay_in_orbit_00_050_80c.hist", hist_dir); 
//-----------------------------------------------------------------------------
// different settings (colors etc) for talk, note and paper modes
//-----------------------------------------------------------------------------
  // if (PlotMode == kNoteMode) {
  // }
  // else if (PlotMode == kTalkMode) {
  // }
  // else if (PlotMode == kPaperMode) {
  // }

}

//_____________________________________________________________________________
TMu2e3189::~TMu2e3189() {
}

//-----------------------------------------------------------------------------
const char* TMu2e3189::GetFiguresDir() {
  return fFiguresDir.Data();
}

void TMu2e3189::get_filename(int Figure, char* Filename) {

  if      (Figure ==  -1) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_dt");
  else if (Figure ==  -2) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_ep");
  else if (Figure ==  -3) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_dy");
  else if (Figure ==  -4) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_dz");
  else if (Figure ==  -5) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_d0");
  else if (Figure ==  -6) sprintf(Filename,"fig_%03i_%s",Figure,"trk_1_z0");
  else if (Figure == 100) sprintf(Filename,"fig_%03i_%s",Figure,"EoverP");
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
void TMu2e3189::remake_plots() {

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
void TMu2e3189::plot(Int_t Figure, const char* CanvasName) {
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
//  Fig.1: not yet defined
//-----------------------------------------------------------------------------
  if (Figure == 1) {
  
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(dio_00_050_80c,"TCalm002", "cls_4/r");
    
    fP1 ->cd(1);

    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("Radius of Clusters in 30cm Calorimeter");
    fH1[0]->GetXaxis()->SetTitle("mm");
    fH1[0]->Draw();

    TLegend* leg = new TLegend(0.2,0.7,0.4,0.8,"");
    leg->AddEntry(fH1[0],"electrons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }

//-----------------------------------------------------------------------------
//  Fig.2: not yet defined
//-----------------------------------------------------------------------------
  if (Figure == 2) {
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(dio_00_050_80c,"TCalm002","cls_4/energy");
    
    fP1 ->cd(1);

    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("Cluster Energy 30cm with 10 MeV cut");
    fH1[0]->GetXaxis()->SetTitle("MeV");
    fH1[0]->Draw();

    TLegend* leg = new TLegend(0.2,0.7,0.4,0.8,"");
    leg->AddEntry(fH1[0],"electrons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }
 
//-----------------------------------------------------------------------------
//  Fig.3: not yet defined
//-----------------------------------------------------------------------------
  if (Figure == 3) {
  
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(dio_00_050_80c,"TCalm002", "cls_6/energy");
    
    fP1 ->cd(1);

    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("Cluster Energy 33cm with 10 MeV cut");
    fH1[0]->GetXaxis()->SetTitle("MeV");
    fH1[0]->Draw();

    TLegend* leg = new TLegend(0.2,0.7,0.4,0.8,"");
    leg->AddEntry(fH1[0],"electrons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }

//-----------------------------------------------------------------------------
//  Fig.4: not yet defined
//-----------------------------------------------------------------------------
  if (Figure == 4) {
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(dio_00_050_80c,"TCalm002","cls_10/energy");
    
    fP1 ->cd(1);

    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("Cluster E 34cm w/ cut");
    fH1[0]->GetXaxis()->SetTitle("MeV");
    fH1[0]->Draw();

    TLegend* leg = new TLegend(0.2,0.7,0.4,0.8,"");
    leg->AddEntry(fH1[0],"electrons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }
 
//-----------------------------------------------------------------------------
//  Fig.5: not yet defined
//-----------------------------------------------------------------------------
  if (Figure == 5) {
  
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(dio_00_050_80c,"TCalm002", "cls_12/energy");
    
    fP1 ->cd(1);

    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("Cluster E hist 35cm w/ cut");
    fH1[0]->GetXaxis()->SetTitle("MeV");
    fH1[0]->Draw();

    TLegend* leg = new TLegend(0.2,0.7,0.4,0.8,"");
    leg->AddEntry(fH1[0],"electrons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }

//-----------------------------------------------------------------------------
//  Fig.6: not yet defined
//-----------------------------------------------------------------------------
  if (Figure == 6) {
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fH1[0] = gh1(dio_00_050_80c,"TCalm002","cls_8/energy");
    
    fP1 ->cd(1);

    //    fH1[0]->SetOptStat(0);
    fH1[0]->SetTitle("Cluster E 36cm w/ cut");
    fH1[0]->GetXaxis()->SetTitle("MeV");
    fH1[0]->Draw();

    TLegend* leg = new TLegend(0.2,0.7,0.4,0.8,"");
    leg->AddEntry(fH1[0],"electrons","f");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  } 
//-----------------------------------------------------------------------------
//  Fig.7: N(clusters) vs radius 
//-----------------------------------------------------------------------------
  if (Figure == 7) {
    TLegend* leg ;
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fP1->cd(1);
    fH1[0] = gh1(dio_00_050_80c,"TCalm002","cls_4/r");
    //    double q1 = fH1[0]->Integral();
    
    double sf1 = (19931-1754.95)/160000.;
    fH1[0]->Scale(sf1);
    
    fH1[0]->GetXaxis()->SetRangeUser(300.,500.);
    fH1[0]->Draw("ep");
    //     fH1[1]->SetLineColor(kBlue-7);
    //     fH1[1]->SetMarkerColor(kBlue-7);
    //     fH1[1]->Draw("same,ep");

    leg = new TLegend(0.15,0.75,0.40,0.85,"");
    leg->AddEntry(fH1[0],"all tracks","h");
    leg->AddEntry(fH1[1],"Set C tracks","h");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
    
    fP1->cd(2);
    fH1[0] = gh1(dio_50_105_50b,"TCalm002","cls_4/r");
    //    double q2 = fH1[0]->Integral();
    
    double sf2 = (1754.95)/75000.;
    fH1[0]->Scale(sf2);
    
    fH1[0]->GetXaxis()->SetRangeUser(300.,500.);
    fH1[0]->SetLineColor(kRed-7);
    fH1[0]->SetMarkerColor(kRed-7);
    fH1[0]->Draw("same,ep");

    leg = new TLegend(0.15,0.75,0.40,0.85,"");
    leg->AddEntry(fH1[0],"all tracks","h");
    leg->AddEntry(fH1[1],"Set C tracks","h");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();
  }

//-----------------------------------------------------------------------------
//  Fig.8: Cluster Energy 
//-----------------------------------------------------------------------------
  if (Figure == 8) {
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fP1->cd(1);
    TH1F* h1 = gh1(dio_00_050_80c,"TCalm002","cls_4/energy");
    TH1F* h3 = (TH1F*) h1->Clone("h3");
    h3->Reset();
    double sf1 = (19931-1754.95)/160000.;
    h1->Scale(sf1);
    TH1F* h2 = gh1(dio_50_105_50b,"TCalm002","cls_4/energy");
    double sf2 = 1754.95/75000.;
    h2->Scale(sf2);
    h3->Add(h1,h2);
    h3->GetXaxis()->SetRangeUser(0.,100.);
    h3->GetYaxis()->SetRangeUser(0.,15.);
    h3->SetLineColor(kGreen-5);
    h3->Draw();

    h1->SetLineColor(kBlue-7);
    h1->SetFillStyle(3004);
    h1->SetFillColor(kBlue-7);
    h1->Draw("same");

    h2->SetLineColor(kRed-7);
    h2->SetFillStyle(3005);
    h2->SetFillColor(kRed-7);
    h2->Draw("same");
  }

//-----------------------------------------------------------------------------
//  Fig.9: Cluster Size 33cm
//-----------------------------------------------------------------------------
  if (Figure == 9) {
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fP1->cd(1);
    TH1F* h1 = gh1(dio_00_050_80c,"TCalm002","cls_6/ncr0");
    TH1F* h3 = (TH1F*) h1->Clone("h3");
    h3->Reset();
    double sf1 = (19931-1754.95)/160000.;
    TH1F* h2 = gh1(dio_50_105_50b,"TCalm002","cls_6/ncr0");
    double sf2 = (1754.95)/75000.;
    h3->Add(h1,h2,sf1,sf2);
    h3->GetXaxis()->SetRangeUser(0.,16.);
    h3->GetYaxis()->SetRangeUser(0.,12.);
    h3->SetLineColor(kGreen-2);
    h3->Draw();
   
    fP1->cd(2);
    TH1F* h4 = gh1(dio_00_050_80c,"TCalm002","cls_6/ncr1");
    TH1F* h6 = (TH1F*) h4->Clone("h6");
    h6->Reset();
    TH1F* h5 = gh1(dio_50_105_50b,"TCalm002","cls_6/ncr1");
    h6->Add(h4,h5,sf1,sf2);
    h6->GetXaxis()->SetRangeUser(0.,16.);
    h6->GetYaxis()->SetRangeUser(0.,12.);
    h6->SetLineColor(kCyan-3);
    h6->Draw("same");
  }

//-----------------------------------------------------------------------------
//  Fig.10: Cluster Size 30cm
//-----------------------------------------------------------------------------
  if (Figure == 10) {
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fP1->cd(1);
    TH1F* h1 = gh1(dio_00_050_80c,"TCalm002","cls_4/ncr0");
    TH1F* h3 = (TH1F*) h1->Clone("h3");
    h3->Reset();
    double sf1 = (19931-1754.95)/160000.;
    TH1F* h2 = gh1(dio_50_105_50b,"TCalm002","cls_4/ncr0");
    double sf2 = (1754.95)/75000.;
    h3->Add(h1,h2,sf1,sf2);
    h3->GetXaxis()->SetRangeUser(0.,16.);
    h3->GetYaxis()->SetRangeUser(0.,65.);
    h3->SetLineColor(kGreen-2);
    h3->Draw();
   
    fP1->cd(2);
    TH1F* h4 = gh1(dio_00_050_80c,"TCalm002","cls_4/ncr1");
    TH1F* h6 = (TH1F*) h4->Clone("h6");
    h6->Reset();
    TH1F* h5 = gh1(dio_50_105_50b,"TCalm002","cls_4/ncr1");
    h6->Add(h4,h5,sf1,sf2);
    h6->GetXaxis()->SetRangeUser(0.,16.);
    h6->GetYaxis()->SetRangeUser(0.,65.);
    h6->SetLineColor(kCyan-3);
    h6->Draw("same");
  }
//-----------------------------------------------------------------------------
//  Fig.11: Cluster time 
//-----------------------------------------------------------------------------
  if (Figure == 11) {
    
    fCanvas  = new_slide(name,title,1,1,1000,700);
    fP1      = (TPad*) fCanvas->GetPrimitive("p1");
    
    fP1->cd(1);
    TH1F* h1 = gh1(dio_00_050_80c,"TCalm002","cls_4/t0");
    h1->SetTitle("");
    h1->GetXaxis()->SetTitle("Cluster T0, ns");
    h1->Draw();
  }

}


/*

TH1F  *h1, *h2;

TH1F* h3 = (TH1F*) h1->Clone("h3");
h3->Reset();

h3->Add(h1,h2,w1,w2);

*/
