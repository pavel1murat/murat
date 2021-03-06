//

#include "murat/scripts/datasets.hh"
#include "murat/scripts/plot_utilities.hh"

namespace local {

  const char* e41s5721_track_comp = "~/hist/mu2e/v5_7_0/e41s5721.track_comp.hist"; // HEE+BGR , matcorr in CalPatRec
  const char* e42s5721_track_comp = "~/hist/mu2e/v5_7_0/e42s5721.track_comp.hist"; // CE+BGR  , matcorr in CalPatRec

  const char* figures_dir =  "~/work/mu2e/2016-04-26-trkqual/figures";
};

//-----------------------------------------------------------------------------
// use TrackComp histograms
// e42s5721: CE + MIXCD3-cut-v2 
// e41s5721: Hee 95-105 MeV + MIXCD3-cut-v2 
//-----------------------------------------------------------------------------

class make_money_plot {
public:
  
  char fFiguresDir[200];

  const aa_dataset_t* fDatasetDio;
  const aa_dataset_t* fDatasetCnv;

public:
  make_money_plot() {
    fFiguresDir[0] = 0;
    fDatasetDio    = &e41s5721;
    fDatasetCnv    = &e42s5721;
  }
  
  void set_figures_dir(const char* Dir) { strcpy(fFiguresDir,Dir); }
  
  void plot           (const char* TrkFolder = "trk_3",
		       int         IMin      = 478,
		       int         IMax      = 510,
		       int         Print     =   0);

  void plot           (const aa_dataset_t*  DatasetCnv,
		       const aa_dataset_t*  DatasetDio,
		       const char*       TrkFolder = "trk_3",
		       int               IMin      = 478,
		       int               IMax      = 510,
		       int               Print     =   0);
};


//-----------------------------------------------------------------------------
void make_money_plot::plot(const char* TrkFolder, int IMin, int IMax, int Print) {

  // dataset_t* dio_dataset = &e41s5721;
  // dataset_t* cnv_dataset = &e42s5721;

  // 1. get DIO histogram and normalize it

  char h_dio_name[100], h_cnv_name[100], title[200];

  // sprintf(h_dio_name,"%s/%s",TrkFolder,"pdio");
  // sprintf(h_cnv_name,"%s/%s",TrkFolder,"p");

  sprintf(h_dio_name,"%s/%s",TrkFolder,"pdiolw");
  sprintf(h_cnv_name,"%s/%s",TrkFolder,"plw");

  if      (strcmp(TrkFolder,"trk_0"  ) == 0) strcpy(title,"TrkPatRec, TrkQual > 0.4");
  else if (strcmp(TrkFolder,"trk_101") == 0) strcpy(title,"TrkPatRec, TrkQual > 0.4");
  else if (strcmp(TrkFolder,"trk_201") == 0) strcpy(title,"CalPatRec, TrkQual > 0.4");
  else if (strcmp(TrkFolder,"trk_3"  ) == 0) strcpy(title,"Optimized track selection");

  TCanvas* c = new TCanvas("c_money_plot","Money Plot",1500,800);

  TH1F* h_dio = (TH1F*) gh1(fDatasetDio->fn_track_comp,"TrackComp",h_dio_name)->Clone("h_dio");

					// upper limit of the plot
  float ymax = 0.6;
					// normalization
//-----------------------------------------------------------------------------
// now - normalization:
//-----------------------------------------------------------------------------
					// assumed conversion fraction
  double R(1.e-16);
					// capture probability in Al
  double mu_capture_prob(0.609);
					// roughly speaking, fraction of captures within the 
					// active window [700,1700] ns  (0.53 = 2594./4889.)
					// is accounted for by already applied cut T0>700ns
  double live_fraction(1.); // (0.53);
					// total number of 8 GeV protons on target in 3 years
  double n_prot      (1.2e20 * 3.);
					// muon (transport+stopping) efficiency : 1.6e-3

  double muon_transport_eff = 1.86e-3;

  double n_stopped_muons = n_prot*muon_transport_eff;

  printf("n_stopped_muons: %12.5e\n",n_stopped_muons);
  
					// fraction DIO events in [95,105] MeV/c
  double f_95_105(3.9821e-11);
					// track reco eff for DIO's is accounted for later

  double n_dio = n_stopped_muons*(1-mu_capture_prob)*live_fraction;

  double bin = h_dio->GetBinWidth(1) ; // MeV

  double emax(105.), emin(95.);  // in MeV

  double norm_dio = n_dio*(emax-emin)/fDatasetDio->ngen;

  printf("norm_dio: %12.5e\n",norm_dio);
//-----------------------------------------------------------------------------
// calculate signal and background
//-----------------------------------------------------------------------------
//  int imin(478), imax(510);
  
  TH1F* h_cnv = (TH1F*) gh1(fDatasetCnv->fn_track_comp,"TrackComp",h_cnv_name)->Clone("h_cnv");

  double emin1 = h_cnv->GetBinLowEdge(IMin);
  double emax1 = h_cnv->GetBinLowEdge(IMax)+h_cnv->GetBinWidth(IMax);

  h_dio->SetStats(0);
  h_dio->SetTitle("");
  h_dio->GetXaxis()->SetTitle("P, MeV/c");
  h_dio->GetXaxis()->SetRangeUser(102.0,107.0);
  h_dio->GetYaxis()->SetTitle(Form("N / %5.2f MeV/c", bin));
  h_dio->Scale(norm_dio);
  h_dio->GetYaxis()->SetRangeUser(0.,ymax);
  h_dio->SetMarkerStyle(20);
  h_dio->SetMarkerSize (0.8);
  h_dio->SetMarkerColor(kBlue+3);
  h_dio->Draw();
//-----------------------------------------------------------------------------
// signal window
//-----------------------------------------------------------------------------
  TBox* b = new TBox(emin1,0.,emax1,ymax);
  b->Draw();
					// draw on top...
  h_dio->Draw("same");

  printf("done drawing\n");
//-----------------------------------------------------------------------------
// 2. get normalized conversion histogram for R = 1.e-16
//-----------------------------------------------------------------------------
  double n_conv = n_stopped_muons*mu_capture_prob*R*live_fraction;

  double norm_cnv = n_conv/fDatasetCnv->ngen;

  h_cnv->Scale(norm_cnv);

  h_cnv->SetStats(0);
  h_cnv->SetLineColor(kRed-3);
  h_cnv->SetLineWidth(2);
  h_cnv->SetFillColor(kRed-3);
  h_cnv->SetFillStyle(3002);
  h_cnv->Draw("sames");
//-----------------------------------------------------------------------------
// calculate signal and background within the window
//-----------------------------------------------------------------------------  
  double N_sig = h_cnv->Integral(IMin,IMax);
  double N_dio = h_dio->Integral(IMin,IMax);
  double ses   = 1/N_sig*R;
//-----------------------------------------------------------------------------
// main title and numbers
//-----------------------------------------------------------------------------
  int font = 106;
  double x1 = 0.61;
  draw_label_ndc(Form("%s"                  ,title            ),x1,0.85,20,font);
  draw_label_ndc(Form("N(POT)       : %8.1e",n_prot           ),x1,0.75,20,font);
  draw_label_ndc(Form("N(stopped #mu) : %8.1e",n_stopped_muons),x1,0.70,20,font);
  draw_label_ndc(Form("R_{#mu #rightarrow e}        : %8.1e",R),x1,0.65,20,font);

  draw_label_ndc(Form("Window       : [%6.2f,%6.2f] MeV/c",emin1,emax1),x1,0.60,20,font);

  draw_label_ndc(Form("N(CE)        : %9.2e",N_sig            ),x1,0.55,20,font);
  draw_label_ndc(Form("N(DIO)       : %9.2e",N_dio            ),x1,0.50,20,font);
  draw_label_ndc(Form("SES          : %9.2e",ses              ),x1,0.45,20,font);
  
  if (Print == 1) c->Print(Form("%s/money_plot_%s.eps",local::figures_dir,TrkFolder));

}


//-----------------------------------------------------------------------------
void make_money_plot::plot(const aa_dataset_t*  DatasetCnv,
			   const aa_dataset_t*  DatasetDio,
			   const char*       TrkFolder ,
			   int               IMin      ,
			   int               IMax      ,
			   int               Print     ) {

  fDatasetCnv = DatasetCnv;
  fDatasetDio = DatasetDio;

  plot(TrkFolder,IMin,IMax,Print);
}



