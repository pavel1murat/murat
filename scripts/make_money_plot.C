//

namespace {

  const char* e41s5721_track_ana  = "~/hist/mu2e/v5_7_0/e41s5721.track_ana.hist" ; // HEE+BGR , matcorr in CalPatRec
  const char* e41s5721_track_comp = "~/hist/mu2e/v5_7_0/e41s5721.track_comp.hist"; // HEE+BGR , matcorr in CalPatRec

  const char* e42s5721_track_ana  = "~/hist/mu2e/v5_7_0/e42s5721.track_ana.hist" ; // CE+BGR  , matcorr in CalPatRec
  const char* e42s5721_track_comp = "~/hist/mu2e/v5_7_0/e42s5721.track_comp.hist"; // CE+BGR  , matcorr in CalPatRec

  const char* figures_dir =  "~/work/mu2e/2016-04-12-calpatrec-update/figures";

  int N_gen_e41s5721 = 5110000;

  int N_gen_e42s5721 = 880000;
};

//-----------------------------------------------------------------------------
void draw_label_ndc(const char* Text, double X1, double Y1, double TextSize = 0.035, int Font = 52) {

  TLatex* label = new TLatex(X1,Y1,Text); label->SetNDC();
  label->SetTextSize(TextSize);
  label->SetTextFont(Font);
  label->Draw();
}

//-----------------------------------------------------------------------------
// use TrackComp histograms
//-----------------------------------------------------------------------------
void make_money_plot(const char* TrkFolder = "trk_3", int Print = 0) {

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

  TH1F* h_dio = (TH1F*) gh1(e41s5721_track_comp,"TrackComp",h_dio_name)->Clone("h_dio");

					// normalization
//-----------------------------------------------------------------------------
// now - normalization:
//-----------------------------------------------------------------------------
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

  double norm_dio = n_dio*(emax-emin)/N_gen_e41s5721;

  printf("norm_dio: %12.5e\n",norm_dio);

  h_dio->SetStats(0);
  h_dio->SetTitle("");
  h_dio->GetXaxis()->SetTitle("P, MeV/c");
  h_dio->GetXaxis()->SetRangeUser(101.0,107.0);
  h_dio->GetYaxis()->SetTitle(Form("N / %5.2f MeV/c", bin));
  h_dio->Scale(norm_dio);
  h_dio->GetYaxis()->SetRangeUser(0.,1.);
  h_dio->SetMarkerStyle(20);
  h_dio->SetMarkerSize (0.8);
  h_dio->SetMarkerColor(kBlue+3);
  h_dio->Draw();
  
//-----------------------------------------------------------------------------
// 2. get normalized conversion histogram for R = 1.e-16
//-----------------------------------------------------------------------------
  double R(1.e-16);
  double n_conv = n_stopped_muons*mu_capture_prob*R*live_fraction;

  TH1F* h_cnv = (TH1F*) gh1(e42s5721_track_comp,"TrackComp",h_cnv_name)->Clone("h_cnv");

  double norm_cnv = n_conv/N_gen_e42s5721;

  h_cnv->Scale(norm_cnv);

  h_cnv->SetStats(0);
  h_cnv->SetLineColor(kRed-3);
  h_cnv->SetLineWidth(2);
  h_cnv->SetFillColor(kRed-3);
  h_cnv->SetFillStyle(3002);
  h_cnv->Draw("sames");
//-----------------------------------------------------------------------------
// calculate signal and background
//-----------------------------------------------------------------------------
  int imin(477), imax(510);
  
  double emin = h_cnv->GetBinLowEdge(imin);
  double emax = h_cnv->GetBinLowEdge(imax)+h_cnv->GetBinWidth(imax);

  double N_sig = h_cnv->Integral(imin,imax);
  double N_dio = h_dio->Integral(imin,imax);
  double ses   = 1/N_sig*R;
//-----------------------------------------------------------------------------
// main title
//-----------------------------------------------------------------------------
  int font = 106;
  draw_label_ndc(Form("%s",title),0.65,0.85,20,font);
//-----------------------------------------------------------------------------
// print numbers
//-----------------------------------------------------------------------------

  draw_label_ndc(Form("N(prot on target) : %8.1e",n_prot         ),0.14,0.75,18,font);
  draw_label_ndc(Form("N(stopped muons)  : %8.1e",n_stopped_muons),0.14,0.70,18,font);
  draw_label_ndc(Form("R(mu --> e)       : %8.1e",R              ),0.14,0.65,18,font);

  draw_label_ndc(Form("Window            : [%6.2f,%6.2f] MeV",emin,emax),0.14,0.60,18,font);

  draw_label_ndc(Form("N(CE)             : %9.2e",N_sig            ),0.14,0.55,18,font);
  draw_label_ndc(Form("N(DIO)            : %9.2e",N_dio            ),0.14,0.50,18,font);
  draw_label_ndc(Form("SES               : %9.2e",ses              ),0.14,0.45,18,font);
  
  if (Print == 1) c->Print(Form("%s/money_plot_%s.eps",figures_dir,TrkFolder));

}
