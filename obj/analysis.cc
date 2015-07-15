///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TLegend.h"
#include "THStack.h"
#include "TROOT.h"
#include "TColor.h"
#include "TLatex.h"
#include "TMath.h"

#include "obj/analysis.hh"

#include <iostream>

ClassImp(analysis)


//-----------------------------------------------------------------------------
analysis::analysis(const char* Name, const char* Title, const char* Version, int LumiBin): 
  TNamed(Name,Title) {
  // 

  fVersion    = Version;
  fNRunRanges = 0;
  fRR         = 0;
  fDsMetadata = 0;
  fLumiBin    = LumiBin;
  fDatHist    = 0;
  fBgrHist    = 0;
  fBgrStack   = 0;

  fLegendXMin = 0.60;
  fLegendYMin = 0.65;
  fLegendXMax = 0.85;
  fLegendYMax = 0.81;

  for (int i=0; i<kDebugFlags; i++) {
    fDebugFlag[i] = 0;
  }
}


//-----------------------------------------------------------------------------
int analysis::init() {
  
  fNRunRanges = fDsMetadata->NRunRanges();
  fRR         = new run_range_dat_t[fNRunRanges];

  run_range_dat_t* rr;

  for (int i=0; i<fNRunRanges; i++) {
    rr = fRR+i;
    rr->fIntLumi = -1;
    rr->fDat     = 0;
    rr->fBgr     = 0;
    rr->fListOfHistograms    = new TObjArray();
    rr->fListOfDmbHistograms = new TObjArray();

    rr->fQSig    = new double[aprocess::kNBins];
    rr->fQBgr    = new double[aprocess::kNBins];
    rr->fESig    = new double[aprocess::kNBins];
    rr->fEBgr    = new double[aprocess::kNBins];
    rr->fChi2    = new double[aprocess::kNBins];
    rr->fQMcSig  = new double[aprocess::kNBins];
  }

  return 0;
}

//-----------------------------------------------------------------------------
analysis::~analysis() {

  run_range_dat_t* rr;

  if (fNRunRanges > 0) {

    for (int i=0; i<fNRunRanges; i++) {
      rr = fRR+i;

      if (rr->fDat) delete rr->fDat;

      if (rr->fBgr) {
	rr->fBgr->Delete();
	delete rr->fBgr;
      }

      delete rr->fListOfHistograms;
      delete rr->fListOfDmbHistograms;
      
      delete [] rr->fQSig;
      delete [] rr->fQBgr;
      delete [] rr->fESig;
      delete [] rr->fEBgr;
      delete [] rr->fChi2;
      delete [] rr->fQMcSig;
    }

    delete fRR;
  }

  if (fDsMetadata) delete fDsMetadata;
}

//-----------------------------------------------------------------------------
// every time creates new histograms, so it is supposed to leak memory
// McFlag == 0: plot everything
// McFlag == 1: plot data and MC backgrounds only
// McFlag == 3: data + MC backgrounds + MC signal
//-----------------------------------------------------------------------------
void analysis::plot(const char* Module  , 
		    const char* HistSet , 
		    int         RunRange,
		    int         Bin     , 
		    const char* HistName, 
		    int         Rebin   ,
		    int         McFlag  ,
		    double      XMin    ,
		    double      XMax    ,
		    const char* Axislabel,
		    int         Lumi    ,
		    int         YMax    ,
                    const char* YUnits  ,
		    float       AxislabelSize,
		    int         Ndivisions,
                    int         CombineProcesses)
{
  // note, that as overlaying stacked 2D histograms doesn't make much sense, 
  // assume that we're not doing it...
  int         rc(0), nbgr_label;
  double      q_bgr(0), n_bgr(0);
  TH1         *h_dat, /* *h_sig,*/ *h_bgr, *h_label[100]; // , *h_bgr_accumulate;
  aprocess    *bgr; //, *sig;
  char        hist_name[200]; //, title[200];
  const char  *name_label[100];
  double      _total_bgr_integral = 0;

  //  int add_WW_WZ_ZZ = 0;
  //  int add_Wjets_Wgamma = 0;
  // if(CombineProcesses==1){
  //   add_WW_WZ_ZZ = 1;
  //   add_Wjets_Wgamma = 1;
  // }
  // if(CombineProcesses==2){
  //   add_WW_WZ_ZZ = 1;
  // }

  int verbose=1;
//-----------------------------------------------------------------------------
// stack backgrounds and signal... mmm it doesnt' work exactly this way
//-----------------------------------------------------------------------------
  run_range_dat_t* rr = fRR+RunRange;

  THStack* stk = new THStack("stack","test"); // here we intentionally leak memory...
//-----------------------------------------------------------------------------
// get data histogram - need it at this point to choose the normalization
//-----------------------------------------------------------------------------
  rc = GetHistogram(rr->fDat,Module,HistSet,RunRange,Bin,HistName,Rebin,h_dat);
  fBgrHist = (TH1*) h_dat->Clone("bgr_hist");
  fBgrHist->Reset();
//-----------------------------------------------------------------------------
// create legend to be displayed and add data histogram as the first entry
//-----------------------------------------------------------------------------
  TLegend* leg = new TLegend(fLegendXMin,fLegendYMin,fLegendXMax,fLegendYMax);

  leg->SetFillColor(0);
  leg->SetLineWidth(0);

  leg->AddEntry(h_dat,rr->fDat->GetLabel(),"pl");

  nbgr_label = 0;

  for (int i=0; i<rr->fNBgr; i++) {
    bgr = (aprocess*) rr->fBgr->At(i);
    if (fDebugFlag[0] != 0) {
      if(verbose) printf("[analysis::plot]: BGR=%s %s\n",bgr->GetName(),bgr->GetTitle());
    }
//-----------------------------------------------------------------------------
// handle QCD histograms (fXSec = -1) separately in get_h1 - need to define 
// number of events 
//-----------------------------------------------------------------------------
    if (McFlag != 0) {
      if (bgr->GetMcFlag() > McFlag)                        goto NEXT_PROCESS;
      rc = GetHistogram(bgr,Module,HistSet,RunRange,Bin,HistName,Rebin,h_bgr);
      if (i < rr->fNBgr-1) {
	h_bgr->SetFillColor(bgr->GetColor());
	h_bgr->SetFillStyle(bgr->GetStyle());
      }
    }
    else {
//-----------------------------------------------------------------------------
// do not check fXSec and draw QCD
// aprocess:get_h1 returns possibly rebinned clone
//-----------------------------------------------------------------------------
      rc = GetHistogram(bgr,Module,HistSet,RunRange,Bin,HistName,Rebin,h_bgr);
      if (h_bgr == 0) {
	GetHistogramName(HistSet,RunRange,Bin,HistName,hist_name);
	printf("analysis::plot ERROR: hist %s = 0 for %s, SKIP it\n",
	       hist_name,bgr->GetName()); 
	                                                    goto NEXT_PROCESS;
      }

      if (i < rr->fNBgr-1) {
	h_bgr->SetFillColor(bgr->GetColor());
	h_bgr->SetFillStyle(bgr->GetStyle());
      }
    }
    
   //  //Kludge to add up the Alpgen
   //  if( strstr(bgr->GetName(),"alpgen_zee")!=0){
   //    //Assumes they are together in the list with 1p first and 4p last...
   //    if( strcmp(bgr->GetName(),"alpgen_zee_1p")==0){
   // 	h_bgr_accumulate=(TH1F*)h_bgr->Clone(h_bgr->GetTitle());
   //    }
   //    else{
   // 	h_bgr_accumulate->Add(h_bgr);
   //    }
      
   //    if( strcmp(bgr->GetName(),"alpgen_zee_4p")!=0) goto NEXT_PROCESS;
   //    else
   // 	h_bgr = h_bgr_accumulate;
   //  }
   //  if( strstr(bgr->GetName(),"alpgen_zmm")!=0){
   //    //Assumes they are together in the list with 1p first and 4p last...
   //    if( strcmp(bgr->GetName(),"alpgen_zmm_1p")==0)
   // 	h_bgr_accumulate=(TH1F*)h_bgr->Clone(h_bgr->GetTitle());
   //    else
   // 	h_bgr_accumulate->Add(h_bgr);
   //    if( strcmp(bgr->GetName(),"alpgen_zmm_4p")!=0) goto NEXT_PROCESS;
   //    else
   // 	h_bgr = h_bgr_accumulate;
   //  }

   //  //Kludge to add up the WW, WZ, ZZ if flag set
   //  if(add_WW_WZ_ZZ){
   //    if( (strcmp(bgr->GetName(),"ww")==0)||(strcmp(bgr->GetName(),"wz")==0)||(strcmp(bgr->GetName(),"zz")==0)){
   // 	//Assumes they are together in the list in the order ww, wz, zz
   // 	if( strcmp(bgr->GetName(),"ww")==0){
   // 	  h_bgr_accumulate=(TH1F*)h_bgr->Clone(h_bgr->GetTitle());
   // 	}
   // 	else{
   // 	  h_bgr_accumulate->Add(h_bgr);
   // 	}
   // 	if( strcmp(bgr->GetName(),"zz")!=0) goto NEXT_PROCESS;
   // 	else
   // 	  h_bgr = h_bgr_accumulate;
   //    }
   //  }
   // if(add_Wjets_Wgamma){
   //   if( (strcmp(bgr->GetName(),"wenu")==0)||(strcmp(bgr->GetName(),"wenugamma")==0)){
   // 	//Assumes they are together in the list in the order wenugamma, wenu
   // 	if( strcmp(bgr->GetName(),"wenugamma")==0){
   // 	  h_bgr_accumulate=(TH1F*)h_bgr->Clone(h_bgr->GetTitle());
   // 	}
   // 	else{
   // 	  h_bgr_accumulate->Add(h_bgr);
   // 	}
   // 	if( strcmp(bgr->GetName(),"wenu")!=0) goto NEXT_PROCESS;
   // 	else
   // 	  h_bgr = h_bgr_accumulate;
   //    }
    //    }

    if (XMin != -999.) {
      h_bgr->GetXaxis()->SetRangeUser(XMin,XMax);
    }
//-----------------------------------------------------------------------------
// show only histograms with non-zero number of entries
//-----------------------------------------------------------------------------
    q_bgr = h_bgr->Integral();
    n_bgr = h_bgr->GetEntries();

    //---h_bgr->Scale(q_bgr/2.);// for 106, Aidan 110221

    _total_bgr_integral+= q_bgr;

    if (verbose) {
      printf("[analysis::plot] bgr = %s, nev = %10.1f, integral = %10.3f, intlumi=%10.3f\n ",
	     bgr->GetName(),n_bgr,q_bgr,bgr->fIntLumi);
    }

    if (q_bgr != 0) {
      stk->Add(h_bgr,"h");
					// labels in the legend should go in 
					// the same orger as backgrounds in stack
      h_label   [nbgr_label] = h_bgr;
      name_label[nbgr_label] = bgr->GetLabel();
      // if(add_WW_WZ_ZZ && (strcmp(bgr->GetLabel(),"ZZ")==0)) name_label[nbgr_label] = "WW,WZ,ZZ";
      // if(add_Wjets_Wgamma && (strcmp(bgr->GetLabel(),"W+jets")==0)) name_label[nbgr_label] = "W+jets,W#gamma";
      nbgr_label += 1;

      fBgrHist->Add(h_bgr);
      //      h_bgr_accumulate = 0;
    }
  NEXT_PROCESS:;
  }
//-----------------------------------------------------------------------------
// add background descriptions to the legend
//-----------------------------------------------------------------------------
  for (int i=nbgr_label-1; i>=0; i--) {
    leg->AddEntry(h_label[i],name_label[i],"f");
  }
//-----------------------------------------------------------------------------
// draw stacked backgrounds
//-----------------------------------------------------------------------------
  double qmx_stk = stk->GetMaximum();
  double qmx_dat = h_dat->GetMaximum();

  if (qmx_stk < qmx_dat) stk->SetMaximum(qmx_dat);

  stk->SetTitle(h_dat->GetTitle());

  //  const char* tit = h_dat->GetXaxis()->GetTitle();
  //  TAxis* axis = stk->GetXaxis();
  //  axis->SetTitle(tit);  // 'axis' is zero....

  if (XMin != -999.) {
    h_dat->GetXaxis()->SetRangeUser(XMin,XMax);
  }

  //---h_dat->Scale(h_dat->Integral()/2.); // for 106, Aidan 110221 - put it here so scale is funny as reminder

  if (YMax>0) h_dat->SetMaximum(YMax);

  h_dat->SetMarkerStyle(20);
  h_dat->SetMarkerSize (1);

  double bin=h_dat->GetXaxis()->GetBinWidth(1);

  //h_dat->GetYaxis()->SetTitle(Form("Events / %.1f %s",bin,YUnits));
  h_dat->GetYaxis()->SetTitle(Form("Events / %d %s",int(bin),YUnits));
  h_dat->GetYaxis()->SetTitleSize(.05);
  h_dat->GetYaxis()->SetTitleOffset(1.);
  h_dat->GetXaxis()->SetTitle(Axislabel);
  h_dat->GetXaxis()->SetTitleSize(.06);
  h_dat->GetXaxis()->SetTitleOffset(.77);

  h_dat->SetStats(0); //For preblessing 110331
//-----------------------------------------------------------------------------
// finally draw it, with the data distribution on top
//-----------------------------------------------------------------------------
  stk->Draw();
  if(AxislabelSize>0.){
    stk->GetXaxis()->SetLabelSize(AxislabelSize);//must be drawn before accessing axis...
    stk->GetYaxis()->SetLabelSize(AxislabelSize);//must be drawn before accessing axis...
    stk->GetXaxis()->SetLabelFont(42);
    stk->GetXaxis()->SetTitleFont(42);
    stk->GetYaxis()->SetLabelFont(42);
    stk->GetYaxis()->SetTitleFont(42);
    h_dat->GetXaxis()->SetLabelSize(AxislabelSize);//must be drawn before accessing axis...
    h_dat->GetYaxis()->SetLabelSize(AxislabelSize);//must be drawn before accessing axis...
    h_dat->GetXaxis()->SetLabelFont(42);
    h_dat->GetXaxis()->SetTitleFont(42);
    h_dat->GetYaxis()->SetLabelFont(42);
    h_dat->GetYaxis()->SetTitleFont(42);
  }
  if(Ndivisions>0.){
    stk->GetXaxis()->SetNdivisions(Ndivisions);
    h_dat->GetXaxis()->SetNdivisions(Ndivisions);
  }
  h_dat->Draw("ep");
  stk->Draw("same");
  h_dat->Draw("same");

  leg->Draw();

//   if (Lumi==6) {
//     TLatex* tex = new TLatex;
//     tex->SetNDC(kTRUE);
//     tex->SetTextSize(0.03);
//     tex->DrawLatex(0.42,0.845,"CDF Run II Preliminary #int L dt = 6 fb^{-1}");
//   }

  fDatHist  = h_dat;
  fLegend   = leg;

  //print out data and background integrals
  double _total_data = h_dat->Integral();
  double _ratio = _total_bgr_integral>0.?_total_data/_total_bgr_integral:0.;
  if(verbose) printf(" Total bgr integral = %10.3f; data = %10.3f ; ratio = %10.3f\n",_total_bgr_integral,_total_data,_ratio);

  if (rc < 0) printf(" RC = %i\n",rc);
  
}

//-----------------------------------------------------------------------------
// every time creates new histograms, so it is supposed to leak memory
// McFlag == 0: plot everything
// McFlag == 1: plot data and MC backgrounds only
// McFlag == 2: MC backgrounds only
// McFlag == 3: MC backgrounds + MC signal
//-----------------------------------------------------------------------------
void analysis::plot_no_data(const char* Module  , 
			    const char* HistSet , 
			    int         RunRange,
			    int         Bin     , 
			    const char* HistName, 
			    int         Rebin   ,
			    int         McFlag  ,
			    double      XMin    ,
			    double      XMax    ,
			    const char* Axislabel,
			    int         Lumi    )
{
  // note, that as overlaying stacked 2D histograms doesn't make much sense, 
  // assume that we're not doing it...
  int         rc(0), nbgr_label, first;
  double      q_bgr(0);
  TH1         *h_dat(0), /**h_sig(0),*/ *h_bgr(0), *h_label[100], *h_bgr_accumulate(0);
  aprocess    *bgr; //, *sig;
  char        hist_name[200]; //, title[200];
  const char  *name_label[100];
  const char  /* *tit,*/ *hist_tit(NULL); 

//-----------------------------------------------------------------------------
// stack backgrounds and signal... mmm it doesnt' work exactly this way
//-----------------------------------------------------------------------------
  run_range_dat_t* rr = fRR+RunRange;

  THStack* stk = new THStack("stack","test");
//-----------------------------------------------------------------------------
// create legend to be displayed and add data histogram as the first entry
//-----------------------------------------------------------------------------
  //TLegend* leg = new TLegend(0.50,0.50,0.75,0.82);
  TLegend* leg = new TLegend(0.60,0.50,0.85,0.81);
  leg->SetFillColor(0);
  leg->SetLineWidth(0);
  //  leg->AddEntry(h_dat,rr->fDat->GetTitle(),"pl");

  nbgr_label = 0;
  first      = 1;

  for (int i=0; i<rr->fNBgr; i++) {
    bgr = (aprocess*) rr->fBgr->At(i);
    if (fDebugFlag[0] != 0) {
      printf("[analysis::plot]: BGR=%s %s\n",bgr->GetName(),bgr->GetTitle());
    }
//-----------------------------------------------------------------------------
// handle QCD histograms (fXSec = -1) separately in get_h1 - need to define 
// number of events 
//-----------------------------------------------------------------------------
    if (McFlag) {
      if (bgr->GetMcFlag() > McFlag)                        goto NEXT_PROCESS;
      //      if (bgr->fXSec > 0) {
	rc = GetHistogram(bgr,Module,HistSet,RunRange,Bin,HistName,Rebin,h_bgr);
	if (i < rr->fNBgr-1) {
	  h_bgr->SetFillColor(bgr->GetColor());
	  h_bgr->SetFillStyle(bgr->GetStyle());
	  //      sprintf(title,"%s",bgr->GetName());
	  //    leg->AddEntry(h_bgr,title); 
	}
	// }
    }
    else {
//-----------------------------------------------------------------------------
// do not check fXSec and draw QCD
// aprocess:get_h1 returns possibly rebinned clone
//-----------------------------------------------------------------------------
      rc = GetHistogram(bgr,Module,HistSet,RunRange,Bin,HistName,Rebin,h_bgr);
      if (h_bgr == 0) {
	GetHistogramName(HistSet,RunRange,Bin,HistName,hist_name);
	printf("analysis::plot ERROR: hist %s = 0 for %s, SKIP it\n",
	       hist_name,bgr->GetName()); 
	                                                    goto NEXT_PROCESS;
      }
      if (i < rr->fNBgr-1) {
	h_bgr->SetFillColor(bgr->GetColor());
	h_bgr->SetFillStyle(bgr->GetStyle());
	  //      sprintf(title,"%s",bgr->GetName());
	  //    leg->AddEntry(h_bgr,title); 
      }
    }

    //Kludge to add up the Alpgen
    if( strstr(bgr->GetName(),"alpgen_zee")!=0){
      //Assumes they are together in the list with 1p first and 4p last...
      if( strcmp(bgr->GetName(),"alpgen_zee_1p")==0){
	h_bgr_accumulate=(TH1F*)h_bgr->Clone(h_bgr->GetTitle());
      }
      else{
	h_bgr_accumulate->Add(h_bgr);
      }
      
      if( strcmp(bgr->GetName(),"alpgen_zee_4p")!=0) goto NEXT_PROCESS;
      else
	h_bgr = h_bgr_accumulate;
    }
    if( strstr(bgr->GetName(),"alpgen_zmm")!=0){
      //Assumes they are together in the list with 1p first and 4p last...
      if( strcmp(bgr->GetName(),"alpgen_zmm_1p")==0)
	h_bgr_accumulate=(TH1F*)h_bgr->Clone(h_bgr->GetTitle());
      else
	h_bgr_accumulate->Add(h_bgr);
      if( strcmp(bgr->GetName(),"alpgen_zmm_4p")!=0) goto NEXT_PROCESS;
      else
	h_bgr = h_bgr_accumulate;
    }

    if (first == 1) {
      first = 0;
      h_dat = h_bgr;
      hist_tit = h_bgr->GetTitle();

      fBgrHist = (TH1*) h_bgr->Clone("bgr_hist");
      fBgrHist->Reset();
    }

    if (XMin != -999.) {
      h_bgr->GetXaxis()->SetRangeUser(XMin,XMax);
    }
//-----------------------------------------------------------------------------
// show only histograms with non-zero number of entries
//-----------------------------------------------------------------------------
    q_bgr = h_bgr->Integral();

    printf(" bgr = %s , integral = %10.3f\n ",bgr->GetName(),q_bgr);
    if (q_bgr != 0) {
      stk->Add(h_bgr,"h");
					// labels in the legend should go in 
					// the same orger as backgrounds in stack
      h_label   [nbgr_label] = h_bgr;
      name_label[nbgr_label] = bgr->GetLabel();
      nbgr_label += 1;

      fBgrHist->Add(h_bgr);
      h_bgr_accumulate = 0;
    }
  NEXT_PROCESS:;
  }
//-----------------------------------------------------------------------------
// add background descriptions to the legend
//-----------------------------------------------------------------------------
  for (int i=nbgr_label-1; i>=0; i--) {

//     char labeltext[128];
//     char labellabel[128];
//     sprintf(labellabel,name_label[i]);
//     if(      strcmp(labellabel,"wz"      )==0) sprintf(labeltext,"WZ");
//     else if(      strcmp(labellabel,"bhelzl"      )==0) sprintf(labeltext,"data");
//     else if(      strcmp(labellabel,"gh0325"      )==0) sprintf(labeltext,"G*, m=325GeV");
//     else if(      strcmp(labellabel,"gh0600"      )==0) sprintf(labeltext,"G*, m=600GeV");
//     else if(      strcmp(labellabel,"gpt325"      )==0) sprintf(labeltext,"boosted G* 325GeV");
//     else if( strcmp(labellabel,"alpgen_wen_4p")==0) sprintf(labeltext,"W#rightarrow e#nu + 4p");
//     else if( strcmp(labellabel,"alpgen_wen_3p")==0) sprintf(labeltext,"W#rightarrow e#nu + 3p");
//     else if( strcmp(labellabel,"alpgen_wen_2p")==0) sprintf(labeltext,"W#rightarrow e#nu + 2p");
//     else if( strcmp(labellabel,"alpgen_wen_1p")==0) sprintf(labeltext,"W#rightarrow e#nu + 1p");
//     else if( strcmp(labellabel,"alpgen_wmn_4p")==0) sprintf(labeltext,"W#rightarrow #mu #nu + 4p");
//     else if( strcmp(labellabel,"alpgen_wmn_3p")==0) sprintf(labeltext,"W#rightarrow #mu #nu + 3p");
//     else if( strcmp(labellabel,"alpgen_wmn_2p")==0) sprintf(labeltext,"W#rightarrow #mu #nu + 2p");
//     else if( strcmp(labellabel,"alpgen_wmn_1p")==0) sprintf(labeltext,"W#rightarrow #mu #nu + 1p");
//     else if( strcmp(labellabel,"alpgen_zee_4p")==0) sprintf(labeltext,"Z#rightarrow ee + jets");
    //else if( strcmp(labellabel,"alpgen_zee_4p")==0) sprintf(labeltext,"Z#rightarrow ee + 4p");
//     else if( strcmp(labellabel,"alpgen_zee_3p")==0) sprintf(labeltext,"Z#rightarrow ee + 3p");
//     else if( strcmp(labellabel,"alpgen_zee_2p")==0) sprintf(labeltext,"Z#rightarrow ee + 2p");
//     else if( strcmp(labellabel,"alpgen_zee_1p")==0) sprintf(labeltext,"Z#rightarrow ee + 1p");
//     else if( strcmp(labellabel,"alpgen_zmm_4p")==0) sprintf(labeltext,"Z#rightarrow #mu #mu + jets");
    //else if( strcmp(labellabel,"alpgen_zmm_4p")==0) sprintf(labeltext,"Z#rightarrow #mu #mu + 4p");
//     else if( strcmp(labellabel,"alpgen_zmm_3p")==0) sprintf(labeltext,"Z#rightarrow #mu #mu + 3p");
//     else if( strcmp(labellabel,"alpgen_zmm_2p")==0) sprintf(labeltext,"Z#rightarrow #mu #mu + 2p");
//     else if( strcmp(labellabel,"alpgen_zmm_1p")==0) sprintf(labeltext,"Z#rightarrow #mu #mu + 1p");
//     else if( strcmp(labellabel,"zee")==0) sprintf(labeltext,"Z#rightarrow ee");
//     else if( strcmp(labellabel,"zmm")==0) sprintf(labeltext,"Z#rightarrow #mu#mu");
//     else if( strcmp(labellabel,"wmunu")==0) sprintf(labeltext,"W#rightarrow #mu#nu");
//     else if( strcmp(labellabel,"zee_jets")==0) sprintf(labeltext,"Z#rightarrow ee + jets");
//     else if( strcmp(labellabel,"zmm_jets")==0) sprintf(labeltext,"Z#rightarrow #mu #mu + jets");
//     else if( strcmp(labellabel,"zeegamma")==0) sprintf(labeltext,"Z#rightarrow ee + #gamma");
//     else if( strcmp(labellabel,"zmmgamma")==0) sprintf(labeltext,"Z#rightarrow #mu #mu + #gamma");
//     else if( strcmp(labellabel,"wenugamma")==0) sprintf(labeltext,"W#gamma#rightarrow e#nu#gamma");
    //    else if( strcmp(labellabel,"fakes"    )==0) sprintf(labeltext,"W+jets");
//     else if( strcmp(labellabel,"zz"      )==0) sprintf(labeltext,"ZZ");
//     else if( strcmp(labellabel,"ww"      )==0) sprintf(labeltext,"WW");
//     else if( strcmp(labellabel,"ttbar"   )==0) sprintf(labeltext,"t#bar{t}");

//     else sprintf(labeltext,labellabel);
    
    //    leg->AddEntry(h_label[i],labeltext,"f");
    leg->AddEntry(h_label[i],name_label[i],"f");
  }
//-----------------------------------------------------------------------------
// draw stacked backgrounds
//-----------------------------------------------------------------------------
  double qmx_stk = stk->GetMaximum();
  double qmx_dat = h_dat->GetMaximum();

  //  stk->SetMinimum(0);
  if (qmx_stk > qmx_dat) h_dat->SetMaximum(qmx_stk*1.1);

  stk->SetTitle(hist_tit);

  //  tit = h_dat->GetXaxis()->GetTitle();
  TAxis* xa = h_dat->GetXaxis();
  //  axis->SetTitle(tit);  // 'axis' is zero....

  if (XMin != -999.) {
    //    h_dat->GetXaxis()->SetRangeUser(XMin,XMax);
    xa->SetRangeUser(XMin,XMax);
  }
//-----------------------------------------------------------------------------
// finally draw data distribution on top
//-----------------------------------------------------------------------------
//   h_dat->SetMarkerStyle(20);
//   h_dat->SetMarkerSize (1);

  double bin=h_dat->GetXaxis()->GetBinWidth(1);

  //h_dat->GetYaxis()->SetTitle(Form("Events/%f",bin));

  h_dat->GetYaxis()->SetTitle(Form("Events/%.1f",bin));
  h_dat->GetXaxis()->SetTitleSize(.05);
  h_dat->GetYaxis()->SetTitleOffset(1.);
  h_dat->GetXaxis()->SetTitle(Axislabel);
  h_dat->GetXaxis()->SetTitleSize(.06);
  h_dat->GetXaxis()->SetTitleOffset(.77);
  
  h_dat->SetStats(0); //For preblessing 110331

  h_dat->Draw("");

  stk->Draw("ep,same");

  leg->Draw();

  if(Lumi==6){
    TLatex* tex = new TLatex;
    tex->SetNDC(kTRUE);
    tex->SetTextSize(0.03);
    tex->DrawLatex(0.42,0.845,"CDF Run II Preliminary #int L dt = 6 fb^{-1}");
  }

  fDatHist  = h_dat;

  if (rc < 0) printf(" RC = %i\n",rc);
}

//-----------------------------------------------------------------------------
// do not plot QCD background
//-----------------------------------------------------------------------------
void analysis::plot_mc_only(const char* Module  , 
			    const char* HistSet , 
			    int         RunRange,
			    int         Bin     , 
			    const char* HistName, 
			    int         Rebin   )
{
  plot(Module,HistSet,RunRange,Bin,HistName,Rebin,1);
}


///-----------------------------------------------------------------------------
// data minus sum of all the backgrounds
// skip backgrounds defined from the data - in our case this is only QCD
// remember that signal is the last 'background' in the list....
// creates a clone histogram Hist provided by a caller
// Tight applies only to the signal...
//-----------------------------------------------------------------------------
void analysis::plot_dmboverb(const char* Module  , 
			     const char* HistSet , 
			     int         RunRange,
			     int         Bin     ,
			     const char* HistName,
			     int         Rebin   ,
			     int         McFlag  ,
			     double      XMin    ,
			     double      XMax    ,
			     const char* Axislabel,
			     int         Lumi    ,
			     int    SubtractSignal){

  TH1        *Hist(0);
  TH1        *h_dat, *h_bgr; //, *hist, *h2;
  aprocess   *bgr; //, *sig;
  TObjArray  list_of_bgr_hist;

  // char    name[500];
  int     nbins;
  double  q, err, err2, qbgr, n_proc, np;

  //  sprintf(name,"%s_%s_%s_%i_dmb_rebin_%i",Module,HistSet,HistName,Bin,Rebin);

  run_range_dat_t* rr = fRR+RunRange;

  if (SubtractSignal == 0) n_proc = rr->fNBgr-1;
  else                     n_proc = rr->fNBgr;

  rr->fDat->get_h1(Module,HistSet,Bin,HistName,Rebin,h_dat);

  //if (Hist) delete Hist;
  Hist = (TH1*) h_dat->Clone();
    
  for (int i=0; i<n_proc; i++) {
    bgr   = (aprocess*) rr->fBgr->At(i);
//-----------------------------------------------------------------------------
// skip QCD ... and, may be, other data-based histograms (fXSec = -1)
// for the time being it is only QCD ....
//-----------------------------------------------------------------------------
    if (bgr->fXSec > 0) {
      bgr->get_h1(Module,HistSet,Bin,HistName,Rebin,h_bgr);
//-----------------------------------------------------------------------------
// skip NULL histograms - skip undefined backgrounds....
//-----------------------------------------------------------------------------
      if (h_bgr != NULL) {
	list_of_bgr_hist.Add(h_bgr);
      }
      else {
	printf(" >>> analysis::get_dmb: ERROR SKIP UNDEFINED %s\n",
	       bgr->GetName());
      }
    }
  }
//-----------------------------------------------------------------------------
// form DMB histogram (assume no underflows and overflows)
//-----------------------------------------------------------------------------
  nbins = Hist->GetNbinsX();
//-----------------------------------------------------------------------------
// number fo histograms may be less then the number of processes - 
// do not subtract QCD
//-----------------------------------------------------------------------------
  np    = list_of_bgr_hist.GetEntries();
  for (int i=1; i<=nbins; i++) {
    q    = h_dat->GetBinContent(i);
    err  = h_dat->GetBinError  (i);
    err2 = err*err;
    qbgr = 0;
    for (int j=0; j<np; j++) {
      h_bgr = (TH1*) list_of_bgr_hist.At(j);
      qbgr += h_bgr->GetBinContent(i);
      //err   = h_bgr->GetBinError  (i); //Aidan 110309 just want stat error on data.  
                                         // show bck error as band on ratio plot
      //err2 += err*err;
    }
    q   -= qbgr;
    err  = TMath::Sqrt(err2);

    // (dt-bg)/bg
    if(qbgr>0){
      q /= qbgr;
      err /= qbgr;
    }
    else{
      q = 0.;
      err = 0.;
    }
    
    Hist->SetBinContent(i,q);
    Hist->SetBinError  (i,err);
  }

  TH1* uncband_plus  = (TH1*)Hist->Clone("uncband_plus");
  TH1* uncband_minus = (TH1*)Hist->Clone("uncband_minus");
  for (int i=1; i<=nbins; i++) {
    if(fabs(Hist->GetBinContent(i))>0.){
      uncband_plus->SetBinContent(i,0.1);// flat 10% uncertainty on bck, dominated by DY
      uncband_plus->SetBinError(i,0.);// flat 10% uncertainty on bck, dominated by DY
      uncband_minus->SetBinContent(i,-0.1);// flat 10% uncertainty on bck, dominated by DY
      uncband_minus->SetBinError(i,-0.);// flat 10% uncertainty on bck, dominated by DY
    }
  }      
  uncband_plus->SetLineColor(10);
  uncband_plus->SetFillColor(5);
  uncband_minus->SetLineColor(10);
  uncband_minus->SetFillColor(5);

  if (XMin != -999.) {
    Hist->GetXaxis()->SetRangeUser(XMin,XMax);
  }

  Hist->SetStats(0);
  Hist->GetYaxis()->SetTitle("(Dt-Exp)/Exp");
  Hist->GetYaxis()->SetTitleSize(.15);
  Hist->GetYaxis()->SetLabelSize(.2);
  Hist->GetYaxis()->SetTitleOffset(.25);
  Hist->GetYaxis()->SetNdivisions(503);
  //Hist->GetXaxis()->SetTitle(Axislabel);
  Hist->GetXaxis()->SetTitleSize(.2);
  Hist->GetXaxis()->SetLabelSize(.2);
  Hist->GetXaxis()->SetTitleOffset(.3);

  Hist->SetMaximum(.5);
  Hist->SetMinimum(-.5);

  Hist->Draw("ep");

  uncband_plus->Draw("same");
  uncband_minus->Draw("same");
  Hist->Draw("ep,same");

  list_of_bgr_hist.Clear();

  return;
}



//-----------------------------------------------------------------------------
// data minus sum of all the backgrounds
// skip backgrounds defined from the data - in our case this is only QCD
// remember that signal is the last 'background' in the list....
// creates a clone histogram Hist provided by a caller
// Tight applies only to the signal...
//-----------------------------------------------------------------------------
int analysis::get_dmb(const char* Module        , 
		      const char* HistSet       , 
		      int         RunRange      ,
		      int         Bin           ,
		      const char* HistName      ,
		      int         SubtractSignal, 
		      int         Rebin         ,
		      TH1*&       Hist          ) {
  //

  TH1        *h_dat, *h_bgr; // , *hist, *h2;
  aprocess   *bgr; // , *sig;
  TObjArray  list_of_bgr_hist;

  // char    name[500];
  int     nbins;
  double  q, err, err2, qbgr, n_proc, np;

  //  sprintf(name,"%s_%s_%s_%i_dmb_rebin_%i",Module,HistSet,HistName,Bin,Rebin);

  run_range_dat_t* rr = fRR+RunRange;

  if (SubtractSignal == 0) n_proc = rr->fNBgr-1;
  else                     n_proc = rr->fNBgr;

  rr->fDat->get_h1(Module,HistSet,Bin,HistName,Rebin,h_dat);

//   hist = (TH1*) fListOfDmbHistograms->FindObject(name);
//   if (hist != 0) {
//     fListOfDmbHistograms->Remove(hist);
//     delete hist;
//   }

  if (Hist) delete Hist;
  Hist = (TH1*) h_dat->Clone();
    
//  fListOfDmbHistograms->Add(hist);
//-----------------------------------------------------------------------------
// store (possibly rebinned) background histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<n_proc; i++) {
    bgr   = (aprocess*) rr->fBgr->At(i);
//-----------------------------------------------------------------------------
// skip QCD ... and, may be, other data-based histograms (fXSec = -1)
// for the time being it is only QCD ....
//-----------------------------------------------------------------------------
    if (bgr->fXSec > 0) {
      bgr->get_h1(Module,HistSet,Bin,HistName,Rebin,h_bgr);
//-----------------------------------------------------------------------------
// skip NULL histograms - skip undefined backgrounds....
//-----------------------------------------------------------------------------
      if (h_bgr != NULL) {
	list_of_bgr_hist.Add(h_bgr);
      }
      else {
	printf(" >>> analysis::get_dmb: ERROR SKIP UNDEFINED %s\n",
	       bgr->GetName());
      }
    }
  }
//-----------------------------------------------------------------------------
// form DMB histogram (assume no underflows and overflows)
//-----------------------------------------------------------------------------
  nbins = Hist->GetNbinsX();
//-----------------------------------------------------------------------------
// number fo histograms may be less then the number of processes - 
// do not subtract QCD
//-----------------------------------------------------------------------------
  np    = list_of_bgr_hist.GetEntries();
  for (int i=1; i<=nbins; i++) {
    q    = h_dat->GetBinContent(i);
    err  = h_dat->GetBinError  (i);
    err2 = err*err;
    qbgr = 0;
    for (int j=0; j<np; j++) {
      h_bgr = (TH1*) list_of_bgr_hist.At(j);
      qbgr += h_bgr->GetBinContent(i);
      err   = h_bgr->GetBinError  (i);
      err2 += err*err;
    }
    q   -= qbgr;
    err  = TMath::Sqrt(err2);
    
    Hist->SetBinContent(i,q);
    Hist->SetBinError  (i,err);
  }

  list_of_bgr_hist.Clear();

  return 0;
}


//-----------------------------------------------------------------------------
// data minus sum of all the backgrounds
// skip backgrounds defined from the data - in our case this is only QCD
// remember that signal is the last 'background' in the list....
// creates a clone histogram Hist provided by a caller
// Tight applies only to the signal...
//-----------------------------------------------------------------------------
int analysis::get_dmb_1(const char* Module        , 
			const char* HistSet       , 
			int         RunRange      ,
			int         Bin           ,
			const char* HistName      ,
			double      XMin          ,
			double      XMax          ,
			int         Rebin         ,
			TH1*&       Hist          ) {
  //

  TH1        *h_dat, *h_bgr; // , *hist, *h2;
  aprocess   *bgr; // , *sig;
  TObjArray  list_of_bgr_hist;

  // char    name[500];
  int     nbins, np;
  double  q, err, err2, qbgr, n_proc, e1, e2, sum_dat, sum_sig, x, sf;

  //  sprintf(name,"%s_%s_%s_%i_dmb_rebin_%i",Module,HistSet,HistName,Bin,Rebin);

  run_range_dat_t* rr = fRR+RunRange;

  n_proc = rr->fNBgr-1;

  rr->fDat->get_h1(Module,HistSet,Bin,HistName,Rebin,h_dat);

//   hist = (TH1*) fListOfDmbHistograms->FindObject(name);
//   if (hist != 0) {
//     fListOfDmbHistograms->Remove(hist);
//     delete hist;
//   }

  if (Hist) delete Hist;
  Hist = (TH1*) h_dat->Clone();
    
//  fListOfDmbHistograms->Add(hist);
//-----------------------------------------------------------------------------
// store (possibly rebinned) background histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<n_proc; i++) {
    bgr   = (aprocess*) rr->fBgr->At(i);
//-----------------------------------------------------------------------------
// skip QCD ... and, may be, other data-based histograms (fXSec = -1)
// for the time being it is only QCD ....
//-----------------------------------------------------------------------------
    if (bgr->fXSec > 0) {
      bgr->get_h1(Module,HistSet,Bin,HistName,Rebin,h_bgr);
//-----------------------------------------------------------------------------
// skip NULL histograms - skip undefined backgrounds....
//-----------------------------------------------------------------------------
      if (h_bgr != NULL) {
	list_of_bgr_hist.Add(h_bgr);
      }
      else {
	printf(" >>> analysis::get_dmb: ERROR SKIP UNDEFINED %s\n",
	       bgr->GetName());
      }
    }
  }
//-----------------------------------------------------------------------------
// form DMB histogram (assume no underflows and overflows)
//-----------------------------------------------------------------------------
  nbins = Hist->GetNbinsX();
//-----------------------------------------------------------------------------
// number fo histograms may be less then the number of processes - 
// do not subtract QCD
//-----------------------------------------------------------------------------
  np    = list_of_bgr_hist.GetEntries();
  for (int i=1; i<=nbins; i++) {
    q    = h_dat->GetBinContent(i);
    err  = h_dat->GetBinError  (i);
    err2 = err*err;
    qbgr = 0;
    for (int j=0; j<np; j++) {
      h_bgr = (TH1*) list_of_bgr_hist.At(j);
      qbgr += h_bgr->GetBinContent(i);
      err   = h_bgr->GetBinError  (i);
      err2 += err*err;
    }
    q   -= qbgr;
    err  = TMath::Sqrt(err2);
    
    Hist->SetBinContent(i,q);
    Hist->SetBinError  (i,err);
  }
//-----------------------------------------------------------------------------
// now scale the MC signal histogram
//-----------------------------------------------------------------------------
  sum_dat = 0;
  sum_sig = 0;
  rr->fSig->get_h1(Module,HistSet,Bin,HistName,Rebin,h_bgr);
  for (int i=1; i<=nbins; i++) {
    x    = Hist->GetXaxis()->GetBinCenter(i);
    if ((x >= XMin) && (x <= XMax)) {
      sum_dat += Hist->GetBinContent(i);
      sum_sig += h_bgr->GetBinContent(i);
    }
  }
  sf = sum_dat/sum_sig;
//-----------------------------------------------------------------------------
// normalization coefficient is defined, perform final scaled subtraction
//-----------------------------------------------------------------------------
  for (int i=1; i<=nbins; i++) {
    x    = Hist->GetXaxis()->GetBinCenter(i);
    if ((x >= XMin) && (x <= XMax)) {
      Hist->SetBinContent(i,0);
      Hist->SetBinError  (i,0);
    }
    else {
      q   = Hist->GetBinContent(i)-h_bgr->GetBinContent(i)*sf;
      if (q < 0) q = 0;

      e1  = Hist->GetBinError(i);
      e2  = h_bgr->GetBinError(i);
      err = TMath::Sqrt(e1*e1+sf*sf*e2*e2);
      Hist->SetBinContent(i,q);
      Hist->SetBinError  (i,err);
    }
  }

  list_of_bgr_hist.Clear();

  return 0;
}


//-----------------------------------------------------------------------------
// sum of all the backgrounds
// remember that signal is the last 'background' in the list,
// but it has McFlag=3
// creates a clone histogram Hist provided by a caller
// what is a background estimate is data-driven ? - still count it in!
//-----------------------------------------------------------------------------
int analysis::GetBackgroundHistogram(const char* Module        , 
				     const char* HistSet       , 
				     int         RunRange      ,
				     int         Bin           ,
				     const char* HistName      ,
				     int         Rebin         ,
				     TH1*&       Hist          ) {

  TH1        *h_dat, *h_bgr; // , *hist, *h2;
  aprocess   *bgr ; //, *sig;
  TObjArray  list_of_bgr_hist;

  int     nbins;
  double  /*q,*/ err, err2, qbgr, n_proc, np;

  run_range_dat_t* rr = fRR+RunRange;

  n_proc = rr->fNBgr;

  rr->fDat->get_h1(Module,HistSet,Bin,HistName,Rebin,h_dat);

  if (Hist) delete Hist;
  Hist = (TH1*) h_dat->Clone();
//-----------------------------------------------------------------------------
// data are not in rr->fBgr, they are in fDat
// ignore MC signal (MC_FLAG == 3) histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<n_proc; i++) {
    bgr   = (aprocess*) rr->fBgr->At(i);
    if (bgr->GetMcFlag() < 3) {
      bgr->get_h1(Module,HistSet,Bin,HistName,Rebin,h_bgr);
//-----------------------------------------------------------------------------
// skip NULL histograms - skip undefined backgrounds....
//-----------------------------------------------------------------------------
      if (h_bgr != NULL) {
	list_of_bgr_hist.Add(h_bgr);
      }
      else {
	printf(" >>> analysis::get_dmb: ERROR SKIP UNDEFINED %s\n",bgr->GetName());
      }
    }
  }
//-----------------------------------------------------------------------------
// form DMB histogram - save also overflows and underflows
//-----------------------------------------------------------------------------
  nbins = Hist->GetNbinsX();
//-----------------------------------------------------------------------------
// number of histograms may be less then the number of processes - 
// do not subtract QCD
//-----------------------------------------------------------------------------
  np    = list_of_bgr_hist.GetEntries();
  for (int i=0; i<=nbins+1; i++) {
    err  = 0;
    err2 = err*err;
    qbgr = 0;
    for (int j=0; j<np; j++) {
      h_bgr = (TH1*) list_of_bgr_hist.At(j);
      qbgr += h_bgr->GetBinContent(i);
      err   = h_bgr->GetBinError  (i);
      err2 += err*err;
    }
    err  = TMath::Sqrt(err2);
    
    Hist->SetBinContent(i,qbgr);
    Hist->SetBinError  (i,err);
  }

  list_of_bgr_hist.Clear();

  return 0;
}


//-----------------------------------------------------------------------------
// default: for backward compatibility, 'Process' defines the run range
//-----------------------------------------------------------------------------
double analysis::GetFitQEvents(const aprocess* Bgr    , 
			       const char*     Module ,
			       const char*     HistSet,
			       int             RunRange,
			       int             Bin    ,
			       const char*     HistName) const {

  double qev = Bgr->GetFitQEvents(Module,HistSet,Bin,HistName);

  return qev;
}

//-----------------------------------------------------------------------------
// default: for backward compatibility, 'Process' defines the run range
//-----------------------------------------------------------------------------
double analysis::GetIDEff     (const aprocess* Proc    , 
			       const char*     Module  ,
			       const char*     HistSet ,
			       int             RunRange,
			       int             Bin     ,
			       const char*     HistName) const {

  double id_eff = 1;

  return id_eff;
}

//-----------------------------------------------------------------------------
// default: for backward compatibility, 'Process' defines the run range
//-----------------------------------------------------------------------------
double analysis::GetIDEffSF   (const aprocess* Proc    , 
			       const char*     Module  ,
			       const char*     HistSet ,
			       int             RunRange,
			       int             Bin     ,
			       const char*     HistName) const {

  double id_eff_sf = Proc->GetIDEffSF(Module,HistSet,Bin,HistName);

  return id_eff_sf;
}

//-----------------------------------------------------------------------------
// default: for backward compatibility, 'Process' defines the run range
//-----------------------------------------------------------------------------
double analysis::GetTrEff     (const aprocess* Proc    , 
			       const char*     Module  ,
			       const char*     HistSet ,
			       int             RunRange,
			       int             Bin     ,
			       const char*     HistName) const {

  double tr_eff = Proc->GetTrEff(Module,HistSet,Bin,HistName);

  return tr_eff;
}

//-----------------------------------------------------------------------------
// default: for backward compatibility, 'Process' defines the run range
//-----------------------------------------------------------------------------
double analysis::GetTrEffSF   (const aprocess* Proc    , 
			       const char*     Module  ,
			       const char*     HistSet ,
			       int             RunRange,
			       int             Bin     ,
			       const char*     HistName) const {

  double tr_eff_sf = Proc->GetIDEffSF(Module,HistSet,Bin,HistName);

  return tr_eff_sf;
}


//-----------------------------------------------------------------------------
// default: for backward compatibility, 'Process' defines the run range
//-----------------------------------------------------------------------------
double analysis::GetAcceptance(const aprocess* Process , 
			       const char*     Module  ,
			       const char*     HistSet ,
			       int             RunRange,
			       int             Bin     ,
			       const char*     HistName) const {

  double acc = Process->GetAcceptance(Module,HistSet,Bin,HistName);

  return acc;
}

//-----------------------------------------------------------------------------
// assume all processes in fBgr and that no two processes have the same name
//-----------------------------------------------------------------------------
aprocess* analysis::GetProcess(const char* Name, int RunRange) const {
  aprocess  *p, *process(0);

  run_range_dat_t* rr = fRR+RunRange;
				        // first check data 
  p = rr->fDat;
  if (strcmp(p->GetName(),Name) == 0) {
    process = p;
  }
  else {
				        // check backgrounds
    int n = rr->fBgr->GetEntriesFast();
    for (int i=0; i<n; i++) {
      p = (aprocess*) rr->fBgr->UncheckedAt(i);
      if (strcmp(p->GetName(),Name) == 0) {
	process = p;
	break;
      }
    }
  }

  return process;
}


//-----------------------------------------------------------------------------
// I < 0: data, i >= 0: background
//-----------------------------------------------------------------------------
aprocess* analysis::GetProcess(int I, int RunRange) const {
  aprocess  *process(0);

  run_range_dat_t* rr = fRR+RunRange;

  if (I < 0) {
    process = rr->fDat;
  }
  else {
    process = (aprocess*) rr->fBgr->UncheckedAt(I);
  }

  return process;
}

//-----------------------------------------------------------------------------
double analysis::GetCrossSection(const char*   Module  ,
				 const char*   HistSet ,
				 int           RunRange, 	
				 int           Bin     ,
				 const char*   HistName) const {

  Error("GetCrossSection","NOT IMPLEMENTED YET");
  return -1;
}


//-----------------------------------------------------------------------------
void analysis::Print(const char* Module  ,
		     const char* HistSet ,
		     int         Bin     ,
		     const char* HistName,
		     const char* Opt     ) const {

  int        n_proc;
  aprocess*  bgr; 
  run_range_dat_t* rr;

  for (int ir=0; ir<fNRunRanges; ir++) {
    rr = fRR+ir;
    rr->fDat->Print(Module,HistSet,Bin,HistName,"banner+data");
    n_proc = rr->fBgr->GetEntries();
    for (int i=0; i<n_proc; i++) {
      bgr   = (aprocess*) rr->fBgr->At(i);
      bgr->Print(Module,HistSet,Bin,HistName,"data");
    }
  }
}

//-----------------------------------------------------------------------------
void analysis::Print(const char* Opt     ) const {

  int        n_proc;
  aprocess*  bgr;
  run_range_dat_t* rr;

  printf(" analysis channel : Name= %s  Title=%s\n",GetName(),GetTitle());

  for (int ir=0; ir<fNRunRanges; ir++) {
    rr = fRR+ir;
    n_proc = rr->fBgr->GetEntries();
    for (int i=0; i<n_proc; i++) {
      bgr   = (aprocess*) rr->fBgr->At(i);
      bgr->Print(Opt);
    }
  }
}
