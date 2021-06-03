///////////////////////////////////////////////////////////////////////////////
// use numbers from 
// to be compiled

#include "Stntuple/val/stntuple_val_functions.hh"
#include "TEnv.h"
#include "TF2.h"

#include "murat/limits/channel.hh"

namespace murat {

double f_bgr_cosmics(double* X, double* P) {
//-----------------------------------------------------------------------------
// to begin with, assume cosmics uniform
// missing scale factor to CORSIKA - set to 1, add a 5% uncertainty
// assume P = X[0] ; T = X[1]
//-----------------------------------------------------------------------------
  float p              = X[0];
  float t              = X[1];
  //  float rho_per_mev    = (15.15+0.4037*(p-105.))/5./110.73 ; // background per MeV, integrated over time 700-1695 = 995 ns
  float rho_per_mev    = (14.57+0.3199*(p-105.))/5./110.73 ;     // 2025lo background per MeV, integrated over time 700-1695 = 995 ns
  rho_per_mev         += 0.5/50/3.7;                             // half event from 2025hi (positrons+electrons = 1 event)/2
					                         // account for T>1650... - correction on the tail "eats up" 2.5%,
					                         // renormalize the integral
  if (t > 1650) rho_per_mev = rho_per_mev*(1.-(t-1650)/40.)*1.025; 
  if (t > 1690) rho_per_mev = 0;
  
  float time_window    = 995.;
  float rho_per_mev_ns = rho_per_mev/time_window;

  return rho_per_mev_ns;
}

//-----------------------------------------------------------------------------
void plot_bgr_cosmics() {

  TF2* bgr_cosmics_2D = new TF2("bgr_cosmics_2D",f_bgr_cosmics,100,110,500,1700);

  bgr_cosmics_2D->Draw();
}

channel::channel(const char* ChannelName, float NPOT_1B, float NPOT_2B, float ExtraSF) : TNamed(ChannelName,ChannelName) {
    // HistName: generic: "time", "mom", "time_vs_mom"
    // allow additional scaling of the histograms, by default, ExtraSF = 1

  const char* mu2e_hist_dir = "/projects/mu2e/hist";
    
  fHist      = nullptr;
  fCapture   = 0.609;	                // probability of muon capture in Al
  fMusr      = 1.59e-3;	                // muon stopping rate (stopped muons/POT)
    
  kLumiSF_1B = NPOT_1B;	                // SU2020 1B-mode scale factor, see su2020/analysis/dio.org
  kLumiSF_2B = NPOT_2B;                 // SU2020 2B-mode scale factor, see su2020/analysis/dio.org

  const char* hist_dir = gEnv->GetValue("mu2e.HistDir",mu2e_hist_dir);
  fHistDir = Form("%s/su2020",hist_dir);

  printf("ChannelName, SF1B, SF2B: %-10s  %9.3e %9.3e ",ChannelName,kLumiSF_1B,kLumiSF_2B);
    
  TString channel_name = ChannelName;
    
  if (channel_name == "DIO") {
    const char* dsid    = "su2020.fele2s51b1";
    const char* ana_job = "su2020_track_ana.1011";
    double      ngen (1.e7), erange(10.);            // properties of the dataset 

    double sf1b = NPOT_1B*fMusr*(1-fCapture)*erange/ngen*ExtraSF;
    double sf2b = NPOT_2B*fMusr*(1-fCapture)*erange/ngen*ExtraSF;

    printf("DIO : sf1b, sf2b = %12.5e %12.5e\n",sf1b,sf2b);
    
    fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_TrackAna","trk_2006/p_vs_time")->Clone("DIO_t_vs_p");
    fTimeVsMom->Scale(sf1b);
    TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_TrackAna","trk_2007/p_vs_time")->Clone("tmp");
    fTimeVsMom->Add(h,sf2b);
    delete h;
  }
  else if (channel_name == "CE") {
    const char* dsid    = "su2020.cele0s61b1";
    const char* ana_job = "su2020_track_ana.1011";
    double      ngen (1.e6); 
    
    double sf1b = NPOT_1B*fMusr*fCapture/ngen*ExtraSF;
    double sf2b = NPOT_2B*fMusr*fCapture/ngen*ExtraSF;
    
    printf("CE : sf1b, sf2b = %12.5e %12.5e\n",sf1b,sf2b);
    
    fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_TrackAna","trk_2004/p_vs_time")->Clone("CE_t_vs_p");
    fTimeVsMom->Scale(sf1b);
    TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_TrackAna","trk_2005/p_vs_time")->Clone("tmp");
    fTimeVsMom->Add(h,sf2b);
    delete h;
  }
  else if (channel_name == "PbarAnni") {
    const char* dsid    = "su2020.pbar0s61b0";
    const char* ana_job = "su2020_pbar_ana.1010";
    
    double ngen       (2.e7);         // number of pbars generated at ST
    double nsa_per_pot(4.84e-18);     // number of pbars stopped in ST per POT

    double sf1b = NPOT_1B*nsa_per_pot/ngen*ExtraSF;
    double sf2b = NPOT_2B*nsa_per_pot/ngen*ExtraSF;
    
    printf("PbarAnni : NPOT_1B, NPOT_2B, sf1b, sf2b = %12.5e %12.5e %12.5e %12.5e\n",NPOT_1B,NPOT_2B,sf1b,sf2b);
     
    fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_PbarAna","trk_2000/p_vs_time")->Clone("PbarAnni_t_vs_p");
    fTimeVsMom->Scale(sf1b);
    TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_PbarAna","trk_2000/p_vs_time")->Clone("tmp");
    fTimeVsMom->Add(h,sf2b);
    delete h;
  }
  else if (channel_name == "PbarRPCe") {
    const char* dsid    = "su2020.rpce3s41b0";
    const char* ana_job = "su2020_pbar_ana.1010";

    double ngen(5.e7);                 // number of generated RPCe events 
    double nsp_per_pot(1.87e-16);      // number of stopped pi^- produced by pbars, per POT
    double br_rpc     (2.15e-2);       // BR(RPC)

    double sf1b = NPOT_1B*nsp_per_pot/ngen*br_rpc*ExtraSF;
    double sf2b = NPOT_2B*nsp_per_pot/ngen*br_rpc*ExtraSF;

    printf("PbarRPCe : sf1b, sf2b = %12.5e %12.5e\n",sf1b,sf2b);
    
    fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_PbarAna","trk_2000/p_vs_time")->Clone("PbarRPCe_t_vs_p");
    fTimeVsMom->Scale(sf1b);
    TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_PbarAna","trk_2000/p_vs_time")->Clone("tmp");
    fTimeVsMom->Add(h,sf2b);
    delete h;
  }
  else if (channel_name == "PbarRPCi") {
    const char* dsid    = "su2020.rpci3s41b0";
    const char* ana_job = "su2020_pbar_ana.1010";

    double ngen(5.e5);                 // number of generated RPCe events 
    double nsp_per_pot(1.87e-16);      // number of stopped pi^- produced by pbars, per POT
    double br_rpc     (2.15e-2);       // BR(RPC)
    double rho        (6.94e-3);       // internal conversion fraction, RPC

    double sf1b = NPOT_1B*nsp_per_pot/ngen*br_rpc*rho*ExtraSF;
    double sf2b = NPOT_2B*nsp_per_pot/ngen*br_rpc*rho*ExtraSF;
     
    printf("PbarRPCi : sf1b, sf2b = %12.5e %12.5e\n",sf1b,sf2b);

    fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_PbarAna","trk_2000/p_vs_time")->Clone("PbarRPCi_t_vs_p");
    fTimeVsMom->Scale(sf1b);
    TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_PbarAna","trk_2000/p_vs_time")->Clone("tmp");
    fTimeVsMom->Add(h,sf2b);
    delete h;
  }
  else if (channel_name == "PbarRPC") {
//-----------------------------------------------------------------------------
// sum of internal and external , just assume gPbarRPCe and gPbarRPCi are already initialized
//-----------------------------------------------------------------------------
    printf("\n");
    channel pbar_rpci("PbarRPCi",NPOT_1B,NPOT_2B,ExtraSF);
    channel pbar_rpce("PbarRPCe",NPOT_1B,NPOT_2B,ExtraSF);
    fTimeVsMom = (TH2F*) pbar_rpci.fTimeVsMom->Clone("PbarRPC_t_vs_p");
    fTimeVsMom->Add(pbar_rpce.fTimeVsMom);
  }
  else if (channel_name == "PbarTOT") {
//-----------------------------------------------------------------------------
// sum of internal and external , just assume gPbarRPCe and gPbarRPCi are already initialized
//-----------------------------------------------------------------------------
    printf("\n");
    channel pbar_anni("PbarAnni",NPOT_1B,NPOT_2B,ExtraSF);
    channel pbar_rpc ("PbarRPC" ,NPOT_1B,NPOT_2B,ExtraSF);
    fTimeVsMom = (TH2F*) pbar_anni.fTimeVsMom->Clone("PbarTOT_t_vs_p");
    fTimeVsMom->Add(pbar_rpc.fTimeVsMom);
  }
  else if (channel_name == "RPCe") {
//-----------------------------------------------------------------------------
// external RPC, generation parameters
// ($npot/$npot_gen)*($nstops/$eff_ngen)*$br_rpc*$fr
//-----------------------------------------------------------------------------
    const char* dsid    = "su2020.rpce0s51b1";
    const char* ana_job = "su2020_rpc_ana.1011";

    double npot_gen(1.e8);          // N(POT) for pi^- beam transport
    double eff_450 (2.423e-2);      // efficiency of the T>450 ns cut-off

    double nstops  (210551);        // N(stopped pi^-)

    double ngen_ph(1.e8);           // N(generated photons), from ST
    double br_rpc (2.15e-2);        // BR(RPC)

    double eff_ngen = ngen_ph/eff_450; // effective number of generated photons

    double sf1b = (NPOT_1B/npot_gen)*(nstops/eff_ngen)*br_rpc*ExtraSF;
    double sf2b = (NPOT_2B/npot_gen)*(nstops/eff_ngen)*br_rpc*ExtraSF;
    
    printf("RPCe : sf1b, sf2b = %12.5e %12.5e\n",sf1b,sf2b);

    fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_RPCAna","trk_2002/p_vs_time")->Clone("RPCe_t_vs_p");
    double sum1 = GetIntegral(103.85,105.1,700,1700);
    fTimeVsMom->Scale(sf1b);
    double sum2 = GetIntegral(103.85,105.1,700,1700);
    TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_RPCAna","trk_2002/p_vs_time")->Clone("tmp");
    fTimeVsMom->Add(h,sf2b);
    double sum3 = GetIntegral(103.85,105.1,700,1700);

    printf("sum1, sum2, sum3: %12.5e %12.5e %12.5e\n",sum1,sum2,sum3);
    
    delete h;
  }
  else if (channel_name == "RPCi") {
//-----------------------------------------------------------------------------
// internal RPC, generation parameters 
//-----------------------------------------------------------------------------
    const char* dsid    = "su2020.rpci0s51b1";
    const char* ana_job = "su2020_rpc_ana.1011";
    
    double npot_gen(1.e8);          // N(POT) for pi^- beam transport
    double eff_450 (2.423e-2);      // efficiency of the T>450 ns cut-off
    double nstops  (210551);        // N(stopped pi^-)

    double ngen_ph (1.e6);          // N(generated photons), from ST
    double rho     (6.9e-3);        // internal conversion fraction
    double br_rpc  (2.15e-2);       // BR(RPC)

    double eff_ngen = ngen_ph/eff_450; // effective number of generated photons
    
    double sf1b = (NPOT_1B/npot_gen)*(nstops/eff_ngen)*br_rpc*rho*ExtraSF;
    double sf2b = (NPOT_2B/npot_gen)*(nstops/eff_ngen)*br_rpc*rho*ExtraSF;
    
    printf("RPCi : sf1b, sf2b = %12.5e %12.5e\n",sf1b,sf2b);

    fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_RPCAna","trk_2002/p_vs_time")->Clone("RPCi_t_vs_p");
    fTimeVsMom->Scale(sf1b);
    TH2F* h    = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_RPCAna","trk_2002/p_vs_time")->Clone("tmp");
    fTimeVsMom->Add(h,sf2b);
    delete h;
  }
  else if (channel_name == "RPC") {
//-----------------------------------------------------------------------------
// sum of internal and external , just assume gRPCe and gRPCi are already initialized
//-----------------------------------------------------------------------------
    printf("\n");
    channel rpci("RPCi",NPOT_1B,NPOT_2B,ExtraSF);
    channel rpce("RPCe",NPOT_1B,NPOT_2B,ExtraSF);
    fTimeVsMom = (TH2F*) rpci.fTimeVsMom->Clone("RPC_t_vs_p");
    fTimeVsMom->Add(rpce.fTimeVsMom);
  }
  else if (channel_name == "Cosmics") {
//-----------------------------------------------------------------------------
// kludge - initialize using a parameterization
// initlaize "RPCi" just to clone a histogram - that's the safest
//-----------------------------------------------------------------------------
    double x[2], p[2];
      
    const char* dsid    = "su2020.cry33s51b0";
    const char* ana_job = "su2020_cosmic_ana.1010";
    
    fTimeVsMom = (TH2F*) gh2(Form("%s/%s.%s.hist",GetHistDir(),dsid,ana_job),"su2020_CosmicAna","trk_2000/p_vs_time")->Clone("Cosmics_t_vs_p");
    fTimeVsMom->Reset();
    
    double xmin = fTimeVsMom->GetXaxis()->GetXmin();
    double binx = fTimeVsMom->GetXaxis()->GetBinWidth(1);

    double ymin = fTimeVsMom->GetYaxis()->GetXmin();
    double biny = fTimeVsMom->GetYaxis()->GetBinWidth(1);

    printf("Cosmics: binx = %10.4f biny = %10.4f\n",binx,biny);
    
    int nbx = fTimeVsMom->GetNbinsX();
    int nby = fTimeVsMom->GetNbinsY();
    for(int ix=1; ix<=nbx; ix++) {
      x[0] = xmin+(ix-0.5)*binx;
      for(int iy=1; iy<=nby; iy++) {
	x[1] = ymin+(iy-0.5)*biny;
	double bgr = f_bgr_cosmics(x,p)*binx*biny*ExtraSF;
	fTimeVsMom->SetBinContent(ix,iy,bgr);
	fTimeVsMom->SetBinError  (ix,iy,bgr*0.1);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// create momentum histogram with given  timing cutoff
//-----------------------------------------------------------------------------
TH1D* channel::CreateMomHist(double TMin, double TMax) {

    double tmin = fTimeVsMom->GetYaxis()->GetXmin();
    double tmax = fTimeVsMom->GetYaxis()->GetXmax();
    double bin  = fTimeVsMom->GetYaxis()->GetBinWidth(1);

    int iy1(0), iy2(-1); // defaults

    if (TMin > tmin) iy1 = (TMin+1.e-6*bin - tmin)/bin + 1;
    if (TMax < tmax) iy2 = (TMax - tmin)/bin + 1;
    
    printf("iy1, iy2 = %5i %5i \n",iy1,iy2);

    TH1D* htemp   = fTimeVsMom->ProjectionX("htemp",iy1,iy2);

    TString name  = Form("h_%s_mom_hist",GetName());
    TH1D* hist    = (TH1D*) htemp->Clone(name.Data());
    delete htemp;
    
    return hist;
  }

//-----------------------------------------------------------------------------
// create momentum histogram with given  timing cutoff
//-----------------------------------------------------------------------------
TH1D* channel::CreateTimeHist(double PMin, double PMax) {

  double pmin = fTimeVsMom->GetXaxis()->GetXmin();
  double pmax = fTimeVsMom->GetXaxis()->GetXmax();
  double bin  = fTimeVsMom->GetXaxis()->GetBinWidth(1);

  int ix1(0), ix2(-1); // defaults

  printf ("PMin, PMax, pmin, pmax, bin = %10.3f  %10.3f  %10.3f  %10.3f %10.4f\n",
	  PMin, PMax, pmin, pmax,bin);

  if (PMin > pmin) ix1 = (PMin+1.e-6*bin - pmin)/bin + 1;
  if (PMax < pmax) ix2 = (PMax - pmin)/bin + 1;

  printf("ix1, ix2 = %5i %5i \n",ix1,ix2);
    
  TH1D* htemp  = fTimeVsMom->ProjectionY("htemp",ix1,ix2);
  TString name = Form("h_%s_mom_hist",GetName());
  TH1D* hist   = (TH1D*) htemp->Clone(name);
  delete htemp;
    
  return hist;
}

//-----------------------------------------------------------------------------
// create running integral (T_i, TMax) histogram with given timing cutoff for
// the time distribution within given (PMin,PMax) momentum band
//-----------------------------------------------------------------------------
TH1D* channel::CreateTimeIntegralHist(double PMin, double PMax, double TMin, double TMax) {

  double pmin = fTimeVsMom->GetXaxis()->GetXmin();
  double pmax = fTimeVsMom->GetXaxis()->GetXmax();
  double binp = fTimeVsMom->GetXaxis()->GetBinWidth(1);
  
  double tmin = fTimeVsMom->GetYaxis()->GetXmin();
  double tmax = fTimeVsMom->GetYaxis()->GetXmax();
  double bint = fTimeVsMom->GetYaxis()->GetBinWidth(1);
  int    nbt  = fTimeVsMom->GetYaxis()->GetNbins();

  int ix1(0), ix2(-1), iy1(0), iy2(nbt); // defaults

  printf ("PMin, PMax, pmin, pmax, binp = %10.3f  %10.3f  %10.3f  %10.3f %10.4f\n",
	  PMin, PMax, pmin, pmax,binp);

  if (PMin > pmin) ix1 = (PMin+1.e-6*binp - pmin)/binp + 1;
  if (PMax < pmax) ix2 = (PMax - pmin)/binp + 1;
  
  if (TMin > tmin) iy1 = (TMin+1.e-6*bint - tmin)/bint + 1;
  if (TMax < tmax) iy2 = (TMax - tmin)/bint + 1;

  printf("ix1, ix2 = %5i %5i \n",ix1,ix2);
  printf("iy1, iy2 = %5i %5i \n",iy1,iy2);
    
  TH1D* htemp  = fTimeVsMom->ProjectionY("htemp",ix1,ix2);
  TString name = Form("h_%s_time_integral_hist",GetName());
  TH1D* hist   = (TH1D*) htemp->Clone(name);

  hist->Reset();

  for (int iy = iy1; iy<=iy2; iy++) {
    double sy  = 0;
    double sw = 0;
    for (int iyy = iy; iyy<=iy2; iyy++) {
      sy += htemp->GetBinContent(iyy);
      sw += htemp->GetBinError(iyy)*fTimeVsMom->GetBinError(iyy);
      //	printf("iyy, y, sy, err, sw : %5i %12.5e  %12.5e  %12.5e  %12.5e\n",
      //     iyy, fTimeVsMom->GetBinContent(iyy),sy,fTimeVsMom->GetBinError(iyy), sw);
    }
    hist->SetBinContent(iy,sy);
    hist->SetBinError  (iy,sqrt(sw));
    // printf("--------- iy, y, sw : %5i %12.5e  %12.5e\n",iy,sy,sw);
  }
    
  // delete htemp;
  return hist;
}
  
//-----------------------------------------------------------------------------
// create running integral (T_i, TMax) histogram with given timing cutoff for
// the time distribution within given (PMin,PMax) momentum band
//-----------------------------------------------------------------------------
TH1D* channel::CreateMomIntegralHist(double PMin, double PMax, double TMin, double TMax) {

  double pmin = fTimeVsMom->GetXaxis()->GetXmin();
  double pmax = fTimeVsMom->GetXaxis()->GetXmax();
  double binp = fTimeVsMom->GetXaxis()->GetBinWidth(1);
  //  int    nbp  = fTimeVsMom->GetXaxis()->GetNbins();
  
  double tmin = fTimeVsMom->GetYaxis()->GetXmin();
  double tmax = fTimeVsMom->GetYaxis()->GetXmax();
  double bint = fTimeVsMom->GetYaxis()->GetBinWidth(1);
  int    nbt  = fTimeVsMom->GetYaxis()->GetNbins();

  int ix1(0), ix2(-1), iy1(0), iy2(nbt); // defaults

  printf ("PMin, PMax, pmin, pmax, binp = %10.3f  %10.3f  %10.3f  %10.3f %10.4f\n",
	  PMin, PMax, pmin, pmax,binp);

  if (PMin > pmin) ix1 = (PMin+1.e-6*binp - pmin)/binp + 1;
  if (PMax < pmax) ix2 = (PMax - pmin)/binp + 1;
  
  if (TMin > tmin) iy1 = (TMin+1.e-6*bint - tmin)/bint + 1;
  if (TMax < tmax) iy2 = (TMax - tmin)/bint + 1;

  printf("ix1, ix2 = %5i %5i \n",ix1,ix2);
  printf("iy1, iy2 = %5i %5i \n",iy1,iy2);
  
  TH1D*   htemp = fTimeVsMom->ProjectionX("htemp",iy1,iy2);
  TString name  = Form("h_%s_mom_integral_hist",GetName());
  TH1D*   hist  = (TH1D*) htemp->Clone(name);
  hist->Reset();
  
  for (int ix = ix1; ix<=ix2; ix++) {
    double sx  = 0;
    double sw = 0;
    for (int ixx = ix; ixx<=ix2; ixx++) {
      sx += htemp->GetBinContent(ixx);
      sw += htemp->GetBinError(ixx)*fTimeVsMom->GetBinError(ixx);
      //	printf("iyy, y, sy, err, sw : %5i %12.5e  %12.5e  %12.5e  %12.5e\n",
      //     iyy, fTimeVsMom->GetBinContent(iyy),sy,fTimeVsMom->GetBinError(iyy), sw);
    }
    hist->SetBinContent(ix,sx);
    hist->SetBinError  (ix,sqrt(sw));
    // printf("--------- iy, y, sw : %5i %12.5e  %12.5e\n",iy,sy,sw);
  }
  
  delete htemp;
  return hist;
}

//-----------------------------------------------------------------------------
double channel::GetIntegral(float PMin, float PMax, float TMin, float TMax) {

  TAxis* ax = fTimeVsMom->GetXaxis();
  TAxis* ay = fTimeVsMom->GetYaxis();
  
  int nbp   = ax->GetNbins();
  int nbt   = ay->GetNbins();

  float sw  = 0;
  
  for (int ix=0; ix<nbp; ix++) {
    float p = ax->GetBinCenter(ix);
    if ((p < PMin) or (p > PMax))   continue;
    for (int iy=0; iy<nbt; iy++) {
      float t = ay->GetBinCenter(iy);
      if ((t < TMin) or (t > TMax))   continue;
      float w = fTimeVsMom->GetBinContent(ix,iy);
      sw += w;
    }
  }

  return sw;
}

  
//-----------------------------------------------------------------------------
double channel::FluctuateBackground(float PMin, float PMax, float TMin, float TMax) {
  // in a simplest form, dont fluctuate, just return the mean
  // to fluctuate, need uncertainties...
  
  double bgr = GetIntegral(PMin, PMax, TMin, TMax);
  // not finished

  return bgr;
}
};
