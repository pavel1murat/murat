///////////////////////////////////////////////////////////////////////////////
// 10 numbers per line are assumed ?
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
TH1* get_mu2e_hist(const char* Filename, const char* Histname) {
  TFile  *f;
  TH1    *hist(0);

  f = gROOT->GetFile(Filename);
  if (f == 0) f = TFile::Open(Filename);

  if (f) {
    hist = (TH1*) f->Get(Histname);
  }

  return hist;
}

//-----------------------------------------------------------------------------
void generate_dt_templates(const char* HistFile, const char* Module, const char* HistName, const char* Fout) {

  TH1F* he = gh1(HistFile,Module,HistName);
  he->Draw();

  int nb = he->GetXaxis()->GetNbins();

  FILE* f = fopen(Fout,"w");

  fprintf(f,"title:          %-s/%s/%s\n",HistFile,Module,HistName);
  fprintf(f,"name:           %-s\n",HistName);
  fprintf(f,"nbx,xmin,xmax: %5i %10.4f %10.4f\n",
	  nb,he->GetXaxis()->GetXmin(),he->GetXaxis()->GetXmax());

  int nmax(10), ip(0);
  
  for (int ib=1; ib<=nb; ib++) {
    float x = he->GetBinContent(ib);
    fprintf(f,"%10.1f",x);
    ip++;
    if (ip >= nmax) {
      fprintf(f,"\n");
      ip = 0;
    }
  }

  if (ip >= 0) {
    fprintf(f,"\n");
    ip = 0;
  }

  fclose(f);
  
}


//-----------------------------------------------------------------------------
void generate_xs_templates(const char* HistFile, const char* Module, const char* HistName, const char* Fout) {

  TH1F* he = gh1(HistFile,Module,HistName);
  he->Draw();

  int nb = he->GetXaxis()->GetNbins();

  FILE* f = fopen(Fout,"w");

  fprintf(f,"title:          %-s/%s/%s\n",HistFile,Module,HistName);
  fprintf(f,"name:           %-s\n",HistName);
  fprintf(f,"nbx,xmin,xmax: %5i %10.4f %10.4f\n",
	  nb,he->GetXaxis()->GetXmin(),he->GetXaxis()->GetXmax());

  int nmax(10), ip(0);
  
  for (int ib=1; ib<=nb; ib++) {
    float x = he->GetBinContent(ib);
    fprintf(f,"%10.1f",x);
    ip++;
    if (ip >= nmax) {
      fprintf(f,"\n");
      ip = 0;
    }
  }

  if (ip >= 0) {
    fprintf(f,"\n");
    ip = 0;
  }

  fclose(f);
  
}


//-----------------------------------------------------------------------------
// E/P vs Path: a 2D histogram
//-----------------------------------------------------------------------------
void generate_ep_templates(const char* HistFile, const char* Module, const char* HistName, const char* Fout) {

  TH2F* h = gh2(HistFile,Module,HistName);
  h->Draw();

  TAxis  *ax, *ay;

  ax = h->GetXaxis();
  ay = h->GetYaxis();

  int nbx = ax->GetNbins();
  int nby = ay->GetNbins();

  FILE* f = fopen(Fout,"w");

  fprintf(f,"title:          %-s/%s/%s\n",HistFile,Module,HistName);
  fprintf(f,"name:           %-s\n",HistName);
  fprintf(f,"nbx,xmin,xmax,nby,ymin,ymax: %5i %10.4f %10.4f %5i %10.4f %10.4f\n",
	  nbx,ax->GetXmin(),ax->GetXmax(),
	  nby,ay->GetXmin(),ay->GetXmax());

  int nmax(10), ip(0);
  
  for (int iy=1; iy<=nby; iy++) {
    for (int ix=1; ix<=nbx; ix++) {
      float val = h->GetBinContent(ix,iy);
      fprintf(f,"%10.1f",val);
      ip++;
      if (ip >= nmax) {
	fprintf(f,"\n");
	ip = 0;
      }
    }
  }

  if (ip >= 0) {
    fprintf(f,"\n");
    ip = 0;
  }

  fclose(f);
  
}


//-----------------------------------------------------------------------------
// E/P vs Path: a 2D histogram
// number of X bins: 200, 0-20 mm
// the 2D histogram is 200(0-20mm):1000(0-0.005keV)
// each X-slice has 10 bins (1mm) and needs to be rebinned x10 to end up with 100 bins
// Particle: "e" or "m"
//-----------------------------------------------------------------------------
void generate_dedx_templates(const char* HistFile, const char* HistName, const char* Particle, const char* Fout) {

  char   pname[200];
  
  TH2F* h = (TH2F*) get_mu2e_hist(HistFile,HistName);
  h->Draw();

  TAxis  *ax, *ay;

  ax = h->GetXaxis();
  ay = h->GetYaxis();

  int nbx = ax->GetNbins();
  int nby = ay->GetNbins();

  TFile* fout = TFile::Open(Fout,"recreate");

  int nxslices(11);

  for (int i=0; i<nxslices; i++) {
    TH1D* h_proj = (TH1D*) h->ProjectionY("py",10*i,10*i+9)->Rebin(10)->Clone(Form("htemp%s%i",Particle,i));
  }

  fout->Write();
  fout->Close();

  delete fout;
}


//-----------------------------------------------------------------------------
// generate the PID templates for the offline version v5_7
//-----------------------------------------------------------------------------
int generate_pid_templates_v5_7() {
  
  const char* ele_track_ana_fn  = "~/hist/mu2e/v5_7_0/e00s5700.track_ana.hist";
  const char* muo_track_ana_fn  = "~/hist/mu2e/v5_7_0/m00s5700.track_ana.hist";

  const char* ele_pid_ana_fn    = "~/hist/mu2e/v5_7_0/e00s5700.pid_ana.hist";
  const char* muo_pid_ana_fn    = "~/hist/mu2e/v5_7_0/m00s5700.pid_ana.hist";

  const char* fn_ele_dedx       = "~/hist/mu2e/v5_7_0/e00s5700.egun_stnmaker.hist";
  const char* fn_muo_dedx       = "~/hist/mu2e/v5_7_0/m00s5700.mgun_stnmaker.hist";

//-----------------------------------------------------------------------------
// DT templates
//-----------------------------------------------------------------------------
  generate_dt_templates(ele_track_ana_fn,"TrackAna","trk_19/dt","pid_ele_dt.tab");
  generate_dt_templates(muo_track_ana_fn,"TrackAna","trk_19/dt","pid_muo_dt.tab");
  
//-----------------------------------------------------------------------------
// E/P templates
//-----------------------------------------------------------------------------
  generate_ep_templates(ele_track_ana_fn,"TrackAna","trk_19/ep_vs_path","pid_ele_ep_vs_path.tab");
  generate_ep_templates(muo_track_ana_fn,"TrackAna","trk_19/ep_vs_path","pid_muo_ep_vs_path.tab");
  
//-----------------------------------------------------------------------------
// X(dR/dS) templates
//-----------------------------------------------------------------------------
  generate_xs_templates(ele_pid_ana_fn,"PidAna","pid_1/xdrds_vadim_ele","pid_ele_xdrds.tab");
  generate_xs_templates(muo_pid_ana_fn,"PidAna","pid_1/xdrds_vadim_ele","pid_muo_xdrds.tab");
  
//-----------------------------------------------------------------------------
// de/dx templates: read 2D histograms, store 10 slices in 1mm step
//-----------------------------------------------------------------------------
  generate_dedx_templates(fn_ele_dedx,"TrackRecoCheck/ehit_vs_path","e","pid_ele_dedx.root");
  generate_dedx_templates(fn_muo_dedx,"TrackRecoCheck/ehit_vs_path","m","pid_muo_dedx.root");

  return 0;
}

//-----------------------------------------------------------------------------
// generate the calorimeter PID templates offline branch cd3_pion_branch
// do not update dE/dX and tracker-only timing templates
//-----------------------------------------------------------------------------
int generate_pid_templates_cd3_pion() {
  
  const char* ele_track_ana_fn  = "~/hist/mu2e/cd3-pion/e00scd30.track_ana.hist";
  const char* muo_track_ana_fn  = "~/hist/mu2e/cd3-pion/m00scd30.track_ana.hist";

  // const char* ele_pid_ana_fn    = "~/hist/mu2e/cd3_pion/e00scd30.pid_ana.hist";
  // const char* muo_pid_ana_fn    = "~/hist/mu2e/cd3_pion/m00scd30.pid_ana.hist";

  // const char* fn_ele_dedx       = "~/hist/mu2e/v5_7_0/e00s5700.egun_stnmaker.hist";
  // const char* fn_muo_dedx       = "~/hist/mu2e/v5_7_0/m00s5700.mgun_stnmaker.hist";

//-----------------------------------------------------------------------------
// DT templates
//-----------------------------------------------------------------------------
  generate_dt_templates(ele_track_ana_fn,"TrackAna","trk_19/dt","pid_ele_dt.tab");
  generate_dt_templates(muo_track_ana_fn,"TrackAna","trk_19/dt","pid_muo_dt.tab");
  
//-----------------------------------------------------------------------------
// E/P templates
//-----------------------------------------------------------------------------
  generate_ep_templates(ele_track_ana_fn,"TrackAna","trk_19/ep_vs_path","pid_ele_ep_vs_path.tab");
  generate_ep_templates(muo_track_ana_fn,"TrackAna","trk_19/ep_vs_path","pid_muo_ep_vs_path.tab");
  
//-----------------------------------------------------------------------------
// X(dR/dS) templates
//-----------------------------------------------------------------------------
  // generate_xs_templates(ele_pid_ana_fn,"PidAna","pid_1/xdrds_vadim_ele","pid_ele_xdrds.tab");
  // generate_xs_templates(muo_pid_ana_fn,"PidAna","pid_1/xdrds_vadim_ele","pid_muo_xdrds.tab");
  
//-----------------------------------------------------------------------------
// de/dx templates: read 2D histograms, store 10 slices in 1mm step
//-----------------------------------------------------------------------------
  // generate_dedx_templates(fn_ele_dedx,"TrackRecoCheck/ehit_vs_path","e","pid_ele_dedx.root");
  // generate_dedx_templates(fn_muo_dedx,"TrackRecoCheck/ehit_vs_path","m","pid_muo_dedx.root");

  return 0;
}
