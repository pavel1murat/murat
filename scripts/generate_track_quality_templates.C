///////////////////////////////////////////////////////////////////////////////
// 10 numbers per line are assumed ?
///////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
void generate_1D_templates(const char* HistFile, const char* Module, const char* HistName, const char* Fout) {

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
void generate_2D_templates(const char* HistFile, const char* Module, const char* HistName, const char* Fout) {

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
// generate track quality templates for TrkPatRec and CalPatRec
// offline version v5_7_7
//-----------------------------------------------------------------------------
int generate_track_quality_templates() {
  
  const char* fn  = "~/hist/mu2e/v5_7_0/e11s5731.track_comp_use_mva_001.hist";
//-----------------------------------------------------------------------------
// track quality templates:
// ------------------------
// default: linear weight training
// TrackComp/trk_100 : TrkPatRec , all reconstructed tracks
// TrackComp/trk_200 : CalPatRec , all reconstructed tracks
//-----------------------------------------------------------------------------
  generate_1D_templates(fn,"TrackComp","trk_100/mvaout","tpr_qual.tab");
  generate_1D_templates(fn,"TrackComp","trk_200/mvaout","cpr_qual.tab");
}
