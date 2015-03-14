///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
void plot_ehit_vs_path_slices(const char* fn_e = "hist/egun_stnmaker_v4_2_7.trackRecoCheck.hist",
			      const char* fn_m = "hist/mgun_stnmaker_v4_2_7.trackRecoCheck.hist",
			      int Print        = 0                                              ) {

  TFile  *fe(0), *fm(0);
					// assume taht the first file always defined
  fe = new TFile(fn_e);

  if (fn_m[0] != 0) fm = new TFile(fn_m);

  int    nx(2), ny(2);
  int    ib1, ib2, ipad, npads(nx*ny);

  TH2F   *he(0), *hm(0);

  TCanvas* c[100];

  char   name_e[100], name_m[100], name_c[100];

  he = (TH2F*) fe->Get("//TrackRecoCheck/ehit_vs_path");

  if (fm) hm = (TH2F*) fm->Get("//TrackRecoCheck/ehit_vs_path");


  TH1D   *hpe, *hpm ;

  int nc(0);

  int new_canvas = 1;

  for (int i=0; i<10; i++) {

    ib1 = 10*i+1;
    ib2 = 10*(i+1);

    sprintf(name_e,"hpe_%i",i);
    sprintf(name_m,"hpm_%i",i);

    hpe = he->ProjectionY(name_e,ib1,ib2);

    if (new_canvas) {
      sprintf(name_c,"canvas_%i",nc);
      c[nc] = new TCanvas(name_c,name_c,0,0,1200,800);
      c[nc]->Divide(nx,ny);
      new_canvas = 0;
      ipad       = 1;
      nc         += 1;
    }

    c[nc-1]->cd(ipad);

    hpe->SetLineColor(2);
    hpe->Rebin(10);
    hpe->DrawNormalized("");

    if (fm) {
      hpm = hm->ProjectionY(name_m,ib1,ib2);
      hpm->Rebin(10);
      hpm->DrawNormalized("sames");
    }

    ipad += 1;

    if (ipad > npads) {
      new_canvas = 1;
    }
    
  }

  char fn_png[100], date[100];

  TString s = gSystem->GetFromPipe("date +%Y-%m-%d");

  //  fgets(date,100,f);
  //  close(f);

  if (Print) {
    for (int i=0; i<nc; i++) {
      sprintf(fn_png,"%s-canvas_%i.png",s.Data(),i);
      c[i]->Print(fn_png);
    }
  }
}


//-----------------------------------------------------------------------------
// following Vadim: 'Particle' = 'e' or 'm'
//-----------------------------------------------------------------------------
void make_ehit_vs_path_templates(const char* fn_e = "hist/egun_stnmaker_v4_2_7.trackRecoCheck.hist",
				 const char* Particle = "e",
				 const char* fn_o = "electrontemplates_v427.hist") {

  TFile  *fe(0), *fo(0);

  fe = new TFile(fn_e);

  TH2F   *he(0), *hm(0);

  char   name_e[100], name_m[100], name_c[100];

  he = (TH2F*) fe->Get("//TrackRecoCheck/ehit_vs_path");

  TH1D   *hpe ;

  fo = new TFile(fn_o,"recreate");

  int nc(0), ib1, ib2;

  for (int i=0; i<11; i++) {
    printf(" -- i = %i\n",i);
    ib1 = 10*i+1;
    ib2 = 10*(i+1);

    sprintf(name_e,"htemp%s%i",Particle,i);

    hpe = (TH1D*) he->ProjectionY(name_e,ib1,ib2)->Rebin(10);
    hpe->Scale(1./hpe->Integral());

    //    hpe->Write();
  }

  fo->Write();
  

  fo->Close();

  delete fo;
}
