//

TCanvas* c_plot_hits_xy(0);
TCanvas* c_plot_hits_zphi(0);

void muhit_display_plot_hits(const char* Fn) {

  float  x[1000], y[1000], rho[1000], z[1000];
  int    flags[1000], used, done(0);

  int    parameters_read(0);

  float   X0, Y0, R, chi2;

  int np(0), ngood(0);

  char c[1000], name[100], s[100];

  FILE* f  = fopen(Fn,"r");

  if (f == 0) {
    Error("plot_hits",Form("missing file %s\n",Fn));
    return -2;
  }

  while ( ((c[0]=getc(f)) != EOF) && !done) {

					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);
//-----------------------------------------------------------------------------
// first line - parameters
//-----------------------------------------------------------------------------
      if (parameters_read == 0) {
	parameters_read = 1;
	fscanf(f,"%s"  ,name  );

	fscanf(f,"%s ="  ,name  );
	//	printf(" ......... name =%s\n",name);
	fscanf(f,"%f",&X0);

	fscanf(f,"%s =",name  );
	//	printf(" ......... name =%s\n",name);
	fscanf(f,"%f",&Y0);

	fscanf(f,"%s =",name);
	//	printf(" ......... name =%s\n",name);
	fscanf(f,"%f",&R);

	fscanf(f,"%s =",name);
	//	printf(" ......... name =%s\n",name);
	fscanf(f,"%f",&chi2);

// 	printf("X0 = %12.5f Y0 = %12.5f R = %12.5f  chi2 = %12.5e \n",
// 	       X0, Y0, R, chi2);

      }
      else {
					// read points
	fscanf(f,"%s"  ,name  );
	fscanf(f,"%08x",&flags[np]);
	fscanf(f,"%i"  ,&used );
	fscanf(f,"%f"  ,&x[np]);
	fscanf(f,"%f"  ,&y[np]);
	fscanf(f,"%f"  ,&z[np]);
	
// 	printf("name = %s np=%3i flags=%08x %4i x[np], y[np], z[np] = %10.3f %10.3f %10.3f \n",
// 	       name,np,flags[np], flags[np],x[np],y[np],z[np]);

	if (flags[np] < 256) ngood++;
	np++;
      }

    }
					// skip line
    fgets(c,100,f);
  }

  fclose(f);

  printf(">> np = %i, ngood = %3i\n",np,ngood);

  TGraph* gr_xy = new TGraph(np,x,y);
  TGraph* gr_yz = new TGraph(np,z,y);

  if (c_plot_hits_xy != 0) delete c_plot_hits_xy;

  c_plot_hits_xy = new TCanvas(Form("c_plot_hits_xy_%s",name),Form("c_xy %s",name),1200,600);

  c_plot_hits_xy->Divide(2,1);

  c_plot_hits_xy->cd(1);

  TH2F* h2_xy = new TH2F("h2_xy","XY View",140,-700,700,140,-700,700);
  h2_xy->SetStats(0);
  h2_xy->Draw();

  TMarker* m;

  int color;

  for (int i=0; i<np; i++) {
    m = new TMarker(x[i],y[i],2);
    if (flags[i] >= 256) color = kBlack;
    else                color = kRed;
    m->SetMarkerSize(0.7);
    m->SetMarkerColor(color);
    m->Draw();
  }

  TEllipse* e = new TEllipse(X0,Y0,R);

  e->SetFillStyle(0);
  e->Draw();


  c_plot_hits_xy->cd(2);

  TH2F* h2_yz = new TH2F("h2_yz","YZ VIEW",1600,-1600,1600,140,-700,700);
  h2_yz->SetStats(0);
  h2_yz->Draw();


  for (int i=0; i<np; i++) {
    m = new TMarker(z[i],y[i],2);
    if (flags[i] >= 256) color = kBlack;
    else                color = kRed;
    m->SetMarkerColor(color);
    m->SetMarkerSize(0.7);
    m->Draw();
  }


  c_plot_hits_zphi = new TCanvas(Form("c_plot_hits_zphi_%s",name),Form("c_xy %s",name),1200,600);

  //  c_plot_hits_zphi->Divide(2,1);

  c_plot_hits_zphi->cd(1);

  TH2F* h2_zphi = new TH2F("h2_zphi","ZPHI View",1600,-1600,1600,160,-3.2,3.2);
  h2_zphi->SetStats(0);
  h2_zphi->Draw();

  TMarker* m;

  int color;

  for (int i=0; i<np; i++) {
    double phi = TMath::ATan2(y[i],x[i]);
    m = new TMarker(z[i],phi,2);
    if (flags[i] >= 256) color = kBlack;
    else                color = kRed;
    m->SetMarkerSize(0.7);
    m->SetMarkerColor(color);
    m->Draw();
  }

}
