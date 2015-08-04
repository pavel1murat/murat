//

namespace {
  double emu(105.194), mAl(25133.), mmu(105.6583715) ;
  double z(13.), alpha(1./137.036);
};


double dio_weight(double E) {

  double a5(8.6434), a6(1.16874), a7(-1.87828e-2), a8(9.16327e-3);
  double de, de5, w;

  double eb = mmu*(z*alpha)*(z*alpha)/2.;
  double er = E*E/(2*mAl);

  //  de  = mmu-eb-E-E*E/(2*mAl);
  de  = emu-er-E;
  //  de  = emu-E;

  de5 = de*de*de*de*de;

  w   = 1.e-17*de5*(a5 + de*(a6+de*(a7+a8*de)));

  if (de < 0) w = 0;

  return w;
}

//-----------------------------------------------------------------------------
// it is not clear whether emu = mmu*(1-(Z*alpha)**2/2) is calculated correctly
// seems that it should be 
//-----------------------------------------------------------------------------
double dio_weight_rc(double E) {

  double de, de5, w, xe, xe5;

  //  de  = emu-E-E*E/(2*mAl);

  double eb = mmu*(z*alpha)*(z*alpha)/2.;

  double er = E*E/(2*mAl);

  // printf("eb = %12.5f\n",eb);
  // printf("er = %12.5f\n",er);
  
  //  de  = mmu-eb-er-E;
  de  = emu-er-E;
  //  de  = emu-E;

  xe  = de/mmu;

  xe5 = xe*xe*xe*xe*xe;

  w   = (1.44*pow(xe,0.023)-0.22)*1e-4*xe5/mmu;

  if (de < 0) w = 0;

  return w;

}


//-----------------------------------------------------------------------------
double fun2(double* X, double* Par) {
  return dio_weight_rc(X[0]);
}

//-----------------------------------------------------------------------------
double fun1(double* X, double* Par) {
  return dio_weight(X[0]);
}

//-----------------------------------------------------------------------------
void plot_dio_weight_rc() {

  int nb = 1000;

  TH1F* h1 = new TH1F("h1","h1",nb,95,105);
  TH1F* h2 = new TH1F("h2","h2",nb,95,105);

  double bin = h1->GetBinWidth(1);
  
  for (int i=1; i<=nb; i++) {
    x = 95.+(i-0.5)*bin;

    double w    = dio_weight   (x);
    double w_rc = dio_weight_rc(x);


    h1->SetBinContent(i,w);
    h2->SetBinContent(i,w_rc);
  }


  TCanvas* c = new TCanvas("c","c",0,0,1000,800);

  h1->SetLineColor(2);
  h1->Draw();
  h2->Draw("same");


  TF1* f1 = new TF1("f1",fun1,95,105,0);
  TF1* f2 = new TF1("f2",fun2,95,105,0);

  double f1int = f1->Integral(104.0,105);
  double f2int = f2->Integral(104.0,105);

  printf("f1->integral(104,105) = %12.5e\n",f1int);
  printf("f2->integral(104,105) = %12.5e\n",f2int);

  printf("ratio = %12.5f\n",f1int/f2int);
}




//-----------------------------------------------------------------------------
// read the DIO parameterization for Al from Kyle's table 
// ConditionsService/data/czarnecki_Al.tbl
// and compare my DIO parameterization with the tabulated data
//-----------------------------------------------------------------------------
void compare_dio() {

  float x[10000], y[10000];

  float xmin, xmax;

  char   c[1000];
  int    done = 0;

  FILE* f = fopen("ConditionsService/data/czarnecki_Al.tbl","r");

  fscanf(f,"%f",&xmax);
  fscanf(f,"%f",&xmin);
  
  printf (" xmin = %12.3f xmax = %12.3f\n",xmin,xmax); 

  int np = 0;
  
  while ( ((c[0]=getc(f)) != EOF) && !done) {
    ungetc(c[0],f);

    fscanf(f,"%f" ,&x[np]   );
    fscanf(f,"%e" ,&y[np]   );

    np++;

    if (x[np-1] == xmin) done = 1;

    fgets(c,200,f);
  }

  printf (" npt = %i\n",np);


  float bin = 0.1;
  TH1F* h = new TH1F("h","h",1050,0-bin/2,105.-bin/2);

  for (int i=0; i<np; i++) {

    h->SetBinContent(1050-i,y[i]);
  }


  h->GetXaxis()->SetRangeUser(95,105);
  h->Draw();
  
  TF1* f1 = new TF1("f1",fun1,95,105,0);
  f1->SetNpx(1000);

  f1->Draw("same");
}
