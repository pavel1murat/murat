//

///////////////////////////////////////////////////////////////////////////////
// read DIO data from Mu2e table file and make a histogram
///////////////////////////////////////////////////////////////////////////////
void read_dio_spectrum(const char* Filename = "ConditionsService/data/czarnecki_Al.tbl") {
  TTree* t = new TTree("t_dio","DIO spectrum");

  t->ReadFile(Filename,"e/F:w/F");

  int n = t->GetEntries();

  float e, w;

  t->SetBranchAddress("e",&e);
  t->SetBranchAddress("w",&w);

  TH1F* hist = new TH1F("h_dio","DIO spectrum",1051,-0.05,105.05);

  for (int i=0; i<n; i++) {
    t->GetEntry(i);
    int bin = (e+0.05)/0.1 + 1;
    hist->SetBinContent(bin,w);
    hist->SetBinError  (bin,0);
  }

  hist->Draw();
}

//-----------------------------------------------------------------------------
int read_lo_dio_spectrum(const char* Fn, TH1D** Hist) {
  FILE  *f;
  int    done = 0, nbx, loc(0), ix, line(0);
  char   c[1000], title[200], name[200];

  float  xmin, wt;

  if ((*Hist) != NULL) delete (*Hist);

  *Hist = new TH1D("h_dio_lo","LO DIO spctrum",1051,0,105.1);
  
  f = fopen(Fn,"r");
  if (f == 0) {
    printf("EROOR in read_lo_dio_spectrum: missing file %s\n",Fn);
    return -2;
  }

  while ( ((c[0]=getc(f)) != EOF) && !done) {
					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);

      fscanf(f," %g %g"  ,&xmin,&wt);

      int bin = xmin*10+1.1;

      printf("xmin,wt, bin: %10.3f %10.3f %12.5g %5i\n",xmin,xmin*10, wt,bin);

      (*Hist)->SetBinContent(bin,wt);
    }
					// skip the rest of the line
    fgets(c,100,f);
  }

  return 0;
}

//-----------------------------------------------------------------------------
void init_h_dio_wt(TH1D** Hist) {

  if ((*Hist) != NULL) delete (*Hist);
  *Hist = new TH1D("h_dio_wt","LO DIO spctrum",1051,0,105.1);

  double step = 0.1;
  
  for (int i=1; i<=1051; i++) {

    double x = (i-1)*step;
    double q = TStntuple::DioWeightAl(x); 
    (*Hist)->SetBinContent(i,q);
  }
}



//-----------------------------------------------------------------------------
// answer: fraction(95,105) =   3.9821e-11, approx 4.0e-11 ....
// this is 100 bins, ==> 100 uniformly distributed weighted events events
// correspond to the integral of 4e-10
// I have 6e^6, that gives the integral of 6e^6*4e-10         = 2.4e-3
// for 3 years of data taking the integral is n_dio(95,105)   = 4.7532e+06
// so each event in the weighted DIO histogram needs to be rescaled by approximately 2e+09
//-----------------------------------------------------------------------------
int dio_main() {

  TH1D* h_dio(0), *h_dio_wt(0);

  read_lo_dio_spectrum("/projects/mu2e/users/murat/dev/ConditionsService/data/czarnecki_Al.tbl",&h_dio);

  h_dio->SetLineColor(2);
  h_dio->SetFillColor(2);
  h_dio->SetFillStyle(3005);
  h_dio->Draw();
  
  // calculate integral fraction of DIO events in (95,105) MeV/c
  
  double n_95_105 = h_dio->Integral(951,1050);
  double n_total  = h_dio->Integral();

  double fraction = n_95_105/n_total;
  
  printf("fraction(95,105) = %12.5g\n",fraction);

  init_h_dio_wt(&h_dio_wt);

  h_dio_wt->Draw("same");

  return 0;
}
