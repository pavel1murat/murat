///////////////////////////////////////////////////////////////////////////////
// Len - length code: 12 = 12 cm etc
// fit evt_1/emax distributions - max cluster energy for events with track
// or  evt_6/emax distributions - max cluster energy for events with "Set C" track
///////////////////////////////////////////////////////////////////////////////

TPolFitSpectrum* pf;


double dat_12 [12] = {
  1, 70, 92 , 
  1, 91, 98 , 
  1, 97,101 ,
  -1, 0 , 0
};

double dat_14 [12] = {
  1, 70, 92 , 
  1, 91, 98 , 
  1, 97,101 ,
  -1, 0 , 0
};

double dat_16 [12] = {
  1, 70, 92 , 
  1, 91, 98 , 
  1, 97,101 ,
  -1, 0 , 0
};

double dat_18 [12] = {
  1, 70, 92 , 
  1, 91, 99 , 
  1, 98,103 ,
  -1, 0 , 0
};

double dat_20 [12] = {
  1, 70, 92 , 
  1, 91, 98 , 
  1, 97, 101 ,
  -1, 0 , 0
};

double dat_22 [12] = {
  1, 70,  92 , 
  1, 91,  99 , 
  1, 98, 102 ,
  -1, 0 , 0
};


const char* fn_12 = "hist/scan_crystal_length/e20s5600_e12s5600.track_ana.hist";
const char* fn_14 = "hist/scan_crystal_length/e20s5600_e14s5600.track_ana.hist";
const char* fn_16 = "hist/scan_crystal_length/e20s5600_e16s5600.track_ana.hist";
const char* fn_18 = "hist/scan_crystal_length/e20s5600_e18s5600.track_ana.hist";
const char* fn_20 = "hist/scan_crystal_length/e20s5600_e20s5600.track_ana.hist";
const char* fn_22 = "hist/scan_crystal_length/e20s5600_e22s5600.track_ana.hist";

//-----------------------------------------------------------------------------
// CINT doesn't seem to handle correctly passing T*& , so use T**
//-----------------------------------------------------------------------------
void init(int Len, TPolFitSpectrum::Range_t** R, TH1F*& Hist) {
  double*                   dat;
  char*                     fn ;

  if (Len == 12) {
    dat = dat_12;
    fn  = fn_12;
  }
  if (Len == 14) { 
    dat = dat_14;
    fn  = fn_14;
  }
  if (Len == 16) {
    dat = dat_16;
    fn  = fn_16;
  }
  if (Len == 18) {
    dat = dat_18;
    fn  = fn_18;
  }
  if (Len == 20) {
    dat = dat_20;
    fn  = fn_20;
  }
  if (Len == 22) {
    dat = dat_22;
    fn  = fn_22;
  }

  int n = 0; while (dat[3*n] > 0) n++;

  printf(" Len = %5i n = %5i\n",Len,n);

  (*R) = new TPolFitSpectrum::Range_t[n+1];

  TPolFitSpectrum::Range_t* r;
  for (int i=0; i<n; i++) {
    r = (*R)+i;
    printf(" r = 0x%08x \n",r);
    r->fMode = dat[3*i  ];
    r->fXMin = dat[3*i+1];
    r->fXMax = dat[3*i+2];
  }

  r = (*R)+n;
  printf(" last: r = 0x%08x \n",r);

  r->fMode = -1;

  //  Hist = gh1(fn,"TrackAna","trk_1/ep");
  //  Hist = gh1(fn,"TrackAna","evt_1/emax");
  Hist = gh1(fn,"TrackAna","evt_6/emax");

  Hist->Draw();


}


//-----------------------------------------------------------------------------
void test_pol_cluster_e(int Len) {

  double                    x, fwhm; 
  TH1F*                     h;
  TPolFitSpectrum::Range_t* r;

  init(Len,&r,h);

  for (int i=0; ; i++) {
    printf(" i, loc, fMode, fXMin, fXMax = %3i 0x%08x %3i %10.3f %10.3f \n",i, &r[i], r[i].fMode,r[i].fXMin,r[i].fXMax);
    if ((r+i)->fMode <= 0) break;
  }

  pf = new TPolFitSpectrum("test",h,r);

  pf->main_pol();

  pf->GetFWHM(&x,&fwhm);

  printf("x = %10.3f fwhm = %10.3f\n",x,fwhm);
}
