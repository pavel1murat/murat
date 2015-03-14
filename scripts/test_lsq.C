//-----------------------------------------------------------------------------
//  x_0 = 250
//  y_0 =  70
//  r   = 250
//-----------------------------------------------------------------------------
void test_lsq(const char* Fn, double offset = 0.) {


  float  xx[1000], yy[1000], zz[1000], phi[1000];
  int    flags[1000];
  int    loc;
  double x_0, y_0, r, x, y, chi2;

  LsqSums4 s;

  double data1 [] = {
    10., 0., 1.,
    0., 10., 0., 
    -10., 0., 0.,
    0., -10., 0.,
    -1.
  };

  double data2 [] = {
    250.000,      0.000,    504.876,
    1.76776695296636859e+02, 1.76776695296636859e+02, 504.,
      0.000,    250.000,    504.876,
   -1.76776695296636859e+02, 1.76776695296636859e+02, 504.,
   -250.000,      0.000,    504.876,
    -1.
  };

  double data [] = {
    500.000,     70.000,    504.876,
    499.688,     82.495,    506.451,
    498.751,     94.958,    507.710,
    497.193,    107.360,    508.652,
    495.017,    119.667,    509.276,
    492.228,    131.851,    509.581,
    488.834,    143.880,    509.569,
    484.843,    155.724,    509.238,
    480.265,    167.355,    508.589,
    475.112,    178.741,    507.622,
    469.396,    189.856,    506.338,
    463.131,    200.672,    504.737,
    456.334,    211.161,    502.821,
    449.021,    221.297,    500.592,
    441.211,    231.054,    498.049,
    432.922,    240.410,    495.195,
    424.177,    249.339,    492.032,
    414.996,    257.820,    488.562,
    405.402,    265.832,    484.786,
    395.421,    273.354,    480.708,
    385.076,    280.368,    476.329,
    374.393,    286.856,    471.653,
    363.399,    292.802,    466.682,
    352.122,    298.191,    461.419,
    340.589,    303.010,    455.868,
    328.831,    307.246,    450.033,
    316.875,    310.890,    443.917,
    304.752,    313.931,    437.523,
    292.492,    316.362,    430.856,
    280.126,    318.178,    423.920,
    267.684,    319.374,    416.719,
    255.199,    319.946,    409.258,
    242.700,    319.893,    401.541,
    -1
  };

  char c[1000];
  int done(0), np(0), ngood(0);
//-----------------------------------------------------------------------------
// if file is undefined, take inputs from the data
//-----------------------------------------------------------------------------
  if (Fn == 0) {
    for (int i=0; data[3*i]!=-1; i++) {
      xx[np] = data[3*i  ];
      yy[np] = data[3*i+1] + offset;
      printf(" test_lsq: x = %10.3f y = %10.3f \n",x,y);
      np++;
    }
  }
  else {
//-----------------------------------------------------------------------------
// file is specified, read it
//-----------------------------------------------------------------------------
    FILE* f  = fopen(Fn,"r");
  
    if (f == 0) {
      Error("[test_lsq.C]",Form("missing file %s\n",Fn));
      return;
    }

    while ( ((c[0]=getc(f)) != EOF) && !done) {
					// check if it is a comment line
      if (c[0] != '#') {
	ungetc(c[0],f);
	// read channel number
	fscanf(f,"%f"  ,&xx [np]);
	fscanf(f,"%f"  ,&yy [np]);
	fscanf(f,"%f"  ,&phi[np]);

	printf("[test_lsq]: np=%3i flags=%08x x[np], y[np], phi[np] = %10.3f %10.3f %10.3f \n",
	       np,flags[np],xx[np],yy[np],phi[np]);

	if (flags[np] < 256) ngood++;
	np++;
      }
					// skip line
      fgets(c,100,f);
    }

    fclose(f);
  }
//-----------------------------------------------------------------------------
// calculate averages and o
//-----------------------------------------------------------------------------
  double xm(0), ym(0);
  for (int i=0; i<np; i++) {
    xm += xx[i];
    ym += yy[i];
  }

  xm = xm/np;
  ym = ym/np;

  s.fX0 = xm;
  s.fY0 = ym;

  for (int i=0; i<np; i++) {
    s.addPoint(xx[i],yy[i]);
  }
//-----------------------------------------------------------------------------
// retrieve the answer
//-----------------------------------------------------------------------------
  x_0 = s.x0();
  y_0 = s.y0();
  r   = s.radius();

  chi2 = s.chi2DofCircle();

  double sum(0), sum_dr2(0), sum_01(0), dr2, dx;

  for (int i=0; data[3*i]>=0; i++) {
    x = data[3*i  ];
    y = data[3*i+1] + offset;
    dr2 = (x-x_0)*(x-x_0) + (y-y_0)*(y-y_0) - r*r;

    sum += dr2*dr2/(4*r*r);

    sum_dr2 += dr2;

    dx = (x*x-s.x2Mean())+(y*y-s.y2Mean())-2*x_0*(x-s.xMean())-2*y_0*(y-s.yMean());

    sum_01 += dx*dx;
  }

  double uu = s.sigX2X2()+2*s.sigX2Y2()+s.sigY2Y2();

  printf (" test_lsq: X_0 = %10.3f Y_0 = %10.3f  R = %10.3f  chi2 = %12.5e sum = %12.5e sum_dr2 = %12.5e\n",
	  x_0, y_0, r, chi2, sum/s.qn(),sum_dr2/s.qn());

  printf(" test_lsq: sum_01 = %12.5e\n",sum_01);
  printf(" test_lsq: uu     = %12.5e\n",uu);

  printf(" test_lsq: xmean  = %12.5e\n",s.xMean());
  printf(" test_lsq: ymean  = %12.5e\n",s.yMean());
  printf(" test_lsq: x2mean = %12.5e\n",s.x2Mean());
  printf(" test_lsq: y2mean = %12.5e\n",s.y2Mean());

  printf(" test_lsq: sigxx  = %12.5e\n",s.sigXX());
  printf(" test_lsq: sigxy  = %12.5e\n",s.sigXY());
  printf(" test_lsq: sigyy  = %12.5e\n",s.sigYY());

  printf(" test_lsq: sigx2x  = %12.5e\n",s.sigX2X());
  printf(" test_lsq: sigx2y  = %12.5e\n",s.sigX2Y());
  printf(" test_lsq: sigxy2  = %12.5e\n",s.sigXY2());
  printf(" test_lsq: sigyy2  = %12.5e\n",s.sigYY2());

  printf(" test_lsq: sigxx3  = %12.5e\n",s.sigXX3());
  printf(" test_lsq: sigx2x2 = %12.5e\n",s.sigX2X2());
  printf(" test_lsq: sigx3y  = %12.5e\n",s.sigX3Y());
  printf(" test_lsq: sigx2y2 = %12.5e\n",s.sigX2Y2());
  printf(" test_lsq: sigxy3  = %12.5e\n",s.sigXY3());
  printf(" test_lsq: sigy2y2 = %12.5e\n",s.sigY2Y2());
  printf(" test_lsq: sigyy3  = %12.5e\n",s.sigYY3());

  TH2F* h2 = new TH2F("h2","h2",120,0,600,120,0,600);

  h2->Draw();


  TGraph* gr = new TGraph(np,xx,yy);
  gr->SetMarkerStyle(2);
  gr->SetMarkerSize(0.7);

  gr->Draw("P");

  TEllipse* e = new TEllipse(x_0,y_0,r,r);
  e->SetFillStyle(0);

  e->Draw();
}
