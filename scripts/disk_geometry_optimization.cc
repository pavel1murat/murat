///////////////////////////////////////////////////////////////////////////////
//
// .L disk_geometry_optimization.cc+

// scan the geometry (staggered):
//
// scan_rmin(340,380,660,33,0.2,1,0)
// scan_rmax(360,640,680,34,0.2,1,0)
//
// plot  the best layout:
//
// crystals(360,660,34.,0.160,1,1)
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TH2.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TEllipse.h"
//-----------------------------------------------------------------------------
// Ic runs  from 0
//-----------------------------------------------------------------------------
void get_indices(int Ic, int Layer, int&Ix, int& Iy) {


  int delta;

  int nside = 2*Layer+1;

  if (Ic <= nside-1) {
    Ix = -Layer+Ic;
    Iy = Layer;
  }
  else if (Ic < 2*nside-2) {
    delta = Ic-(nside-1);
    Ix = Layer;
    Iy = Layer-delta;
  }
  else if (Ic < 3*nside-3) {
    delta = Ic-(2*nside-2);
    Ix = Layer-delta;
    Iy = -Layer;
  }
  else if (Ic < 4*nside - 4) {
    delta = Ic-(3*nside-3);
    Ix = -Layer;
    Iy = -Layer+delta;
  }
  else {
    printf(" Ic, Layer, nside: %5i %5i %5i : we're in trouble\n",
	   Ic,Layer,nside);
  }
}

//-----------------------------------------------------------------------------
int crystal_is_inside(double X0, double Y0, double Side, double Rin, double Rout) {

  int n0(0), n1(0), n2(0);

  double x[4], y[4], r;


  x[0] = X0-Side/2;
  y[0] = Y0+Side/2;

  x[1] = X0+Side/2;
  y[1] = Y0+Side/2;

  x[2] = X0+Side/2;
  y[2] = Y0-Side/2;

  x[3] = X0-Side/2;
  y[3] = Y0-Side/2;


  for (int i=0; i<4; i++) {
    r = sqrt(x[i]*x[i] + y[i]*y[i]);

    if      (r < Rin ) n0++;
    else if (r < Rout) n1++;
    else               n2++;
  }

  if      (n0 > 0) {
    if (n0 < 4) return 1;  // partially inside the inner ring
    else        return 2;  // fully inside the inner ring
  }

  else if (n2 > 0) {
    if (n2 < 4) return 21;  // partially outside the outer ring
    else        return 22;  // fully outside teh outer ring
  }
  else {
    // all 4 vertices are inside the annulus, check the center

    r = sqrt(X0*X0+Y0*Y0);
    if ((r > Rin) && (r < Rout)) return 11;  // fully inside the annulus
    else return -1;
  }
}



//-----------------------------------------------------------------------------
int crystals(double RMin = 340, double RMax = 660, double CrystalSize = 34., double Wrap=0.160, int Staggered = 1, int Plot=0) {

  //   double   size = 34.; 
  //   double   wrap = 0.16; 

//   double   size = 10.; 
//   double   wrap = 0.; 

  double   side = CrystalSize+2*Wrap;

  int    n_inside, n_below, n_above;  // for each layer or "ring"
  int    ncr ;             // n(crystals) in a ring
  int    flag;

  double    x0, y0;

  int nside;  // number of crystals on a side

  int ncr_disk = 0;
  
  char   title[200];

  TCanvas* c;
  TH2F* h2;

  if (Plot != 0) {
    c = new TCanvas("c","c",1000,1000);

    h2 = new TH2F("h2","h2",280,-700,700,280,-700,700);
    h2->SetStats(0);

    sprintf(title," RMin=%5.0f RMax=%5.0f Crystal Size=%5.0f wrap=%5.3f",
	  RMin,RMax,CrystalSize,Wrap);

    h2->SetTitle(title);
    h2->Draw();
  }

  for (int il=0; ; il++) {

    //    printf( " ring : %5i\n", il);

    n_inside = 0;
    n_below  = 0;
    n_above  = 0;

    if (il == 0) ncr = 1;
    else         ncr = 8*il;

    nside = 2*il+1;

    for (int ic=0; ic<ncr; ic++) {
//-----------------------------------------------------------------------------
// check if a crystal is fully inside the annulus
//-----------------------------------------------------------------------------
      int ix, iy;

      get_indices(ic,il,ix,iy);

      x0 = side*ix;
      if (Staggered != 0) {
	if (iy%2 != 0) x0 += 0.;      // side/4;
	else           x0 -= side/2.; // side/4; 
      }
      y0 = side*iy;

      flag = crystal_is_inside(x0,y0,side,RMin,RMax);

      if (flag < 10) {
	n_below += 1;
      }
      else if (flag == 11) {
	// printf("crystal inside : ix, iy = %5i %5i\n", ix,iy);
	n_inside += 1;
      }
      else if (flag > 20) {   // outside
	n_above += 1;
      }


      if (Plot != 0) {
	if ((flag == 11) || (flag == 21) || (flag == 1)) {
	  // draw the crystal
	  
	  int color, style;
	  
	  if      (flag == 11) {
	    color = kGreen+3;
	    style = 1;
	  }
	  else if (flag == 21) {
	    color = 2;
	    style = 3001;
	  }
	  else if (flag ==  1) {
	    color = 2;
	    style = 3001;
	  }
	  
	  TBox* b = new TBox(x0-side/2,y0-side/2,x0+side/2,y0+side/2);
	  
	  b->SetFillColor(color);
	  b->SetFillStyle(color);
	  b->SetLineColor(color);
	  b->SetLineWidth(1);
	  b->SetFillStyle(1);
	  
	  b->Draw();
	  
	  TBox* b1 = new TBox(x0-side/2,y0-side/2,x0+side/2,y0+side/2);
	  
	  b1->SetFillColor(color);
	  b1->SetFillStyle(color);
	  b1->SetLineColor(color);
	  b1->SetLineWidth(1);
	  b1->SetFillStyle(style);
	  
	  b1->Draw();
	}
      }
    }

    ncr_disk += n_inside;

    if (ncr == n_above) break;
  }

  //   printf ("ncr_disk = %5i\n",ncr_disk);


  if (Plot != 0) {

    TEllipse* e1 = new TEllipse(0.,0.,RMin,RMin,0.,360.,0);
    e1->SetLineColor(kGreen);
    e1->SetFillStyle(0);

    e1->Draw();
    
    TEllipse* e2 = new TEllipse(0.,0.,RMax,RMax,0.,360.,0);
    e2->SetLineColor(kGreen);
    e2->SetFillStyle(0);
    
    e2->Draw();
  }
    
  double sann = M_PI*(RMax*RMax-RMin*RMin);

  double scr  = ncr_disk*side*side;

  double eff = scr/sann;

  printf ("CrystalSize, rmin, rmax, ncr_disk, sann, scr, 1-eff = %5.1f %5.1f %5.1f %5i %10.5f %10.5f %10.5f\n",
	  CrystalSize,RMin, RMax,ncr_disk,sann/1e6,scr/1e6,1-eff);

  return 0;
}


//-----------------------------------------------------------------------------
void scan_rmin(double RMin1, double RMin2, double RMax, double CrystalSize, double Wrap, int Staggered) {

  for (double r=RMin1; r<=RMin2; r+=1.) {
    crystals(r, RMax,CrystalSize,Wrap,Staggered,0);
  }
}

//-----------------------------------------------------------------------------
void scan_rmax(double RMin, double RMax1, double RMax2, double CrystalSize, double Wrap, int Staggered) {

  for (double r=RMax1; r<=RMax2; r+=1.) {
    crystals(RMin, r,CrystalSize,Wrap,Staggered,0);
  }
}

//-----------------------------------------------------------------------------
void scan_crystal_size(double RMin, double RMax, double Wrap, int Staggered) {

  for (double size=30.; size<35; size+=0.1) {
    crystals(RMin,RMax,size,Wrap,Staggered,0);
  }
}
