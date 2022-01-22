//
// Excised from Federica's "test_grid.C"
//

#include <TVector2.h>

// #include "federica.h"  - don't need any more

const int nXPS =  89, nXTSu  = 201, nXTSd = 205, nXDS =  97;
const int nZPS = 281, nZTSu  = 149, nZTSd = 157, nZDS = 521;
const int nY   =  97, nYhalf =  48;

const double epsilon=1e-7;

double step =    25;
double R    =  2929;         // radius of the TSu /TSd C's
double L0   =  7000;
double L1   =     0;
double L3   =  1950;         // total length of the TS3 collimator
double L5   =     0;
double L6   = 13142;

double X0PS = 2804, X0TSu = 4,   X0TSd =-5096, X0DS =-5096;
double Y0PS =-1200, Y0TSu =-1200,Y0TSd =-1200, Y0DS =-1200;
double Z0PS =-9929, Z0TSu =-2929,Z0TSd =- 829, Z0DS = 3071;

//-----------------------------------------------------------------------------
void xyz2uys(double Mu2eX, double Mu2eY, double Mu2eZ)  {

  double X1PS  = X0PS +nXPS *step;
  double X1TSu = X0TSu+nXTSu*step;
  double X1TSd = X0TSd+nXTSd*step;
  double X1DS  = X0DS +nXDS *step;

  double Y1PS=Y0PS+nY*step;
  double Y1TSu=Y0TSu+nY*step;
  double Y1TSd=Y0TSd+nY*step;
  double Y1DS=Y0DS+nY*step;

  double Z1PS=Z0PS+nZPS*step;
  double Z1TSu=Z0TSu+nZTSu*step;
  double Z1TSd=Z0TSd+nZTSd*step;
  double Z1DS=Z0DS+nZDS*step;

  double   u,y,s;
  TVector3 r0,r1,r2;
//-----------------------------------------------------------------------------
// PS
//-----------------------------------------------------------------------------
  if (Mu2eX>=X0PS && Mu2eX<=X1PS && Mu2eY>=Y0PS && Mu2eY<=Y1PS && Mu2eZ>=Z0PS && Mu2eZ<Z1PS) {
    // Determine the point on the magnetic axis closest to the input point
    // and get coordinates of the point on a system attached to the axis

    r0.SetXYZ(R+L3/2,0,Mu2eZ);
    r1.SetXYZ(Mu2eX,Mu2eY,Mu2eZ);
    Mu2eToLocal(r0,r1,r2);
    u=r2.X();
    y=r2.Y();
    s=Mu2eZ-Z0PS;
  }

  // TSu

  if (Mu2eX>=X0TSu && Mu2eX<=X1TSu && Mu2eY>=Y0TSu && Mu2eY<=Y1TSu && Mu2eZ>=Z0TSu && Mu2eZ<Z1TSu) {
    // Determine the point on the magnetic axis closest to the input point
    // and get coordinates of the point on a system attached to the axis

    if (Mu2eX>=0 && Mu2eX<=L3/2) {
      r0.SetXYZ(Mu2eX,0,0);  // TS3u
    }
    else  {                                            // TS2
      x3=Mu2eX-L3/2;
      z3=Mu2eZ+R;
      d=TMath::Sqrt(x3*x3+z3*z3);
      if(d<=0) r0.SetXYZ(L3/2,0,-R);
      else r0.SetXYZ(R*x3/d+L3/2,0,R*z3/d-R);
    }

    r1.SetXYZ(Mu2eX,Mu2eY,Mu2eZ);
    Mu2eToLocal(r0,r1,r2);
    u=r2.X();
    y=r2.Y();
    
    if (Mu2eX>=0 && Mu2eX<=L3/2) s=L0+L1+TMath::Pi()*R/2+L3/2-Mu2eX;
    else s=L0+L1+R*TMath::ASin(fabs(x3)/d);
  }
//-----------------------------------------------------------------------------
// TSd
//-----------------------------------------------------------------------------
  if (Mu2eX>=X0TSd && Mu2eX<=X1TSd && Mu2eY>=Y0TSd && Mu2eY<=Y1TSd && Mu2eZ>=Z0TSd && Mu2eZ<Z1TSd) {
    // Determine the point on the magnetic axis closest to the input point
    // and get coordinates of the point on a system attached to the axis

    if(Mu2eX>=-L3/2 && Mu2eX<=0) r0.SetXYZ(Mu2eX,0,0); // TS3d
    else  {                                                          // TS4
      x3=Mu2eX+L3/2;
      z3=Mu2eZ-R;
      d=TMath::Sqrt(x3*x3+z3*z3);
      if(d<=0) r0.SetXYZ(-L3/2,0,R);
      else r0.SetXYZ(R*x3/d-L3/2,0,R*z3/d+R);
    }

    r1.SetXYZ(Mu2eX,Mu2eY,Mu2eZ);
    Mu2eToLocal(r0,r1,r2);
    u=r2.X();
    y=r2.Y();

    if(Mu2eX>=-L3/2 && Mu2eX<=0) s=L0+L1+TMath::Pi()*R/2+L3/2-Mu2eX;
    else s=L0+L1+TMath::Pi()*R/2+L3+R*TMath::ASin(fabs(x3)/d);
  }
//-----------------------------------------------------------------------------
// DS
//-----------------------------------------------------------------------------
  if (Mu2eX>=X0DS && Mu2eX<=X1DS && Mu2eY>=Y0DS && Mu2eY<=Y1DS && Mu2eZ>=Z0DS && Mu2eZ<=Z1DS) {
    // Determine the point on the magnetic axis closest to the input point
    // and get coordinates of the point on a system attached to the axis
    
    r0.SetXYZ(-(R+L3/2),0,Mu2eZ);   // R+L3/2 = 3904
    r1.SetXYZ(Mu2eX,Mu2eY,Mu2eZ);
    Mu2eToLocal(r0,r1,r2);
    u=r2.X();
    y=r2.Y();
    s=L0+L1+L3+L5+TMath::Pi()*R+Mu2eZ-Z0DS;
  }

  cout << "u = " << u << endl;
  cout << "y = " << y << endl;
  cout << "s = " << s << endl;
}

//-----------------------------------------------------------------------------
void Mu2eToLocal(TVector3 r0, TVector3 r1, TVector3& r2) {
  double x0,y0,z0,x1,y1,z1,x2,y2,z2,theta;

  // Load the radius vectors for the on-axis point (local origin)
  // and the arbitrary space point

  x0=r0.X();
  y0=r0.Y();
  z0=r0.Z();

  x1=r1.X();
  y1=r1.Y();
  z1=r1.Z();

  // Determine where the local origin lies on the axis

  if (y0 != 0) {
    cout<<"Local system must have origin on the global XZ plane"<<endl;
    break;
  }
  
  if      (x0==R+L3/2 && z0>-(R+L0+L1) && z0<=-(R+L1)                        ) theta = TMath::Pi()/2;            // PS
  else if (x0==R+L3/2 && z0>-(R+L1) && z0<=-R                                ) theta = TMath::Pi()/2;            // TS1
  else if (TMath::Abs((x0-L3/2)*(x0-L3/2)+(z0+R)*(z0+R)-R*R)<epsilon && z0<=0) theta = TMath::ASin((x0-L3/2)/R); // TS2
  else if (x0<L3/2 && x0>=-L3/2 && z0==0                                     ) theta = 0;                        // TS3
  else if (TMath::Abs((x0+L3/2)*(x0+L3/2)+(z0-R)*(z0-R)-R*R)<epsilon && z0>0 ) theta = TMath::ACos((R-z0)/R);    // TS4
  else if (x0==-(R+L3/2) && z0>R && z0<=R+L5                                 ) theta = TMath::Pi()/2;            // TS5
  else if (x0==-(R+L3/2) && z0>R+L5 && z0<=R+L5+L6                           ) theta = TMath::Pi()/2;            // DS
  else {
    cout<<"Local system must have origin on the Mu2e magnetic axis"<<" x0="<<x0<<" z0="<<z0<<endl;
    break;
  }

  // Move to the local origin and rotate clockwise around the Y axis

  x2=(z1-z0)*TMath::Cos(theta)+(x1-x0)*TMath::Sin(theta);
  y2=y1-y0;
  z2=(z1-z0)*TMath::Sin(theta)-(x1-x0)*TMath::Cos(theta);

  // Return the radius vector in the local coordinates

  r2.SetXYZ(x2,y2,z2);
}

