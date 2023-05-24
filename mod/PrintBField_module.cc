///////////////////////////////////////////////////////////////////////////////
// print B-field on the axis of the Mu2e magnetic channel
///////////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"

#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <cstdio>
#include <iostream>
#include <vector>

namespace mu2e {

  class PrintBField : public art::EDAnalyzer {

    int    _npoints;
    double _step ; // mm

    double _x [100000],  _y[100000],  _z[100000];
    double _bx[100000], _by[100000], _bz[100000];
    double _s [100000];
    
  public:
    explicit PrintBField(const fhicl::ParameterSet& pset) : art::EDAnalyzer(pset) {
      _npoints = 0;
      _step = 25.;
    }
    
    void beginRun(const art::Run& run);
    void analyze(const art::Event&){};

    void  define_points();
  };
  
//-----------------------------------------------------------------------------
  void PrintBField::beginRun(const art::Run& run) {
    GeomHandle<BFieldManager> bfmgr;

    define_points();

    for (int i=0; i<_npoints; i++) {
      CLHEP::Hep3Vector field = bfmgr->getBField(CLHEP::Hep3Vector(_x[i], _y[i], _z[i]));
      _bx[i] = field.x();
      _by[i] = field.y();
      _bz[i] = field.z();
      _s[i]  = _step*i;

      printf("%5i %9.1f %9.1f %9.1f %9.1f %10.3f %10.3f %10.3f\n",
             i,_s[i],_x[i],_y[i],_z[i],_bx[i],_by[i],_bz[i]);
    }
  }

//-----------------------------------------------------------------------------
void PrintBField::define_points() {

  double R    =  2929;         // radius of the TSu /TSd C's
  //  double L0   =  7000;
  // double L1   =     0;
  double L3   =  1950;         // total length of the TS3 collimator
  // double L5   =     0;
  double L6   = 13142;

  //  double X0PS = 2804, X0TSu = 4,   X0TSd =-5096, X0DS =-5096;
  //  double Y0PS =-1200, Y0TSu =-1200,Y0TSd =-1200, Y0DS =-1200;
  double Z0PS =-9929; // , Z0TSu =-2929,Z0TSd =- 829, Z0DS = 3071;
//-----------------------------------------------------------------------------
// calculate coordinates on the magnet axis
//-----------------------------------------------------------------------------
  int off(0), loc(0);
  //  double s(0);

  int nps = 281;

  for (int i=0; i<nps; i++) {
    // s      = _step*i;
    loc    = off+i;
    _x[loc] = (R+L3/2);
    _y[loc] = 0;
    _z[loc] = Z0PS+i*_step;
  }
//-----------------------------------------------------------------------------
// first point out of PS
//-----------------------------------------------------------------------------
  double x0tsu = L3/2 + R;
  double z0tsu = Z0PS+(nps-1)*_step;

  off       = nps;
  double len = R*M_PI/2;
  int ntsu  = len/_step + 1;

  for (int i=0; i<ntsu; i++) {
    double phi = i*_step/len*(M_PI/2);

    loc    = off+i;
    _x[loc] = x0tsu-R*(1-cos(phi));
    _y[loc] = 0;
    _z[loc] = z0tsu+R*sin(phi);
  }
//-----------------------------------------------------------------------------
// TS3 collimator - straight X-section through the TS3 collimator
//-----------------------------------------------------------------------------
  double x0ts3 = L3/2;

  off      = off+ntsu;
  int nts3 = L3/_step + 1;
  for (int i=0; i<nts3; i++) {
    loc = off+i;

    _x[loc] = x0ts3-i*_step;
    _y[loc] = 0;
    _z[loc] = 0;
  }
//-----------------------------------------------------------------------------
// TSd
//-----------------------------------------------------------------------------
  double x0tsd = -L3/2;
  double z0tsd = 0;

  off      = off+nts3;

  int ntsd = len/_step + 1;
  for (int i=0; i<ntsd; i++) {
    loc = off+i;

    double phi = i*_step/len*(M_PI/2);

    _x[loc] = x0tsd-R*sin(phi);
    _y[loc] = 0;
    _z[loc] = z0tsd+R*(1-cos(phi));
  }
//-----------------------------------------------------------------------------
// finally, straight line till the end
//-----------------------------------------------------------------------------
  double z0ds = _z[loc];
  double x0ds = _x[loc];

  off = off+ntsd;

  for (int i=0; 1 ; i++) {
    loc = i+off;
    
    _x[loc] = x0ds;
    _y[loc] = 0;
    _z[loc] = z0ds+i*_step;
    if (_z[loc] > L6) break;
  }

  _npoints = loc+1;
  printf ("npoints = %i\n",_npoints);

}

}  // namespace mu2e

DEFINE_ART_MODULE(mu2e::PrintBField);
