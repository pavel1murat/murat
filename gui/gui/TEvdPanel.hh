//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdPanel__
#define __murat_gui_TEvdPanel__

#include "murat/gui/TEvdNumerology.hh"
#include "TEveGeoShape.h"

#include "murat/gui/TEvdStraw.hh"


namespace murat {

class TEvdPlane;

class TEvdPanel : public TEveGeoShape {
public:
  int         fNumber;
  murat::TEvdPlane*  fPlane; // backward pointer to the mother plane
  murat::TEvdStraw*  fStraw[kNStraws];
  double      fNx;    // as fas as the event display is concerned, all wires are parallel
  double      fNy;
  double      fPhi;

  TEvdPanel(int I = -1);
  ~TEvdPanel();

  TEvdStraw*  Straw(int I) { return fStraw[I]; }

  void        InitGeometry();
  void        InitStraw(int straw, int StrawID, int Plane, int Panel, int Layer,
			double rho, double z,
			double nx, double ny, double half_length);

  double      Nx () { return fNx;  }
  double      Ny () { return fNy;  }
  double      Phi() { return fPhi; }
  double      Z  () const { return (fStraw[0]->Z()+fStraw[1]->Z())/2.; }

  void        Print(Option_t* Opt) const;   // *MENU* 

  ClassDef(murat::TEvdPanel,0)
};
}
#endif
