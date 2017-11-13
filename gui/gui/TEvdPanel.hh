//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdPanel__
#define __murat_gui_TEvdPanel__

#include "murat/gui/TEvdNumerology.hh"
#include "TEveGeoShape.h"

class TEvdStraw;

class TEvdPanel : public TEveGeoShape {
public:
  int         fNumber;
  TEvdStraw*  fStraw[kNStraws];
  double      fNx;  // as fas as event display is concerned, all wires are parallel
  double      fNy;

  TEvdPanel(int I = -1);
  ~TEvdPanel();

  TEvdStraw*  Straw(int I) { return fStraw[I]; }

  void        InitGeometry();
  void        InitStraw(int straw, int StrawIndex, double rho, double z,
			double nx, double ny, double half_length);

  ClassDef(TEvdPanel,0)
};

#endif
