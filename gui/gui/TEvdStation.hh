//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdStation__
#define __murat_gui_TEvdStation__

#include "TEveElement.h"

#include "murat/gui/TEvdNumerology.hh"
#include "murat/gui/TEvdPlane.hh"

namespace murat {
//-----------------------------------------------------------------------------
class TEvdStation : public TEveElementList {
public:
  int         fNumber;
  TEvdPlane*  fPlane[2];

  TEvdStation(int I = -1);

  TEvdPlane*  Plane(int I) { return fPlane[I]; }

  ClassDef(murat::TEvdStation,0)
};

}
#endif
