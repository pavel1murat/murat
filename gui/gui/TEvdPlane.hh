//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdPlane__
#define __murat_gui_TEvdPlane__

#include "TEveElement.h"

#include "murat/gui/TEvdNumerology.hh"
#include "murat/gui/TEvdPanel.hh"

//-----------------------------------------------------------------------------
class TEvdPlane : public TEveElementList {
public:
  int         fNumber;
  TEvdPanel*  fPanel[6];

  TEvdPlane(int I = -1);

  TEvdPanel*  Panel(int I)   { return fPanel[I]; }

  int         Number() const { return fNumber; }

  ClassDef(TEvdPlane,0);
};

#endif
