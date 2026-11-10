//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdPanel__
#define __murat_gui_TEvdPanel__

#include "murat/gui/TEvdNumerology.hh"
#include "murat/gui/TEvdStraw.hh"
#include "murat/gui/TEvdSubdetector.hh"

namespace mu2e { class Panel; }

namespace murat {

class TEvdPlane;

class TEvdPanel : public TEvdSubdetector {
public:
  int         fNumber;
  murat::TEvdPlane*  fPlane; // backward pointer to the mother plane
  murat::TEvdStraw*  fStraw[kNStraws];
  double      fNx;    // as fas as the event display is concerned, all wires are parallel
  double      fNy;
  double      fPhi;

  // panel volume: sector of a disk of radius Rout and half-length Dz, from Phi1 to Phi2 [deg],
  // minus the triangle with the origin and the two sector ends as vertices (R = Rout)
  TEvdPanel(const char* Name, double Rout, double Dz, double Phi1, double Phi2,
            TGeoMedium* Medium, int Color, int Transparency);

  // straws as daughter volumes, positioned in the panel frame (z centre ZCenter subtracted)
  // visibility is off for now
  void        AddStraws(const mu2e::Panel& Panel, double StrawInnerRadius, double StrawOuterRadius,
                        double ZCenter, TGeoMedium* Medium);
  ~TEvdPanel();

  TEvdStraw*  Straw(int I) { return fStraw[I]; }

  virtual int  InitGeometry(const char* Fn) override;
  
  void        InitStraw(int straw, int StrawID, int Plane, int Panel, int Layer,
			double rho, double z,
			double nx, double ny, double half_length);

  double      Nx () { return fNx;  }
  double      Ny () { return fNy;  }
  double      Phi() { return fPhi; }
  double      Z  () const { return (fStraw[0]->Z()+fStraw[1]->Z())/2.; }

  void        Print(Option_t* Opt = "") const override;   // *MENU* 

  ClassDefOverride(murat::TEvdPanel,0)
};
}
#endif
