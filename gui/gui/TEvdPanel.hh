//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdPanel__
#define __murat_gui_TEvdPanel__

#include "murat/gui/TEvdNumerology.hh"
#include "murat/gui/TEvdStraw.hh"
#include "murat/gui/TEvdSubdetector.hh"

namespace mu2e { class Panel; }

namespace murat {

class TEvdPlane;
class TEvdTracker;

// geometry of a panel volume: sector from Phi1 to Phi2 [deg] of a disk of radius Rout and
// half-length Dz, centred at z = Zc in the tracker frame
struct PanelGeometry {
  double rin;      // inner radius: the innermost straw's inner edge
  double rout;
  double dz;
  double phi1;
  double phi2;
  double zc;
};

// geometry of the tracker panel: the panel volume and its straws' z centre
PanelGeometry panelGeometry(const TEvdTracker* Tracker, int Plane, int Panel);

class TEvdPanel : public TEvdSubdetector {
public:
  int         fNumber;
  murat::TEvdPlane*  fPlane; // backward pointer to the mother plane
  murat::TEvdStraw*  fStraw[kNStraws];
  double      fNx;    // as fas as the event display is concerned, all wires are parallel
  double      fNy;
  double      fPhi;

  // the panel volume of the tracker panel, placed in the tracker with its straws
  TEvdPanel(TEvdTracker* Tracker, int Plane, int Panel, bool CreateStraws = false);

  // panel volume: the sector of PanelGeometry minus the triangle with the origin and the two
  // sector ends as vertices
  TEvdPanel(const char* Name, const PanelGeometry& Geometry,
            TGeoMedium* Medium, int Color);

  // transparency of the panel volume (ROOT scale: 0 opaque, 100 invisible), set by the owner
  void        SetPanelTransparency(int Transparency);

  // straws as daughter volumes, positioned in the panel frame (z centre ZCenter subtracted)
  // visibility is off for now
  void        AddStraws(const mu2e::Panel& Panel, double StrawInnerRadius, double StrawOuterRadius,
                        double ZCenter, TGeoMedium* Medium, bool Visible = false, int Transparency = 100);
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
