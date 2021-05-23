//-----------------------------------------------------------------------------
#ifndef __murat_gui_View3D_hh__
#define __murat_gui_View3D_hh__

#include "TAtt3D.h"
#include "TObjArray.h"

#include "TGeoMedium.h"
#include "Stntuple/gui/TStnGeoManager.hh"
#include "murat/gui/TEvdTrack.hh"
#include "murat/gui/TEvdTrajectory.hh"

class View3D : public TObject, public TAtt3D {
public:

  struct Materials_t {
    TGeoMaterial*  fVacuum;
    TGeoMaterial*  fVac1;
    TGeoMaterial*  fAl;
  };
  
  struct Media_t {
    TGeoMedium*  fVac1;
    TGeoMedium*  fVacuum;
    TGeoMedium*  fAl;
  };
  
  Materials_t  fMaterials;
  Media_t      fMedia;
  
  TStnGeoManager* fGeoManager;
  TGeoVolume*     fMu2eGeom;

  
  TEvdTrack*      fTrack;

  TObjArray*      fListOfTrajectories;
					// materials
  TGeoMaterial*   matVacuum;
public:
  View3D();
  ~View3D();

  int InitGeometry();
    
  int InitCalorimeterGeometry();
  int InitCrvGeometry        ();
  int InitTrackerGeometry    ();
  
  int ReadTrajectories(const char* DirName);
//-----------------------------------------------------------------------------
// overloaded virtual functions of TObject
//-----------------------------------------------------------------------------
  virtual void Draw (Option_t *option);
  virtual void Paint(Option_t *option);

   ClassDef(View3D,0);
};

#endif
