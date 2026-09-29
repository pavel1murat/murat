///////////////////////////////////////////////////////////////////////////////
#include "TGeoMatrix.h"
#include "murat/gui/TEvdSubdetector.hh"

ClassImp(murat::TEvdSubdetector)

namespace murat {
  //-----------------------------------------------------------------------------
  TEvdSubdetector::TEvdSubdetector() : TGeoVolumeAssembly() {
    fListOfSubdetectors = new TObjArray();
  }

  //-----------------------------------------------------------------------------
  TEvdSubdetector::TEvdSubdetector(const char* Name)
    : TGeoVolumeAssembly(Name),
      fName(Name),
      fListOfSubdetectors(new TObjArray())
  {
  }

  //-----------------------------------------------------------------------------
  TEvdSubdetector::TEvdSubdetector(const char* Name, TGeoShape* Shape, TGeoMedium* Medium,
                                   int Color, int Transparency)
    : TGeoVolumeAssembly(Name),
      fListOfSubdetectors(new TObjArray())
  {
    TGeoVolume* vol = new TGeoVolume(Name,Shape,Medium);
    vol->SetLineColor(Color);
    vol->SetFillColor(Color);
    vol->SetTransparency(Transparency);
    AddNode(vol,1,new TGeoRotation());
  }

  //-----------------------------------------------------------------------------
  TEvdSubdetector::~TEvdSubdetector() {
    fListOfSubdetectors->Delete();
    delete fListOfSubdetectors;
  }

  //-----------------------------------------------------------------------------
  int TEvdSubdetector::InitGeometry(const char* Fn) {
    return 0;
  }
  
  int TEvdSubdetector::InitEvent() {
    return 0;
  }

  //-----------------------------------------------------------------------------
  void TEvdSubdetector::AddSubdetector(TEvdSubdetector* sd) {
    if (!sd) return;
    
    fListOfSubdetectors->Add(sd);

    // upon creation, a subdetector has to have a copy number...
    AddNode(sd, sd->CopyNumber(), new TGeoTranslation());
  }
}
