///////////////////////////////////////////////////////////////////////////////
#include "TGeoMatrix.h"
#include "murat/gui/TEvdSubdetector.hh"
#include "murat/gui/TEvdManager.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"

ClassImp(murat::TEvdSubdetector)

namespace murat {
  //-----------------------------------------------------------------------------
  TEvdSubdetector::TEvdSubdetector() : TGeoVolumeAssembly() {
    fListOfSubdetectors = new TObjArray();
  }

  //-----------------------------------------------------------------------------
  TEvdSubdetector::TEvdSubdetector(const char* Name, int CopyNumber)
    : TGeoVolumeAssembly(Name),
      fName(Name),
      fCopyNumber(CopyNumber),
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

 
  //-----------------------------------------------------------------------------
  bool TEvdSubdetector::Initialized() {

    auto vm = TEvdManager::Instance();
    
    TStnHeaderBlock* hb = (TStnHeaderBlock*) vm->GetDataBlock("HeaderBlock");
    int evn = hb->EventNumber ();
    int srn = hb->SubrunNumber();
    int run = hb->RunNumber   ();

    bool initialized = false;
    if ((evn == f_EventNumber) and (srn == f_SubrunNumber) and (run == f_RunNumber)) {
      initialized = true;
    }
    return initialized;
  }

  //-----------------------------------------------------------------------------
  int TEvdSubdetector::InitEvent() {
    return 0;
  }

  //-----------------------------------------------------------------------------
  void TEvdSubdetector::AddSubdetector(TEvdSubdetector* sd, TGeoMatrix* Matrix) {
    if (!sd) return;
    
    fListOfSubdetectors->Add(sd);

    // upon creation, a subdetector has to have a copy number...
    if (Matrix == nullptr) {
      GetVolume()->AddNode(sd, sd->CopyNumber(), new TGeoTranslation());
    }
    else {
      GetVolume()->AddNode(sd, sd->CopyNumber(), Matrix);
    }
  }
}
