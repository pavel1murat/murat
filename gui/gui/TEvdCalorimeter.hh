//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdCalorimeter__
#define __murat_gui_TEvdCalorimeter__

#include "murat/gui/TEvdDisk.hh"
#include "murat/gui/TEvdSubdetector.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

namespace murat {
  
//-----------------------------------------------------------------------------
// fMu2eDisk is owned by the Mu2e calorimeter
//-----------------------------------------------------------------------------
class TEvdCalorimeter: public TEvdSubdetector {
public:
  mu2e::Disk*  fMu2eDisk;
  
  std::unique_ptr<mu2e::DiskCalorimeter> fCaloPtr;
  TEvdDisk*                              fDisk[2];

  TEvdCalorimeter();
  
  TEvdCalorimeter(const char* Fn = "murat/fcl/geom_common_extracted_v04.txt");

  virtual int InitEvent() override;
  virtual int InitGeometry(const char* Fn) override;
  
  virtual void Print (Option_t* Opt = "") const override;  // *MENU*
  
  virtual void Draw(Option_t* Opt = "") override;

  ClassDefOverride(murat::TEvdCalorimeter,0)
};

}
#endif
