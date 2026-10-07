///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <format>

#include "murat/gui/TEvdCrystal.hh"

ClassImp(murat::TEvdCrystal)

namespace murat {
//-----------------------------------------------------------------------------
  TEvdCrystal::TEvdCrystal(const char* Name, TGeoShape* Shape, TGeoMedium* Medium):
    TGeoVolume(Name,Shape,Medium) {
    fListOfHits = new TObjArray();
    fEDep = 0;
  }
  
//-----------------------------------------------------------------------------
  void TEvdCrystal::AddHit(TCaloHit* Hit) {
    fListOfHits->Add(Hit);
    fEDep += Hit->EDep();

    SetLineColor(kRed + 1);
    SetFillColor(kRed + 1);
  }

//-----------------------------------------------------------------------------
void TEvdCrystal::Clear(Option_t* Opt) {
  fListOfHits->Clear();
  fEDep = 0;
                                        // crystal color will depend on whether the crystal has hits
  SetLineColor(kOrange + 1);
  SetFillColor(kOrange + 1);
  SetTransparency(0);
}

//-----------------------------------------------------------------------------
void TEvdCrystal::Print(Option_t* Opt) const {
  std::cout << std::format("emoe\n");
}

//-----------------------------------------------------------------------------
void TEvdCrystal::PrintA() const {
  std::cout << std::format("emoe AAAAAAA\n");
}

}
