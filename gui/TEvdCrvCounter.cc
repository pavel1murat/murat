///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <format>

#include "murat/gui/TEvdCrvCounter.hh"

ClassImp(murat::TEvdCrvCounter)

namespace murat {
//-----------------------------------------------------------------------------
TEvdCrvCounter::TEvdCrvCounter(const char* Name, TGeoShape* Shape, TGeoMedium* Medium):
  TGeoVolume(Name,Shape,Medium) {
  fListOfHits = new TObjArray();
          // the counter color will depend on whether the bar has hits
  SetLineColor(kCyan + 1);
  SetFillColor(kCyan + 1);
  SetTransparency(0);
}

//-----------------------------------------------------------------------------
  void TEvdCrvCounter::AddRecoPulse(TCrvRecoPulse* Hit) {
    fListOfHits->Add(Hit);
    
    SetLineColor(kRed + 1);
    SetFillColor(kRed + 1);
  }

//-----------------------------------------------------------------------------
void TEvdCrvCounter::Clear(Option_t* Opt) {
  fListOfHits->Clear();
  SetLineColor(kCyan + 1);
  SetFillColor(kCyan + 1);
}

//-----------------------------------------------------------------------------
void TEvdCrvCounter::Print(Option_t* Opt) const {
  std::cout << std::format("TEvdCrvCounter::Print emoe\n");
}

//-----------------------------------------------------------------------------
void TEvdCrvCounter::PrintA() const {
  std::cout << std::format("TEvdCrvCounter::PrintA emoe AAAAAAA\n");
}

}
