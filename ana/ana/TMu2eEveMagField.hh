//-----------------------------------------------------------------------------
// apparenty , the field is in Tesla, not in Gauss
//-----------------------------------------------------------------------------
#ifndef murat_ana_TMu2eEveMagField
#define murat_ana_TMu2eEveMagField

#include "TObject.h"
#include "TEveTrackPropagator.h"
#include "TEveVector.h"

namespace mu2e {
  class BFieldManagerMaker;
  class BFieldConfig;
  class Beamline;
  class BFieldManager;
}

class TMu2eEveMagField: public TEveMagField {
protected:
  int    fN;
  float* fZ;
  float* fBz;

  mu2e::Beamline*           fBeamline;
  mu2e::BFieldConfig*       fBfc;
  mu2e::BFieldManagerMaker* fBfmm;
  mu2e::BFieldManager*      fBmgr;
  
public:
  TMu2eEveMagField();
  ~TMu2eEveMagField();
  
  virtual Float_t    GetMaxFieldMag() const { return 5. ; }
  virtual TEveVector GetField(float X, float Y, float Z) const ;

  ClassDef(TMu2eEveMagField, 0);
};

#endif
