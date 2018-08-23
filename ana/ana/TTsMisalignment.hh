//
#ifndef __murat_ana_TTsMisalignment
#define __murat_ana_TTsMisalignment

#include "TObject.h"
#include "TPolyLine.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TVector3.h"

#include <math.h>

class TTsMisalignment: public TObject {

  struct TTs_t {
    TPolyLine* fTsuL;
    TPolyLine* fTsuLMeas;
  
    TPolyLine* fTsuR;
    TPolyLine* fTsuRMeas;

    TPolyLine* fTsdL;
    TPolyLine* fTsdLMeas;

    TPolyLine* fTsdR;
    TPolyLine* fTsdRMeas;
  };

public:

  TTs_t  fTs;
  TTs_t  fTsRot;    // TS rotated to have fTsuRMeas and fTsdRMeas positioned around zero

  TPolyLine* fTsuL;
  TPolyLine* fTsuLMeas;
  
  TPolyLine* fTsuR;
  TPolyLine* fTsuRMeas;

  TPolyLine* fTsdL;
  TPolyLine* fTsdLMeas;

  TPolyLine* fTsdR;
  TPolyLine* fTsdRMeas;

  double     fR;   // inner radius;

  double     fScale;

  TTsMisalignment();
  ~TTsMisalignment();
  
  void  Init(double Scale = 1);
  void  InitPolyLine(TVector3* Pos, TPolyLine*& Pl);
  void  Transform(double X0, double Y0, double Phi, TVector3* X, TVector3* Xr);
  
  virtual void Paint(Option_t* Opt);
  
};

#endif
