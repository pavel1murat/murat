#ifndef TTrkTZView_hh
#define TTrkTZView_hh


#include "TMarker.h"
#include "TNamed.h"

#include "Stntuple/gui/TStnView.hh"

class TTrkTZView: public TStnView {
protected:
public:
  TTrkTZView();
  virtual ~TTrkTZView();

  virtual void  Paint              (Option_t* option = "");
  virtual void  ExecuteEvent       (Int_t event, Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitive(Int_t px, Int_t py);

  void    SetStations(int I1, int I2);   // *MENU* 
  void    SetTimeCluster(int I);         // *MENU* 

  ClassDef(TTrkTZView,0)
};

#endif
