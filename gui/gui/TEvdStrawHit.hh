//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdStrawHit__
#define __murat_gui_TEvdStrawHit__

#include "TEveGeoShape.h"

class TEvdStrawHit : public TEveGeoShape {
public:
  int        fIndex;			// in the data collection
  int        fStrawIndex;		// straw index - channel ID
  int        fPlane;
  int        fPanel;
  int        fLayer;
  int        fStraw;
  int        fDeltaID;			// reconstructed delta electron ID or -1
  float      fTime;
  float      fDt;
  float      fEDep;
  float      fWDist;
  float      fWRes;
  int        fSimID;
  int        fPdgID;
  float      fX;
  float      fY;
  float      fZ;
  float      fMcMom;			// momentum of the associated MC particle

  TEveGeoShape*  fErrorBars;

  TEvdStrawHit();
  ~TEvdStrawHit();
  
  void Init(int Index, int StrawIndex, int Plane, int Panel,
	    int Layer, int Straw, int DeltaID,
	    float Time, float Dt, float EDep,
	    float WDist, float WRes,
	    float X, float Y, float Z);

  void SetSimID(int ID)  { fSimID = ID; }
  void SetPdgID(int ID)  { fPdgID = ID; }
  void SetMcMom(float P) { fMcMom = P ; }

  void     Print(Option_t* Opt) const; // *MENU*

  ClassDef(TEvdStrawHit,0);
};

#endif
