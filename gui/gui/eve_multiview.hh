#ifndef __murat_gui_eve_multiview_hh__
#define __murat_gui_eve_multiview_hh__

#include <TEveManager.h>
#include <TEveViewer.h>
#include <TGLViewer.h>
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>
#include <TEveBrowser.h>
#include <TEveWindow.h>

namespace murat {
class Mu2eMultiView {
public:
  
  TEveProjectionManager *fRPhiMgr;
  TEveProjectionManager *fRhoZMgr;

  TEveViewer            *f3DView;
  TEveViewer            *fRPhiView;
  TEveViewer            *fRhoZView;

  TEveScene             *fRPhiGeomScene;
  TEveScene             *fRhoZGeomScene;
  TEveScene             *fRPhiEventScene;
  TEveScene             *fRhoZEventScene;

  Mu2eMultiView();
  ~Mu2eMultiView();

  void SetDepth(float Depth);

  void ImportGeomRPhi(TEveElement* el);
  void ImportGeomRhoZ(TEveElement* el);

  void ImportEventRPhi(TEveElement* el);
  void ImportEventRhoZ(TEveElement* el);
  
  void DestroyEventRPhi();
  void DestroyEventRhoZ();
  
};
}
#endif
