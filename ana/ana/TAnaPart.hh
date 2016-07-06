//
#ifndef __murat_TAnaPart__
#define __murat_TAnaPart__

#include "TObject.h"
#include "murat/ana/AnaDefs.hh"

class TStnTrack;
class TTrackStrawHitBlock;
class TStnPid;

class TAnaPart: public TObject {
public:

  TStnTrack*            fTrack    [8];
  TTrackStrawHitBlock*  fTrackHits[8];
  TStnPid*              fPid      [8];

  int                   fIndex    [8];      

  int                   fBestHyp;

  TStnTrack*            fBestTrack;  // fTrack[fBestHyp]
};

#endif
