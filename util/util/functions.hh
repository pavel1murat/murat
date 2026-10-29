#ifndef  __murat_util_hh__
#define  __murat_util_hh__

#include "TH1.h"

int fit_edep      (TH1* Hist, const char* FitOpt = "sql", double XMin = 0, double XMax = 1);
int fit_panel_edep(int IMin = 0, int IMax = 216,
                   const char* FileName = "hst.mu2e.kpp00s000r104n013.make_det_003.125470_000001.root");
#endif
