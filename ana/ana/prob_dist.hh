//  include file murat/ana/prob_dist.hh

#ifndef __murat_ana_prob_dist_hh__
#define __murat_ana_prob_dist_hh__

#include "TH1.h"

class prob_dist {
public:

  static int fgIndex;

  TH1*  h;
  TH1F* hprob;

  prob_dist();
  prob_dist(TH1F* Hist);
  prob_dist(const char* Fn, const char* Folder, const char* Hist);

  ~prob_dist();

  double prob(double X);

};

#endif
