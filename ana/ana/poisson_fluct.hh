//
#ifndef __murat_ana_poisson_fluct_hh__
#define __murat_ana_poisson_fluct_hh__

// histograms for E/Mu + MIXP-x1
// histset: 1412 for offline v4_2_1
//          0041 for offline v4_2_4 (new dataset naming conventions)

class poisson_fluct: public TObject {
public:

  int fNTries;


  poisson_fluct();
  
  static double fun    (double*x, double* p);
  
  int    poi1   (double N, double x   );
  int    poi    (double N, double Mean);

  int    l95    (double Signal);
  
  int    limit95(double Expected, double Observed);

  int    fluct_prob(double Expected, double Observed);


  ClassDef(poisson_fluct, 0)
};

#endif
