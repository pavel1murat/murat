// plot radiation length X0 of a GaAs-lead mix vs x - relative thickness of lead

#include "TCanvas.h"
#include "TF1.h"

void plot_gaas_lead_x0() {

  TCanvas* c  = new TCanvas("c","c",1000,700);

  TF1* f = new TF1("f","1/(x/2.23+(1-x)/0.56)",0,1);

  f->Draw();
}
