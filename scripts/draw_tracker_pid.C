//


//-----------------------------------------------------------------------------
double f_trajectory(double* X, double* P) {
  double x, y;

  x = X[0];
  y = P[0]+P[1]*TMath::Sin(P[2]*x+P[3]);
  return y;
}

//-----------------------------------------------------------------------------
void draw_tracker_pid() {

  TCanvas* c = new TCanvas("c_pid","c_pid",1400,700);
//-----------------------------------------------------------------------------
// draw tracker cartoon
//-----------------------------------------------------------------------------
  TLine* l;
  TBox*  b;

  double  x1, x2, y1, y2;

  for (int i=0; i<22; i++) {
    x1 = 0.1+0.03*i;
    x2 = 0.1+0.03*i;
    y1 = 0.2;
    y2 = 0.4;
    l = new TLine(x1,y1,x2,y2);
    l->Draw();

    x1 = 0.1+0.03*i;
    x2 = 0.1+0.03*i;
    y1 = 0.6;
    y2 = 0.8;
    l = new TLine(x1,y1,x2,y2);
    l->Draw();
  }
//-----------------------------------------------------------------------------
// draw calorimeter
//-----------------------------------------------------------------------------
  x1 = 0.8;
  x2 = 0.83;
  y1 = 0.2;
  y2 = 0.4;

  b = new TBox(x1,y1,x2,y2);
  b->SetLineColor(1);
  b->SetFillColor(601);
  b->SetFillStyle(3001);
  b->Draw();

  x1 = 0.8;
  x2 = 0.83;
  y1 = 0.6;
  y2 = 0.8;

  b = new TBox(x1,y1,x2,y2);
  b->SetFillColor(601);
  b->SetFillStyle(3001);
  b->Draw();

  x1 = 0.9;
  x2 = 0.93;
  y1 = 0.2;
  y2 = 0.4;

  b = new TBox(x1,y1,x2,y2);
  b->SetFillColor(601);
  b->SetFillStyle(3001);
  b->Draw();

  x1 = 0.9;
  x2 = 0.93;
  y1 = 0.6;
  y2 = 0.8;

  b = new TBox(x1,y1,x2,y2);
  b->SetFillColor(601);
  b->SetFillStyle(3001);
  b->Draw();

  c->Draw();
//-----------------------------------------------------------------------------
// draw tracks
//-----------------------------------------------------------------------------

  TF1* f1 = new TF1("f1",f_trajectory,0.1,0.9,4);
  f1->SetParameters(0.5-0.15,0.15,16.,1.5);
  f1->Draw("same");

  TF1* f2 = new TF1("f2",f_trajectory,0.05,0.8,4);
  f2->SetParameters(0.5+0.15,0.15,16.,2.5);
  f2->Draw("same");

}
