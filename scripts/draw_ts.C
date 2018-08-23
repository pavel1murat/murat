// to be interpreted
#include "TLine.h"
#include "TArrow.h"

TTsMisalignment*  gTs;




//-----------------------------------------------------------------------------
void init(double Scale) {

  if (gTs == NULL) {
    gTs = new TTsMisalignment();
    gTs->Init(Scale);
  }
}

//-----------------------------------------------------------------------------
void draw_all(double Scale = 1.) {

  init(Scale);

  TCanvas* c = new TCanvas("a","a",900,900);

  c->Range(-600,-600,600,600);


  TArrow* ax = new TArrow(0,0,200,0,0.015,">");
  TArrow* ay = new TArrow(0,0,0,200,0.015,">");

  ax->SetLineColor(kRed+2);
  ay->SetLineColor(kRed+2);

  ax->Draw();
  ay->Draw();
					// invert x1, y1
  double x1 = -417.50;
  double y1 = +405.00;
  
  TArrow* ax1 = new TArrow(x1,y1,x1+200,y1    ,0.015,">");
  TArrow* ay1 = new TArrow(x1,y1,x1    ,y1-200,0.015,">");

  ax1->SetLineColor(kRed+2);
  ay1->SetLineColor(kRed+2);

  ax1->Draw();
  ay1->Draw();
					// X and Y inverted
  double x2 = +435.00;
  double y2 = -405.00;
  
  TArrow* ax2 = new TArrow(x2,y2,x2    ,y2+200,0.015,">");  // x2'
  TArrow* ay2 = new TArrow(x2,y2,x2-200,y2    ,0.015,">");  // y2'

  ax2->SetLineColor(kRed+2);
  ay2->SetLineColor(kRed+2);

  ax2->Draw();
  ay2->Draw();

  gTs->Draw();
}


//-----------------------------------------------------------------------------
void draw_tsul(double Scale = 1.) {

  init(Scale);

  TCanvas* c = new TCanvas("c_tsul","c_tsul",1400,800);

  c->Range(-430,380,-345,425);


  TArrow* ax = new TArrow(0,0,200,0,0.015,">");
  TArrow* ay = new TArrow(0,0,0,200,0.015,">");

  ax->SetLineColor(kRed+2);
  ay->SetLineColor(kRed+2);

  ax->Draw();
  ay->Draw();
					// invert x1, y1
  double x1 = -417.50;
  double y1 = +405.00;
  
  TArrow* ax1 = new TArrow(x1,y1,x1+70 ,y1   ,0.015,">");
  TArrow* ay1 = new TArrow(x1,y1,x1    ,y1-20,0.015,">");

  ax1->SetLineColor(kRed+2);
  ay1->SetLineColor(kRed+2);

  ax1->Draw();
  ay1->Draw();
					// X and Y inverted
  double x2 = +435.00;
  double y2 = -405.00;
  
  TArrow* ax2 = new TArrow(x2,y2,x2    ,y2+200,0.015,">");  // x2'
  TArrow* ay2 = new TArrow(x2,y2,x2-200,y2    ,0.015,">");  // y2'

  ax2->SetLineColor(kRed+2);
  ay2->SetLineColor(kRed+2);

  ax2->Draw();
  ay2->Draw();

  gTs->Draw();
}


//-----------------------------------------------------------------------------
void draw_tsur(double Scale = 1.) {

  init(Scale);

  TCanvas* c = new TCanvas("c_tsur","c_tsur",800,900);

  c->Range(-25,-10,25,110);


  TArrow* ax = new TArrow(0,0,200,0,0.015,">");
  TArrow* ay = new TArrow(0,0,0,200,0.015,">");

  ax->SetLineColor(kRed+2);
  ay->SetLineColor(kRed+2);

  ax->Draw();
  ay->Draw();
					// invert x1, y1 to position the axes
  double x1 = -417.50;
  double y1 = +405.00;
  
  TArrow* ax1 = new TArrow(x1,y1,x1+70 ,y1   ,0.015,">");
  TArrow* ay1 = new TArrow(x1,y1,x1    ,y1-20,0.015,">");

  ax1->SetLineColor(kRed+2);
  ay1->SetLineColor(kRed+2);

  ax1->Draw();
  ay1->Draw();
					// X and Y inverted
  double x2 = +435.00;
  double y2 = -405.00;
  
  TArrow* ax2 = new TArrow(x2,y2,x2    ,y2+200,0.015,">");  // x2'
  TArrow* ay2 = new TArrow(x2,y2,x2-200,y2    ,0.015,">");  // y2'

  ax2->SetLineColor(kRed+2);
  ay2->SetLineColor(kRed+2);

  ax2->Draw();
  ay2->Draw();

  gTs->Draw();
}


//-----------------------------------------------------------------------------
void draw_tsdr(double Scale = 1.) {

  init(Scale);

  TCanvas* c = new TCanvas("c_tsdr","c_tsdr",800,900);

  c->Range(-25, -90,25, 10);

  // common axes (at 0,0)
  TArrow* ax = new TArrow(0,0,200,0  ,0.015,">");
  TArrow* ay = new TArrow(0,0,0  ,200,0.015,">");

  ax->SetLineColor(kRed+2);
  ay->SetLineColor(kRed+2);

  ax->Draw();
  ay->Draw();
					// invert x1, y1 to position the axes
  double x1 = -417.50;
  double y1 = +405.00;
  
  TArrow* ax1 = new TArrow(x1,y1,x1+70 ,y1   ,0.015,">");
  TArrow* ay1 = new TArrow(x1,y1,x1    ,y1-20,0.015,">");

  ax1->SetLineColor(kRed+2);
  ay1->SetLineColor(kRed+2);

  ax1->Draw();
  ay1->Draw();
					// X and Y inverted
  double x2 = +435.00;
  double y2 = -405.00;
  
  TArrow* ax2 = new TArrow(x2,y2,x2    ,y2+200,0.015,">");  // x2'
  TArrow* ay2 = new TArrow(x2,y2,x2-200,y2    ,0.015,">");  // y2'

  ax2->SetLineColor(kRed+2);
  ay2->SetLineColor(kRed+2);

  ax2->Draw();
  ay2->Draw();

  gTs->Draw();
}

//-----------------------------------------------------------------------------
void draw_tsdl(double Scale = 1.) {

  init(Scale);

  TCanvas* c = new TCanvas("c_tsdl","c_tsdl",1400,800);

  c->Range(370,-430, 450, -380);

  // common axes (at 0,0)
  TArrow* ax = new TArrow(0,0,200,0  ,0.015,">");
  TArrow* ay = new TArrow(0,0,0  ,200,0.015,">");

  ax->SetLineColor(kRed+2);
  ay->SetLineColor(kRed+2);

  ax->Draw();
  ay->Draw();
					// invert x1, y1 to position the axes
  double x1 = -417.50;
  double y1 = +405.00;
  
  TArrow* ax1 = new TArrow(x1,y1,x1+70 ,y1   ,0.015,">");
  TArrow* ay1 = new TArrow(x1,y1,x1    ,y1-20,0.015,">");

  ax1->SetLineColor(kRed+2);
  ay1->SetLineColor(kRed+2);

  ax1->Draw();
  ay1->Draw();
					// X and Y inverted
  double x2 = +435.00;
  double y2 = -405.00;
  
  TArrow* ax2 = new TArrow(x2,y2,x2    ,y2+200,0.015,">");  // x2'
  TArrow* ay2 = new TArrow(x2,y2,x2-200,y2    ,0.015,">");  // y2'

  ax2->SetLineColor(kRed+2);
  ay2->SetLineColor(kRed+2);

  ax2->Draw();
  ay2->Draw();

  gTs->Draw();
}
