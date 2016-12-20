// Helper script for showing of extracted / simplified geometries.
// By default shows a simplified ALICE geometry.

#if ! defined(__CINT__) || defined(__MAKECINT__)

#include "TEveTrackPropagator.h"
#include "TEveTrack.h"
#include "TEveVSDStructs.h"
#include "TEveManager.h"
#include "TEveViewer.h"
#include "TSystem.h"
#include "TGLViewer.h"
#include "TMath.h"

#include "TEveViewer.h"
#include "TEvePointSet.h"
#include "TGeoTube.h"
#include "TEveGeoShape.h"

#include <iostream>
#endif

TEveTrackPropagator* g_prop = 0;

TEveTrack* g_track;

namespace {
  double gPos[3] = {0.0, 0.0, -150.};

  double gMom[3] = {-0.08, -0.045, 0.05} ;

  double gMom2[3] = {-0.04, -0.023, 0.025} ;

  double phi0, p, pt, r, Nz, X0, Y0;

  double Rin(35.);

  int Sign (-1);

  double Step = 1.;
  
}

//-----------------------------------------------------------------------------
TEveTrack* make_track(TEveTrackPropagator* prop, Int_t sign, double* Pos, double* Mom) {
  // Make track with given propagator.
  // Add to math-marks to test fit.

  TEveRecTrackD *rc = new TEveRecTrackD();
  
  rc->fV.Set(Pos[0], Pos[1], Pos[2]);
  rc->fP.Set(Mom[0], Mom[1], Mom[2]);
  rc->fSign  = sign;

  TEveTrack* track = new TEveTrack(rc, prop);
  return track;
}



//-----------------------------------------------------------------------------
void print_numbers() {

  p  = sqrt(gMom[0]*gMom[0]+gMom[1]*gMom[1]+gMom[2]*gMom[2]);

  pt = sqrt(gMom[0]*gMom[0]+gMom[1]*gMom[1]);

  r  = pt/2.997*1000.; // in cm
  Nz        = gMom[2]/p;
					// normalized in 2D
  double nx = gMom[0]/pt;
  double ny = gMom[1]/pt;

  X0   =  gPos[0]-ny*r*(p/pt)*Sign;
  Y0   =  gPos[1]+nx*r*(p/pt)*Sign;

  phi0 = TMath::ATan2(Y0,X0) - TMath::Pi();

  printf("x0,y0,phi0 = %10.3f %10.3f %10.3f\n",X0,Y0,phi0);

  printf(" r = %10.3f\n",r);

  double s  = 0;
  double z  = gPos[2];

  double x1, y1, z1, phi;

  z1 = z;
  
  while  (z1 <  151.) {

    phi = phi0+Sign*s*(pt/p)/r;
    x1 = X0 + r*cos(phi);
    y1 = Y0 + r*sin(phi);
    z1 = gPos[2]+ Nz*s;

    printf(" x,y,z,s,nx,ny, r = %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n", x1,y1, z1, s, -Sign*sin(phi), Sign*cos(phi), sqrt(x1*x1+y1*y1));

    s += Step;
  }
}

//-----------------------------------------------------------------------------
void murat() {
  
  TEveManager::Create();

  TEveGeoShape* gs = new TEveGeoShape("emoe");

  TGeoTube* tube = new TGeoTube(35,70,150);

  gs->SetShape(tube);
  //  gs->SetMainTransparency(1);
  gs->SetMainColor(40);
  gs->SetMainAlpha(0.5);
  //  gs->SetMainAlpha(1);

  gEve->AddGlobalElement(gs);
//------------------------------------------------------------------------------
// print numbers
//-----------------------------------------------------------------------------
  print_numbers();
//-----------------------------------------------------------------------------
// make track
//-----------------------------------------------------------------------------
  TEveTrackPropagator* prop = new TEveTrackPropagator("prop");

  printf("emoe 000\n");

  TEveTrackList* list;

  list = new TEveTrackList("Heix Propagator",prop);

  printf("emoe 001\n");
  
  if (prop == 0) printf("Prop=0\n");
  
  prop->SetFitDaughters(kFALSE);

  printf("emoe 002\n");

  prop->SetMaxZ(150.);

  prop->SetMagFieldObj(new TEveMagFieldConst(0., 0., 1.));
  prop->SetMaxOrbs(20);


  printf("emoe 003\n");

  //  char name[100];
  // sprintf(name,"%s", list->GetElementName());
  
  // printf("emoe 003\n");
  
  // list->SetElementName(name);

  printf("emoe 0035\n");
  
  TEveTrack* track = make_track(prop, -1, gPos, gMom);
  
  printf("emoe 0036\n");

  if (track == 0) printf("emoe 0036 track == 0\n");

  //  track->SetName("track1");
  printf("emoe 004\n");

  list->AddElement(track);
  
  list->SetLineColor(kRed+3);
  track->SetLineColor(list->GetLineColor());
  track->SetLineWidth(1);
  track->SetLineStyle(2);
  
  printf("emoe 005\n");
  track->MakeTrack();
  
  gEve->AddElement(list);
//-----------------------------------------------------------------------------
// make another track
//-----------------------------------------------------------------------------
  // list = new TEveTrackList();
  // list->SetName("Heix Propagator");

  // prop = list->GetPropagator();
  // prop->SetFitDaughters(kFALSE);
  // prop->SetMaxZ(150.);

  // prop->SetMagFieldObj(new TEveMagFieldConst(0., 0., 1.));
  // prop->SetMaxOrbs(20);
  // list->SetElementName(Form("%s, constB", list->GetElementName()));
  
  // track = make_track(prop, -1, gPos, gMom2);
  // track->SetName("track11");

  // list->AddElement(track);
  
  // list->SetLineColor(kRed);
  // track->SetLineColor(list->GetLineColor());
  // track->SetLineWidth(1);
  // track->MakeTrack();
  
  //  gEve->AddElement(list);
//-----------------------------------------------------------------------------
// make first highlighted part
//-----------------------------------------------------------------------------
  int const nw = 7;

  double xx1[14] = {  // x,y,z,s,nx,ny, r
    -34.53925,    5.22155, -131.34377,   39.00000,   -0.82201,    0.56947,   34.93171, 
    -35.25384,    5.73192, -130.86541,   40.00000,   -0.80535,    0.59281,   35.71678 
  };

  double xx2[14] = {  // x,y,z,s,nx,ny, r
    13.46639,   32.34354,  -63.89432,  180.00000,    0.06356,   -0.99798,   35.03496, 
    13.50963,   31.46648,  -63.41596,  181.00000,    0.03492,   -0.99939,   34.24397 
  };

  double ds1 =  (Rin-xx1[6])/(xx1[6+nw]-xx1[6]);
  double ds2 =  (xx2[6]-Rin)/(xx2[6]-xx2[6+nw]);

  double s1 = xx1[3]+ds1*Step;
  double s2 = xx2[3]+ds2*Step;

  double pos[3], mom[3];


  double phi = phi0+Sign*s1*(pt/p)/r;

  dy = 3.5;
  
  pos[0] = X0 + r*cos(phi)+1.5;
  pos[1] = Y0 + r*sin(phi)-dy;
  pos[2] = gPos[2]+ Nz*s1;

  double nx  = -Sign*sin(phi);
  double ny  =  Sign*cos(phi);

  mom[0] = pt*nx;
  mom[1] = pt*ny;
  mom[2] = p*Nz;

  double z2 = gPos[2]+s2*Nz;
  printf(" pos, z2 = %10.4f %10.4f %10.4f %10.4f \n",pos[0],pos[1],pos[2],z2);
  printf(" mom     = %10.4f %10.4f %10.4f \n",mom[0],mom[1],mom[2]);

  TEveTrackPropagator* p2 = new TEveTrackPropagator("p2");
  p2->SetFitDaughters(kFALSE);
  p2->SetMagFieldObj(new TEveMagFieldConst(0., 0., 1.));

  double orbs = 0.61;
  
  p2->SetMaxOrbs(orbs);
  
  //  p2->SetMaxOrbs(1.01);

  TEveTrackList *list2 = new TEveTrackList("Helix Propagator P2",p2);
  list2->SetElementName(Form("%s_constB", list2->GetElementName()));

  TEveTrack *t2 = make_track(p2, Sign, pos, mom);

  t2->SetName("track2");

  list2->AddElement(t2);

  list2->SetLineColor(kRed+3);
  t2->SetLineColor(list2->GetLineColor());
  t2->SetLineWidth(5);
  t2->MakeTrack();
  
  gEve->AddElement(list2);

//-----------------------------------------------------------------------------
// 2nd segment
//========================================================================================
  phi = phi+2*TMath::Pi();

  pos[0] = X0 + r*cos(phi)+1.5;
  pos[1] = Y0 + r*sin(phi)-dy;

  s1 = s1+ 2*TMath::Pi()*r*p/pt;
    
  pos[2] = gPos[2]+ Nz*s1;

  nx  = -Sign*sin(phi);
  ny  =  Sign*cos(phi);

  mom[0] = pt*nx;
  mom[1] = pt*ny;
  mom[2] = p*Nz;

  z2 = gPos[2]+s2*Nz;
  printf(" pos, z2 = %10.4f %10.4f %10.4f %10.4f \n",pos[0],pos[1],pos[2],z2);
  printf(" mom     = %10.4f %10.4f %10.4f \n",mom[0],mom[1],mom[2]);

  p2 = new TEveTrackPropagator("p3");
  p2->SetFitDaughters(kFALSE);
  p2->SetMagFieldObj(new TEveMagFieldConst(0., 0., 1.));
  p2->SetMaxOrbs(orbs);
  //  p2->SetMaxOrbs(1.01);

  list2 = new TEveTrackList("Helix Propagator P3",p2);
  list2->SetElementName(Form("%s_constB", list2->GetElementName()));

  t2 = make_track(p2, Sign, pos, mom);
  t2->SetName("track2");

  list2->AddElement(t2);

  list2->SetLineColor(kRed+3);
  t2->SetLineColor(list2->GetLineColor());
  t2->SetLineWidth(5);
  t2->MakeTrack();
  
  gEve->AddElement(list2);

//-----------------------------------------------------------------------------
// 2nd segment
//========================================================================================
  phi = phi+2*TMath::Pi();

  pos[0] = X0 + r*cos(phi)+1.5;
  pos[1] = Y0 + r*sin(phi)-dy;

  s1 = s1+ 2*TMath::Pi()*r*p/pt;
    
  pos[2] = gPos[2]+ Nz*s1;

  nx  = -Sign*sin(phi);
  ny  =  Sign*cos(phi);

  mom[0] = pt*nx;
  mom[1] = pt*ny;
  mom[2] = p*Nz;

  z2 = gPos[2]+s2*Nz;
  printf(" pos, z2 = %10.4f %10.4f %10.4f %10.4f \n",pos[0],pos[1],pos[2],z2);
  printf(" mom     = %10.4f %10.4f %10.4f \n",mom[0],mom[1],mom[2]);

  p2 = new TEveTrackPropagator("p4");
  p2->SetFitDaughters(kFALSE);
  p2->SetMagFieldObj(new TEveMagFieldConst(0., 0., 1.));
  p2->SetMaxOrbs(orbs);
  //  p2->SetMaxOrbs(1.01);

  list2 = new TEveTrackList("Helix Propagator P4",p2);
  list2->SetElementName(Form("%s_constB", list2->GetElementName()));

  t2 = make_track(p2, Sign, pos, mom);
  t2->SetName("track2");

  list2->AddElement(t2);

  list2->SetLineColor(kRed+3);
  t2->SetLineColor(list2->GetLineColor());
  t2->SetLineWidth(5);
  t2->MakeTrack();
  
  gEve->AddElement(list2);


  gEve->Redraw3D(kTRUE);
}
