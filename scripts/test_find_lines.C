///////////////////////////////////////////////////////////////////////////////
// determine parameters of 4 lines tangent to two given circles
//
// direction cosine convention: nx*(xb-xa)+ny*(yb-ya) > 0
// (from the first to second point)
///////////////////////////////////////////////////////////////////////////////
void find_lines(double xa, double ya, double ra,
		double xb, double yb, double rb,
		double *Nx, double* Ny, int* Invert, int *ambStrawA, int *ambStrawB) {
  
  double            cosVec[4], sinVec[4], dx, dy, dr2;

//   CLHEP::Hep3Vector dir1;              // traj direction
//   CLHEP::Hep3Vector dir2(0., 0., -1.); // wire direction
  
//   CLHEP::Hep3Vector pos1, pos2, delta, between;

  double alpha, dr, dr2, lx, ly;

  dx  = xb-xa;
  dy  = yb-ya;
  dr2 = dx*dx+dy*dy;
  dr  = sqrt(dr2);

  lx  = dx/dr;
  ly  = dy/dr;
//-----------------------------------------------------------------------------
// sign ordering convention: (++) , (+-), (--), (-+)
//-----------------------------------------------------------------------------
  double sa[4] = {+1, +1, -1, -1};
  double sb[4] = {+1, -1, -1, +1};

//   alpha = (rb*sb[0]-ra*sa[0])/dr;
//   Nx[0] = -sa[0]*ly*alpha+sa[0]*lx*sqrt(1-alpha*alpha); // (+ +)
//   Ny[0] =  sa[0]*lx*alpha+sa[0]*ly*sqrt(1-alpha*alpha);
    
//   alpha = (rb*sb[1]-ra*sa[1])/dr;
//   Nx[1] = -sa[1]*ly*alpha+sa[1]*lx*sqrt(1-alpha*alpha); // (+ -)
//   Ny[1] =  sa[1]*lx*alpha+sa[1]*ly*sqrt(1-alpha*alpha);

//   alpha = (rb*sb[2]-ra*sa[2])/dr;
//   Nx[2] = -sa[2]*ly*alpha+sa[2]*lx*sqrt(1-alpha*alpha); // (- -)
//   Ny[2] =  sa[2]*lx*alpha+sa[2]*ly*sqrt(1-alpha*alpha);
  
//   alpha = (rb*sb[3]-ra*sa[3])/dr;

//   Nx[3] = -sa[3]*ly*alpha+sa[3]*lx*sqrt(1-alpha*alpha); // (- +)
//   Ny[3] =  sa[3]*lx*alpha+sa[3]*ly*sqrt(1-alpha*alpha);

  for (int i=0; i<4; i++) {
    alpha     = (rb*sb[i]-ra*sa[i])/dr;
    Nx[i]     = -sa[i]*ly*alpha+sa[i]*lx*sqrt(1-alpha*alpha); // (- +)
    Ny[i]     =  sa[i]*lx*alpha+sa[i]*ly*sqrt(1-alpha*alpha);
    Invert[i] =  1;

    if (Nx[i]*dx+Ny[i]*dy < 0) {
     Nx[i]     = -Nx[i];
     Ny[i]     = -Ny[i];
     Invert[i] = -1;
   }
  }


//     for (int i=0; i<4; i++) {
//       drho      = sa[i]*ra+sb[i]*rb;
    
//       cosVec[i] = (dx*drho + dy*sqrt(dr2 - drho*drho))/dr2;
//       sinVec[i] = (dy*drho - dx*sqrt(dr2 - drho*drho))/dr2;

//       slopes[i] = -cosVec[i]/sinVec[i]

//       dir1.setX(-cosVec[i]);
//       dir1.setY( sinVec[i]);

//       pos1.setX(xa - ra*signA[i]*cosVec[i]);
//       pos1.setY(ya - ra*signA[i]*sinVec[i]);
    
//       pos2.setX(xa);
//       pos2.setY(ya);
//       pos2.setZ(0);

//       delta.set(pos2.x()-pos1.x(),pos2.y()-pos1.y(),pos2.z()-pos1.z());

//       between = dir1.cross(dir2).unit();
//       ambStrawA[i] = delta.dot(between) > 0.0 ? 1 : -1;
  
//       pos1.setX(xb + rb*signB[i]*cosVec[i]);
//       pos1.setY(yb + rb*signB[i]*sinVec[i]);
    
//       pos2.set(xb,yb,0);
//       //      delta   = pos2 - pos1;
//       delta.set(pos2.x()-pos1.x(),pos2.y()-pos1.y(),pos2.z()-pos1.z());
  
//       ambStrawB[i] = delta.dot(between) > 0.0 ? 1 : -1;
//     }
}


//-----------------------------------------------------------------------------
void test_find_lines() {

  double x1, y1, x2, y2, r1, r2;
  double nx[4], ny[4];

  int    iamb1[4], iamb2[4], invert[4];

  x1 = 1.;
  y1 = 2.;

  x2 = -10.;
  y2 = -8.;

  r1 = 2.;
  r2 = 3.;

  TH2F* h2 = new TH2F("h2","h2",1000,-25,25,1000,-25,25);
  h2->SetStats(0);
  h2->Draw();

  TEllipse* e1 = new TEllipse(x1,y1,r1);
  TEllipse* e2 = new TEllipse(x2,y2,r2);


  e1->Draw();
  e2->Draw();

  double s = 20;

  double ux, uy, x0, y0;

  find_lines(x1,y1,r1,x2,y2,r2,nx,ny,invert,iamb1,iamb2);

  for (int i=0; i<4; i+=1) {
    printf("i, nx, ny = %3i %10.4f %10.4f\n",i,nx[i],ny[i]);

    ux = -ny[i]*invert[i];
    uy =  nx[i]*invert[i];

    x0  = x1+r1*ux;
    y0  = y1+r1*uy;

    line = new TArrow(x0-nx[i]*s,y0-ny[i]*s,x0+nx[i]*s,y0+ny[i]*s,0.015);

    line->Draw();
  }
}
