//

void print_circle(int NPoints = 50, double Dr = 0) {

  double r0 = 250.;

  double x0 = 250.;
  double y0 = 70.;

  double phi;

  double x, y, rho, r, phi;
  TRandom3   rn;

  for (int i=0; i<NPoints; i++) {

    phi = 0.02*i;

    if (i%2 == 0) r = r0+Dr*rn.Gaus();
    else          r = r0-Dr*rn.Gaus();

    x = x0+r*cos(phi);
    y = y0+r*sin(phi);

    rho = sqrt(x*x+y*y);

    printf("  %10.3f %10.3f %10.3f\n",x,y,rho);
  }
}
