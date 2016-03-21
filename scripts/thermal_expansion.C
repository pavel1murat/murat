//

//-----------------------------------------------------------------------------
// thermal expansion
// tungsten (W) : from  IntJournThermophysics v11, 4, p619 (1990)
//
// in the formula, T is in K
//
// (L(T)-L(T0))/L(T0) = l.3896e-3 -8.2797e-7*T+4.0557e-9*T^2 -1.2164e-12*T^3+ 1.7034e-16T^4
//
// calculate dL(T)/dT = k(T)*L(T), where k(T) is the thermal expansion coefficient
//
// input for the function: T in C
//-----------------------------------------------------------------------------
double thermal_expansion_tungsten(double* T, double* P) {

  double c[5] = {
    1.3896e-3, -8.2797e-7, 4.0557e-9, -1.2164e-12, 1.7034e-16
  };


  double f, dfdt;

  double t  = T[0]+273;  // convert C --> K
  double t2 = t*t;
  
  f = c[0] + c[1]*t+c[2]*t2+c[3]*t*t2+c[4]*t2*t2;

  dfdt = c[1]+2*c[2]*t+3*c[3]*t2+4*c[4]*t*t2;

  double x = dfdt/(1+f);

  return x;
}


//-----------------------------------------------------------------------------
void thermal_expansion(const char* Material) {
  // so far: W only

  TString material = Material;
  
  if (material == "W") {
    TF1* f = new TF1("f_W",thermal_expansion_tungsten,0,2500,0);

    f->SetTitle("Thermal expansion, tungsten");
    f->GetXaxis()->SetTitle("T, C");

    f->Draw();
  }
}
