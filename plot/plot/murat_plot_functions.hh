#ifndef murat_plot_murat_plot_functions_hh
#define murat_plot_murat_plot_functions_hh

					// from plot_ejected_proton_spectrum.cc

void  plot_ejected_proton_spectrum(const char* Variable          = "e", 
				   double      ProtonsPerCapture = 0.05,
				   const char* Opt               = "");

                                        // from plot_pbar_yield.cc

void  plot_pbar_yield     (float P0, float Theta);
void  plot_pbar_kinematics(float P0);
float get_pbar_yield      (float P0, float P, float Theta);
void  pbar_bgr            (double PBeam = 8.9, long int NEvents = 1000000);
void  pbar_xsec_norm      (double PBeam = 8.9, long int NEvents = 1000000);

					// from Stntuple/alg/pbar.f
// pa2pbarx_ returns E*d^3Sigma/dP^3
extern "C" float pa2pbarx_(float* P0, float* P, float* Theta);

#endif
