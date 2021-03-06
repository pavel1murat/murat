///////////////////////////////////////////////////////////////////////////////
// digitized plots from J.A.Bistirlich et el, Phys. Rev. C 5, 1867–1883 (1972)
///////////////////////////////////////////////////////////////////////////////
// 87
// e,n
//

#include "string.h"
#include "TGraphErrors.h"
#include "Stntuple/alg/smooth.hh"
#include "TF1.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TCanvas.h"

using std::string;

//-----------------------------------------------------------------------------
// digitized figure 7b, efficiency folded in, v1
//-----------------------------------------------------------------------------
double data_7b_v1[] = { 
  50.9062, 5.59259e-11,
  51.7342, 9.83393e-13,
  52.4444,-1.27969e-10,
  53.2710, 2.04145e-11,
  54.2748,2.05875e-10,
  55.2217,2.16164e-11,
  56.1084,2.21628e-11,
  57.2316,2.28548e-11,
  58.3546,4.20310e-11,
  59.3018,-1.79196e-10,
  60.2462, 4.31965e-11,
  61.0724, 2.65516e-10,
  62.5   , 2.57685e-11,
  63.0826, 2.11302e-10,
  64.0902,-2.13214e-10,
  65.0338, 1.38567e-10,
  66.1568, 1.57744e-10,
  67.0982, 8.60724e-10,
  68.0431, 1.00918e-09,
  69.1081, 8.61963e-10,
  70.2332, 5.48424e-10,
  71.0625, 2.71672e-10,
  72.1805, 1.08567e-09,
  73.2459, 8.82996e-10,
  74.1895, 1.23478e-09,
  75.1959, 9.95103e-10,
  76.1991, 1.27298e-09,
  77.0262, 1.34743e-09,
  78.0924, 9.96888e-10,
  79.1541, 1.36723e-09,
  80.1017, 1.09055e-09,
  81.5   , 2.21863e-09,
  82.1061, 1.96054e-09,
  83.1686, 2.20149e-09,
  84.1763, 1.75849e-09,
  85.5   , 2.55378e-09,
  86.1802, 2.70242e-09,
  87.1907, 1.81580e-09,
  88.1226, 4.03448e-09,
  89.1395, 2.13123e-09,
  90.0818, 2.68634e-09,
  91.0176, 4.27656e-09,
  92.0903, 2.89090e-09,
  93.5   , 4.25928e-09,
  94.0934, 3.96422e-09,
  95.1036, 3.13305e-09,
  96.0463, 3.63271e-09,
  97.5   , 5.64803e-09,
  98.1686, 4.53974e-09,
  99.1019, 6.53662e-09,
  100.114, 5.37273e-09,
  101.5  , 5.87232e-09,
  102.237, 6.26128e-09,
  103.5  , 7.27835e-09,
  104.241, 7.11279e-09,
  105.243, 7.66794e-09,
  106.068, 7.96419e-09,
  107.070, 8.42692e-09,
  108.130, 9.11149e-09,
  109.259, 8.16949e-09,
  110.898, 9.92588e-09,
  111.081, 9.90813e-09,
  112.5  , 1.07959e-08,
  113.034, 9.41026e-09,
  114.5  , 1.01872e-08,
  115.099, 1.01879e-08,
  116.097, 1.12421e-08,
  117.047, 1.05772e-08,
  118.812, 1.19462e-08,
  119.941, 1.09117e-08,
  120.777, 9.63684e-09,
  121.777, 9.63684e-09,
  122.018, 9.65609e-09,
  123.019, 1.02297e-08,
  124.153, 8.64077e-09,
  125.095, 9.25133e-09,
  126.231, 7.20028e-09,
  127.182, 6.33211e-09,
  128.375, 4.61381e-09,
  129.377, 5.11350e-09,
  130.264, 4.96617e-09,
  131.161, 3.34011e-09,
  132.104, 3.93219e-09,
  133.067, 1.16015e-09,
  134.134, 6.43249e-10,
  135.379, 5.25205e-11,
  -1
};
//-----------------------------------------------------------------------------
// digitized figure 7b: efficiency folded in, second digitization attempt
//-----------------------------------------------------------------------------
double data_7b[] = { 
50.5133, 1.64222e-13,
51.5794, 1.89768e-11,
52.3499,-1.47143e-10,
53.1786, 1.94695e-11,
54.1851, 1.67653e-10,
55.0743,-7.23674e-11,
56.2585, 5.73866e-11,
57.8579, 2.42684e-12,
58.3313, 1.31962e-10,
59.1617,-2.37466e-10,
60.8187, 2.06665e-10,
61.7079,-1.48713e-11,
62.0102, 2.07340e-10,
63.0102, 2.07340e-10,
64.0776,-1.98983e-10,
65.0836, 9.70737e-11,
66.1497, 1.15886e-10,
67.0166, 8.37024e-10,
68.0131, 9.85207e-10,
69.0105, 8.37644e-10,
70.0575, 5.60710e-10,
71.0063, 1.91319e-10,
72.0697, 1.06040e-09,
73.1957, 8.57424e-10,
74.1424, 1.17195e-09,
75.2093, 9.31981e-10,
76.0374, 1.28344e-09,
77.1034, 1.32073e-09,
78.1930, 9.69807e-10,
79.1173, 1.32135e-09,
80.0065, 1.06285e-09,
81.0098, 2.20918e-09,
82.0769, 1.91376e-09,
83.3199, 2.17292e-09,
84.1505, 1.72955e-09,
85.1772, 2.52463e-09,
86.2205, 2.69137e-09,
87.1711, 1.76745e-09,
88.1118, 3.98584e-09,
89.1840, 2.06382e-09,
90.0706, 2.65559e-09,
91.0132, 4.24552e-09,
92.0245, 2.85952e-09,
93.0863, 4.24615e-09,
94.1759, 3.89523e-09,
95.1036, 3.13773e-09,
96.1090, 3.61862e-09,
97.1688, 5.63373e-09,
98.1794, 4.46953e-09,
99.0020, 6.54001e-09,
100.190, 5.33891e-09,
101.196, 5.76436e-09,
102.164, 6.28215e-09,
103.168, 7.24364e-09,
104.272, 7.09616e-09,
105.218, 7.65098e-09,
106.224, 7.91007e-09,
107.170, 8.33550e-09,
108.175, 9.09366e-09,
109.303, 8.15131e-09,
110.008, 9.88904e-09,
111.134, 9.87090e-09,
112.197, 1.07770e-08,
113.090, 9.35395e-09,
114.176, 1.01675e-08,
115.160, 1.01679e-08,
116.105, 1.12587e-08,
117.114, 1.05012e-08,
118.000, 1.11854e-08,
119.005, 1.19066e-08,
120.015, 1.08903e-08,
121.148, 9.61512e-09,
122.151, 9.63400e-09,
123.216, 1.01704e-08,
124.227, 8.61801e-09,
125.232, 9.22830e-09,
126.187, 7.13988e-09,
127.433, 6.30848e-09,
128.505, 4.58978e-09,
129.451, 5.08915e-09,
130.399, 4.94156e-09,
131.352, 3.33373e-09,
132.298, 3.87007e-09,
133.076, 1.17162e-09,
134.263, 6.17458e-10,
135.508, 4.48327e-11,
  -1
};
//-----------------------------------------------------------------------------
// digitized figure 7a, efficiency divided out, first attempt
// turns out, the background is subtracted here, this is the plot to use
//-----------------------------------------------------------------------------
double data_7a_v1[] = {
  60.6693, 9.39130e-05,
  61.5551, 3.13043e-05,
  62.6181, 6.95652e-05,
  63.6220,-2.95652e-05,
  64.7441, 1.65217e-05,
  65.8071, 1.13043e-05,
  66.6339, 0.000140870,
  67.5197, 0.000155652,
  68.4646, 0.000119130,
  69.4094, 7.04348e-05,
  70.5906, 3.21739e-05,
  71.5354, 0.000117391,
  72.5984, 0.000102609,
  73.6614, 0.000109565,
  74.6063, 0.000118261,
  75.4921, 8.95652e-05,
  76.3780, 9.91304e-05,
  77.3228, 0.000121739,
  78.5039, 9.39130e-05,
  79.5079, 0.000108696,
  80.5709, 0.000176522,
  81.6929, 0.000168696,
  82.5787, 0.000160000,
  83.4055, 0.000144348,
  84.4685, 0.000202609,
  85.4724, 0.000202609,
  86.4764, 0.000133913,
  87.5394, 0.000306087,
  88.4843, 0.000161739,
  89.5472, 0.000197391,
  90.5512, 0.000306087,
  91.4961, 0.000200870,
  92.4409, 0.000289565,
  93.5630, 0.000252174,
  94.5669, 0.000206957,
  95.4528, 0.000225217,
  96.3386, 0.000353043,
  97.5197, 0.000276522,
  98.5827, 0.000384348,
  99.8819, 0.000310435,
  100.591, 0.000326087,
  101.476, 0.000354783,
  102.480, 0.000401739,
  103.543, 0.000391304,
  104.488, 0.000401739,
  105.374, 0.000408696,
  106.496, 0.000433913,
  107.441, 0.000465217,
  108.327, 0.000405217,
  109.272, 0.000489565,
  110.335, 0.000478261,
  111.457, 0.000516522,
  112.579, 0.000440000,
  113.524, 0.000476522,
  114.528, 0.000476522,
  115.354, 0.000521739,
  116.240, 0.000470435,
  117.244, 0.000503478,
  118.366, 0.000531304,
  119.606, 0.000478261,
  120.728, 0.000417391,
  121.728, 0.000417391,
  122.264, 0.000449565,
  123.327, 0.000383478,
  124.449, 0.000407826,
  125.394, 0.000318261,
  126.457, 0.000281739,
  127.520, 0.000205217,
  128.524, 0.000231304,
  129.469, 0.000220000,
  130.472, 0.000146957,
  131.417, 0.000202609,
  132.421, 4.17391e-05,
  133.543, 1.30435e-05,
  134.606, 0.0,
  -1
};
//-----------------------------------------------------------------------------
// digitized figure 7a, efficiency divided out, second attempt
// background is subtracted, and this is the plot to use
//-----------------------------------------------------------------------------
double data_7a[] = {
60.6034, 0.000100222,
61.6116, 3.86105e-05,
62.5537, 7.88520e-05,
63.6824,-1.95622e-05,
64.8605, 2.58181e-05,
65.8058, 1.89843e-05,
66.5649, 0.000146527,
67.5087, 0.000161947,
68.3970, 0.000125155,
69.4632, 7.80947e-05,
70.5288, 4.04496e-05,
71.4680, 0.000123487,
72.5910, 0.000109808,
73.5943, 0.000118382,
74.6570, 0.000123532,
75.4855, 9.78667e-05,
76.3707, 0.000106438,
77.3140, 0.000130418,
78.3790, 0.000101332,
79.5590, 0.000116755,
80.5584, 0.000183531,
81.7992, 0.000174133,
82.5674, 0.000167297,
83.4544, 0.000149335,
84.3362, 0.000209262,
85.4582, 0.000210134,
86.4078, 0.000140818,
87.4002, 0.000310304,
88.4728, 0.000169093,
89.4744, 0.000203344,
90.4709, 0.000312059,
91.4819, 0.000208508,
92.3620, 0.000293256,
93.5455, 0.000258181,
94.4933, 0.000214542,
95.3780, 0.000230817,
96.2553, 0.000356650,
97.5005, 0.000283059,
98.6154, 0.000387496,
99.6833, 0.000315615,
100.568, 0.000330178,
101.629, 0.000361006,
102.453, 0.000404669,
103.457, 0.000395269,
104.402, 0.000405553,
105.405, 0.000412415,
106.466, 0.000438107,
107.409, 0.000468078,
108.299, 0.000409032,
109.238, 0.000492926,
110.302, 0.000481814,
111.421, 0.000518634,
112.608, 0.000445042,
113.491, 0.000479292,
114.436, 0.000479305,
115.260, 0.000525536,
116.208, 0.000474195,
117.210, 0.000505022,
118.389, 0.000533284,
119.396, 0.000481944,
120.523, 0.000422045,
121.290, 0.000422912,
122.233, 0.000452027,
123.301, 0.000386992,
124.480, 0.000410975,
125.372, 0.000323684,
126.496, 0.000287752,
127.446, 0.000212445,
128.508, 0.000237281,
129.571, 0.000226170,
130.344, 0.000155996,
131.462, 0.000210790,
132.536, 5.07492e-05,
133.601, 2.33751e-05,
134.843, 2.83059e-07,
  -1
};
//-----------------------------------------------------------------------------
// efficiency: digitized figure 3 from the paper
//-----------------------------------------------------------------------------
double data_eff[] = {
  58.8721,2.31616e-06,
  62.2951,4.01606e-06,
  65.9811,5.78126e-06,
  69.8423,7.54634e-06,
  73.3527,9.21345e-06,
  76.5992,1.06189e-05,
  79.3201,1.19592e-05,
  82.3896,1.29719e-05,
  86.4236,1.42459e-05,
  90.7203,1.55197e-05,
  95.7181,1.69239e-05,
  100.540,1.82628e-05,
  105.100,1.96019e-05,
  109.834,2.06463e-05,
  115.092,2.14939e-05,
  120.435,2.19487e-05,
  124.901,2.19785e-05,
  128.578,2.17142e-05,
  130.852,2.12217e-05,
  133.564,2.05980e-05,
  136.536,1.96141e-05,
  140.034,1.85643e-05,
  143.707,1.74490e-05,
  148.429,1.61365e-05,
  151.840,1.52832e-05,
  155.864,1.44295e-05,
  158.838,1.38384e-05,
  161.376,1.34112e-05,
 
  // 60.3203,2.75507e-06,
  // 70.2947,7.57709e-06,
  // 80.1892,1.15990e-05,
  // 90.2580,1.29546e-05,
  // 100.153,1.69765e-05,
  // 110.050,2.00651e-05,
  // 119.867,2.25535e-05,
  // 129.862,2.15757e-05,
  // 139.949,1.75983e-05,
  // 150.114,1.48209e-05,
  -1
};

//-----------------------------------------------------------------------------
TH1F*         h_fig7b;
TH1F*         h_fig7a;
TH1F*         h_665;
TH1F*         h_fig7b_eff;
TH1F*         h_fig7a_fit(0);

TGraphErrors* gr_eff;
smooth*       sgr_eff;
//------------------------------------------------------------------------------
// Ivano's parameterization from mu2e-665 ? from PiCaptureEffects
//-----------------------------------------------------------------------------
double rpc_ivano(double E) {

  constexpr double emax  = 138.2;
  constexpr double alpha =   2.691;
  constexpr double gamma =   1.051;
  constexpr double tau   =   8.043;
  constexpr double c0    =   2.741;
  constexpr double c1    =  -0.005;

  double de = emax-E;

  if (de < 0) return 0;

  return pow(emax-E,alpha) * exp(-(emax-gamma*E)/tau) * (c0 + c1*E)/10000.;
}

//-----------------------------------------------------------------------------
void plot_ivano() {
  
  TH1F* h = new TH1F("h1","RPC mu2e-665 ?",200,0,200);

  for (int i=0; i<200; i++) {
    double e = 0.5+i;

    double w = rpc_ivano(e);

    h->Fill(e,w);
  }

  h->Draw("h");
}

//-----------------------------------------------------------------------------
void make_eff() {

  double e[1000], eff[1000], ee[1000], eeff[1000];

  int np(0);
  
  for (int i=0; data_eff[2*i]>0; i++) {
    e   [np] = data_eff[2*i];
    eff [np] = data_eff[2*i+1];
    ee  [np] = 0;
    eeff[np] = eff[np]/10.;
    np++;
  }

  gr_eff = new TGraphErrors(np,e,eff,ee,eeff);

  TCanvas* c_eff = new TCanvas("c_eff","e_eff",1000,600);

  //  gr_eff->Fit("pol3","w");
  gr_eff->Draw();

  sgr_eff = new smooth(gr_eff,60,160);
  sgr_eff->fFunc->Draw("same");
}


//-----------------------------------------------------------------------------
// main function to call
//-----------------------------------------------------------------------------
void rpc_bistirlich() {

  double e[1000], w[1000], ee[1000], ew[1000];

  /* efficiency fit from Phys. Rev. Lett. 25, 689 (1970)
    p0                        =  8.65145e-06   +/-   2.17037e-05 
    p1                        = -7.30012e-07   +/-   6.69052e-07 
    p2                        =   1.4327e-08   +/-   6.58486e-09 
    p3                        = -6.13523e-11   +/-   2.07986e-11 
  */

  double fp[4] = {8.65145e-06, -7.30012e-07, 1.4327e-08, -6.13523e-11};

  int n(0);

  h_fig7b  = new TH1F("h_fig7b","RPC photon spectrum for ^{24}Mg",150,0,150);
  h_fig7a  = new TH1F("h_fig7a","RPC photon spectrum for ^{24}Mg, corrected",150,0,150);
  h_665  = new TH1F("h_665","Ivano\'s fit from Mu2e Offline",150,0,150);
  h_fig7b_eff  = new TH1F("h_fig7b_eff","RPC photon spectrum for ^{24}Mg, eff div off",150,0,150);


  make_eff();
  
  for (int i=0; data_7b[2*i]>0; i++) {
    e[n] = data_7b[2*i];
    w[n] = data_7b[2*i+1]*1.e9;

    ee[n] = 0.;
    ew[n] = 0.;

    
    int bin = h_fig7b->GetXaxis()->FindBin(e[n]);

    h_fig7b->SetBinContent(bin,w[n]);

    double w_ivano = rpc_ivano(e[n]);

    h_665->SetBinContent(bin,w_ivano);

    double x = e[n];

    //    double eff = (fp[0]+fp[1]*x+fp[2]*x*x+fp[3]*x*x*x);

    double eff = sgr_eff->GetFunc()->Eval(x);

    if (x < 60) eff = 1.e10;
    
    printf(" e = %10.4f w_ivano = %12.4f eff = %12.5e\n",e[n],w_ivano,eff);

    h_fig7b_eff->SetBinContent(bin,w[n]/eff);

    n += 1;
    
  }

  for (int i=0; data_7a[2*i]>0; i++) {
    double energy = data_7a[2*i];
    double weight = data_7a[2*i+1]*1.e4;

    int bin = h_fig7a->GetXaxis()->FindBin(energy);

    h_fig7a->SetBinContent(bin,weight);
  }

  //  TGraphErrors* gr = new TGraphErrors(n,e,w,ee,ew);

  //  gr->Draw("ALP");

  TCanvas* c_dat = new TCanvas("c_dat","e_dat",1100,700);

  h_fig7b->SetStats(0);
  h_fig7b->GetXaxis()->SetTitle("E(photon), MeV");
  h_fig7b->DrawNormalized();

  h_fig7a->SetStats(0);
  h_fig7a->SetLineColor(kBlue+2);
  h_fig7a->SetFillStyle(3002);
  h_fig7a->SetFillColor(kBlue+2);
  h_fig7a->DrawNormalized("same",1);

  h_665->SetStats(0);
  h_665->SetLineColor(2);
  h_665->DrawNormalized("same",1);

  h_fig7b_eff->SetStats(0);
  h_fig7b_eff->SetLineColor(kRed-2);
  h_fig7b_eff->SetFillColor(kRed-2);
  h_fig7b_eff->SetFillStyle(3003);
  h_fig7b_eff->DrawNormalized("same",1);

  TLegend* leg = new TLegend(0.15,0.7,0.6,0.8);
  leg->SetLineWidth(0);

  leg->AddEntry(h_fig7b    ,"Phys Rev C5 1867 figure 7b");
  leg->AddEntry(h_fig7a    ,"Phys Rev C5 1867 figure 7a");
  leg->AddEntry(h_665      ,"fit from Mu2e Offline, mu2e-665");
  leg->AddEntry(h_fig7b_eff,"Phys Rev C5 1867 figure 7b/eff");

  leg->Draw();
  
}

//-----------------------------------------------------------------------------
// fit function
//-----------------------------------------------------------------------------
double f_rpc(double* X, double* P) {

  double f(0.);

  double emax  = P[0];
  double alpha = P[1];
  double gamma = P[2];
  double tau   = P[3];
  double c0    = P[4];
  double c1    = P[5];

  double e  = X[0];
  double de = emax-e;

  if (de > 0) {
    f = pow(de,alpha) * exp(-(emax-gamma*e)/tau) * (c0 + c1*e);
  }

  return f;
}

//-----------------------------------------------------------------------------
// fit results, function normalized to 1
//
// FCN=8.82901 FROM HESSE     STATUS=NOT POSDEF     40 CALLS         993 TOTAL
//                     EDM=5.49549e-08    STRATEGY= 1      ERR MATRIX NOT POS-DEF
//  EXT PARAMETER                APPROXIMATE        STEP         FIRST   
//  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
//   1  p0           1.34530e+02   4.25845e-01   1.17475e-04  -8.14601e-04
//   2  p1           1.29931e+00   2.39423e-01   1.38819e-06   2.41928e-05
//   3  p2           9.28705e-01   6.53886e-02   1.96419e-06  -5.20820e-03
//   4  p3           9.39676e+00   1.01635e+00   5.10886e-05   5.65493e-05
//   5  p4           4.77611e-02   8.86851e-03   2.33297e-07  -6.48587e-02
//   6  p5          -3.26349e-04   6.48155e-05   1.98579e-09  -8.57054e+00
//
//-----------------------------------------------------------------------------
void fit_h1_fig7a() {

  TF1* f = new TF1("f1",f_rpc,60,140,6);

  f->SetParameter(0,138.2);
  f->SetParameter(1,2.691);
  f->SetParameter(2,1.051);
  f->SetParameter(3,8.043);
  f->SetParameter(4,0.016);
  f->SetParameter(5,-1.1e-4);


  TCanvas* c_fit_h_fig7a = new TCanvas("c_fit_h_fig7a","fit h_fig7a",1100,700);

  if (h_fig7a_fit) delete h_fig7a_fit;
  
  h_fig7a_fit = (TH1F*) h_fig7a->Clone("h_fig7a_fit");

  int nb      = h_fig7a_fit->GetNbinsX();
  double norm = h_fig7a->Integral();
  
  for (int i=1; i<=nb; i++) {
    h_fig7a_fit->SetBinContent(i,h_fig7a->GetBinContent(i)/norm);
    h_fig7a_fit->SetBinError  (i,h_fig7a->GetBinError  (i)/norm);
  }

  h_fig7a_fit->Fit(f);
}


//-----------------------------------------------------------------------------
void plot_rpc_spectrum(int NEvents = 100000) {

  double fint  {0.98445};  // integral

  double emax  {  134.530  };
  double alpha {  1.29931  };
  double gamma {  0.928705 };
  double tau   {  9.39676  };
  double c0    {  4.77611e-02 /fint };
  double c1    { -3.26349e-04 /fint };
  
  TF1* f = new TF1("f1",f_rpc, 0,140,6);

  f->SetParameter(0, emax);
  f->SetParameter(1,alpha);
  f->SetParameter(2,gamma);
  f->SetParameter(3,  tau);
  f->SetParameter(4,   c0);
  f->SetParameter(5,   c1);

  TH1F* h = new TH1F("h1","h1",1400,0,140);

  TRandom3 rn3;

  for (int i=0; i<NEvents; i++) {
    double e = rn3.Rndm(i)*emax;
    double w = f->Eval(e)*emax;  // this gives correct normalization
    h->Fill(e,w);
  }

  h->Draw();

  printf(" integral: %12.5e\n",f->Integral(0,emax));
}
