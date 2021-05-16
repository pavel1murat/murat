///////////////////////////////////////////////////////////////////////////////
// use:
// mumem* ana = new mumem(1);
// ana->
///////////////////////////////////////////////////////////////////////////////
const char* FiguresDir         = "/projects/mu2e/talks/2021-03-11-su2020-backgrounds-collab/figures";

#include "murat/limits/mumem.hh"

namespace murat {
//-----------------------------------------------------------------------------
void mumem::plot(int Figure, int Print = 0) {

  if (Figure == 1) {
//-----------------------------------------------------------------------------
// momentum distribution : 
//-----------------------------------------------------------------------------
    int nbgr = GetNBgrChannels();

    plot_data_t  p(nbgr+1);

    su2020::channel* ce  = GetSignal();
    
    TH1D* h1             = ce->CreateMomHist(fTMin,fTMax);
    h1->SetName("CE");
    
    h1->Scale(fRmue);
    p.hd[0]              = hist_data_t(h1);
    p.hd[0].fLabel       = "CE";
    p.hd[0].fLineColor   = kRed+1;
    p.hd[0].fMarkerColor = kRed+1;
    p.hd[0].fMarkerStyle = 20;
    p.hd[0].fMarkerSize  = 0.8;

    int col[4] = { EColor::kBlue-1, EColor::kGreen+2, EColor::kMagenta, EColor::kRed+2 };
    
    for (int i=0; i<nbgr; i++) {
      int ih = i+1;			// signal first

      su2020::channel* bgr   = GetBgrChannel(i);
      TH1D* h2               = bgr->CreateMomHist(fTMin,fTMax);
      h2->SetName(bgr->GetName());
      
      p.hd[ih]                = hist_data_t(h2);

      p.hd[ih].fLabel         = bgr->GetName();
      p.hd[ih].fLineColor     = col[i];
      p.hd[ih].fMarkerColor   = col[i];
      p.hd[ih].fMarkerStyle   = 20;
      p.hd[ih].fMarkerSize    = 0.8;
    }
//-----------------------------------------------------------------------------
// common for the plot opt stat: 1110001
// by default, assume that all histograms have the same binning,
// so the rebinning factor can be specified just once per plot
//-----------------------------------------------------------------------------    
					// Y-scale and rebinning
    p.fYLogScale     = 0;
    p.fRebin         = 1;

    p.fXMin          = 102.5;
    p.fXMax          = 105.199;
    p.fXAxisTitle    = "e^{-} momentum, MeV/c";

    p.fYMin          = 1e-5;
    p.fYMax          = 0.30;
    
    p.fStatBoxYMin   = 0.83;         // up to 0.90, dY=0.07
    p.fOptStat       = 1000001;

    p.fLabel         = Form("SU2020: R = %12.5e, T0 in [%5.0f,%5.0f] ns",fRmue,fTMin,fTMax);
    p.fLabelFontSize = 0.04;

    p.fLegendXMin    = 0.15; p.fLegendXMax  = 0.35; p.fLegendYMin  = 0.50; p.fLegendYMax  = 0.65;

    p.fCanvasName    = Form("Figure_%04i",Figure);
    p.fName          = Form("figure_%05i_su2020_total_p",Figure);

    plot_hist_1d(&p,-1);

    if (Print == 1) p.print();
  }
  
  if (Figure == 2) {
//-----------------------------------------------------------------------------
// timing distribution : 103.5-105.1 MeV
//-----------------------------------------------------------------------------
    int nbgr = GetNBgrChannels();

    plot_data_t  p(nbgr+1);

    su2020::channel* ce  = GetSignal();
    
    TH1D* h1             = ce->CreateTimeHist(fPMin,fPMax);
    h1->SetName("CE");
    
    h1->Scale(fRmue);
    p.hd[0]              = hist_data_t(h1);
    p.hd[0].fLabel       = "CE";
    p.hd[0].fLineColor   = kRed+1;
    p.hd[0].fMarkerColor = kRed+1;
    p.hd[0].fMarkerStyle = 20;
    p.hd[0].fMarkerSize  = 0.8;

    int col[4] = { EColor::kBlue-1, EColor::kGreen+2, EColor::kMagenta, EColor::kRed+2 };
    
    for (int i=0; i<nbgr; i++) {
      int ih = i+1;			// signal first

      su2020::channel* bgr   = GetBgrChannel(i);
      TH1D* h2               = bgr->CreateTimeHist(fPMin,fPMax);
      h2->SetName(bgr->GetName());
      
      p.hd[ih]                = hist_data_t(h2);

      p.hd[ih].fLabel         = bgr->GetName();
      p.hd[ih].fLineColor     = col[i];
      p.hd[ih].fMarkerColor   = col[i];
      p.hd[ih].fMarkerStyle   = 20;
      p.hd[ih].fMarkerSize    = 0.8;
    }
//-----------------------------------------------------------------------------
// common for the plot opt stat: 1110001
// by default, assume that all histograms have the same binning,
// so the rebinning factor can be specified just once per plot
//-----------------------------------------------------------------------------    
					// Y-scale and rebinning
    p.fYLogScale     = 1;
    p.fRebin         = 1;

    p.fXMin          = 500;
    p.fXMax          = 1700;
    p.fXAxisTitle    = "e^{-} time, ns";

    p.fYMin          = 1e-5;
    p.fYMax          = 1e3;
    
    p.fStatBoxYMin   = 0.83;         // up to 0.90, dY=0.07
    p.fOptStat       = 1000001;

    p.fLabel         = Form("SU2020: R = %12.5e, P in [%7.2f,%7.2f] ns",fRmue,fPMin,fPMax);
    p.fLabelFontSize = 0.04;

    p.fLegendXMin    = 0.15; p.fLegendXMax  = 0.35; p.fLegendYMin  = 0.50; p.fLegendYMax  = 0.65;

    p.fCanvasName    = Form("Figure_%04i",Figure);
    p.fName          = Form("figure_%05i_su2020_total_time",Figure);

    plot_hist_1d(&p,-1);

    if (Print == 1) p.print();
  }
}
}
