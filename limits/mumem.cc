///////////////////////////////////////////////////////////////////////////////
// use:
// mumem* ana = new mumem(1);
// ana->
///////////////////////////////////////////////////////////////////////////////
const char* FiguresDir         = "/projects/mu2e/talks/2021-03-11-su2020-backgrounds-collab/figures";

#include "murat/limits/mumem.hh"

namespace murat {
//-----------------------------------------------------------------------------
mumem::mumem(int Mode, const char* Name, const char* Title) : su2020::analysis(Name,Title) {
  init_channels(Mode) ;
  fPMin = 103.85;
  fPMax = 105.10;
  fTMin =  700;
  fTMax = 1700;
}

//-----------------------------------------------------------------------------
// channel names - predefined
// 1. "CE",
// 2. "Cosmics",
// 3. "DIO",
// 4. "PbarTOT" = "PbarAnni" + "PbarRPCe" + "PbarPRCi",
// 4. "PbarRPC" = "PbarRPCe" + "PbarPRCi",
// 5. "RPC"     = "RPCe" + "RPCi"
//
// any change in the names needs to be propagated to the channel initialization code
//-----------------------------------------------------------------------------
int mumem::init_channels(int Mode) {
					// N(POT) for Mu2e Run-I 1-batch and 2-batch modes
  float npot_1b = 2.86e19; 
  float npot_2b = 9.03e18;

  if (Mode == 1) {			// tan(dip) < 1.00;
    fSFSignal  = 1.00;
    fSFCosmics = 1.00;
  }
  else if (Mode == 2) {			// 30% degradation of the light yield
     fSFSignal  = 1.;
     fSFCosmics = 0.21/0.036;
  }
  // else if (Mode == 3) {			// tan(dip) < 0.90; approx
  //   fSFSignal  = 0.90;
  //   fSFCosmics = 0.33;
  // }
//-----------------------------------------------------------------------------
// 1. CE signal. For CE, specify signal in units of acceptance ?
//-----------------------------------------------------------------------------
  // Rmue comes independently
  float Rmue = 1.e-16;
  su2020::channel* ce = new su2020::channel("CE",npot_1b,npot_2b,fSFSignal);
  SetSignal(ce,Rmue);
//-----------------------------------------------------------------------------
// 2. cosmics: fit , NPOT numbers are not used
//-----------------------------------------------------------------------------
  su2020::channel* cosmics = new su2020::channel("Cosmics",npot_1b,npot_2b,fSFCosmics);
  AddBgrChannel(cosmics);
//-----------------------------------------------------------------------------
// 3. DIO
//-----------------------------------------------------------------------------
  su2020::channel* dio = new su2020::channel("DIO",npot_1b,npot_2b,fSFSignal);
  AddBgrChannel(dio);
//-----------------------------------------------------------------------------
// 4. pbars
//-----------------------------------------------------------------------------
  su2020::channel* pbar = new su2020::channel("PbarTOT",npot_1b,npot_2b,fSFSignal);
  AddBgrChannel(pbar);
//-----------------------------------------------------------------------------
// 5. RPC
//-----------------------------------------------------------------------------
  su2020::channel* rpc = new su2020::channel("RPC",npot_1b,npot_2b,fSFSignal);
  AddBgrChannel(rpc);

  return 0;
}

}
