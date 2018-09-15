///////////////////////////////////////////////////////////////////////////////
// using:
// ------
// .L murat/scripts/print_genid_enum.C+
// print_genid_enum()
//------------------------------------------------------------------------------
#include "MCDataProducts/inc/GenId.hh"


namespace mu2e {
  
//-----------------------------------------------------------------------------
int print_genid_enum() {

  for (int i=GenId::unknown; i<GenId::lastEnum; i++) {

    GenId gi(i);
    
    cout << gi << endl;
  
  }
  return 0;
};

}
