//

#include "TMatrixD.h"
#include "TArrayD.h"

//-----------------------------------------------------------------------------
void test_matrix_001() {
  TMatrixD h(3,3);
  TArrayD data(9);

  h(0,0) = 1;
  h(0,1) = 2;
  h(0,2) = 3;
  h(1,0) = 10;
  h(1,1) = 10;
  h(1,2) = 11;
  h(2,0) = 20;
  h(2,1) = 20;
  h(2,2) = 21;

  TMatrixD b(TMatrixD::kInverted,h);

  h.Print("a");

  b.Print();

  TMatrixD m3 = h*b;

  v(0) = 2;
  v(1) = 2;
  v(2) = 2;


  v2 = m3*v;

  v2.Print();

  //  h.SetMatrixArray(data.GetArray());
  
}
