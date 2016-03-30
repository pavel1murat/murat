//

//

void lk(int thisIndex, int& ix, int& iy) {         

  ix = 0;
  iy = 0;

  if (thisIndex==0) return ; 

  int nRing = int(0.5*sqrt(thisIndex) + 0.5);
  
  int nSeg  = (thisIndex - (2*nRing-1)*(2*nRing-1)) / (2*nRing);
  int nPos  = (thisIndex - (2*nRing-1)*(2*nRing-1)) % (2*nRing);
  
  if (nSeg==0) {ix = -nRing+nPos ; iy =  nRing     ; }
  if (nSeg==1) {ix =  nRing      ; iy =  nRing-nPos; }
  if (nSeg==2) {ix =  nRing-nPos ; iy = -nRing     ; }
  else         {ix = -nRing      ; iy = -nRing+nPos; }
} 


int test_square_mapper(int N) {

  int ix, iy;
  
  for (int i=0; i<N; i++) {
    lk(i,ix,iy);
    printf(" i, ix, iy = %5i %5i %5i\n",i, ix, iy);

  }
}
