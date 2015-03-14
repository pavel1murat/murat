      integer function pbar_common_address()
c-----------------------------------------------------------------------
c *0002 Aug 03 1998 P.Murat: use /hepevt/ definition from STDHEP
c-----------------------------------------------------------------------
      include "plot/pbar_common.inc"
c      integer       pbar_common_address
      pointer :: pbar_common_address
c      real*8  x, var_address

      pbar_common_address = 1
      pbar_common_address = ifirst_time

c      x = var_address(ifirst_time)
c      pbar_common_address = x;
      
      return
      end
