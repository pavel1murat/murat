c-----------------------------------------------------------------------
c pbar.f  - provides dN/dp/d(tet)   - #pbar/(GeV/c)/rad   in proton gold collision at 8.9 GeV/c .
c It is ready to simulate Mu2e events.
c 
c To compare with figure
c  1) put plab =10. 
c 2) to get ed3g/dp3  multiply result on 1.539e6 mub *E/p**2/(2 pi sin(tet)
c 
c Good luck.
c Sergei
c

c-----------------------------------------------------------------------
      block data
      include "plot/pbar_common.inc"
c
      data ifirst_time/0/
      data p2max/1.3/
      data p2c  /1.1/
c-----------------------------------------------------------------------
c was used for Mu2e simulation - 2011-2012
c      p2max=1.3
c new data shows that better use p2max=2 (p > 1.2 GeV/c only)
c      p2max=2.
c      p2c=1.1
c-----------------------------------------------------------------------
      end

      subroutine set_p2max(p)
      include "plot/pbar_common.inc"
      p2max = p
      end
c-----------------------------------------------------------------------
      function pbar_yield(p0,p,tet)

      include "plot/pbar_common.inc"

c     Probability of pbar production with momentum P (GeV/c)
c     and angle TET (rad) in proton-nucleus collision.

c     Parametrization is valid for tantalum and nucleus with
c     close atomic weight. Atomic weight dependece is about A^(1/3).
      
c     Pbar yieled from Mu2e target was simulated for
c     tantalum target and 8.89 GeV/c proton beam  

c     Tantalum inelstic cross section = 1.539e6 mubarn

      PARAMETER (PI=3.141592)
      PARAMETER (PM=0.938272)

c      print *, "privet"

      ss=sin(tet)
c      plab=8.89
      plab=p0
      e=sqrt(p**2+pm**2)
      conv=2.*pi*ss*p**2/e/1.539e6
      val =pa2pbarx(p,tet,plab)*conv
      pbar_yield=val;
      return
      end

      FUNCTION PA2PBARX(P,T,PLAB)
C...........................................
C   INPUT: 
C          PLAB - PRIMARY PROTON MOMENTUM (GEV/C)
C          P    - SECONDARY ANTIPROTON MOMENTUM (GEV/C)
C          T    - SECONDARY ANTIPROTON ANGLE (RAD)
C   OUTPUT:
C     Ed^3Sigma/dP^3 - INVARIANT CROSS SECTION (mub/GeV^2/c^3)
C...........................................
      include "plot/pbar_common.inc"

      PARAMETER (JP=7)
      PARAMETER (PM=0.938272)
   
      ss=sin(t)
      xlq=(sqrt(p**2+pm**2)-p*cos(t))/pm
      PT=P*SS
      ELAB=SQRT(PLAB**2+PM**2)
      SQS=SQRT(2.*PM*(ELAB+PM))
      xfa=xfam(plab,p,t,jp)
      if(plab.le.12.) then
         f1=TANNG1(JP,XFA,PT,SQS)
      else
         f1=fsmlan1(plab,t,p)
      endif
c-----------------------------------------------------------------------
c was used for Mu2e simulation - 2011-2012
c      p2max=1.3
c new data shows that better use p2max=2 (p > 1.2 GeV/c only)
c      p2max=2.
c      p2c=1.1
c-----------------------------------------------------------------------
      p0max=plmax(plab,t,jp) 
      p1max=max(p2max,p0max)
      if(p.ge.p2c) then
         xlq0=(sqrt(p2c**2+pm**2)-p2c*cos(t))/pm
         f20=flrgan1(plab,t,xlq0)
         f2=(1.-(p-p2c)/(p1max-p2c))**5*f20
      else
         f2=flrgan1(plab,t,xlq)
      endif
      pa2pbarx=f1+f2
      if(pa2pbarx.le.0.) pa2pbarx=0.
      return
      end

      function flrgan1(plab,t,xlq)
C     large angle approximation
      PARAMETER (PM=0.938272)
      ss=sin(t)
      if(xlq.ge.ss) then 
         dsqs=2.*pm*(plab-7.*pm)
         if(plab.le.12.9) then
            x01=0.023*(dsqs-2.93)
            x02=0.5598*exp(-.664*dsqs)
         else
            x01=0.15
            x02=0.0002
            x01=0.2059
            x02=0.00021
         endif
         c2=4.5
         c1=x02
         x1=xlq-x01
         if(x1.ge.0) then
            fc=x1**c2/(x1**c2+c1)
         else
            fc=0.
         endif
         flrgan1=9674.7*exp(-10.127*xlq)*(1.-xlq/181.)**9*fc
      else
      flrgan1=0.
      endif
      return
      end
      function fsmlan1(plab,t,p)
C     Parametrization of  antiproton production cross section > 12 GeV/c
C     R.P. Duperray et al, Phys. Rev. D68, 094017 (2003)   
      PARAMETER (PM=0.938272)
      PARAMETER (JP=7)
      pm2=pm**2
      e0=sqrt(plab**2+pm**2)
      S=2.*PM*E0+2.*PM2
      SQS=SQRT(S)
      s0=45.*(1.+0.016*sin(5.3))
      csin=s0*(1.-0.62*exp(-e0/200.)*sin(10.9/e0**0.28))
      XR=XRAD(plab,P,T,JP)
      ptt=p*sin(t)
      fsmlan1=fpbar1(xr,ptt,sqs)*csin*1.e3*29.
      return
      end



      function fpbar1(xr,pt,sqs)
C     Parametrization of  antiproton production cross section > 12 GeV/c
C     R.P. Duperray et al, Phys. Rev. D68, 094017 (2003) 
      if(xr.le.1) then
      c1=3.461
      c2=4.340
      c3=0.007855
      c4=0.5121
      c5=3.6620
c  Duperray fit
c          c31=c3*sqs**c4
c  corrected fit
      if(sqs.le.201.) then
         c31=c3*6.1509**(c4-0.35)*sqs**0.36
         c5=c5*(6.1509/sqs)**0.095
      else 
         c31=0.07162
         c5=2.6003
      endif
      c6=0.02307
      c7=3.2540
      f1=(1.-xr)**c1*exp(-c2*xr)
      f2=c31*exp(-c5*pt)+c6*exp(-c7*pt*pt)
      fpbar1=f1*f2
      else
      fpbar1=0.
      endif
      return
      end  
      FUNCTION XRAD(PO,P,TET,JP)
C     P&TET TO RADIAL SCALING XR
      DIMENSION PM(7)
      DIMENSION SMX(7)
      DATA PM/0.938272,0.9396,2*0.1396,2*0.493677,0.938272/
      DATA SMX/0.93827231,1.079,1.878,2.016,2.053,2.370,2.814/
      
      EO=SQRT(PO**2+PM(1)**2)
      S=2.*PM(1)*EO+PM(1)**2+PM(1)**2
      SQS=SQRT(S)
      emax=(s+pm(jp)**2-smx(jp)**2)/2./sqs
      GACM=(EO+PM(1))/SQS
      GABCM=PO/SQS
      EPI=SQRT(P**2+pm(jp)**2 )
      PMP=P*COS(TET)
      ecm=gacm*epi-gabcm*pmp      
      XRAD=ecm/emax
      RETURN
      END
      FUNCTION XFAM(PO,P,TET,JP)
C     P&TET TO FEINMAN SCALING XF (2. PZCM/SQS !!!, NOT PMAX !!!)
      DIMENSION PM(7)
      DIMENSION SMX(7)
      DATA PM /0.938272,0.9396,2*0.1396,2*0.493677,0.938272/
      DATA SMX/0.93827231,1.079,1.878,2.016,2.053,2.370,2.814/
      
      EO=SQRT(PO**2+PM(1)**2)
      S=2.*PM(1)*EO+PM(1)**2+PM(1)**2
      SQS=SQRT(S)
      emax=(s+pm(jp)**2-smx(jp)**2)/2./sqs
      GACM=(EO+PM(1))/SQS
      GABCM=PO/SQS
      EPI=SQRT(P**2+pm(jp)**2 )
      PMP=P*COS(TET)
      pzcm=-gabcm*epi+gacm*pmp
      XFAM=2.*pzcm/sqs
      RETURN
      END
      FUNCTION TAN0(JT,XF,PT,SQS)
C     Tan and Ng parametrization of antiproton production cross section < 12 GeV/c
      DIMENSION A(12,5),G00(5),SMX(5),PM(5),ST(5)
      DATA A/
     * 143.,9.54,.7,2.72,3.53,3.48,1.88,.212,1.64E-4,.779,4.70,7.82,
     * 150.,13.8,.5,3.38,4.20,4.82,1.07,.984,.873,.455,1.71,2.01,
     * 2.61,5.36,.7,2.49,1.96,1.03E-3,1.67,.263,8.24E3,6.41,4.72,6.09,
     * 6.97,11.8,.5,3.27,2.27,4.32,4.89E-2,6.03,22.6,.849,2.81,.375,
     * 1.05E-4,10.1,.5,7.90,.465,3.7E-2,2.31,1.4E-2,3.02E-2,3.19,.399,
     * 8.39/
      DATA PM/0.1396,0.1396,.4936,.4936,.9383/
      DATA SMX/1.878,2.016,2.053,2.370,2.814/
      DATA G00/163.,163.,7.33,7.33,3.15/
      DATA ST/2.018,2.156,2.547,2.864,3.752/
        J=JT-2
        IF(SQS-ST(J)) 5,5,8
  8     S=SQS*SQS
        EMAX=(S-SMX(J)**2+PM(J)**2)/2./SQS
        PMAX=SQRT(EMAX**2-PM(J)**2)
        XM=PM(J)/EMAX
C ------  IN THIS PROGRAMM INPUT XF IS DEFINED AS 2*P_PAR/SQRT(S)
C        PCM=PMAX*XF
        PCM=XF*SQS/2.
        ECM=SQRT(PCM**2+PM(J)**2+PT**2)
        X=ECM/EMAX
        IF(X-1.) 4,4,5
  4     XT=(X-XM)/(1.-XM)
        IF(A(3,J)-X) 1,1,2
  2     F1=A(1,J)*EXP(-A(2,J)*X)
        GO TO 3
  1     F1=0.
  3     F=(G00(J)-A(1,J))*(1.-X)**A(4,J)+F1
        AH=A(5,J)*EXP(-A(6,J)*X)+A(7,J)*EXP(A(8,J)*X)
        BH=A(9,J)*EXP(-A(10,J)*(X+A(11,J)))*(X+A(11,J))**A(12,J)
        IF(SQS-15.) 6,6,7
  6     RLE=RTAN(J,XT,SQS)
        TAN0=F*EXP(-(AH+BH*PT)*PT)*RLE
        RETURN
  7     TAN0=F*EXP(-(AH+BH*PT)*PT)
        RETURN
  5     TAN0=0.
        RETURN
        END
C**************************************************
        FUNCTION RTAN(J,XT,SQS)
        DIMENSION B(9,5),ST(5)
        DATA B/
     * 1.53,3.20,1.60,0.,.229,-.198,.068,.649,.147,
     * 3.00,3.20,1.14,0.,.264,-.679,1.08,.259,2.74,
     * 12.40,5.65,.522,0.521,.482,-.613,.616,.659,0.,
     * 2.99,2.09,1.70,1.090,.495,-.946,1.85,.505,2.88,
     * .3060,.120,.0552,2.720,.758,-.680,1.54,.594,2.87/
        DATA ST/2.018,2.156,2.547,2.864,3.752/
        AH=B(1,J)*EXP(-B(2,J)*XT)
        BH=B(3,J)*EXP(B(4,J)*XT)
        CH=B(5,J)+B(6,J)*XT+B(7,J)*XT**2
        DH=B(8,J)*EXP(B(9,J)*XT)
        Q=SQS-ST(J)
        R1=1.-EXP(-AH*Q**BH)
        R2=EXP(CH*Q-DH)
        RTAN=1./(1.-EXP(-R1*R2))
        RETURN
        END
        FUNCTION TANNG1(JT,XF,PT,SQS)
C       XF=2*PL/SQS
        DATA P0/.2/

        IF(JT-7) 2,1,2
  1     IF(PT-P0) 3,3,2
  3     A=TAN0(JT,XF,P0,SQS)*EXP(1.4*P0)
        TANNG1=A*EXP(-1.4*PT)*1.e3*29
        RETURN
  2     TANNG1=TAN0(JT,XF,PT,SQS)*1.e3*29
        RETURN
        END

	FUNCTION PLMAX(PO,TET,JP)
	DIMENSION PM(7)
	DIMENSION SMX(7)
	DATA PM/0.938272,0.9396,2*0.1396,2*0.493677,0.938272/
	DATA SMX/0.93827231,1.079,1.878,2.016,2.053,2.370,2.814/

	EO=SQRT(PO**2+PM(1)**2)
	S=2.*PM(1)*EO+PM(1)**2+PM(1)**2
	SQS=SQRT(S)
	emax=(s+pm(jp)**2-smx(jp)**2)/2./sqs
	pmax=sqrt(emax**2-pm(jp)**2)
	GACM=(EO+PM(1))/SQS
	GABCM=PO/SQS
        c2=-gabcm**2+(gacm/cos(tet))**2
        c1=gabcm*emax
        disc=gacm**2*(gabcm**2*pm(jp)**2+(emax**2-gacm**2*pm(jp)**2)
     &   /(cos(tet))**2)
        if(disc.ge.0.) then
           plmax=(c1+sqrt(disc))/c2
        else
           plmax=0.
        endif
	RETURN
      	END

