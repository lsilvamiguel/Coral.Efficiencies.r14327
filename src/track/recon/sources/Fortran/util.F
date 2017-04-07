

      SUBROUTINE STFIT5(NPLANE,COSN,SINN,FD,XPLANE,THETA,SGPLN,NHIT,ANS)

*********************************************************************
*                                                                   *
*    General Least Squares fit to get slopes and intercepts in 2    *
*    projections plus their variances                               *
*                                                                   *
*********************************************************************
*    Input :
*
*    ANS(1)    : X coordinate at which Fit results are given
*                 This is expected to be the position of first plane
*                 (FD(1)=0)
*
*    NPLANE    : Number of planes involved in the Fit
*
*    COSN,SINN : Angle of each plane involved in the Fit
*
*    FD        : double integral of main component of field from first plane
*                 to present plane
*
*    XPLANE    : X position of planes involved in the Fit
*
*    THETA     : HIT (U) coordiante on each plane
*                HIT=Y sin - Z cos
*
*    Output :
*
*    ANS(2) : Y intercept of the fitted track at ANS(1)
*
*    ANS(3) : Z intercept of the fitted track at ANS(1)
*
*    ANS(4) : Y Slope (dY/dX) of the fitted track at ANS(1)
*
*    ANS(5) : Z Slope (dZ/dX) of the fitted track at ANS(1)
*
*    ANS(6) : variance on Y intercept of the fitted track at ANS(1)
*
*    ANS(7) : variance on Z intercept of the fitted track at ANS(1)
*
*    ANS(8) : variance on Y Slope of the fitted track at ANS(1)
*
*    ANS(9) : variance on Z Slope of the fitted track at ANS(1)
*
*    ANS(10): CHISquared by degree of freedom
*
*    ANS(11): pb (=1./momentum) with sign = charge
*
*    ANS(12): variance on pb
*
*    NHIT   : Number of hit used for the fit
*
*
*    Assumed relation is :
*
*    THETA(I)= SINN(I)*(ANS(2)+(XPLANE(I)-ANS(1))*ANS(4)+FD(I)*ANS(11))
*             -COSN(I)*(ANS(3)+(XPLANE(I)-ANS(1))*ANS(5))
*
*========================================================================
 
      DOUBLE PRECISION A(5,5),B(5,5),TT(5),                     
     1 S1,S2,S3,S4,AA,BB,TH,THX,XX,SSS1,SSS2,D1,D2
 
      DIMENSION COSN(1),SINN(1),FD(*),THETA(1),XPLANE(1),SGPLN(1),
     1 ANS(12),wk(5)

       logical pout

       pout=.false.
  
      if (abs(ans(1)-xplane(1)).gt.0.005) then
         print *, "stfit5---WARNING---ANS(1) reset to
     1  pos. of plane 1--check FD(1)=0.  !!!!!"
          ans(1) = xplane(1)
      endif                    
    
      XG = ANS(1)
C
C *** SET DEFAULT VALUES
C
      NHIT=0

      do i=1,12
      ans(i)=88.
      enddo

          SSS1=0.0
          SSS2=0.0
*-----require at least 5 usable hits....................
      do 20 i=1,nplane
      if (theta(i).gt.80.0)  go to 20
      nhit = nhit + 1
       WEIGT=1.0/(SGPLN(I)**2)
       SSS1=SSS1+XPLANE(I)*WEIGT
       SSS2=SSS2+WEIGT
20    continue
      if (nhit.le.5)  return
      ANS(1)=SSS1/SSS2

C
C *** SET UP NORMAL EQUATIONS
C
      do i=1,5
      do j=1,5
      a(i,j) = 0.d0
      b(i,j) = 0.d0
      enddo
      tt(i) = 0.d0
      enddo
C
      ANS(10)=1.0E5
      DO 3 K=1,NPLANE
          AA=-COSN(K)
          BB=SINN(K)
          XX=XPLANE(K)-XG
          TH=THETA(K)
          FFD = FD(K)
C
C *** CHECK THERE IS A DATA POINT FOR THIS PLANE
C
          IF(TH.GT.80.) GO TO 3
          SIG2=1.0/SGPLN(K)**2
          if(pout.and.sgpln(k).eq.0)print*,k,sig2
          THX=TH*XX
          S1=AA*AA*SIG2
          S2=BB*BB*SIG2
          S3=AA*BB*SIG2
          AA=AA*SIG2
          BB=BB*SIG2
          S4=XX*XX

          A(1,1)=A(1,1)+S1
          A(1,2)=A(1,2)+S1*XX
          A(1,3)=A(1,3)+S3
          A(1,4)=A(1,4)+S3*XX
          A(1,5)=A(1,5)+S3*FFD
          A(2,2)=A(2,2)+S4*S1
          A(2,3)=A(1,4)
          A(2,4)=A(2,4)+S3*S4
          A(2,5)=A(2,5)+S3*XX*FFD
          A(3,3)=A(3,3)+S2
          A(3,4)=A(3,4)+S2*XX
          A(3,5)=A(3,5)+S2*FFD
          A(4,4)=A(4,4)+S2*S4
          A(4,5)=A(4,5)+S2*XX*FFD
          A(5,5)=A(5,5)+S2*(FFD**2)
          TT(1)=TT(1)+AA*TH
          TT(2)=TT(2)+AA*THX
          TT(3)=TT(3)+BB*TH
          TT(4)=TT(4)+THX*BB
          TT(5)=TT(5)+BB*FFD*TH
3     CONTINUE
C    
*     Symmetrize......................
          a(2,1) = a(1,2)
          a(3,1) = a(1,3)
          a(4,1) = a(1,4)
          a(5,1) = a(1,5)
          a(3,2) = a(2,3)
          a(4,2) = a(2,4)
          a(5,2) = a(2,5)
          a(4,3) = a(3,4)
          a(5,3) = a(3,5)
          a(5,4) = a(4,5)
C
C *** INVERSE OF 5*5 STRICTLY POSITIVE SYMMETRIC MATRIX
C

      if(pout)then
       print *,"FD:",(fd(i),i=1,nplane)
       do k=1,5
       print *,"A: ",(a(k,j),j=1,5)
       enddo
      endif

      CALL DINV(5,A,5,wk,IFAIL)
       if (ifail.ne.0) then
         print *,'======  ERROR in STFIT5===STOP'
         stop
       endif
    
      do i=1,5
      do j=1,5
      b(i,j) = a(i,j)
      enddo
      enddo
C
C *** FITTED VARIABLES
C
      ANS(6)=B(3,3)
      ANS(7)=B(1,1)
      ANS(8)=B(4,4)
      ANS(9)=B(2,2)
      ANS(12)=B(5,5)
      ANS(3)=B(1,1)*TT(1)+B(1,2)*TT(2)+B(1,3)*TT(3)+B(1,4)*TT(4)
     1  +B(1,5)*TT(5)
      ANS(5)=B(2,1)*TT(1)+B(2,2)*TT(2)+B(2,3)*TT(3)+B(2,4)*TT(4)
     1  +B(2,5)*TT(5)
      ANS(2)=B(3,1)*TT(1)+B(3,2)*TT(2)+B(3,3)*TT(3)+B(3,4)*TT(4)
     1  +B(3,5)*TT(5)
      ANS(4)=B(4,1)*TT(1)+B(4,2)*TT(2)+B(4,3)*TT(3)+B(4,4)*TT(4)
     1  +B(4,5)*TT(5)
      ANS(11)=B(5,1)*TT(1)+B(5,2)*TT(2)+B(5,3)*TT(3)+B(5,4)*TT(4)
     1  +B(5,5)*TT(5)
C
C *** FOR NO D.O.F SET CHISQ = 0.0
C
      ANS(10)=0.0
      IF(NHIT.LE.5) THEN
*--- kk -December 2001 - small correction just to be sure
       return
      ELSE
C
C *** WORK OUT CHISQUARED
C
          DO 5 J=1,NPLANE
              IF(THETA(J).GT.80.) GOTO 5
              AA=-COSN(J)
              BB=SINN(J)
              FFD = FD(J)
              XX=XPLANE(J)-XG
              THFIT=AA*(ANS(3)+ANS(5)*XX)+BB*(ANS(2)+ANS(4)*XX
     1         +ANS(11)*FFD)
              ANS(10)=ANS(10)+((THETA(J)-THFIT)/SGPLN(J))**2
*           print*,'stfit5 - j,ans(10) ',j,ans(10)
5         CONTINUE
          ANS(10)=ANS(10)/(NHIT-5)
      ENDIF

      RETURN
      END



      SUBROUTINE STFIT1(NPLANE,COSN,SINN,XPLANE,THETA,SGPLN,NHIT,ANS)

*********************************************************************
*                                                                   *
*    General Least Squares fit to get slopes and intercepts in 2    *
*    projections plus their variances                               *
*                                                                   *
*********************************************************************
*    Input :
*
*    ANS(1)    : X coordinate at which Fit results are given
*                If ANS(1)<-1000 Then the mean value of the x positions
*                of the planes used for the fit is used (weighted by plane
*                the plane RMS)
*    NPLANE    : Number of planes involved in the Fit
*
*    COSN,SINN : Angle of each plane involved in the Fit
*
*    XPLANE    : X position of planes involved in the Fit
*
*    THETA     : HIT (U) coordiante on each plane
*                HIT=Y sin - Z cos
*
*    Output :
*
*    ANS(2) : Y intercept of the fitted track at ANS(1)
*
*    ANS(3) : Z intercept of the fitted track at ANS(1)
*
*    ANS(4) : Y Slope (dY/dX) of the fitted track at ANS(1)
*
*    ANS(5) : Z Slope (dZ/dX) of the fitted track at ANS(1)
*
*    ANS(6) : variance on Y intercept of the fitted track at ANS(1)
*
*    ANS(7) : variance on Z intercept of the fitted track at ANS(1)
*
*    ANS(8) : variance on Y Slope of the fitted track at ANS(1)
*
*    ANS(9) : variance on Z Slope of the fitted track at ANS(1)
*
*    ANS(10): CHISquared by degree of freedom
*
*    NHIT   : Number of hit used for the fit
*
*
*    Assumed relation is :
*
*    THETA(I)= SINN(I)*(ANS(2)+(XPLANE(I)-ANS(1))*ANS(4))
*             -COSN(I)*(ANS(3)+(XPLANE(I)-ANS(1))*ANS(5))
*
*========================================================================
 
      DOUBLE PRECISION A11,A12,A13,A14,A22,A23,A24,A33,A34,A44,
     1 B11,B12,B13,B14,B22,B23,B24,B33,B34,B44,TT1,TT2,TT3,TT4,
     2 S1,S2,S3,S4,AA,BB,TH,THX,XX,SSS1,SSS2,D1,D2,R311,R312,R321,
     3 R322,R511,R512,R522
 
      DIMENSION COSN(1),SINN(1),THETA(1),XPLANE(1),SGPLN(1),ANS(10)

      dimension ERROR4(10)

  
      XSAVE=ANS(1)
C
C *** SET DEFAULT VALUES
C
      NHIT=0

      do i=1,9
      ans(i)=88.
      enddo

*      CALL VFILL(ANS,9,88.)
C
C *** TAKE X COORDINATE FROM INPUT ??
C
      IF(XSAVE.GT.-1000.0) THEN
 
          DO 20 I=1,NPLANE
              IF(THETA(I).GT.80.0) GO TO 20
              NHIT=NHIT+1
  20      CONTINUE
          IF(NHIT.LT.4) RETURN
	  ANS(1) = XSAVE
      ELSE
 
C
C *** GET X COORD FOR RESULTS
C
          SSS1=0.0
          SSS2=0.0
          DO 30 I=1,NPLANE
              IF(THETA(I).GT.80.0) GO TO 30
              NHIT=NHIT+1
              WEIGT=1.0/(SGPLN(I)**2)
              SSS1=SSS1+XPLANE(I)*WEIGT
              SSS2=SSS2+WEIGT
  30      CONTINUE
          IF(NHIT.LT.4) RETURN
          ANS(1)=SSS1/SSS2
      ENDIF
C
C *** PROTECT AGAINST STUPID VALUES
C
      IF(ABS(ANS(1)).GT.20.0) ANS(1)=0.0
      XG=ANS(1)
C
C *** SET UP NORMAL EQUATIONS
C
      A11=0.d0
      A12=0.d0
      A13=0.d0
      A14=0.d0
      A22=0.d0
      A23=0.d0
      A24=0.d0
      A33=0.d0
      A34=0.d0
      A44=0.d0
      TT1=0.d0
      TT2=0.d0
      TT3=0.d0
      TT4=0.d0
C
      ANS(10)=1.0E5
      DO 3 K=1,NPLANE
          AA=-COSN(K)
          BB=SINN(K)
          XX=XPLANE(K)-XG
          TH=THETA(K)
C
C *** CHECK THERE IS A DATA POINT FOR THIS PLANE
C
          IF(TH.GT.80.) GO TO 3
          SIG2=1.0/SGPLN(K)**2
          THX=TH*XX
          S1=AA*AA*SIG2
          S2=BB*BB*SIG2
          S3=AA*BB*SIG2
          AA=AA*SIG2
          BB=BB*SIG2
          S4=XX*XX
          A11=A11+S1
          A12=A12+S1*XX
          A13=A13+S3
          A14=A14+S3*XX
          A22=A22+S4*S1
          A24=A24+S3*S4
          A33=A33+S2
          A34=A34+S2*XX
          A44=A44+S2*S4
          TT1=TT1+AA*TH
          TT2=TT2+AA*THX
          TT3=TT3+BB*TH
          TT4=TT4+THX*BB
3     CONTINUE
      A23=A14
C
C
C *** INVERSE OF 4*4 STRICTLY POSITVE SYMMETRIC MATRIX
C
      D1   = 1.d0 / (A12*A12 - A11*A22)
 
      R311 = (A22*A13 - A12*A23) * D1
      R312 = (A22*A14 - A12*A24) * D1
      R321 = (A11*A23 - A12*A13) * D1
      R322 = (A11*A24 - A12*A14) * D1
 
      R511 = A13*R311 + A23*R321 + A33
      R512 = A13*R312 + A23*R322 + A34
      R522 = A14*R312 + A24*R322 + A44
 
      D2   = 1.d0 / ( R511*R522 - R512*R512 )
 
      B33 =   R522 * D2
      B34 = - R512 * D2
      B44 =   R511 * D2
 
      B13 = R311 * B33 + R312 * B34
      B14 = R311 * B34 + R312 * B44
      B23 = R321 * B33 + R322 * B34
      B24 = R321 * B34 + R322 * B44
 
      B11 = R311*B13 + R312*B14 - A22 * D1
      B12 = R311*B23 + R312*B24 + A12 * D1
      B22 = R321*B23 + R322*B24 - A11 * D1
 
C
C *** FILL PACKED ERROR MATRIX ON Y',Z',Y,Z (IN THAT ORDER)
C
      ERROR4(1)=B44
      ERROR4(2)=B24
      ERROR4(3)=B34
      ERROR4(4)=B14
      ERROR4(5)=B22
      ERROR4(6)=B23
      ERROR4(7)=B12
      ERROR4(8)=B33
      ERROR4(9)=B13
      ERROR4(10)=B11
C
C *** FITTED VARIABLES
C
      ANS(6)=B33
      ANS(7)=B11
      ANS(8)=B44
      ANS(9)=B22
      ANS(3)=B11*TT1+B12*TT2+B13*TT3+B14*TT4
      ANS(5)=B12*TT1+B22*TT2+B23*TT3+B24*TT4
      ANS(2)=B13*TT1+B23*TT2+B33*TT3+B34*TT4
      ANS(4)=B14*TT1+B24*TT2+B34*TT3+B44*TT4
C
C *** FOR NO D.O.F SET CHISQ = 0.0
C
      ANS(10)=0.0
      IF(NHIT.EQ.4) THEN
*--- kk -December 2001 - small correction just to be sure
      ELSEif(nhit.gt.4)then
C
C *** WORK OUT CHISQUARED
C
          DO 5 J=1,NPLANE
              IF(THETA(J).GT.80.) GOTO 5
              AA=-COSN(J)
              BB=SINN(J)
              XX=XPLANE(J)-XG
              THFIT=AA*(ANS(3)+ANS(5)*XX)+BB*(ANS(2)+ANS(4)*XX)
              ANS(10)=ANS(10)+((THETA(J)-THFIT)/SGPLN(J))**2
5         CONTINUE
          ANS(10)=ANS(10)/(NHIT-4)
      ENDIF

      RETURN
      END


      
      subroutine matmu2(m,n,k,a,b,c,nla,nlb,nlc)
************************************************
*                                              *
*     new one                                  *
*                                              *
************************************************

      dimension a(nla,*),b(nlb,*),c(nlc,*)
      if(m.gt.nlc) then
         write(6,*) 'wrong input param for matmu2'
         write(6,*) 'm,n,k..',m,n,k, '  nla,nlb,nlc ',nla,nlb,nlc
         stop 34
      endif

      do 10 j=1,m
      do 20 l=1,k
      c(j,l) = 0.
 20   continue
 10   continue
 
       do 11 j=1,m
       do 12 l=1,k
       do 13 ls=1,n
       c(j,l) = c(j,l) + a(j,ls) * b(ls,l)
 13     continue
 12     continue
 11     continue
      end

      subroutine trans2(m,n,a,at,nla,nlat)
**********************************************
*                                            *
*      new one                               *
*                                            *
**********************************************

      dimension a(nla,*),at(nlat,*)
      do 10 jl = 1,n
      do 20 jc = 1,m
      at(jl,jc) = a(jc,jl)
20    continue
10    continue

      return
      end



      SUBROUTINE SPNPOL(NN,XX,YY,AX,BX,XIN,YT) 

*****************************************************
*                                                   *
*                                                   *
*                                                   *
*****************************************************
              
      DIMENSION XX(NN),YY(NN),AX(NN),BX(NN)
                 
     
*----- GT replaced by GE to fix problem of last physics plane...RW..Nov.23,1999
      IF (XIN.LE.XX(1)) THEN               
       YT = YY(1) + ((XIN-XX(1))/(XX(2)-XX(1))) * (YY(2) - YY(1))    
      ELSEIF (XIN.GE.XX(NN)) THEN              
       YT = YY(NN) + ((XIN-XX(NN))/(XX(NN-1)-XX(NN)))*(YY(NN-1)-YY(NN)) 
      ELSE
                                                             
       DO 21 K=2,NN                                                
        IF (XIN.LT.XX(K)) GO TO 22                                 
21      CONTINUE 
                                                 
22      IR = K - 1                                               
        DT = XIN - XX(IR)                                       
        DD = XX(IR+1) - XX(IR)                                 
        YT = YY(IR) + DT*AX(IR) + 0.5 * (DT**2) * BX(IR)      
        YT = YT + (DT**3)*(BX(IR+1)- BX(IR))/(6*DD)         
       ENDIF

      return                
      END
                     

      SUBROUTINE SPNFIT (NN,XX,YY,AX,BX,ZIN1,ZIN2,IRET) 

*****************************************************
*                                                   *
*                  Spline                           *
*                                                   *
*****************************************************
         
      DIMENSION XX(NN),YY(NN),AX(NN),BX(NN),ZIN1(NN),ZIN2(NN)

      iret = 0

      IF (NN.LT.3)  THEN                                    
**    IF (NN.LT.4)  THEN                                    
*         PRINT *,' NN TOO LOW - NO SPLINE!' 
        iret = -90               
         RETURN                                          
      ENDIF 
                                           
      AX(1) = 1.                                      
      BX(1) = 0.                                     
      NM = NN - 1  
                                 
      DO 10 K= 2,NM                                
      DEN = XX(K+1) - XX(K-1)                     
      TA = (XX(K) - XX(K-1))/(2*DEN)                              
      TB = (XX(K+1) - XX(K))/(2*DEN)                             
      AX(K) = - TB/(1+TA*AX(K-1))                               
      RR = 3*((YY(K+1)-YY(K))/(XX(K+1)-XX(K)) -                
     1 (YY(K) - YY(K-1))/(XX(K)-XX(K-1)))/DEN                 
10    BX(K) = (RR - TA*BX(K-1))/(1.+TA*AX(K-1)) 
             
      BX(NN) = BX(NM)/(1.-AX(NM))                           
      BX(NM) = BX(NN)                                      
      NMM = NN - 2    
                                   
      DO 20  K=2,NMM                                    
      KA = NMM + 2 - K                                              
20    BX(KA) = AX(KA)*BX(KA+1)+BX(KA) 
                             
      BX(1) = BX(2)  
                                             
      DO 30 K=1,NM                                               
      DD = XX(K+1) - XX(K)                                      
30    AX(K) = (YY(K+1)-YY(K))/DD -DD*BX(K)/3 -DD*BX(K+1)/6 
    
      ax(nn) = ax(nm) +(xx(nn)-xx(nm)) * bx(nm)
      zin1(1) = 0.

      do 40 k=2,nn
      dd = xx(k) - xx(k-1)
      zin1(k) = zin1(k-1) + 0.5 * dd*(yy(k)+yy(k-1)) -(1./24.0)*
     1  (dd**3) * (bx(k) + bx(k-1))
40    continue

      zin2(1) = 0.

      do 50 k=2,nn
      dd = xx(k) - xx(k-1)
      zin2(k) = zin2(k-1) + zin1(k-1) * dd + (1./6.0) * (dd**2) *
     1 (2.*yy(k-1)+yy(k)) -(1./36.) * (dd**4)*(0.8*bx(k-1)+0.7*bx(k))
50    continue

      return
      END 

      subroutine track_fd(par_in,z_start,z_end,par_out)
*
*---- tracking in solenoid field.......
*
      dimension par_in(*),par_out(*)
       double precision xx(3),h(3)   !  IMPORTANT >>>>>>>>>>>
          step = 0.018
          xx(1) = par_in(1)
          xx(2) = par_in(2)
          xx(3) = z_start   
          xp = par_in(3)
          yp = par_in(4)
          jfin = 0
19        continue
          call mag_field(xx,h)
          facr = sqrt(1. + xp**2 + yp**2)
          ff1 = h(3) * yp + h(1) * xp * yp - h(2)*(1.+xp**2)
          ff2 = -h(3) * xp - h(2) * xp * yp + h(1)*(1.+yp**2)
          x2 =  0.3 * facr * ff1 /par_in(5)
          y2 =  0.3 * facr * ff2 /par_in(5)
          if (xx(3)+step.gt.z_end) then
              jfin = 1
              dz = z_end - xx(3)
              xx(3) = z_end     
          else
              dz =  step
              xx(3) = xx(3) + step
          endif
          xx(1) = xx(1) + xp * dz + 0.5 * x2 * dz**2 
          xx(2) = xx(2) + yp * dz + 0.5 * y2 * dz**2
          xp = xp + dz * x2 
          yp = yp + dz * y2
          if (jfin.eq.0) go to 19
*
          par_out(1) = xx(1)
          par_out(2) = xx(2)
          par_out(3) = xp
          par_out(4) = yp
      end



      subroutine get_par2(jpl,htt,xcen,th_par,ifail)
#include "PARAM.h"
#include "PLANES.h"
      dimension th_par(4),htt(5),jpl(5)  
      dimension cnn(5),snn(5),xpp(5),sig(5),ans(10)
*                                    
*-----------------------------------------------------------------------------------
*........ Find space point at X=XCEN from 5 hits htt(1...5) in planes jpl(1...5)
*-----------------------------------------------------------------------------------
*
      do kp = 1,5
      cnn(kp) = cosn2(jpl(kp))      
      snn(kp) = sinn2(jpl(kp))
      xpp(kp) = xplan2(jpl(kp))
      sig(kp) = sgpln2(kp)
      enddo
*
*------------------------------------
*
      ans(1) = xcen
      call stfit1(5,cnn,snn,xpp,htt,sig,nhit,ans)
*
*------------------------------------
*
      do k = 1,4
      th_par(k) = ans(k+1)
      enddo
      print '("space-pt:",4f10.4)',(th_par(k),k=1,4)
      print '("chi2:",e13.4)',ans(10)
*
*----- reject phony space points.......................
*
*
      xpmin = -0.5 + 0.4 * th_par(1)
      xpmax = xpmin + 1.0
      ypmin = -0.25 + 0.3 * th_par(2)
      ypmax = ypmin + 0.5
      if (th_par(3).lt.xpmin.or.th_par(3).gt.xpmax) then
         ifail = -98
         print '("Get-par2,pr3,min,max",3f10.4)',th_par(3),xpmin,xpmax
      endif

      if (th_par(4).lt.ypmin.or.th_par(4).gt.ypmax) then
          ifail = -99
         print '("Get-par2,pr4,min,max",3f10.4)',th_par(4),ypmin,ypmax
      endif
*
      end



      real function erscal(itel,iplane)
#include "PARAM.h"
#include "PLANES.h"
*


*---TEMP ..................
*      do k=1,5
*      vrscal(k) = 2.0
*      enddo
*----end TEMP.................
      iv = 0
      if (itel.eq.1)  then
       if (ityp1(iplane).eq.111) then
*--------------MMegas....................
            iv = 1
       elseif (ityp1(iplane).eq.222) then
*--------------DC.......................
            iv = 2
       elseif (ityp1(iplane).eq.333) then
*--------------SciFis...................
            iv = 3
       endif
      elseif (itel.eq.2) then
        if (ityp2(iplane).eq.222)  then
*--------------DC.......................
           iv = 2
        elseif (ityp2(iplane).eq.333) then
*--------------SciFis...................
           iv = 3
        elseif (ityp2(iplane).eq.444) then
*--------------GEMS.....................
           iv = 4
        elseif (ityp2(iplane).eq.555)  then
*--------------STRAWS...................
           iv = 5
        endif
      endif


      if (iv.eq.0) then
       print*,itel,iplane
       print *,' erscal err.......set to 1.0'
       erscal = 1
      else
       erscal = vrscal(iv)
      endif
      end






