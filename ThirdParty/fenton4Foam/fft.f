      SUBROUTINE DSINFT(Y,N)
C Double-precision version of Numerical Recipes fast sine transform.
      implicit none
      integer n, j
      real*8 y(n)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,y1,y2,sum

      THETA=3.141592653589793D0/DBLE(N)
      WR=1.0D0
      WI=0.0D0
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      Y(1)=0.0d0
      DO 11 J=1,N/2
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
        Y1=WI*(Y(J+1)+Y(N-J+1))
        Y2=0.5d0*(Y(J+1)-Y(N-J+1))
        Y(J+1)=Y1+Y2
        Y(N-J+1)=Y1-Y2
11    CONTINUE
      CALL DREALFT(Y,N,+1)
      SUM=0.0d0
      Y(1)=0.5d0*Y(1)
      Y(2)=0.0d0
      DO 12 J=1,N-1,2
        SUM=SUM+Y(J)
        Y(J)=Y(J+1)
        Y(J+1)=SUM
12    CONTINUE
      RETURN
      END

      SUBROUTINE dcosft(Y,N)
C
C This is a double-precision version of a fast cosine transform routine 
C from Numerical Recipes.  Computes the transform of Y(1:N+1) where N is a 
C power of two.  This is also the inverse transform but in that case 
C multiply Y by 2/N.  
C 
      IMPLICIT NONE
      INTEGER n, j
      REAL*8 Y(N+1)
      REAL*8 WR, WI, WPR, WPI, WTEMP, THETA, SUM, y1, y2
      THETA=3.141592653589793D0/DBLE(n)
      WR=1.0D0
      WI=0.0D0
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      SUM=.5d0*(Y(1)-y(n+1))
      y(1)=.5d0*(y(1)+y(n+1))
      DO 11 J=1,n/2 -1
         WTEMP=WR
         WR=WR*WPR-WI*WPI+WR
         WI=WI*WPR+WTEMP*WPI+WI
         Y1=0.5d0*(Y(J+1)+Y(n-J+1))
         Y2=(Y(J+1)-Y(n-J+1))
         Y(J+1)=Y1-WI*Y2
         Y(n-J+1)=Y1+WI*Y2
         SUM=SUM+WR*Y2
 11   CONTINUE
      CALL drealft(Y,n,+1)
      y(n+1)=y(2)
      Y(2)=SUM
      DO 12 J=4,N,2
         SUM=SUM+Y(J)
         Y(J)=SUM
 12   CONTINUE
      RETURN
      END SUBROUTINE DCOSFT


      SUBROUTINE DREALFT(DATA,N,ISIGN)
C
C Double-precision version of Numerical recipes REALFT routine.  
C
      implicit none
      integer n, isign
      real*8 DATA(n)
Clocal variables
      integer i,i1,i2,i3,i4,n2p3
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,c1,c2,h1i,h1r,h2i,h2r,wis,wrs

      THETA=3.141592653589793D0/DBLE(N/2)
      C1=0.5
      IF (ISIGN.EQ.1) THEN
         C2=-0.5d00
         CALL dFOUR1(DATA,N/2,+1)
      ELSE
         C2=0.5d00
         THETA=-THETA
      ENDIF
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.0D0+WPR
      WI=WPI
      N2P3=N+3
      DO 11 I=2,N/4
         I1=2*I-1
         I2=I1+1
         I3=N2P3-I2
         I4=I3+1
         WRS= (WR)
         WIS= (WI)
         H1R=C1*(DATA(I1)+DATA(I3))
         H1I=C1*(DATA(I2)-DATA(I4))
         H2R=-C2*(DATA(I2)+DATA(I4))
         H2I=C2*(DATA(I1)-DATA(I3))
         DATA(I1)=H1R+WRS*H2R-WIS*H2I
         DATA(I2)=H1I+WRS*H2I+WIS*H2R
         DATA(I3)=H1R-WRS*H2R+WIS*H2I
         DATA(I4)=-H1I+WRS*H2I+WIS*H2R
         WTEMP=WR
         WR=WR*WPR-WI*WPI+WR
         WI=WI*WPR+WTEMP*WPI+WI
 11   CONTINUE
      IF (ISIGN.EQ.1) THEN
         H1R=DATA(1)
         DATA(1)=H1R+DATA(2)
         DATA(2)=H1R-DATA(2)
      ELSE
         H1R=DATA(1)
         DATA(1)=C1*(H1R+DATA(2))
         DATA(2)=C1*(H1R-DATA(2))
         CALL dFOUR1(DATA,N/2,-1)
      ENDIF
      RETURN
      END SUBROUTINE DREALFT

      SUBROUTINE dFOUR1(DATA,NN,ISIGN)
C
C Double-precision version of Numerical recipes FOUR1 routine.  
C
      implicit none
      integer nn,isign
      real*8 data(2*nn)
C Local variables
      integer i,istep,j,m,mmax,n
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,tempi,tempr

      N=2*NN
      J=1
      DO 11 I=1,N,2
         IF(J.GT.I)THEN
            TEMPR=DATA(J)
            TEMPI=DATA(J+1)
            DATA(J)=DATA(I)
            DATA(J+1)=DATA(I+1)
            DATA(I)=TEMPR
            DATA(I+1)=TEMPI
         ENDIF
         M=N/2
 1       IF ((M.GE.2).AND.(J.GT.M)) THEN
            J=J-M
            M=M/2
            GO TO 1
         ENDIF
         J=J+M
 11   CONTINUE
      MMAX=2
 2    IF (N.GT.MMAX) THEN
         ISTEP=2*MMAX
         THETA=6.28318530717959D0/(ISIGN*MMAX)
         WPR=-2.D0*DSIN(0.5D0*THETA)**2
         WPI=DSIN(THETA)
         WR=1.D0
         WI=0.D0
         DO 13 M=1,MMAX,2
            DO 12 I=M,N,ISTEP
               J=I+MMAX
               TEMPR= (WR)*DATA(J)- (WI)*DATA(J+1)
               TEMPI= (WR)*DATA(J+1)+ (WI)*DATA(J)
               DATA(J)=DATA(I)-TEMPR
               DATA(J+1)=DATA(I+1)-TEMPI
               DATA(I)=DATA(I)+TEMPR
               DATA(I+1)=DATA(I+1)+TEMPI
 12         CONTINUE
            WTEMP=WR
            WR=WR*WPR-WI*WPI+WR
            WI=WI*WPR+WTEMP*WPI+WI
 13      CONTINUE
         MMAX=ISTEP
         GO TO 2
      ENDIF
      RETURN
      END SUBROUTINE dFOUR1
