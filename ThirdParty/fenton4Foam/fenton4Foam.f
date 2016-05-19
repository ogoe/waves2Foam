      module common_variables

      CHARACTER*10 DEPTH,CASE,CURRNT
      integer :: n, num
      real*8 :: PI,HOVERD,HEIGHT,VALUE
      real*8, allocatable :: Z(:),COSA(:),SINA(:),COEFF(:),SOL(:,:),Y(:)

      end module common_variables

C-----------------------------------------------------------------------
      PROGRAM STEADY
C-----------------------------------------------------------------------
C
C     Calculation of steady waves
C
C-----------------------------------------------------------------------
C
      use common_variables
      IMPLICIT none
      integer :: j, iter, i, number, nstep, ns, info
      real*8 :: h, criter, crit, dho, dhe, sum

      integer, allocatable :: IPVT(:)
      real*8, allocatable :: RHS1(:),RHS2(:),A(:,:),B(:)
C
      OPEN (UNIT=4,FILE='fenton.inp',STATUS='unknown')
      OPEN (UNIT=7,FILE='fenton.out',STATUS='unknown')
      OPEN (UNIT=11,FILE='constant/fenton4Foam',STATUS='unknown')
C
C.....Input data
C
C     "depth" IS EITHER 'deep' OR 'finite'
C     "hoverd" IS WAVE HEIGHT/DEPTH
C
      READ (4,*) DEPTH,HOVERD
C
C     "case" IS EITHER 'period' OR 'wavelength'
C     "height" IS HEIGHT/LENGTH IF "case" IS 'wavelength'
C     "height" IS HEIGHT/(g*T**2) IF "case" IS 'period'
C
      READ (4,*) CASE,HEIGHT
C
C     "currnt" IS EITHER 'Euler' OR 'Stokes'
C     "value" is the magnitude of the mean Eulerian or Stokes velocities
C             non-dimensionalized with respect to wave height
C
      READ (4,*) CURRNT,VALUE
C
C     "n" is the number of terms in the Fourier series and the number
C         of intervals in half a wavelength
C     "nstep" is the number of steps in wave height
C
      READ (4,*) N,NSTEP
      NUM=2*N+10
C allocate the commonly accesed variables
      allocate (Z(NUM),COSA(0:NUM),SINA(0:NUM),COEFF(NUM),
     $     SOL(NUM,2),Y(NUM))
C allocate the locally accesed variables
      allocate (IPVT(NUM),RHS1(NUM),RHS2(NUM),A(NUM,NUM),B(NUM))
C
C     "number" is the number of iterations for each height step
C
      NUMBER=60
C
C     "crit" is the criterion for convergence. If the sum of magnitudes
C            of corrections is smaller than crit, the iteration stops
C      
      CRIT=1.D-6
      WRITE(6,20) DEPTH,HOVERD
      WRITE(6,21) HEIGHT,CASE
      WRITE(6,22) CURRNT,VALUE
      WRITE(6,30) N,NSTEP
      WRITE(7,20) DEPTH,HOVERD
      WRITE(7,21) HEIGHT,CASE
      WRITE(7,22) CURRNT,VALUE
      WRITE(7,30) N,NSTEP
      PI=4.D0*DATAN(1.D0)
      DHE=HEIGHT/NSTEP
      DHO=HOVERD/NSTEP
C
C.....Commence stepping through steps in wave height
C
      DO 1 NS=1,NSTEP
         WRITE(6,23) NS,NSTEP
         WRITE(7,23) NS,NSTEP
         HEIGHT=NS*DHE
         HOVERD=NS*DHO
C
C.....Calculate initial linear solution
C
         IF(NS.LE.1) THEN
            CALL INIT
         ELSE
C
C.....or extrapolate for next wave height, is neccessary
C
            DO 3 I=1,NUM
    3          Z(I)=2.D0*SOL(I,2)-SOL(I,1)
         ENDIF
C
C.....Commence iterative solution
C
         DO 4 ITER=1,NUMBER
c            WRITE(6,24) ITER
            WRITE(7,24) ITER
C
C.....Calculate right sides of equations and differentiate numerically
C.....to obtain Jacobian matrix
C
            CALL EQNS(RHS1)
            DO 5 I=1,NUM
               H=0.01D0*Z(I)
               IF(DABS(Z(I)).LT.1.D-4) H=1.D-5
               Z(I)=Z(I)+H
               CALL EQNS(RHS2)
               Z(I)=Z(I)-H
               B(I)=-RHS1(I)
               DO 6 J=1,NUM
    6             A(J,I)=(RHS2(J)-RHS1(J))/H
    5       CONTINUE
C
C.....Solve matrix equation and correct variables, using LINPACK 
C
            CALL DGEFA(A,NUM,NUM,IPVT,INFO)
            IF (INFO.NE.0) THEN
               WRITE(6,27)
               WRITE(7,27)
               STOP
            ENDIF
            CALL DGESL(A,NUM,NUM,IPVT,B,0)
C
C.....The b(i) are now corrections to each variable
C
            SUM=0.D0
            DO 7 I=1,NUM
               SUM=SUM+DABS(B(I))
    7          Z(I)=Z(I)+B(I)
C            WRITE(6,25) (Z(I),I=1,NUM)
c            WRITE(7,25) (Z(I),I=1,NUM)
            CRITER=CRIT
            IF(NS.EQ.NSTEP) CRITER=0.01D0*CRIT
            IF(SUM.LT.CRITER) GOTO 8
    4    CONTINUE
         WRITE(6,26) NUMBER,NS
         WRITE(7,26) NUMBER,NS
         STOP
    8    IF(NS.EQ.1) THEN
            DO 9 I=1,NUM
    9          SOL(I,2)=Z(I)
         ELSE
            DO 10 I=1,NUM
               SOL(I,1)=SOL(I,2)
   10          SOL(I,2)=Z(I)
         ENDIF
         write(6,28)iter
    1 CONTINUE
C
C.....Output of results
C
      CALL OUTPUT
      print*,'return from output'
C
      print *, ' '
      print *, 'Eta and phi written to spectank.init, more complete '
      print *, 'flow quantities written to fenton.plot.'
      print *, ' '

   20 FORMAT(//,'   Depth: ',A6,', Height/Depth',F7.4)
   21 FORMAT(/,'   Wave height',F9.6,', dimensionless with respect to ',
     .         A10)
   22 FORMAT(/,'   Current criterion ',A6,', Magnitude ', F5.2)
   23 FORMAT(/,'   Height step ',I2,' of ',I3)
   24 FORMAT('   Iteration ',I3)
   25 FORMAT(//,'   Solution vector',10(/,6E13.6))
   26 FORMAT(/,' did not converge after ',I3,' iterations, on height ',
     $ 'step ',i3)
   27 FORMAT(/,'   MATRIX SINGULAR')
 28   FORMAT(I3,'  Iterations ')
 30   FORMAT(/,'   Number of modes', I5, ', Steps of increasing ', 
     $ 'wave height',I5)
      STOP 
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE INIT
C-----------------------------------------------------------------------
C
C     Calculate initial solution from linear wave theory
C
C-----------------------------------------------------------------------
C
      use common_variables
      implicit none

      integer :: i
      real*8 :: t, b, a
C
      IF (DEPTH.EQ.'finite') THEN
         IF (CASE.EQ.'period') THEN
            A=4.D0*PI*PI*HEIGHT/HOVERD
            B=A/DSQRT(DTANH(A))
            T=DTANH(B)
            Z(1)=B+(A-B*T)/(T+B*(1.D0-T*T))
         ELSE
            Z(1)=2.D0*PI*HEIGHT/HOVERD
         ENDIF
         Z(2)=Z(1)*HOVERD
         Z(4)=DSQRT(DTANH(Z(1)))
      ELSE
         Z(1)=-1.D0
         Z(4)=1.D0
         IF (CASE.EQ.'period') THEN
            Z(2)=4.D0*PI*PI*HEIGHT
         ELSE
            Z(2)=2.D0*PI*HEIGHT
         ENDIF
      ENDIF
      Z(3)=2.D0*PI/Z(4)
      IF (CURRNT.EQ.'Euler') THEN
         Z(5)=VALUE*DSQRT(Z(2))
         Z(6)=0.D0
      ELSE
         Z(5)=0.D0
         Z(6)=VALUE*DSQRT(Z(2))
      ENDIF
      Z(7)=Z(4)
      Z(8)=0.D0
      Z(9)=0.5D0*Z(7)**2
      Z(10)=0.5D0*Z(2)
      COSA(0)=1.D0
      SINA(0)=0.D0
      DO 1 I=1,N
         COSA(I)=DCOS(I*PI/N)
         COSA(I+N)=DCOS((I+N)*PI/N)
         SINA(I)=DSIN(I*PI/N)
         SINA(I+N)=DSIN((I+N)*PI/N)
         Z(N+I+10)=0.D0
         Z(I+10)=0.5D0*Z(2)*COSA(I)
    1 CONTINUE
      Z(N+11)=0.5D0*Z(2)/Z(7)
C
C      WRITE(6,2) (Z(I),I=1,NUM)
c      WRITE(7,2) (Z(I),I=1,NUM)
    2 FORMAT(//,'   INITIAL LINEAR SOLUTION',10(/,6E13.6))
      DO 3 I=1,9
    3     SOL(I,1)=Z(I)
      DO 4 I=10,NUM
    4    SOL(I,1)=0.D0
      DO 5 I=1,NUM
    5    SOL(I,2)=0.D0
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE EQNS(RHS)
C-----------------------------------------------------------------------
C
C     Evaluation of equations
C
C-----------------------------------------------------------------------
C
      use common_variables
      implicit none

      integer :: i, it, nm, m, j
      real*8 :: e, c, s, psi, u, v
c      IMPLICIT REAL*8 (A-H,K,L,O-Z)
      real*8 :: RHS(*)
C
      IF(DEPTH.EQ.'finite') THEN
         RHS(1)=Z(2)-Z(1)*HOVERD
      ELSE
         RHS(1)=Z(1)+1.D0
      ENDIF
      IF(CASE.EQ.'wavelength') THEN
         RHS(2)=Z(2)-2.D0*PI*HEIGHT
      ELSE
         RHS(2)=Z(2)-HEIGHT*Z(3)**2
      ENDIF
      RHS(3)=Z(4)*Z(3)-PI-PI
      RHS(4)=Z(5)+Z(7)-Z(4)
      RHS(5)=Z(6)+Z(7)-Z(4)
      IF(DEPTH.EQ.'finite') THEN
         RHS(5)=RHS(5)-Z(8)/Z(1)
         DO 2 I=1,N
    2       COEFF(I)=Z(N+I+10)/DCOSH(I*Z(1))        
      ENDIF
      IT=6
      IF(CURRNT.EQ.'Euler') IT=5
      RHS(6)=Z(IT)-VALUE*DSQRT(Z(2))
      RHS(7)=Z(10)+Z(N+10)
      DO 1 I=1,N-1
    1   RHS(7)=RHS(7)+Z(10+I)+Z(10+I)
      RHS(8)=Z(10)-Z(N+10)-Z(2)
      DO 3 M=0,N
         PSI=0.D0
         U=0.D0
         V=0.D0
         IF(DEPTH.EQ.'finite') THEN
            DO 4 J=1,N
               NM=MOD(M*J,N+N)
               E=DEXP(J*(Z(1)+Z(10+M)))
               S=0.5D0*(E-1.D0/E)
               C=0.5D0*(E+1.D0/E)
               PSI=PSI+COEFF(J)*S*COSA(NM)
               U=U+J*COEFF(J)*C*COSA(NM)
               V=V+J*COEFF(J)*S*SINA(NM)
    4       CONTINUE
         ELSE
            DO 5 J=1,N
               NM=MOD(M*J,N+N)
               E=DEXP(J*Z(10+M))
               PSI=PSI+Z(N+J+10)*E*COSA(NM)
               U=U+J*Z(N+J+10)*E*COSA(NM)
    5          V=V+J*Z(N+J+10)*E*SINA(NM)
         ENDIF
         RHS(M+9)=PSI-Z(8)-Z(7)*Z(M+10)
         RHS(N+M+10)=0.5D0*((-Z(7)+U)**2+V**2)+Z(M+10)-Z(9)
    3 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE OUTPUT
C-----------------------------------------------------------------------
C
C     Output of results
C
C-----------------------------------------------------------------------
C
      use common_variables
      implicit none

      integer :: i, j, m
      real*8 :: pulse, ke, pe, slope, fsbc_dyn, fsbc_kin, q, r, s, 
     $     u, v, ub2, sxx, f, grav, ufactor, dx, sum, press, elevn, 
     $     c, x
C Automatic workspace
      real*8 :: phi(2*n), gradphi(2*n,2), phi_S(2*n), eta(2*n), 
     $     gradeta(2*n), k_fourier, knyq, dk, factor

C Prepare header for OpenFoam-formatted output
      WRITE(11, '(a, a)') '/*--------------------------------*- C++',  
     A      '-*----------------------------------*\'
      WRITE(11, '(a, a)') '| =========                 |           ',
     A      '                                      |'
      WRITE(11, '(a, a)') '| \\      /  F ield         | OpenFOAM: ',
     A      'The Open Source CFD Toolbox           |'
      WRITE(11, '(a, a)') '|  \\    /   O peration     | Version:  ',
     A      'Cross-version compatible waves2Foam   |'
      WRITE(11, '(a, a)') '|   \\  /    A nd           | Web:      ',
     A      'http://www.OpenFOAM.org               |'
      WRITE(11, '(a,a)') '|    \\/     M anipulation  |            ',
     A      '                                     |'
      WRITE(11, '(a, a)') '\*--------------------------------------',
     A      '-------------------------------------*/'
      WRITE(11, '(a)') 'FoamFile'
      WRITE(11, '(a)') '{'
      WRITE(11, '(a)') '    version     2.0;'
      WRITE(11, '(a)') '    format      ascii;'
      WRITE(11, '(a)') '    class       dictionary;'
      WRITE(11, '(a)') '    object      fenton4FoamCoeffs;'
      WRITE(11, '(a)') '}'
      WRITE(11, '(a, a)') '// * * * * * * * * * * * * * * * * * * * ',
     A      '* * * * * * * * * * * * * * * * * * //'
      WRITE(11, '(a)') ''

C
C.....Calculate Fourier coefficients of surface elevation
C
      DO 10 J=1,N
         SUM=0.5D0*(Z(10)+Z(N+10)*(-1.D0)**J)
         DO 11 M=1,N-1
   11       SUM=SUM+Z(10+M)*COSA(MOD(M*J,N+N))
   10    Y(J)=2.D0*SUM/N
      WRITE(6,1)
      WRITE(7,1)
    1 FORMAT(//,'   Solution, non-dimensionalized by wavenumber',/)
      IF (DEPTH.EQ.'finite') WRITE(6,2) Z(1)
      IF (DEPTH.EQ.'finite') WRITE(7,2) Z(1)
    2 FORMAT('   Water depth                       ',F10.6,/)
      WRITE(6,3)(Z(I),I=2,9)
      WRITE(7,3)(Z(I),I=2,9)
    3 FORMAT('   Wave height                       ',F10.6,//,
     .       '   Wave period                       ',F10.6,//,
     .       '   Wave speed                        ',F10.6,//,
     .       '   Mean Eulerian fluid velocity      ',F10.6,//,
     .       '   Mean mass transport velocity      ',F10.6,//,
     .       '   Mean fluid speed relative to wave ',F10.6,//,
     .       '   Volume flux due to waves          ',F10.6,//,
     .       '   Bernoulli constant                ',F10.6,//)
      WRITE(6,4) (Z(I),I=10,N+10)
      WRITE(7,4) (Z(I),I=10,N+10)

C     OpenFoam output format
      WRITE(11, '(a,I5,a)') 'nModes ', N, ';'
      WRITE(11, '(a,E13.6,a)') 'depth   depth  [0 0 0 0 0 0 0] ', 
     A      Z(1), ';'
      WRITE(11, '(a,E13.6,a)') 'period  period [0 0 0 0 0 0 0] ', 
     A      Z(3), ';'
      WRITE(11, '(a,E13.6,a)') 'uBar    uBar   [0 0 0 0 0 0 0] ', 
     A      Z(7), ';'
chbb
C Print some dimensional output.  We take lambda=2pi, (k=1).  
C
c---------------------------------
ccc commented out by HBR
      if(1.eq.0) then

      grav=9.81d00
      ufactor=sqrt(grav)
      c=ufactor*z(4)
      dx=pi/real(n)
C
C Compute the potential on eta from the velocity at eta.  
C
      knyq=pi/dx
      dk=2.d0*knyq/(2.d0*N)
      factor=1.d0/N
      do i=1,n+1
         x=dble(i-1)*dx
         call point(x,z(10+i-1),u,v,press,elevn,slope)
         eta(i)=elevn
         gradeta(i)=slope
         gradphi(i,1)=ufactor*u
         gradphi(i,2)=ufactor*v
C grad^S of phi^S
         phi_S(i)=gradphi(i,1)+gradeta(i)*gradphi(i,2)
      end do
C The 2nd half of the wave-length
      do i=n+2,2*n
         x=dble(i-1)*dx
         call point(x,z(10+i-1),u,v,press,elevn,slope)
         eta(i)=elevn
         gradeta(i)=slope
      end do
      do i=n+2,2*n
         x=dble(i-1)*dx
         call point(x,eta(i),u,v,press,elevn,slope)
         gradphi(i,1)=ufactor*u
         gradphi(i,2)=ufactor*v
C grad^S of phi^S
         phi_S(i)=gradphi(i,1)+gradeta(i)*gradphi(i,2)
      end do
      call drealft(phi_S,2*n,+1)
      phi(1)=0.d00
      phi(2)=0.d0
      do i=2,n
         k_fourier=(i-1)*dk
         phi(2*i-1)=-phi_S(2*i)/k_fourier
         phi(2*i)=phi_S(2*i-1)/k_fourier
      end do
      call drealft(phi,2*n,-1)
      phi=factor*phi

      open(unit=47,file='fenton.plot',status='unknown')
      open(unit=48,file='spectank.init',status='unknown')
      write(47,762)
 762  format('# Steady wave quantities:  x, eta, eta_x, phi, u, w, etat,
     $ phit')
      write(48,761)
 761  format('Steady wave according to stream function theory.')
      write(48,*)z(1),2.d0*pi,2*n,1
      do i=1,2*n
         x=real(i-1)*dx
         fsbc_dyn = -(grav*eta(i) + .5d00*((gradphi(i,1)-c)**2
     $        +gradphi(i,2)**2)) + grav*z(9)
         fsbc_kin = gradphi(i,2)-gradeta(i)*(gradphi(i,1)-c)
         write(47,763)x,eta(i),gradeta(i),phi(i),gradphi(i,1),
     $        gradphi(i,2),fsbc_kin,fsbc_dyn
         write(48,764)eta(i),phi(i)
      end do
 763  format(8e12.4)
 764  format(2e23.16)
      close(47)
      close(48)
chbb
      endif
ccccc-------------------------------
    4 FORMAT('   Surface elevations - crest to trough',//,110(E13.6,/))
      WRITE(6,5) (I,Z(I+N+10),I=1,N)
      WRITE(7,5) (I,Z(I+N+10),I=1,N)
    5 FORMAT(/,'   FOURIER COEFFICIENTS Bj',//,110((I3,2X,E13.6,3X),/))

C OpenFoam format of the Bj-coefficients
      WRITE(11, '(a)') ''
      WRITE(11, '(a,I5)') 'bCoeffs nonuniform List<scalar>', N
      WRITE(11, '(a)') '('
      DO 187, I=1,N
          WRITE(11, '(a,E13.6)') '    ',Z(I+N+10)
 187  CONTINUE
      WRITE(11, '(a)') ');'
      WRITE(11, '(a)') ''

chbr
      WRITE(6,6) (I,Y(I),I=1,N)
      WRITE(7,6) (I,Y(I),I=1,N)
    6 FORMAT(/,'   FOURIER COEFFICIENTS kAj',//,110((I3,2X,E13.6,3X),/))

C OpenFoam format of the A-coefficients
      WRITE(11, '(a)') ''
      WRITE(11, '(a,I5)') 'aCoeffs nonuniform List<scalar>', N
      WRITE(11, '(a)') '('
      DO 188, I=1,N
          WRITE(11, '(a,E13.6)') '    ',Y(I)
 188  CONTINUE
      WRITE(11, '(a)') ');'
      WRITE(11, '(a)') ''

C
      PULSE=Z(8)+Z(1)*Z(5)
      KE=0.5D0*(Z(4)*PULSE+Z(5)*(Z(8)-Z(7)*Z(1)))
      PE=0.5D0*(Z(10)**2+Z(N+10)**2)
      DO 7 I=1,N-1
    7    PE=PE+Z(10+I)**2
      PE=PE/(2.D0*N)
      UB2=2.D0*Z(9)-Z(4)**2
      SXX=4.D0*KE-3.D0*PE+UB2*Z(1)+2.D0*Z(5)*(Z(7)*Z(1)-Z(8))
      F=Z(4)*(3.D0*KE-2.D0*PE)+0.5D0*UB2*(PULSE+Z(4)*Z(1))+
     .        Z(4)*Z(5)*(Z(7)*Z(1)-Z(8))
      WRITE(6,8) PULSE,KE,PE,UB2,SXX,F
      WRITE(7,8) PULSE,KE,PE,UB2,SXX,F
    8 FORMAT(/,'   INTEGRAL QUANTITIES',//,
     .         '   Impulse(I)                  ',E13.6,//,
     .         '   Kinetic energy (T)          ',E13.6,//,
     .         '   Potential energy (V)        ',E13.6,//,
     .         '   Mean square of bed velocity ',E13.6,//,
     .         '   Radiation stress (Sxx)      ',E13.6,//,
     .         '   Wave power (F)              ',E13.6)
      IF(DEPTH.EQ.'finite') THEN
         Q=Z(7)*Z(1)-Z(8)
         R=Z(9)+Z(1)
         S=SXX-2.D0*Z(4)*PULSE+(Z(4)**2+0.5D0*Z(1))*Z(1)
         WRITE(6,9) Q,R,S
         WRITE(7,9) Q,R,S
    9    FORMAT(//,'   INVARIANTS FOR FINITE DEPTH',//,
     .             '   Volume flux (Q)        ',F9.6,//,
     .             '   Bernoulli constant (R) ',F9.6,//,
     .             '   Momentum flux (S)      ',F9.6)
      ENDIF

      CLOSE(11)
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE POINT(KX,KY,U,V,PRESS,ELEVN,slope)
C-----------------------------------------------------------------------
C
C     Calculation of free surface elevation ELEVN at KX and velocity
C     (U,V) and pressure PRESS at point (KX,KY).  The slope of the 
C     free-surface is also calculated.  -hrb
C
C-----------------------------------------------------------------------
C
      use common_variables
      implicit none

      real*8 :: kx, ky, u, v, press, elevn, slope
      integer :: j
      real*8 :: e, b, c, s
C
      ELEVN=0.5D0*Y(N)*DCOS(N*KX)
C The slope of the free-surface.  
      slope=-0.5D0*n*Y(N)*DSIN(N*KX)
      DO  J=1,N-1
         ELEVN=ELEVN+Y(J)*DCOS(J*KX)
         slope=slope-real(j)*y(j)*dsin(j*kx)
      end do
      U=Z(4)-Z(7)
      V=0.D0
      IF(DEPTH.EQ.'finite') THEN
         DO 4 J=1,N
            E=DEXP(J*(Z(1)+KY))
            S=0.5D0*(E-1.D0/E)
            C=0.5D0*(E+1.D0/E)
            B=Z(N+J+10)/DCOSH(J*Z(1))
            U=U+J*B*C*DCOS(J*KX)
    4       V=V+J*B*S*DSIN(J*KX)
      ELSE   
         DO 5 J=1,N
            E=DEXP(J*KY)
            U=U+J*Z(N+J+10)*E*DCOS(J*KX)
    5       V=V+J*Z(N+J+10)*E*DSIN(J*KX)
      ENDIF
      PRESS=Z(9)-KY-0.5D0*((U-Z(4))**2+V*V)
      RETURN
      END
c
      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(n),info
      double precision a(lda,n)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(n),job
      double precision a(lda,n),b(n)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      integer i,incx,m,mp1,n,nincx
      double precision da,dx(n)
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      integer i,incx,ix,n
      double precision dx(n),dmax
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      double precision function dasum(n,dx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      subroutine dgeco(a,lda,n,ipvt,rcond,z)
      integer lda,n,ipvt(1)
      double precision a(lda,1),z(1)
      double precision rcond
c
c     dgeco factors a double precision matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dgefa is slightly faster.
c     to solve  a*x = b , follow dgeco by dgesl.
c     to compute  inverse(a)*c , follow dgeco by dgesl.
c     to compute  determinant(a) , follow dgeco by dgedi.
c     to compute  inverse(a) , follow dgeco by dgedi.
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dgefa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,dsign
c
c     internal variables
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call dgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of w  where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(a(k,k))) go to 30
            s = dabs(a(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (a(k,k) .eq. 0.0d0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + dabs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(a(k,k))) go to 150
            s = dabs(a(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         t = -z(k)
         call daxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
      
