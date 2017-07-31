      SUBROUTINE GTCF(A,B,N,K,TCF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(N), B(N), TCF(K)
Cf2py INTENT(OUT) TCF
C
C     LOOP OVER DELAY TIMES
C
      DO I=1, K
C        COMPUTE AVAILABLE AMOUNT OF DATA POINTS
         M = N - I  + 1
C
C        LOOP OVER DATA POINTS
C
         TSUM = 0.0D+00
         DO J=1, M
            TSUM = TSUM + A(J) * B(J+I-1)
         END DO
C        AVERAGE
         TCF(I) = TSUM / DFLOAT(M)
      END DO
C
      RETURN
      END
C ------------------------------------------------
      SUBROUTINE EXPRES(F,STEP,NDEL,RR,RI,NMAX)
C
C     Compute exponential correlation function
C
C     RES(t) = < exp [ -i \int_t0^t F(t') dt' ] >_t0
C
C     F     - array of length N
C     STEP  - integration step for trapezoidal rule; STEP = ARG(F(1)) - ARG(F(0))
C     NDEL  - number of delay times (exluding first point)
C     RR    - real part of correlation function
C     RI    - imaginary part of correlation function
C     NMAX  - number of samples
C
Cf2py INTENT(OUT) RR, RI
      INTEGER NMAX, NDEL
      DOUBLE PRECISION F(NMAX), RR(NDEL+1), RI(NDEL+1)
      DOUBLE PRECISION STEP

      INTEGER I, J, K, IK
      DOUBLE PRECISION VH, VJR, VJI, VNJ
      DOUBLE PRECISION S(NMAX-1)
      LOGICAL LWRITE
      VH = STEP/2.0D+00
      RR(1) = 1.0D+00
      RI(1) = 0.0D+00
      LWRITE= .FALSE.

      IF(LWRITE) WRITE(*,*) " COMPUTING PARTIAL SUMS..."
      DO I = 1, NMAX-1
         S(I) = (F(I+1) + F(I))*VH
      ENDDO

      IF(LWRITE) WRITE(*,*) " COMPUTING RESPONSE FUNCTION..."
      DO J = 1, NDEL
         VJR = 0.0D+00
         VJI = 0.0D+00
         DO I = J, NMAX-1
            IK  = I - J + 1

            VH = 0.0D+00
            DO K = IK, I
               VH = VH + S(IK)
            ENDDO
            VJR = VJR + DCOS(VH)
            VJI = VJI - DSIN(VH)
         ENDDO
         VNJ= DFLOAT(NMAX - J)
         VJR= VJR / VNJ 
         VJI= VJI / VNJ
         RR(J+1) = VJR
         RI(J+1) = VJI
      ENDDO

      RETURN
      END
C ------------------------------------------------
C
C  SUBROUTINES BELOW ARE OBSOLETE AND SHOULD NOT BE USED
C
C
      SUBROUTINE GENRES(F,N,STEP,RES)
C
C     Compute exponential correlation response function
C
C     RES(t) = < exp [ -i \int_t0^t F(t') dt' ] >_t0
C
C     F     - array of length N
C     STEP  - integration step for Simpson's rule; STEP = ARG(F(1)) - ARG(F(0))
C     RES   - output function
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION F(N), RES((N-1)/2)
Cf2py INTENT(OUT) RES
C
      MMAX = (N-1)/2
      DO I = 1, MMAX
         NI = MMAX - I + 1
         VSUM = 0.0D+00
         DO J = 1, NI
            VSUM = VSUM + RUNINT(F(I), 2*I + 1,STEP)
         END DO
         RES(NI) = VSUM / DFLOAT(NI)
      END DO
C
      RETURN
      END
C-------------------------------------------------
      SUBROUTINE GENREST(F,N,STEP,RES)
C
C     Compute exponential correlation response function
C
C     RES(t) = < exp [ -i \int_t0^t F(t') dt' ] >_t0
C
C     F     - array of length N
C     STEP  - integration step for Simpson's rule; STEP = ARG(F(1)) - ARG(F(0))
C     RES   - output function (length N-1)
C
C     Use Trapezoidal rule
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION F(N), RES(N-1,2)
Cf2py INTENT(OUT) RES
      PREF = STEP/2.0D+00
      MMAX = N-1 
C
C     SCRATCH FILE
C
      OPEN(UNIT=17,STATUS='SCRATCH',ACTION='READWRITE',
     &     FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
C
C     COMPUTE RUNNING INTEGRALS AND STORE THEM IN SCRATCH FILE
C
      DO IP = 1, MMAX
         X  = 0.0D+00
         DO IM = IP, MMAX
C
            X  = X + PREF * (F(IM) + F(IM+1))
C
            WRITE(17) DCOS(X), -DSIN(X)
C
         END DO
      END DO
C
      REWIND(UNIT=17)
C
C     ACCUMULATE RES INTEGRALS
C 
      DO IP = 1, MMAX
         DO IM = IP, MMAX
            NR = IM - IP + 1
            READ(17) VRE, VIM
            RES(NR,1) = RES(NR,1) + VRE
            RES(NR,2) = RES(NR,2) + VIM
         END DO
      END DO
C
C     AVERAGE OUT RES INTEGRALS
C
      DO M=1,MMAX
         VM = DFLOAT(N-M)
         RES(M,1) = RES(M,1) / VM
         RES(M,2) = RES(M,2) / VM
      END DO
C
      RETURN
      END
C ---------------------------------------------------------------
      SUBROUTINE GENRESF(F,N,STEP,RES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION F(N), RES((N-1)/2, 2)
Cf2py INTENT(OUT) RES
      PI2  = 1.00D+00/6.283185307D+00
      PREF = STEP/3.0D+00
      MMAX =(N-1)/2
C
C     SCRATCH FILE
C
      OPEN(UNIT=17,STATUS='SCRATCH',ACTION='READWRITE',
     &     FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
C
C     COMPUTE RUNNING INTEGRALS AND STORE THEM IN SCRATCH FILE
C
      DO IP = 1, MMAX
         X   = 0.0D+00
         DO II = IP, MMAX
            IIA = 2*(II-1) + 1
            IIB = 2* II    + 1
            IIC = IIA
C
            X   = X + PREF * (F(IIA) + F(IIB) + 4.0D+00*F(IIC))
C
            WRITE(17) DCOS(X), -DSIN(X)
         END DO
      END DO
C
      REWIND(UNIT=17)
C
C     ACCUMULATE RES INTEGRALS
C 
      DO IP = 1, MMAX
         DO J = IP, MMAX
            NR = J - IP + 1
            READ(17) VRE, VIM
            RES(NR,1) = RES(NR,1) + VRE
            RES(NR,2) = RES(NR,2) + VIM
         END DO
      END DO
C
C     AVERAGE OUT RES INTEGRALS
C
      DO I=1,MMAX
         VN = DFLOAT(MMAX-I+1)
         RES(I,1) = RES(I,1) / VN
         RES(I,2) = RES(I,2) / VN
      END DO
C
C3000 FORMAT((D15.6))
      RETURN
      END
C ------------------------------------------------
      DOUBLE PRECISION FUNCTION RUNINT(F,N,STEP)
C
C     Compute the following function V(x):
C
C     V(x) = int_0^x F(x') dx'
C
C     by using Simpson's rule.
C
C     F    - array of length N
C     N    - number of points in F; N needs to be odd
C     STEP - interval between arguments of F
C
C ----------------------------------------------------
C                                       21 Aug 2016
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION F(*)
C     INTENT(OUT) V
C.....number of points of integral evaluation
      M    =(N-1)/2
C.....running integral
      RUNINT = 0.0D+00
C.....prefactor
      PREF = STEP/3.0D+00
C
      DO I=1,M
         IA = 2*(I-1) + 1
         IB = 2* I    + 1
         IC = IA      + 1
         RUNINT = RUNINT + PREF * ( F(IA) + F(IB) + 4.0D+00*F(IC) )
      END DO
C
      RETURN
      END
C            IP = ((I-1) * (2*MMAX + 2 - I) )/2 + J

