C-----|--|---------|---------|---------|---------|---------|---------|--|------|
C
C      filename: efprot.f
C
C                              THE ROTATION 
C                OF EFFECTIVE FRAGMENT POTENTIAL PARAMETERS
C
C                               version 0.0a    28 Oct 2013    Bartosz Błasiak
C NOTES:
C
C            R - Rotation matrix
C      
C                XX XY XZ     1  2  3
C                YX YY YZ     4  5  6 
C                ZX ZY ZZ     7  8  9
C
C            Turn of indices in AOs (PyQuante convention):
C
C            P - shell:       X  Y  Z
C            D - shell:       XX YY ZZ XY XZ YZ
C -----------------------------------------------------------------------------
      SUBROUTINE VECROT(NMOS,NBASIS,VEC,ROT,ITYP)
C -----------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VEC(NMOS,NBASIS),ROT(3,3),ITYP(NBASIS)
Cf2py INTENT(IN,OUT) VEC
      RXX = ROT(1,1)
      RYY = ROT(2,2)
      RZZ = ROT(3,3)
      RXY = ROT(1,2)
      RYX = ROT(2,1)
      RXZ = ROT(1,3)
      RZX = ROT(3,1)
      RYZ = ROT(2,3)
      RZY = ROT(3,2)
C
      NP = 0
      ND = 0
      DO 1010 K=1,NBASIS
         ITYPK = ITYP(K)
C
C        P-SHELL
C
         IF (ITYPK.EQ.1) THEN
             NP = NP + 1
             IF ((MOD((NP+2),3)).EQ.0) THEN
                DO 2010 I=1,NMOS
                   XO = VEC(I,K  )
                   YO = VEC(I,K+1)
                   ZO = VEC(I,K+2)
C
                   XN = RXX * XO + RYX * YO + RZX * ZO
                   YN = RXY * XO + RYY * YO + RZY * ZO
                   ZN = RXZ * XO + RYZ * YO + RZZ * ZO
C
                   VEC(I,K  ) = XN
                   VEC(I,K+1) = YN
                   VEC(I,K+2) = ZN
 2010           CONTINUE
             ENDIF
C
C        D-SHELL
C
         ELSE IF (ITYPK.EQ.2) THEN
             ND = ND + 1
             IF ((MOD((NP+5),6)).EQ.0) THEN
                DO 3010 I=1,NMOS
                   XXO = VEC(I,K  )
                   YYO = VEC(I,K+1)
                   ZZO = VEC(I,K+2)
                   XYO = VEC(I,K+3)
                   XZO = VEC(I,K+4)
                   YZO = VEC(I,K+5)
C
                   XXN = RXX * RXX * XXO + 
     &                   RXX * RYX * XYO + 
     &                   RXX * RZX * XZO +
     &                   RYX * RXX * XYO +
     &                   RYX * RYX * YYO +
     &                   RYX * RZX * YZO +
     &                   RZX * RXX * XZO +
     &                   RZX * RYX * YZO +
     &                   RZX * RZX * ZZO
C
                   YYN = RXY * RXY * XXO + 
     &                   RXY * RYY * XYO + 
     &                   RXY * RZY * XZO +
     &                   RYY * RXY * XYO +
     &                   RYY * RYY * YYO +
     &                   RYY * RZY * YZO +
     &                   RZY * RXY * XZO +
     &                   RZY * RYY * YZO +
     &                   RZY * RZY * ZZO
C
                   ZZN = RXZ * RXZ * XXO + 
     &                   RXZ * RYZ * XYO + 
     &                   RXZ * RZZ * XZO +
     &                   RYZ * RXZ * XYO +
     &                   RYZ * RYZ * YYO +
     &                   RYZ * RZZ * YZO +
     &                   RZZ * RXZ * XZO +
     &                   RZZ * RYZ * YZO +
     &                   RZZ * RZZ * ZZO
C
                   XYN = RXX * RXY * XXO + 
     &                   RXX * RYY * XYO + 
     &                   RXX * RZY * XZO +
     &                   RYX * RXY * XYO +
     &                   RYX * RYY * YYO +
     &                   RYX * RZY * YZO +
     &                   RZX * RXY * XZO +
     &                   RZX * RYY * YZO +
     &                   RZX * RZY * ZZO
C
                   XZN = RXX * RXZ * XXO + 
     &                   RXX * RYZ * XYO + 
     &                   RXX * RZZ * XZO +
     &                   RYX * RXZ * XYO +
     &                   RYX * RYZ * YYO +
     &                   RYX * RZZ * YZO +
     &                   RZX * RXZ * XZO +
     &                   RZX * RYZ * YZO +
     &                   RZX * RZZ * ZZO
C
                   YZN = RXY * RXZ * XXO + 
     &                   RXY * RYZ * XYO + 
     &                   RXY * RZZ * XZO +
     &                   RYY * RXZ * XYO +
     &                   RYY * RYZ * YYO +
     &                   RYY * RZZ * YZO +
     &                   RZY * RXZ * XZO +
     &                   RZY * RYZ * YZO +
     &                   RZY * RZZ * ZZO
C
                   VEC(I,K  ) = XXN
                   VEC(I,K+1) = YYN
                   VEC(I,K+2) = ZZN
                   VEC(I,K+3) = XYN
                   VEC(I,K+4) = XZN
                   VEC(I,K+5) = YZN
 3010           CONTINUE
             ENDIF
         ENDIF
C
 1010 CONTINUE
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|
      SUBROUTINE VC1ROT(NMODES,NMOS,NBASIS,VEC1,ROT,ITYP)
C
C          ROTATE VECL1 TENSORS
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VEC1(NMODES,NMOS,NBASIS),ROT(3,3),ITYP(NBASIS)
Cf2py INTENT(IN,OUT) VEC1
      RXX = ROT(1,1)
      RYY = ROT(2,2)
      RZZ = ROT(3,3)
      RXY = ROT(1,2)
      RYX = ROT(2,1)
      RXZ = ROT(1,3)
      RZX = ROT(3,1)
      RYZ = ROT(2,3)
      RZY = ROT(3,2)
C
      NP = 0
      ND = 0
      DO 1010 K=1,NBASIS
         ITYPK = ITYP(K)
C
C        P-SHELL
C
         IF (ITYPK.EQ.1) THEN
             NP = NP + 1
             IF ((MOD((NP+2),3)).EQ.0) THEN
                DO 2010 M=1,NMODES
                DO 2010 I=1,NMOS
                   XO = VEC1(M,I,K  )
                   YO = VEC1(M,I,K+1)
                   ZO = VEC1(M,I,K+2)
C
                   XN = RXX * XO + RYX * YO + RZX * ZO
                   YN = RXY * XO + RYY * YO + RZY * ZO
                   ZN = RXZ * XO + RYZ * YO + RZZ * ZO
C
                   VEC1(M,I,K  ) = XN
                   VEC1(M,I,K+1) = YN
                   VEC1(M,I,K+2) = ZN
 2010           CONTINUE
             ENDIF
C
C        D-SHELL
C
         ELSE IF (ITYPK.EQ.2) THEN
             ND = ND + 1
             IF ((MOD((NP+5),6)).EQ.0) THEN
                DO 3010 M=1,NMODES
                DO 3010 I=1,NMOS
                   XXO = VEC1(M,I,K  )
                   YYO = VEC1(M,I,K+1)
                   ZZO = VEC1(M,I,K+2)
                   XYO = VEC1(M,I,K+3)
                   XZO = VEC1(M,I,K+4)
                   YZO = VEC1(M,I,K+5)
C
                   XXN = RXX * RXX * XXO + 
     &                   RXX * RYX * XYO + 
     &                   RXX * RZX * XZO +
     &                   RYX * RXX * XYO +
     &                   RYX * RYX * YYO +
     &                   RYX * RZX * YZO +
     &                   RZX * RXX * XZO +
     &                   RZX * RYX * YZO +
     &                   RZX * RZX * ZZO
C
                   YYN = RXY * RXY * XXO + 
     &                   RXY * RYY * XYO + 
     &                   RXY * RZY * XZO +
     &                   RYY * RXY * XYO +
     &                   RYY * RYY * YYO +
     &                   RYY * RZY * YZO +
     &                   RZY * RXY * XZO +
     &                   RZY * RYY * YZO +
     &                   RZY * RZY * ZZO
C
                   ZZN = RXZ * RXZ * XXO + 
     &                   RXZ * RYZ * XYO + 
     &                   RXZ * RZZ * XZO +
     &                   RYZ * RXZ * XYO +
     &                   RYZ * RYZ * YYO +
     &                   RYZ * RZZ * YZO +
     &                   RZZ * RXZ * XZO +
     &                   RZZ * RYZ * YZO +
     &                   RZZ * RZZ * ZZO
C
                   XYN = RXX * RXY * XXO + 
     &                   RXX * RYY * XYO + 
     &                   RXX * RZY * XZO +
     &                   RYX * RXY * XYO +
     &                   RYX * RYY * YYO +
     &                   RYX * RZY * YZO +
     &                   RZX * RXY * XZO +
     &                   RZX * RYY * YZO +
     &                   RZX * RZY * ZZO
C
                   XZN = RXX * RXZ * XXO + 
     &                   RXX * RYZ * XYO + 
     &                   RXX * RZZ * XZO +
     &                   RYX * RXZ * XYO +
     &                   RYX * RYZ * YYO +
     &                   RYX * RZZ * YZO +
     &                   RZX * RXZ * XZO +
     &                   RZX * RYZ * YZO +
     &                   RZX * RZZ * ZZO
C
                   YZN = RXY * RXZ * XXO + 
     &                   RXY * RYZ * XYO + 
     &                   RXY * RZZ * XZO +
     &                   RYY * RXZ * XYO +
     &                   RYY * RYZ * YYO +
     &                   RYY * RZZ * YZO +
     &                   RZY * RXZ * XZO +
     &                   RZY * RYZ * YZO +
     &                   RZY * RZZ * ZZO
C
                   VEC1(M,I,K  ) = XXN
                   VEC1(M,I,K+1) = YYN
                   VEC1(M,I,K+2) = ZZN
                   VEC1(M,I,K+3) = XYN
                   VEC1(M,I,K+4) = XZN
                   VEC1(M,I,K+5) = YZN
 3010           CONTINUE
             ENDIF
         ENDIF
C
 1010 CONTINUE
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|
