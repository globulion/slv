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
c     RXX = ROT(1)
c     RYY = ROT(5)
c     RZZ = ROT(9)
c     RXY = ROT(2)
c     RYX = ROT(4)
c     RXZ = ROT(3)
c     RZX = ROT(7)
c     RYZ = ROT(6)
c     RZY = ROT(8)
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
                   XN = RXX * XO + RXY * YO + RXZ * ZO
                   YN = RYX * XO + RYY * YO + RYZ * ZO
                   ZN = RZX * XO + RZY * YO + RZZ * ZO
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
     &                   RXX * RXY * XYO + 
     &                   RXX * RXZ * XZO +
     &                   RXY * RXX * XYO +
     &                   RXY * RXY * YYO +
     &                   RXY * RXZ * YZO +
     &                   RXZ * RXX * XZO +
     &                   RXZ * RXY * YZO +
     &                   RXZ * RXZ * ZZO
C
                   YYN = RYX * RYX * XXO + 
     &                   RYX * RYY * XYO + 
     &                   RYX * RYZ * XZO +
     &                   RYY * RYX * XYO +
     &                   RYY * RYY * YYO +
     &                   RYY * RYZ * YZO +
     &                   RYZ * RYX * XZO +
     &                   RYZ * RYY * YZO +
     &                   RYZ * RYZ * ZZO
C
                   ZZN = RZX * RZX * XXO + 
     &                   RZX * RZY * XYO + 
     &                   RZX * RZZ * XZO +
     &                   RZY * RZX * XYO +
     &                   RZY * RZY * YYO +
     &                   RZY * RZZ * YZO +
     &                   RZZ * RZX * XZO +
     &                   RZZ * RZY * YZO +
     &                   RZZ * RZZ * ZZO
C
                   XYN = RXX * RYX * XXO + 
     &                   RXX * RYY * XYO + 
     &                   RXX * RYZ * XZO +
     &                   RXY * RYX * XYO +
     &                   RXY * RYY * YYO +
     &                   RXY * RYZ * YZO +
     &                   RXZ * RYX * XZO +
     &                   RXZ * RYY * YZO +
     &                   RXZ * RYZ * ZZO
C
                   XZN = RXX * RZX * XXO + 
     &                   RXX * RZY * XYO + 
     &                   RXX * RZZ * XZO +
     &                   RXY * RZX * XYO +
     &                   RXY * RZY * YYO +
     &                   RXY * RZZ * YZO +
     &                   RXZ * RZX * XZO +
     &                   RXZ * RZY * YZO +
     &                   RXZ * RZZ * ZZO
C
                   YZN = RYX * RZX * XXO + 
     &                   RYX * RZY * XYO + 
     &                   RYX * RZZ * XZO +
     &                   RYY * RZX * XYO +
     &                   RYY * RZY * YYO +
     &                   RYY * RZZ * YZO +
     &                   RYZ * RZX * XZO +
     &                   RYZ * RZY * YZO +
     &                   RYZ * RZZ * ZZO
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
