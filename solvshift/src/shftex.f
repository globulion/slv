C-----|--|---------|---------|---------|---------|---------|---------|--|------|
C
C      filename: shftex.f
C
C                       COARSE-GRAIN SOLVATOCHROMIC 
C                EXCHANGE-REPULSION FREQUENCY SHIFT THEORY
C
C                      version 1.0a   29 Sep 2014    Bartosz Błasiak
C -----------------------------------------------------------------------------
C
C -----------------------------------------------------------------------------
      SUBROUTINE SHFTEX(REDMSS,FREQ,GIJJ,LVEC,RIA,RIB,RNA,RNB,RIA1,
     &                  CIKA,CIKB,CIKA1,SKM,TKM,SK1M,TK1M,ZA,ZB,
     &                  NBSA,NBSB,NMOSA,NMOSB,NATA,NATB,NMODES,MLIST,
     &                  FAIJ,FBIJ,FAIJ1,MODEID,SIJ,TIJ,SMIJ,TMIJ,FI,
     &                  SHFTMA,SHFTEA)
C -----------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION REDMSS(NMODES),FREQ(NMODES),GIJJ(NMODES),
     &          RIA(NMOSA*3),RIB(NMOSB*3),RNA(NATA*3),RNB(NATB*3),
     &          RIA1(NMODES*NMOSA*3),CIKA(NMOSA*NBSA),CIKB(NMOSB*NBSB),
     &          CIKA1(NMODES*NMOSA*NBSA),SKM(NBSA*NBSB),TKM(NBSA*NBSB),
     &          SK1M(NBSA*NBSB*3),TK1M(NBSA*NBSB*3),MLIST(NBSA),
     &          FAIJ(NMOSA*NMOSA),FBIJ(NMOSB*NMOSB),
     &          FAIJ1(NMODES*NMOSA*NMOSA),ZA(NATA),ZB(NATB),
     &          SIJ(NMOSA*NMOSB),SMIJ(NMODES*NMOSA*NMOSB),
     &          TIJ(NMOSA*NMOSB),TMIJ(NMODES*NMOSA*NMOSB),
     &          FI(NMODES)
      DOUBLE PRECISION LVEC(NMODES*NATA*3)
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     &           FOUR=4.0D+00,FIVE=5.0D+00)
Cf2py INTENT(OUT) SHFTMA,SHFTEA
C                                                               
C     calculate <I|O|J> integrals and their derivatives dMM     
C     Operator O is 1 for overlap and -HALF*LAPLACIAN for kinetic
C                                                               
      CALL CALCIJ(NMOSA,NMOSB,NMODES,NATA,NBSA,NBSB,SKM,TKM,
     &            MLIST,SK1M,TK1M,LVEC,CIKA,CIKB,CIKA1,
     &            SIJ,TIJ,SMIJ,TMIJ)
C                                                             
C     calculate first derivatives of exchange-repulsion energy
C     with respect to normal coordinates MM                   
C                                                             
      CALL CALCFI(NMOSA,NMOSB,NMODES,NATA,NATB,ZA,ZB,LVEC,
     &            RIA,RIB,RIA1,RNA,RNB,FAIJ,FBIJ,FAIJ1,
     &            SIJ,TIJ,SMIJ,TMIJ,
     &            FI)
C                                                    
C     calculate exchange-repulsion interaction energy
C                                                    
c      CALL CALCEN
C
C     calculate mechanical and electronic frequency shift!!!
C
      SHFTMA = ZERO                                                    
      SHFTEA = ZERO
      DENOM  = TWO*REDMSS(MODEID)*FREQ(MODEID)                         
      DO 999 I=1,NMODES
         FREQI = FREQ(I)
         SHFTMA = SHFTMA + GIJJ(I) * FI(I) / (REDMSS(I)*FREQI*FREQI)
 999  CONTINUE                                                         
      SHFTMA = SHFTMA / (-ONE*DENOM)                                   
      SHFTEA = SHFTEA / DENOM                                          
C                                                                      
      RETURN                                                           
      END                                                              
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE CALCFI(NMOSA,NMOSB,NMODES,NATA,NATB,ZA,ZB,LVEC,
     &                  RIA,RIB,RIA1,RNA,RNB,FAIJ,FBIJ,FAIJ1,
     &                  SIJ,TIJ,SMIJ,TMIJ,
     &                  FI)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RIA(NMOSA*3),RIB(NMOSB*3),RNA(NATA*3),RNB(NATB*3),
     &          RIA1(NMODES*NMOSA*3),
     &          FAIJ(NMOSA*NMOSA),FBIJ(NMOSB*NMOSB),
     &          FAIJ1(NMODES*NMOSA*NMOSA),ZA(NATA),ZB(NATB),
     &          SIJ(NMOSA*NMOSB),SMIJ(NMODES*NMOSA*NMOSB),
     &          TIJ(NMOSA*NMOSB),TMIJ(NMODES*NMOSA*NMOSB),
     &          FI(NMODES)
      DOUBLE PRECISION LVEC(NMODES*NATA*3)
      DOUBLE PRECISION LVECX, LVECY, LVECZ
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     &           FOUR=4.0D+00,FIVE=5.0D+00,EIGHT=8.0D+00,
     &           ONEPI=3.141592653589793D+00,
     &           TWOPI=6.283185307179586D+00)
C
      ONEPIINV = ONE/ONEPI
      TWOPIINV = ONE/TWOPI
C
      DO I=1,NMOSA
C
C     --- UNPACKING
C 
      IX = 3*(I-1) + 1
C
      RIAIX = RIA(IX  )
      RIAIY = RIA(IX+1)
      RIAIZ = RIA(IX+2)
C
      DO J=1,NMOSB
C
         IJ = NMOSB*(I-1) + J
         JX = 3*(J-1) + 1
C
         SIJV = SIJ(IJ)
         TIJV = TIJ(IJ)
C
         SIJVLN = DLOG(DABS(SIJV))
         SIJVLNINV = ONE/SIJVLN
         SIJVSQ = SIJV*SIJV
C
C        DISTANCES
C
         RIBJX = RIB(JX  )
         RIBJY = RIB(JX+1)
         RIBJZ = RIB(JX+2)
C
         RIJX = RIAIX-RIBJX
         RIJY = RIAIY-RIBJY
         RIJZ = RIAIZ-RIBJZ
C
         RIJ = DSQRT(RIJX*RIJX+RIJY*RIJY+RIJZ*RIJZ)
         RIJINV = ONE / RIJ
C
         SDIVR = SIJV*RIJINV
C
C        ITERATE OVER NORMAL COORDINATES
C
         DO MM=1,NMODES                                    
C
            MIJ = NMOSA*NMOSB*(MM-1) + NMOSB*(I-1) + J
            MIX = NMOSA*3*(MM-1) + 3*(I-1) + 1
C                                                         
            SMIJV = SMIJ(MIJ)
            TMIJV = TMIJ(MIJ)
C
            RMIAX = RIA1(MIX  )
            RMIAY = RIA1(MIX+1)
            RMIAZ = RIA1(MIX+2)
C
            RMIJ = RIJINV * (RIJX*RMIAX + RIJY*RMIAY + RIJZ*RMIAZ)
C
C           +++ TS TERMS +++
C
            AUX = TWO*DSQRT(-TWO*SIJVLN*ONEPIINV)
C
            TS1 = FOUR*SDIVR*SMIJV*( DSQRT(-ONE*TWOPIINV*SIJVLNINV) - 
     &                               AUX - ONE )
C
            TS2 = TWO*SDIVR*SDIVR*RMIJ*( AUX + ONE )
C
            TS3 = FOUR*SMIJV*TIJV
C
            TS4 = FOUR*TMIJV*SIJV
C
C           +++ TA TERMS +++
C                                                         
C           SUM OVER C ---> MOS OF X
C                                                         
            TA1 = ZERO
            TA2 = ZERO
            TA3 = ZERO
            TA4 = ZERO
C
            DO K=1,NMOSA
C
C              --- UNPACKING
C 
               IK = NMOSA*(I-1) + K
               KX = 3*(K-1) + 1 
               KJ = NMOSB*(K-1) + J
               MIK = NMOSA*NMOSA*(MM-1) + NMOSA*(I-1) + K
               MKJ = NMOSA*NMOSB*(MM-1) + NMOSB*(K-1) + J
               MKX = NMOSA*3*(MM-1) + 3*(K-1) + 1
C                                  
               SKJV = SIJ(KJ)
               SMKJV = SMIJ(MKJ)
C                       
               RKJX = RIA(KX  )-RIBJX
               RKJY = RIA(KX+1)-RIBJY
               RKJZ = RIA(KX+2)-RIBJZ
C                                                         
               RKJ = DSQRT(RKJX*RKJX+RKJY*RKJY+RKJZ*RKJZ)
               RKJINV = ONE / RKJ
C
               RMKAX = RIA1(MKX  )
               RMKAY = RIA1(MKX+1)
               RMKAZ = RIA1(MKX+2)
C
               RMKJ = RKJINV * (RKJX*RMKAX + RKJY*RMKAY + RKJZ*RMKAZ)
C                                                         
               FAIJIK = FAIJ(IK)
               FAIJ1IK = FAIJ1(MIK)
C
               TA1 = TA1 - TWO*SMIJV*FAIJIK*SKJV + 
     &                   EIGHT*SIJV*RKJINV*SMIJV - 
     &                     TWO*SIJV*(FAIJ1IK*SKJV+FAIJIK*SMKJV) - 
     &                    FOUR*SIJVSQ*RKJINV*RKJINV*RMKJ
            ENDDO
C                                                         
C           SUM OVER D ---> MOS OF Y
C                                                         
            DO L=1,NMOSB
C
C              --- UNPACKING
C 
               JL = NMOSB*(J-1) + L
               LX = 3*(L-1) + 1
               IL = NMOSB*(I-1) + L
               MIL = NMOSA*NMOSB*(MM-1) + NMOSB*(I-1) + L
C
               SILV = SIJ(IL)
               SMILV = SMIJ(MIL)
C                                                         
               RILX = RIAIX-RIB(LX  )
               RILY = RIAIY-RIB(LX+1)
               RILZ = RIAIZ-RIB(LX+2)
C                                                         
               RIL = DSQRT(RILX*RILX+RILY*RILY+RILZ*RILZ)
               RILINV = ONE / RIL
C
               RMIL = RILINV * (RILX*RMIAX + RILY*RMIAY + RILZ*RMIAZ)
C                                                         
               FBIJJL = FBIJ(JL)
C
               TA2 = TA2 - TWO*SMIJV*FBIJJL*SILV + 
     &                   EIGHT*SIJV*RILINV*SMIJV -
     &                     TWO*SIJV*FBIJJL*SMILV -
     &                    FOUR*SIJVSQ*RILINV*RILINV*RMIL
            ENDDO
C                                                         
C           SUM OVER X ---> ATOMS OF X
C                                                         
            DO M=1,NATA
C
C              --- UNPACKING
C 
               MX = 3*(M-1) + 1
               MMX = NATA*3*(MM-1) + 3*(M-1) + 1
C
               ZM = ZA(M)
C
               RMJX = RNA(MX  )-RIBJX
               RMJY = RNA(MX+1)-RIBJY
               RMJZ = RNA(MX+2)-RIBJZ
C
               RMJ = DSQRT(RMJX*RMJX+RMJY*RMJY+RMJZ*RMJZ)
               RMJINV = ONE / RMJ
C
               LVECX = LVEC(MMX  )
               LVECY = LVEC(MMX+1)
               LVECZ = LVEC(MMX+2)
C
               RMMJ = RMJINV * (RMJX*LVECX+RMJY*LVECY+RMJZ*LVECZ)
C
               TA3 = TA3 + TWO*SIJVSQ*ZM*RMJINV*RMJINV*RMMJ - 
     &                    FOUR*SIJV*SMIJV*ZM*RMJINV
            ENDDO
C                                                         
C           SUM OVER Y ---> ATOMS OF Y
C                                                         
            DO N=1,NATB
C
C              --- UNPACKING
C 
               NX = 3*(N-1) + 1
C
               ZN = ZB(N)
C
               RINX = RIAIX-RNB(NX  )
               RINY = RIAIY-RNB(NX+1)
               RINZ = RIAIZ-RNB(NX+2)
C
               RIN = DSQRT(RINX*RINX+RINY*RINY+RINZ*RINZ)
               RININV = ONE / RIN
C
               RMIN = RININV * (RINX*RMIAX+RINY*RMIAY+RINZ*RMIAZ)
C
               TA4 = TA4 + TWO*SIJVSQ*ZN*RININV*RININV*RMIN -
     &                    FOUR*SIJV*SMIJV*ZN*RININV
            ENDDO
C
C           SAVE FORCE
C
            FI(MM) = FI(MM) + TS1+TS2+TS3+TS4 + TA1+TA2+TA3+TA4
C
         ENDDO
C
      ENDDO
      ENDDO
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|

      SUBROUTINE CALCIJ(NMOSA,NMOSB,NMODES,NATA,NBSA,NBSB,SKM,TKM,
     &                  MLIST,SK1M,TK1M,LVEC,CIKA,CIKB,CIKA1,
     &                  SIJ,TIJ,SMIJ,TMIJ)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SKM(NBSA*NBSB),TKM(NBSA*NBSB),MLIST(NBSA),
     &          SK1M(NBSA*NBSB*3),TK1M(NBSA*NBSB*3),
     &          CIKA(NMOSA*NBSA),CIKB(NMOSB*NBSB),
     &          CIKA1(NMODES*NMOSA*NBSA),
     &          SIJ(NMOSA*NMOSB),SMIJ(NMODES*NMOSA*NMOSB),
     &          TIJ(NMOSA*NMOSB),TMIJ(NMODES*NMOSA*NMOSB)
      DOUBLE PRECISION LVEC(NMODES*NATA*3)
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     &           FOUR=4.0D+00,FIVE=5.0D+00)
C
C     SIJ and TIJ 
C
C      DO I=1,NMOSA
C      DO J=1,NMOSB
C         IJ = NMOSB*(I-1) + J
C
C         SIJV = ZERO
C         TIJV = ZERO
C
C         DO K=1,NBSA
C         MLISTK = MLIST(K)
C         IK = NBSA*(I-1) + K
C
C         CIKAIK = CIKA(IK)
C         DO L=1,NBSB
C            KL = NBSB*(K-1) + L
C            JL = NBSB*(J-1) + L
C
C            CIKBJL = CIKB(JL)
C            COEFFS = CIKAIK * CIKBJL
C
C            SKMKL = SKM(KL)
C            TKMKL = TKM(KL)
C
C            SIJV = SIJV + COEFFS * SKMKL
C            TIJV = TIJV + COEFFS * TKMKL
C         ENDDO
C         ENDDO
C         SIJ(IJ) = SIJV
C         TIJ(IJ) = TIJV
C      ENDDO
C      ENDDO
C
C     SMIJ and TMIJ
C
      DO I=1,NMOSA                                     
      DO J=1,NMOSB
         IJ = NMOSB*(I-1) + J
C
         SIJV = ZERO
         TIJV = ZERO
C
         DO K=1,NBSA
         MLISTK = MLIST(K)
         IK = NBSA*(I-1) + K
C
         CIKAIK = CIKA(IK)
         DO L=1,NBSB
            KL = NBSB*(K-1) + L
            JL = NBSB*(J-1) + L
C
            CIKBJL = CIKB(JL)
            COEFFS = CIKAIK * CIKBJL
C
            SKMKL = SKM(KL)
            TKMKL = TKM(KL)
C
            SIJV = SIJV + COEFFS * SKMKL
            TIJV = TIJV + COEFFS * TKMKL
         ENDDO
         ENDDO
         SIJ(IJ) = SIJV
         TIJ(IJ) = TIJV
C
         DO M=1,NMODES
C
            MIJ = NMOSA*NMOSB*(M-1) + NMOSB*(I-1) + J
C                                                         
            SMIJV = ZERO
            TMIJV = ZERO
C                                                         
            DO K=1,NBSA
            MLISTK = MLIST(K)
            IK = NBSA*(I-1) + K
C
            CIKAIK = CIKA(IK)
C
            MX1 = NATA*3*(M-1) + 3*(MLISTK-1) + 1
C
            RL1 = LVEC(MX1  )
            RL2 = LVEC(MX1+1)
            RL3 = LVEC(MX1+2)
C
            DO L=1,NBSB
               KL = NBSB*(K-1) + L
               IK = NBSA*(I-1) + K
               JL = NBSB*(J-1) + L
C
               KL1 = NBSB*3*(K-1) + 3*(L-1) + 1
               KL2 = KL1 + 1
               KL3 = KL2 + 1
C
               MIK = NMOSA*NBSA*(M-1) + NBSA*(I-1) + K
C                                                         
               SKMKL = SKM(KL)
               TKMKL = TKM(KL)
C
               CIKBJL = CIKB(JL)
               COEFFS = CIKAIK * CIKBJL
C
               SK1MKL1 = SK1M(KL1)
               SK1MKL2 = SK1M(KL2)
               SK1MKL3 = SK1M(KL3)
               TK1MKL1 = TK1M(KL1)
               TK1MKL2 = TK1M(KL2)
               TK1MKL3 = TK1M(KL3)
C
               CIKA1M = CIKA1(MIK)
C                                                        
               S1 = RL1*SK1MKL1 + RL2*SK1MKL2 + RL3*SK1MKL3
               S2 = RL1*TK1MKL1 + RL2*TK1MKL2 + RL3*TK1MKL3
C
               SMIJV = SMIJV + COEFFS * S1 + CIKA1M * CIKBJL * SKMKL
               TMIJV = TMIJV + COEFFS * S2 + CIKA1M * CIKBJL * TKMKL
C                                                         
            ENDDO
            ENDDO
            SMIJ(MIJ) = SMIJV
            TMIJ(MIJ) = TMIJV
         ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C-----|--|---------|---------|---------|---------|---------|---------|--|------|
