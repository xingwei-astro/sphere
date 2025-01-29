      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
*       
      COMMON /FIELD/ GC(KB2, LB, 0:MB) ,  GS(KB2, LB, 0:MB)
     &            ,  HC(KB2, LB, 0:MB) ,  HS(KB2, LB, 0:MB)
*
      COMMON /FLOW/  EC(KU2, LU, 0:MU) ,  ES(KU2, LU, 0:MU)
     &            ,  FC(KU4, LU, 0:MU) ,  FS(KU4, LU, 0:MU)
*
      DIMENSION FP(LB+LU, 181) , FT(KB2+KU4, 201)
     &         ,FTE1(KB2), FTE2(KB2)
*
      PI=DACOS(-1.D0)
*
      READ(3) EC
      READ(3) ES
      READ(3) FC
      READ(3) FS
      CLOSE(3)
      READ(4) GC
      READ(4) GS
      READ(4) HC
      READ(4) HS
      CLOSE(4)
*
      DO J=1,181
       TH=DFLOAT(J-1)
       IF ( J .EQ.   1 ) TH=000.001D0
       IF ( J .EQ. 181 ) TH=179.999D0
       C=DCOS( TH * PI/180.D0 )
       DO L=1,LB+LU
        FP(L,J)=-DFLOAT(L*(L+1)) / DFLOAT(2*L+1)
     &         * (PPP(L+1,0,C) - PPP(L-1,0,C))
       ENDDO
      ENDDO
      DO I=1,201
       X=-1.D0 + 2.D0*DFLOAT(I-1)/200.D0
       DO K=1,KB2+KU4
        FT(K,I)=TTT(0,K-1,X)
       ENDDO
      ENDDO
*
*     This loop does the calculation in real space.
      GMAX=0.D0
      HMAX=0.D0
      EMAX=0.D0
      FMAX=0.D0
      DO J=1,181
       TH=DFLOAT(J-1)
       IF ( J .EQ.   1 ) TH=000.001D0
       IF ( J .EQ. 181 ) TH=179.999D0
       S=DSIN( TH * PI/180.D0 )
       C=DCOS( TH * PI/180.D0 )
       DO I=1,201
        R=R1  +  (R2 - R1)*DFLOAT(I-1)/200.D0
        XX=R * S
        YY=R * C
        GG=0.D0
        HH=0.D0
        EE=0.D0
        FF=0.D0
*
        DO L=1,LB
         DO K=1,KB2
          GG=GG + GC(K,L,0) * FP(L,J)*FT(K,I)
          HH=HH - HC(K,L,0) * FP(L,J)*FT(K,I)
         ENDDO
        ENDDO
        GG=GG/R/S
*
        DO L=1,LU
         DO K=1,KU2
          EE=EE + EC(K,L,0) * FP(L,J)*FT(K,I)
         ENDDO
         DO K=1,KU4
          FF=FF - FC(K,L,0) * FP(L,J)*FT(K,I)
         ENDDO
        ENDDO
        EE=EE/R**2/S**2
*
        IF ( DABS(GG) .GT. GMAX ) GMAX=DABS(GG)
        IF ( DABS(HH) .GT. HMAX ) HMAX=DABS(HH)
        IF ( DABS(EE) .GT. EMAX ) EMAX=DABS(EE)
        IF ( DABS(FF) .GT. FMAX ) FMAX=DABS(FF)
        WRITE(11,10) XX, YY
        WRITE(12,10) GG, HH
        WRITE(13,10) EE, FF
       ENDDO
      ENDDO
      WRITE(6,*) GMAX, HMAX, EMAX, FMAX
*
*     This computes stream function of potential field in the exterior
      R3=0.0D0
      R4=2.0D0
      DO K=1,KB2
       FTE1(K)=TTT(0,K-1,-1.D0)
       FTE2(K)=TTT(0,K-1, 1.D0)
      ENDDO
      DO J=1,181
       TH=DFLOAT(J-1)
       IF ( J .EQ.   1 ) TH=000.001D0
       IF ( J .EQ. 181 ) TH=179.999D0
       S=DSIN( TH * PI/180.D0 )
       C=DCOS( TH * PI/180.D0 )
       DO I=1,101
        RE1=R3  +  (R1 - R3)*DFLOAT(I-1)/100.D0
        RE2=R2  +  (R4 - R2)*DFLOAT(I-1)/100.D0
        XXE1=RE1 * S
        YYE1=RE1 * C
        XXE2=RE2 * S
        YYE2=RE2 * C
        HHE1=0.D0
        HHE2=0.D0
        DO L=1,LB
         DO K=1,KB2
         HHE1=HHE1 - HC(K,L,0)*FTE1(K)*(RE1/R1)**DFLOAT(L+1)*FP(L,J)
         HHE2=HHE2 - HC(K,L,0)*FTE2(K)*(R2/RE2)**DFLOAT(L)  *FP(L,J)
         ENDDO
        ENDDO
        WRITE(14,20) XXE1, YYE1, HHE1
        WRITE(15,20) XXE2, YYE2, HHE2
       ENDDO
      ENDDO
*
      IF(MOD(LN,2).EQ.1) THEN
       CALL EQU
      ELSE
       WRITE(6,*) 'LN MUST BE ODD'
      ENDIF
*
      IF(MOD(LN,2).EQ.1) THEN
       CALL HELICITY
      ELSE
       WRITE(6,*) 'LN MUST BE ODD'
      ENDIF
*
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)
10    FORMAT(2E16.8)
20    FORMAT(3E16.8)
      STOP
      END
*******************************************************************************
*     PPP(N,M,C) is the N,M-th associated Legendre function evaluated at
*     C, using the recursion relation below from Abramowitz & Stegun.
*
      DOUBLE PRECISION FUNCTION PPP(N,M,C)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION P(0:2000)
      IF ( N .LT. M ) THEN
       PPP=0.D0
       RETURN
      END IF
      S=DSQRT(1.D0 - C**2)
      P(M)=1.D0
      DO 10 I=1,M
       P(M)=P(M) * (-0.5D0) * S * DFLOAT(I+M)
10    CONTINUE
      P(M+1)=DFLOAT(2*M+1) * C * P(M)
      DO 20 I=M+2,N
       P(I)=(DFLOAT(2*I-1)*C*P(I-1) - DFLOAT(I+M-1)*P(I-2))/DFLOAT(I-M)
20    CONTINUE
      PPP=P(N)
      RETURN
      END
*******************************************************************************
      DOUBLE PRECISION FUNCTION TTT(K,M,X)
*
*     TTT(K,M,X) = the K-th derivative of Tm(X),
*     the M-th Chebyshev polynomial evaluated at X.
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(0:2000), B(0:2000)
      DO 10 I=0,M-1
       A(I)=0.D0
10    CONTINUE
      A(M)=1.D0
      DO 20 J=1,K
       CALL DIFF(A,M)
20    CONTINUE
      B(M+2)=0.D0
      B(M+1)=0.D0
      DO 30 I=M,0,-1
       B(I)=2.D0 * X * B(I+1)  -  B(I+2)  +  A(I)
30    CONTINUE
      TTT=(A(0) + B(0) - B(2))/2.D0
      RETURN
      END
      SUBROUTINE DIFF(A,M)
      DOUBLE PRECISION A(0:2000), C(0:2000)
      C(M+1)=0.D0
      C(M)=0.D0
      DO 10 I=M-1,0,-1
       C(I)=C(I+2)  +  DFLOAT(2*I+2) * A(I+1)
10    CONTINUE
      C(0)=C(0)/2.D0
      DO 20 I=0,M
       A(I)=C(I)
20    CONTINUE
      RETURN
      END
*******************************************************************************
      SUBROUTINE EQU
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
*
      COMMON /FIELD/ GC(KB2, LB, 0:MB) ,  GS(KB2, LB, 0:MB)
     &            ,  HC(KB2, LB, 0:MB) ,  HS(KB2, LB, 0:MB)
*
      COMMON /FLOW/  EC(KU2, LU, 0:MU) ,  ES(KU2, LU, 0:MU)
     &            ,  FC(KU4, LU, 0:MU) ,  FS(KU4, LU, 0:MU)
*
      COMMON /NLIN/ TT1(KN, KN) , TT2(KN, KN) , TT3(KN, KN)
     &            , TT4(KN, KN) , TT5(KN, KN) , RNV(KN)
     &            , PP1(0:LL, 0:MM, LN) , PP2(0:LL, 0:MM, LN)
     &            , PP3(0:LL, 0:MM, LN) 
     &            , CTH(LN), STH(LN), CSC(LN) , XLL(LL)
     &            , CPHI(2*MN), SPHI(2*MN)
     &            , CS(0:MM, 2*MN), SN(0:MM, 2*MN)
     &            ,CSM(0:MM, 2*MN),SNM(0:MM, 2*MN)
     &            , CT(2*MN, 0:MM), ST(2*MN, 0:MM)
     &            , CALPHA, SALPHA
*
      COMMON /CURL0/ TTT1(KN, KN) , PPP1(LN, LN, 0:MM)
     &             , TTT3(KN, KN) , PPP2(LN, LN, 0:MM)
     &             , TTT4(KN, KN) , PPP4(LN, LN, 0:MM)
*
      COMMON /DTCS/ TTT(KN, KN) , PPP(LN, 0:LN, 0:MM)
*           
      COMMON /OUT/ B1(2*MN, LN, KN), B2(2*MN, LN, KN), B3(2*MN, LN, KN)
     &           , U1(2*MN, LN, KN), U2(2*MN, LN, KN), U3(2*MN, LN, KN)
*
      DIMENSION WHC1(LB, 0:MB) , WHC2(LB, 0:MB) , WGC3(LB, 0:MB)
      DIMENSION WHS1(LB, 0:MB) , WHS2(LB, 0:MB) , WGS3(LB, 0:MB)
*
      DIMENSION WHC11(0:MB), WHS11(0:MB), WHC22(0:MB), WHS22(0:MB)
      DIMENSION WHC23(0:MB), WHS23(0:MB), WGC32(0:MB), WGS32(0:MB)
      DIMENSION WGC33(0:MB), WGS33(0:MB)
*
      DIMENSION WFC1(LU, 0:MU) , WFC2(LU, 0:MU) , WEC3(LU, 0:MU)
      DIMENSION WFS1(LU, 0:MU) , WFS2(LU, 0:MU) , WES3(LU, 0:MU)
*
      DIMENSION WFC11(0:MU), WFS11(0:MU), WFC22(0:MU), WFS22(0:MU)
      DIMENSION WFC23(0:MU), WFS23(0:MU), WEC32(0:MU), WES32(0:MU)
      DIMENSION WEC33(0:MU), WES33(0:MU)
******
      READ(2) TT1
      READ(2) TT2
      READ(2) TT3
      READ(2) TT4
      READ(2) TT5
      READ(2) RNV
      DO 51 M=0,MM
       DO 51 L=M,LL
        READ(2) (PP1(L,M,IL) , IL=1,LN)
        READ(2) (PP2(L,M,IL) , IL=1,LN)
        READ(2) (PP3(L,M,IL) , IL=1,LN)
51    CONTINUE
      READ(2) CTH
      READ(2) STH
      READ(2) CSC
      READ(2) XLL

      READ(2) CPHI
      READ(2) SPHI
      READ(2) CS
      READ(2) SN
      READ(2) CSM
      READ(2) SNM
      READ(2) CT
      READ(2) ST
*
      READ(2) TTT1
      READ(2) TTT3
      READ(2) TTT4
      DO 61 IL=1,LN
       DO 61 L=1,LN
        READ(2) (PPP1(IL,L,M) , M=0,MM)
        READ(2) (PPP2(IL,L,M) , M=0,MM)
        READ(2) (PPP4(IL,L,M) , M=0,MM)
61    CONTINUE
*
      READ(2) TTT
      DO 71 IL=1,LN
       DO 71 L=0,LN
        READ(2) (PPP(IL,L,M) , M=0,MM)
71    CONTINUE
*
      READ(2) CALPHA
      READ(2) SALPHA
*
*     Calculate B
      DO 20 IK=1,KN
       DO 10 M=0,MB
        LS=MAX(M,1)
        DO 10 L=LS,LB
         WHC1(L,M)=0.D0
         WHS1(L,M)=0.D0
         WHC2(L,M)=0.D0
         WHS2(L,M)=0.D0
         WGC3(L,M)=0.D0
         WGS3(L,M)=0.D0
         DO 10 K=1,KB2
          WHC1(L,M)=WHC1(L,M)  +  HC(K,L,M) * TT1(K,IK)
          WHS1(L,M)=WHS1(L,M)  +  HS(K,L,M) * TT1(K,IK)
          WHC2(L,M)=WHC2(L,M)  +  HC(K,L,M) * TT2(K,IK)
          WHS2(L,M)=WHS2(L,M)  +  HS(K,L,M) * TT2(K,IK)
          WGC3(L,M)=WGC3(L,M)  +  GC(K,L,M) * TT3(K,IK)
          WGS3(L,M)=WGS3(L,M)  +  GS(K,L,M) * TT3(K,IK)
10     CONTINUE
       DO 20 IL=1,LN
        DO 15 M=0,MB
         LS=MAX(M,1)
         WHC11(M)=0.D0
         WHS11(M)=0.D0
         WHC22(M)=0.D0
         WHS22(M)=0.D0
         WHC23(M)=0.D0
         WHS23(M)=0.D0
         WGC32(M)=0.D0
         WGS32(M)=0.D0
         WGC33(M)=0.D0
         WGS33(M)=0.D0
         DO 15 L=LS,LB
          WHC11(M)=WHC11(M)  +  WHC1(L,M) * PP1(L,M,IL)
          WHS11(M)=WHS11(M)  +  WHS1(L,M) * PP1(L,M,IL)
          WHC22(M)=WHC22(M)  +  WHC2(L,M) * PP2(L,M,IL)
          WHS22(M)=WHS22(M)  +  WHS2(L,M) * PP2(L,M,IL)
          WHC23(M)=WHC23(M)  +  WHC2(L,M) * PP3(L,M,IL)
          WHS23(M)=WHS23(M)  +  WHS2(L,M) * PP3(L,M,IL)
          WGC32(M)=WGC32(M)  +  WGC3(L,M) * PP2(L,M,IL)
          WGS32(M)=WGS32(M)  +  WGS3(L,M) * PP2(L,M,IL)
          WGC33(M)=WGC33(M)  +  WGC3(L,M) * PP3(L,M,IL)
          WGS33(M)=WGS33(M)  +  WGS3(L,M) * PP3(L,M,IL)
15      CONTINUE
        DO 20 IM=1,2*MN
         B1(IM,IL,IK)=0.D0
         B2(IM,IL,IK)=0.D0
         B3(IM,IL,IK)=0.D0
         DO 20 M=0,MB
          B1(IM,IL,IK)=B1(IM,IL,IK)
     &                + WHC11(M)* CS(M,IM) + WHS11(M)* SN(M,IM)
          B2(IM,IL,IK)=B2(IM,IL,IK)
     &                + WHC22(M)* CS(M,IM) + WHS22(M)* SN(M,IM)
     &                - WGC33(M)*SNM(M,IM) + WGS33(M)*CSM(M,IM)
          B3(IM,IL,IK)=B3(IM,IL,IK)
     &                - WHC23(M)*SNM(M,IM) + WHS23(M)*CSM(M,IM)
     &                - WGC32(M)* CS(M,IM) - WGS32(M)* SN(M,IM)
20    CONTINUE
*
*     Calculate u
      DO 60 IK=1,KN
       DO 50 M=0,MU
        LS=MAX(M,1)
        DO 50 L=LS,LU
         WEC3(L,M)=0.D0
         WES3(L,M)=0.D0
         DO K=1,KU2
          WEC3(L,M)=WEC3(L,M)  +  EC(K,L,M) * TT3(K,IK)
          WES3(L,M)=WES3(L,M)  +  ES(K,L,M) * TT3(K,IK)
         ENDDO
         WFC1(L,M)=0.D0
         WFS1(L,M)=0.D0
         WFC2(L,M)=0.D0
         WFS2(L,M)=0.D0
         DO K=1,KU4
          WFC1(L,M)=WFC1(L,M)  +  FC(K,L,M) * TT1(K,IK)
          WFS1(L,M)=WFS1(L,M)  +  FS(K,L,M) * TT1(K,IK)
          WFC2(L,M)=WFC2(L,M)  +  FC(K,L,M) * TT2(K,IK)
          WFS2(L,M)=WFS2(L,M)  +  FS(K,L,M) * TT2(K,IK)
         ENDDO
50     CONTINUE
       DO 60 IL=1,LN
        DO 55 M=0,MU
         LS=MAX(M,1)
         WFC11(M)=0.D0
         WFS11(M)=0.D0
         WFC22(M)=0.D0
         WFS22(M)=0.D0
         WFC23(M)=0.D0
         WFS23(M)=0.D0
         WEC32(M)=0.D0
         WES32(M)=0.D0
         WEC33(M)=0.D0
         WES33(M)=0.D0
         DO 55 L=LS,LU
          WFC11(M)=WFC11(M)  +  WFC1(L,M) * PP1(L,M,IL)
          WFS11(M)=WFS11(M)  +  WFS1(L,M) * PP1(L,M,IL)
          WFC22(M)=WFC22(M)  +  WFC2(L,M) * PP2(L,M,IL)
          WFS22(M)=WFS22(M)  +  WFS2(L,M) * PP2(L,M,IL)
          WFC23(M)=WFC23(M)  +  WFC2(L,M) * PP3(L,M,IL)
          WFS23(M)=WFS23(M)  +  WFS2(L,M) * PP3(L,M,IL)
          WEC32(M)=WEC32(M)  +  WEC3(L,M) * PP2(L,M,IL)
          WES32(M)=WES32(M)  +  WES3(L,M) * PP2(L,M,IL)
          WEC33(M)=WEC33(M)  +  WEC3(L,M) * PP3(L,M,IL)
          WES33(M)=WES33(M)  +  WES3(L,M) * PP3(L,M,IL)
55      CONTINUE
        DO 60 IM=1,2*MN
         U1(IM,IL,IK)=0.D0
         U2(IM,IL,IK)=0.D0
         U3(IM,IL,IK)=0.D0
         DO 60 M=0,MU
          U1(IM,IL,IK)=U1(IM,IL,IK)
     &                + WFC11(M)* CS(M,IM) + WFS11(M)* SN(M,IM)
          U2(IM,IL,IK)=U2(IM,IL,IK)
     &                + WFC22(M)* CS(M,IM) + WFS22(M)* SN(M,IM)
     &                - WEC33(M)*SNM(M,IM) + WES33(M)*CSM(M,IM)
          U3(IM,IL,IK)=U3(IM,IL,IK)
     &                - WFC23(M)*SNM(M,IM) + WFS23(M)*CSM(M,IM)
     &                - WEC32(M)* CS(M,IM) - WES32(M)* SN(M,IM)
60    CONTINUE
*
      PI=DACOS(-1.D0)
      DO IK=1,KN
       XX=-DCOS(DFLOAT(2*IK - 1)/DFLOAT(2*KN) * PI)
       R=(R2 + R1)/2.D0  +  (R2 - R1)/2.D0*XX
       DO IM=1,2*MN
        PHI=2.D0*PI*DFLOAT(IM)/DFLOAT(2*MN)
        XE=R*DCOS(PHI)
        YE=R*DSIN(PHI)
        WRITE(16,98) XE,YE
       ENDDO
      ENDDO
*
      DO IK=1,KN
       DO IM=1,2*MN
        WRITE(17,99) B1(IM,1,IK), B2(IM,1,IK), B3(IM,1,IK)
        WRITE(18,99) U1(IM,1,IK), U2(IM,1,IK), U3(IM,1,IK)
       ENDDO
      ENDDO
*
      CLOSE(16)
      CLOSE(17)
      CLOSE(18)
98    FORMAT(2E16.8)
99    FORMAT(3E16.8)
      RETURN
      END
*******************************************************************************
      SUBROUTINE HELICITY
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
*
      COMMON /FIELD/  GC(KB2, LB, 0:MB) ,  GS(KB2, LB, 0:MB)
     &             ,  HC(KB2, LB, 0:MB) ,  HS(KB2, LB, 0:MB)
*
      COMMON /FLOW/   EC(KU2, LU, 0:MU) ,  ES(KU2, LU, 0:MU)
     &             ,  FC(KU4, LU, 0:MU) ,  FS(KU4, LU, 0:MU)
*
      COMMON /NLIN/ TT1(KN, KN) , TT2(KN, KN) , TT3(KN, KN)
     &            , TT4(KN, KN) , TT5(KN, KN) , RNV(KN)
     &            , PP1(0:LL, 0:MM, LN) , PP2(0:LL, 0:MM, LN)
     &            , PP3(0:LL, 0:MM, LN) 
     &            , CTH(LN), STH(LN), CSC(LN) , XLL(LL)
     &            , CPHI(2*MN), SPHI(2*MN)
     &            , CS(0:MM, 2*MN), SN(0:MM, 2*MN)
     &            ,CSM(0:MM, 2*MN),SNM(0:MM, 2*MN)
     &            , CT(2*MN, 0:MM), ST(2*MN, 0:MM)
     &            , CALPHA, SALPHA
*
      COMMON /OUT/ B1(2*MN, LN, KN), B2(2*MN, LN, KN), B3(2*MN, LN, KN)
     &           , U1(2*MN, LN, KN), U2(2*MN, LN, KN), U3(2*MN, LN, KN)
*
      DIMENSION WFC1(LU, 0:MU) , WFC2(LU, 0:MU) , WEC3(LU, 0:MU)
      DIMENSION WFS1(LU, 0:MU) , WFS2(LU, 0:MU) , WES3(LU, 0:MU)
*
      DIMENSION WFC11(0:MU), WFS11(0:MU), WFC22(0:MU), WFS22(0:MU)
      DIMENSION WFC23(0:MU), WFS23(0:MU), WEC32(0:MU), WES32(0:MU)
      DIMENSION WEC33(0:MU), WES33(0:MU)
*
      DIMENSION CU1(2*MN, LN, KN), CU2(2*MN, LN, KN), CU3(2*MN, LN, KN)
      DIMENSION HELICT(2*MN, LN, KN)
*
      DIMENSION R(KN), TH(LN)
*
      DIMENSION TMP(LN)
******
*     Calculate curl u and helicity u . curl u
      DO 80 IK=1,KN
       DO 65 M=0,MU
        LS=MAX(M,1)
        DO 65 L=LS,LU
         WFC1(L,M)=0.D0
         WFS1(L,M)=0.D0
         WFC2(L,M)=0.D0
         WFS2(L,M)=0.D0
         DO K=1,KU2
          WFC1(L,M)=WFC1(L,M)  +  EC(K,L,M) * TT1(K,IK)
          WFS1(L,M)=WFS1(L,M)  +  ES(K,L,M) * TT1(K,IK)
          WFC2(L,M)=WFC2(L,M)  +  EC(K,L,M) * TT2(K,IK)
          WFS2(L,M)=WFS2(L,M)  +  ES(K,L,M) * TT2(K,IK)
         ENDDO
         WFC4     =0.D0
         WFS4     =0.D0
         WFC5     =0.D0
         WFS5     =0.D0
         DO K=1,KU4
          WFC4     =WFC4       +  FC(K,L,M) * TT4(K,IK)
          WFS4     =WFS4       +  FS(K,L,M) * TT4(K,IK)
          WFC5     =WFC5       +  FC(K,L,M) * TT5(K,IK)
          WFS5     =WFS5       +  FS(K,L,M) * TT5(K,IK)
         ENDDO
         WEC3(L,M)=XLL(L)*WFC4 - WFC5
         WES3(L,M)=XLL(L)*WFS4 - WFS5
65     CONTINUE
       DO 80 IL=1,LN
        DO 70 M=0,MU
         LS=MAX(M,1)
         WFC11(M)=0.D0
         WFS11(M)=0.D0
         WFC22(M)=0.D0
         WFS22(M)=0.D0
         WFC23(M)=0.D0
         WFS23(M)=0.D0
         WEC32(M)=0.D0
         WES32(M)=0.D0
         WEC33(M)=0.D0
         WES33(M)=0.D0
         DO 70 L=LS,LU
          WFC11(M)=WFC11(M)  +  WFC1(L,M) * PP1(L,M,IL)
          WFS11(M)=WFS11(M)  +  WFS1(L,M) * PP1(L,M,IL)
          WFC22(M)=WFC22(M)  +  WFC2(L,M) * PP2(L,M,IL)
          WFS22(M)=WFS22(M)  +  WFS2(L,M) * PP2(L,M,IL)
          WFC23(M)=WFC23(M)  +  WFC2(L,M) * PP3(L,M,IL)
          WFS23(M)=WFS23(M)  +  WFS2(L,M) * PP3(L,M,IL)
          WEC32(M)=WEC32(M)  +  WEC3(L,M) * PP2(L,M,IL)
          WES32(M)=WES32(M)  +  WES3(L,M) * PP2(L,M,IL)
          WEC33(M)=WEC33(M)  +  WEC3(L,M) * PP3(L,M,IL)
          WES33(M)=WES33(M)  +  WES3(L,M) * PP3(L,M,IL)
70      CONTINUE
        DO 80 IM=1,2*MN
         CU1(IM,IL,IK)=0.D0
         CU2(IM,IL,IK)=0.D0
         CU3(IM,IL,IK)=0.D0
         DO 75 M=0,MU
          CU1(IM,IL,IK)=CU1(IM,IL,IK)
     &            + WFC11(M)* CS(M,IM) + WFS11(M)* SN(M,IM)
          CU2(IM,IL,IK)=CU2(IM,IL,IK)
     &            + WFC22(M)* CS(M,IM) + WFS22(M)* SN(M,IM)
     &            - WEC33(M)*SNM(M,IM) + WES33(M)*CSM(M,IM)
          CU3(IM,IL,IK)=CU3(IM,IL,IK)
     &            - WFC23(M)*SNM(M,IM) + WFS23(M)*CSM(M,IM)
     &            - WEC32(M)* CS(M,IM) - WES32(M)* SN(M,IM)
75       CONTINUE
       HELICT(IM,IL,IK)=CU1(IM,IL,IK)*U1(IM,IL,IK)
     &                + CU2(IM,IL,IK)*U2(IM,IL,IK)
     &                + CU3(IM,IL,IK)*U3(IM,IL,IK) 
80    CONTINUE
*
* sort of helicity
      DO 90 IK=1,KN
       DO 90 IM=1,2*MN
        IF ( MOD(LN,2) .EQ. 0 ) THEN
         DO IL=1,LN/2
          TMP(IL)=HELICT(IM,LN-2*IL+1,IK)
          TMP(LN/2+IL)=HELICT(IM,2*IL,IK)
         ENDDO
        ELSE
         TMP((LN+1)/2)=HELICT(IM,1,IK)
         DO IL=1,(LN-1)/2
          TMP(IL)=HELICT(IM,LN-2*IL+1,IK)
          TMP((LN+1)/2+IL)=HELICT(IM,2*IL+1,IK)
         ENDDO
        ENDIF
        DO IL=1,LN
         HELICT(IM,IL,IK)=TMP(IL)
        ENDDO
90    CONTINUE
*
      PI=DACOS(-1.D0)
      DO IK=1,KN
       XX=-DCOS(DFLOAT(2*IK - 1)/DFLOAT(2*KN) * PI)
       R(IK)=(R2 + R1)/2.D0  +  (R2 - R1)/2.D0*XX
      ENDDO
*
      DO IL=1,LN
       TH(IL)=DACOS(CTH(IL))
      ENDDO
      IF ( MOD(LN,2) .EQ. 0 ) THEN
       DO IL=1,LN/2
        TMP(IL)=TH(LN-2*IL+1)
        TMP(LN/2+IL)=TH(2*IL)
       ENDDO
      ELSE
      TMP((LN+1)/2)=TH(1)
       DO IL=1,(LN-1)/2
        TMP(IL)=TH(LN-2*IL+1)
        TMP((LN+1)/2+IL)=TH(2*IL+1)
       ENDDO
      ENDIF
      DO IL=1,LN
       TH(IL)=TMP(IL)
      ENDDO
*
      DO IK=1,KN
       DO IL=1,LN
        HELCA = 0.D0
        DO IM=1,2*MN
         HELCA = HELCA + HELICT(IM,IL,IK)
        ENDDO
        HELCA = HELCA / DFLOAT(2*MN)
        X=R(IK)*DSIN(TH(IL))
        Y=R(IK)*DCOS(TH(IL))
        WRITE(19,'(3E16.8)') X, Y, HELCA
       ENDDO
      ENDDO
      CLOSE(19)
*
      RETURN
      END
