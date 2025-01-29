      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
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
      DIMENSION TMAT(KN,KN)
      DIMENSION TMT2(KN,KN) , TINV(KN,KN)
      DIMENSION TMT4(KN,KN) , TMT1(KN,KN)
      DIMENSION WK1(KN)
*
      DIMENSION PMT1(LN,LN)
      DIMENSION PMT2(LN,LN) , PMT3(LN,LN)
      DIMENSION PMT4(LN,LN) , PMT5(LN,LN)
      DIMENSION PNV1(LN,LN) , PNV2(LN,LN)
      DIMENSION WK2(LN)
*
      DIMENSION PPP1(LN, LN, 0:MM)
      DIMENSION PPP2(LN, LN, 0:MM)
      DIMENSION PPP4(LN, LN, 0:MM)
*
      DIMENSION PPPP(LN,0:LN,0:MM)
*
      PI=DACOS(-1.D0)
      DO 10 IK=1,KN
       XX=-DCOS(DFLOAT(2*IK - 1)/DFLOAT(2*KN) * PI)
       R=(R2 + R1)/2.D0  +  (R2 - R1)/2.D0*XX
       RNV(IK)=1.D0/R
       DO 10 K=1,KN
        TT1(K,IK)=TTT(0,K-1,XX) / R**2
        TT2(K,IK)=TTT(1,K-1,XX) / R     * DXDR
        TT3(K,IK)=TTT(0,K-1,XX) / R
        TT4(K,IK)=TTT(0,K-1,XX) / R**3
        TT5(K,IK)=TTT(2,K-1,XX) / R     * DXDR**2
10    CONTINUE
      WRITE(2) TT1
      WRITE(2) TT2
      WRITE(2) TT3
      WRITE(2) TT4
      WRITE(2) TT5
      WRITE(2) RNV
*
      CALL CMESH(LN,CTH)
      DO 15 IL=1,LN
       STH(IL)=DSQRT(1.D0 - CTH(IL)**2)
       CSC(IL)=1.D0/STH(IL)
       DO 15 M=0,MM
        DO 15 L=M,LL
         PP1(L,M,IL)=DFLOAT(L*(L+1)) * PPP(L,M,CTH(IL))
         PP2(L,M,IL)=(DFLOAT((L  )*(L-M+1)) * PPP(L+1,M,CTH(IL))
     &               -DFLOAT((L+M)*(L  +1)) * PPP(L-1,M,CTH(IL)))
     &                / DFLOAT(2*L+1) / STH(IL)
         PP3(L,M,IL)=PPP(L,M,CTH(IL)) / STH(IL)
         IF ( L .GT. 0 ) XLL(L)=DFLOAT(L*(L+1))
15    CONTINUE
      DO 20 M=0,MM
       DO 20 L=M,LL
        WRITE(2) (PP1(L,M,IL) , IL=1,LN)
        WRITE(2) (PP2(L,M,IL) , IL=1,LN)
        WRITE(2) (PP3(L,M,IL) , IL=1,LN)
20    CONTINUE
      WRITE(2) CTH
      WRITE(2) STH
      WRITE(2) CSC
      WRITE(2) XLL
*
      DO 25 IM=1,2*MN
       PHI=2.D0*PI*DFLOAT(IM)/DFLOAT(2*MN)
       CPHI(IM)=DCOS(PHI)
       SPHI(IM)=DSIN(PHI)
       DO 25 M=0,MM
        XM=DFLOAT(M)
        CS(M,IM)=DCOS(XM*PHI)
        SN(M,IM)=DSIN(XM*PHI)
        CSM(M,IM)=XM*CS(M,IM)
        SNM(M,IM)=XM*SN(M,IM)
        FAC=1.D0/DFLOAT(MN)
        IF ( M .EQ. 0 ) FAC=FAC/2.D0
        CT(IM,M)=FAC * CS(M,IM)
        ST(IM,M)=FAC * SN(M,IM)
25    CONTINUE
      WRITE(2) CPHI
      WRITE(2) SPHI
      WRITE(2) CS
      WRITE(2) SN
      WRITE(2) CSM
      WRITE(2) SNM
      WRITE(2) CT
      WRITE(2) ST
******
      DO 35 IK=1,KN
       XX=-DCOS(DFLOAT(2*IK - 1)/DFLOAT(2*KN) * PI)
       DO 30 K=1,KN
        TMAT(IK,K)=TTT(0,K-1,XX)
        TMT2(IK,K)=0.D0
30     CONTINUE
       TMT2(IK,IK)=1.D0
35    CONTINUE
      CALL F04AEF(TMAT,KN,TMT2,KN,KN,KN,TINV,
     &             KN,WK1,TMT4,KN,TMT1,KN,IFL)
      DO 40 IK=1,KN
       XX=-DCOS(DFLOAT(2*IK - 1)/DFLOAT(2*KN) * PI)
       R=(R2 + R1)/2.D0  +  (R2 - R1)/2.D0*XX
       DO 40 K=1,KN
        TMT1(IK,K)=TTT(0,K-1,XX) / R
        TMT2(IK,K)=TTT(0,K-1,XX) / R + TTT(1,K-1,XX) * DXDR
40    CONTINUE
      CALL MATMULT(TINV,TMT2,TMAT,KN)
      CALL MATMULT(TMAT,TINV,TMT2,KN)
      CALL MATMULT(TINV,TMT1,TMAT,KN)
      CALL MATMULT(TMAT,TINV,TMT1,KN)
      CALL TRANSPS(TMT1,TMT4,KN)
      WRITE(2) TMT4
      CALL MATMULT(TMAT,TMT1,TINV,KN)
      CALL TRANSPS(TINV,TMT4,KN)
      WRITE(2) TMT4
      CALL MATMULT(TMAT,TMT2,TINV,KN)
      CALL TRANSPS(TINV,TMT4,KN)
      WRITE(2) TMT4
*
      DO 85 M=0,MM
       DO 45 IL=1,LN
        WK2(IL)=( DFLOAT(2*LN+1)/DFLOAT(LN*(LN+1)) )**2
     &        / ( PPP(LN+1,0,CTH(IL)) - PPP(LN-1,0,CTH(IL)) )**2
     &        * ( 1.D0 - CTH(IL)**2 )
45     CONTINUE
       DO 50 L=1,LN
        FAC=DFLOAT(2*L + 2*M - 1)
        DO LF=0,2*M-1
         FAC=FAC / DFLOAT(L+LF)
        ENDDO
        DO 50 IL=1,LN
         PMT3(IL,L)=FAC * WK2(IL) * PPP(L-1+M,M,CTH(IL))
50     CONTINUE
       DO 55 IL=1,LN
        DO 55 L=1,LN
         PNV1(IL,L)=PMT3(IL,L)
         IF ( M .EQ. 0 ) PNV2(IL,L)=PMT3(IL,L) / STH(IL)
         IF ( M .GT. 0 ) PNV2(IL,L)=PMT3(IL,L) * STH(IL)
55     CONTINUE
       DO 60 IL=1,LN
        DO 60 L=1,LN
         LZ=L-1+M
         IF ( M .EQ. 0 )
     &   PMT1(L,IL)=(DFLOAT((LZ+2)*(LZ-M+1)) * PPP(LZ+1,M,CTH(IL))
     &              -DFLOAT((LZ+M)*(LZ  -1)) * PPP(LZ-1,M,CTH(IL)))
     &               / DFLOAT(2*LZ+1)
         IF ( M .GT. 0 )
     &   PMT1(L,IL)=(DFLOAT((LZ  )*(LZ-M+1)) * PPP(LZ+1,M,CTH(IL))
     &              -DFLOAT((LZ+M)*(LZ  +1)) * PPP(LZ-1,M,CTH(IL)))
     &               / DFLOAT(2*LZ+1) / (1.D0 - CTH(IL)**2)
         PMT2(L,IL)=PPP(LZ,M,CTH(IL)) * DFLOAT(M)/(1.D0 - CTH(IL)**2)
         PMT3(L,IL)=(DFLOAT((LZ  )*(LZ-M+1)) * PPP(LZ+1,M,CTH(IL))
     &              -DFLOAT((LZ+M)*(LZ  +1)) * PPP(LZ-1,M,CTH(IL)))
     &               / DFLOAT(2*LZ+1) / STH(IL)
60     CONTINUE
       CALL MATMULT(PNV2,PMT1,PMT4,LN)
       CALL MATMULT(PMT4,PNV1,PMT5,LN)
       DO 65 I=1,LN
        DO 65 J=1,LN
         JJ=J-1+M
         JJ=MOD(JJ-1+LN,LN) + 1
         PPP1(I,JJ,M)=PMT5(I,J)
65     CONTINUE
       CALL MATMULT(PNV2,PMT2,PMT4,LN)
       CALL MATMULT(PMT4,PNV1,PMT5,LN)
       DO 70 I=1,LN
        DO 70 J=1,LN
         JJ=J-1+M
         JJ=MOD(JJ-1+LN,LN) + 1
         PPP2(I,JJ,M)=PMT5(I,J)
70     CONTINUE
       CALL MATMULT(PNV1,PMT3,PMT4,LN)
       CALL MATMULT(PMT4,PNV2,PMT3,LN)
       CALL MATMULT(PNV1,PMT2,PMT4,LN)
       CALL MATMULT(PMT3,PMT1,PMT5,LN)
       DO 75 I=1,LN
        DO 75 J=1,LN
         PNV2(I,J)=DFLOAT(M)*PMT4(I,J) - PMT5(I,J)
75     CONTINUE
       CALL MATMULT(PNV2,PNV1,PMT4,LN)
       DO 80 I=1,LN
        DO 80 J=1,LN
         JJ=J-1+M
         JJ=MOD(JJ-1+LN,LN) + 1
         PPP4(I,JJ,M)=PMT4(I,J)
80     CONTINUE
85    CONTINUE
      DO 90 I=1,LN
       DO 90 J=1,LN
        WRITE(2) (PPP1(I,J,M) , M=0,MM)
        WRITE(2) (PPP2(I,J,M) , M=0,MM)
        WRITE(2) (PPP4(I,J,M) , M=0,MM)
90    CONTINUE
******
      DO 100 IK=1,KN
       XX=-DCOS(DFLOAT(2*IK - 1)/DFLOAT(2*KN) * PI)
       DO 95 K=1,KN
        TMAT(IK,K)=TTT(0,K-1,XX)
        TMT2(IK,K)=0.D0
95     CONTINUE
       TMT2(IK,IK)=1.D0
100   CONTINUE
      CALL F04AEF(TMAT,KN,TMT2,KN,KN,KN,TINV,
     &             KN,WK1,TMT4,KN,TMT1,KN,IFL)
      CALL TRANSPS(TINV,TMAT,KN)
      WRITE(2) TMAT
*
      DO 120 M=0,MM
       DO 105 IL=1,LN
        WK2(IL)=( DFLOAT(2*LN+1)/DFLOAT(LN*(LN+1)) )**2
     &        / ( PPP(LN+1,0,CTH(IL)) - PPP(LN-1,0,CTH(IL)) )**2
     &        * ( 1.D0 - CTH(IL)**2 )
105    CONTINUE
       DO 110 L=1,LN
        FAC=DFLOAT(2*L + 2*M - 1)
        DO LF=0,2*M-1
         FAC=FAC / DFLOAT(L+LF)
        ENDDO
        DO 110 IL=1,LN
         PMT3(IL,L)=FAC * WK2(IL) * PPP(L-1+M,M,CTH(IL))
110    CONTINUE
       DO 115 I=1,LN
        DO 115 J=1,LN
         JJ=J-1+M
         JJ=MOD(JJ+LN,LN)
         PPPP(I,JJ,M)=PMT3(I,J)
115    CONTINUE
120   CONTINUE
      DO 125 I=1,LN
       DO 125 J=0,LN
        WRITE(2) (PPPP(I,J,M) , M=0,MM)
125   CONTINUE
*
      CALPHA=DCOS(ALPHA*PI/180.D0)
      SALPHA=DSIN(ALPHA*PI/180.D0)
      WRITE(2) CALPHA
      WRITE(2) SALPHA
******
      STOP
      END
*******************************************************************************
      SUBROUTINE MATMULT(A,B,C,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(N,N) , B(N,N) , C(N,N)
      DO 10 I=1,N
       DO 10 J=1,N
        C(I,J)=0.D0
        DO 10 K=1,N
         C(I,J)=C(I,J) + A(I,K)*B(K,J)
10    CONTINUE
      RETURN
      END
*******************************************************************************
      SUBROUTINE TRANSPS(A,B,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(N,N) , B(N,N)
      DO 10 I=1,N
       DO 10 J=1,N
        B(I,J)=A(J,I)
10    CONTINUE
      RETURN
      END
*******************************************************************************
*     CMESH computes the angular collocation points, which are the LN zeros
*     of P_LN^0 on (-1,1).  The algorithm is to divide the interval (0,1)
*     into 999 sub-intervals of 1/1000 each.  Within each subinterval, if
*     there is a zero, 30 bisections will give it to an accuracy of 1.E-12.
*
      SUBROUTINE CMESH(LN,C)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION C(LN)
      K=0
      IF ( MOD(LN,2) .EQ. 1 ) THEN
       K=1
       C(K)=0.D0
      END IF
      DO 20 I=1,9999
       C1=(DFLOAT(I) - 0.1D0)/10000.D0
       C2=(DFLOAT(I) + 0.9D0)/10000.D0
       F1=PPP(LN,0,C1)
       F2=PPP(LN,0,C2)
       IF ( F1*F2 .GT. 0.D0 ) GOTO 20
       DO 10 J=1,30
        CM=0.5D0*(C1 + C2)
        FM=PPP(LN,0,CM)
        IF ( F1*FM .GT. 0.D0 ) THEN
         C1=CM
         F1=FM
        ELSE
         C2=CM
         F2=FM
        END IF
10     CONTINUE
       K=K+1
       C(K)=0.5D0*(C1 + C2)
       K=K+1
       C(K)=-C(K-1)
20    CONTINUE
      RETURN
      END
*******************************************************************************
*     PPP(N,M,C) is the N,M-th associated Legendre function evaluated at
*     C, using the recursion relation below from Abramowitz & Stegun.
*
      DOUBLE PRECISION FUNCTION PPP(N,M,C)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION P(0:10000)
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
*     TTT(K,M,X) = the K-th derivative of Tm(X),
*     the M-th Chebyshev polynomial evaluated at X.
*
      DOUBLE PRECISION FUNCTION TTT(K,M,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(0:10000), B(0:10000)
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
      DOUBLE PRECISION A(0:10000), C(0:10000)
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
