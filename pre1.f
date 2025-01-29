      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*      
      INCLUDE 'params.f'
*
      DIMENSION XB(KB) , XU(KU) , XT(KT)
*
      DIMENSION GMAT(KB2,KB2) , HMAT(KB2,KB2)
      DIMENSION BMT0(KB2,KB2) , BMT3(KB2,KB2)
      DIMENSION BMT4(KB2,KB2) , BMT5(KB2,KB2)
      DIMENSION BBBB(KB2,KB2) , GHMT(KB2,KB2)
*
      DIMENSION GMT1(KB2,KB2,LB) , GMT2(KB2,KB2,LB)
      DIMENSION HMT1(KB2,KB2,LB) , HMT2(KB2,KB2,LB)
*
      DIMENSION EMAT(KU2,KU2)
      DIMENSION EMT0(KU2,KU2) , EMT3(KU2,KU2)
      DIMENSION EMT4(KU2,KU2) , EMT5(KU2,KU2)
      DIMENSION EEEE(KU2,KU2) , EEMT(KU2,KU2)
*
      DIMENSION FMAT(KU4,KU4)
      DIMENSION FMT0(KU4,KU4) , FMT3(KU4,KU4)
      DIMENSION FMT4(KU4,KU4) , FMT5(KU4,KU4)
      DIMENSION FFFF(KU4,KU4) , FFMT(KU4,KU4)
*
      DIMENSION EMT1(KU2,KU2,LU) , EMT2(KU2,KU2,LU)
      DIMENSION FMT1(KU4,KU4,LU) , FMT2(KU4,KU4,LU)
*
      DIMENSION TMAT(KT2,KT2)
      DIMENSION TMT0(KT2,KT2) , TMT3(KT2,KT2)
      DIMENSION TMT4(KT2,KT2) , TMT5(KT2,KT2)
      DIMENSION TTTT(KT2,KT2) , TTMT(KT2,KT2)
*
      DIMENSION TMT1(KT2,KT2,0:LT)
      DIMENSION TMT2(KT2,KT2,0:LT)
*
      DIMENSION WK(KB+KU+KT)
******
*     Compute the time-stepping matrices for g and h.
*     Compute the radial collocation points, the KB zeros of T(KB) on (-1,1).
      PI=DACOS(-1.D0)
      DO 10 K=1,KB
       XB(K)=-DCOS(DFLOAT(2*K - 1)/DFLOAT(2*KB) * PI)
10    CONTINUE
      DO 20 I=1,KB2
       DO 15 J=1,KB2
        BBBB(I,J)=0.D0
15     CONTINUE
       BBBB(I,I)=1.D0
20    CONTINUE
      DO 25 I=1,KB
       DO 25 J=1,KB
        BBBB(I,J)=DT * TTT(0,J-1,XB(I))
25    CONTINUE
*
      DO 60 L=1,LB
       DO 30 I=1,KB
        R=(R2 + R1)/2.D0  +  (R2 - R1)/2.D0*XB(I)
        FAC=DFLOAT(L*(L+1))/R**2
        DO 30 J=1,KB2
         T0=TTT(0,J-1,XB(I))
         T2=TTT(2,J-1,XB(I)) * DXDR**2
         GMAT(I,J)=FAC * ( T0 - EKM/PM*0.5D0*DT*(T2 - FAC*T0) )
         HMAT(I,J)=FAC * ( T0 - EKM/PM*0.5D0*DT*(T2 - FAC*T0) )
         GHMT(I,J)=FAC * ( T0 + EKM/PM*0.5D0*DT*(T2 - FAC*T0) )
30     CONTINUE
*      Boundary Conditions on g and h
       DO 35 J=1,KB2
        GMAT(KB+1,J)=(-1.D0)**(J-1)
        GMAT(KB+2,J)=  1.D0
        HMAT(KB+1,J)=(-1.D0)**J*(DXDR*DFLOAT(J-1)**2 + DFLOAT(L+1)/R1)
        HMAT(KB+2,J)=  1.D0    *(DXDR*DFLOAT(J-1)**2 + DFLOAT( L )/R2)
        GHMT(KB+1,J)=  0.D0
        GHMT(KB+2,J)=  0.D0
35     CONTINUE
       DO 45 I=1,KB2
        DO 40 J=1,KB2
         BMT0(I,J)=0.D0
40      CONTINUE
        BMT0(I,I)=1.D0
45     CONTINUE
       CALL F04AEF(GMAT,KB2,BMT0,KB2,KB2,KB2,BMT3,
     &               KB2,WK,BMT4,KB2,BMT5,KB2,IFL)
       CALL MATMULT(BMT3,GHMT,BMT4,KB2)
       CALL MATMULT(BMT3,BBBB,BMT5,KB2)
       DO 50 I=1,KB2
       DO 50 J=1,KB2
        GMT1(I,J,L)=BMT4(J,I)
        GMT2(I,J,L)=BMT5(J,I)
50     CONTINUE
       CALL F04AEF(HMAT,KB2,BMT0,KB2,KB2,KB2,BMT3,
     &               KB2,WK,BMT4,KB2,BMT5,KB2,IFL)
       CALL MATMULT(BMT3,GHMT,BMT4,KB2)
       CALL MATMULT(BMT3,BBBB,BMT5,KB2)
       DO 55 I=1,KB2
       DO 55 J=1,KB2
        HMT1(I,J,L)=BMT4(J,I)
        HMT2(I,J,L)=BMT5(J,I)
55     CONTINUE
60    CONTINUE
      DO 65 I=1,KB2
      DO 65 J=1,KB2
       WRITE(1) (GMT1(I,J,L) , L=1,LB)
       WRITE(1) (GMT2(I,J,L) , L=1,LB)
       WRITE(1) (HMT1(I,J,L) , L=1,LB)
       WRITE(1) (HMT2(I,J,L) , L=1,LB)
65    CONTINUE
******
*     Compute the time-stepping matrices for e.
*     Compute the radial collocation points, the KU zeros of T(KU) on (-1,1).
      PI=DACOS(-1.D0)
      DO 70 K=1,KU
       XU(K)=-DCOS(DFLOAT(2*K - 1)/DFLOAT(2*KU) * PI)
70    CONTINUE
      DO 80 I=1,KU2
       DO 75 J=1,KU2
        EEEE(I,J)=0.D0
75     CONTINUE
       EEEE(I,I)=1.D0
80    CONTINUE
      DO 85 I=1,KU
       DO 85 J=1,KU
        EEEE(I,J)=DT * TTT(0,J-1,XU(I))
85    CONTINUE
*
      DO 115 L=1,LU
       DO 90 I=1,KU
        R=(R2 + R1)/2.D0  +  (R2 - R1)/2.D0*XU(I)
        FAC=DFLOAT(L*(L+1))/R**2
        DO 90 J=1,KU2
         T0=TTT(0,J-1,XU(I))
         T2=TTT(2,J-1,XU(I)) * DXDR**2
         EMAT(I,J)=FAC * ( T0 - EKM*0.8D0*DT*(T2 - FAC*T0) )
         EEMT(I,J)=FAC * ( T0 + EKM*0.2D0*DT*(T2 - FAC*T0) )
90     CONTINUE
*      Boundary Conditions on e
       DO 95 J=1,KU2
*       no-slip
        EMAT(KU+1,J)=(-1.D0)**(J-1)
        EMAT(KU+2,J)=  1.D0
*       stress-free
*        EMAT(KU+1,J)=TTT(1,J-1,-1.D0)*DXDR/R1**2
*     &              -TTT(0,J-1,-1.D0)*2.D0/R1**3
*        EMAT(KU+2,J)=TTT(1,J-1, 1.D0)*DXDR/R2**2
*     &              -TTT(0,J-1, 1.D0)*2.D0/R2**3
        EEMT(KU+1,J)=  0.D0
        EEMT(KU+2,J)=  0.D0
95     CONTINUE
       DO 105 I=1,KU2
        DO 100 J=1,KU2
         EMT0(I,J)=0.D0
100     CONTINUE
        EMT0(I,I)=1.D0
105    CONTINUE
       CALL F04AEF(EMAT,KU2,EMT0,KU2,KU2,KU2,EMT3,
     &               KU2,WK,EMT4,KU2,EMT5,KU2,IFL)
       CALL MATMULT(EMT3,EEMT,EMT4,KU2)
       CALL MATMULT(EMT3,EEEE,EMT5,KU2)
       DO 110 I=1,KU2
       DO 110 J=1,KU2
        EMT1(I,J,L)=EMT4(J,I)
        EMT2(I,J,L)=EMT5(J,I)
110    CONTINUE
115   CONTINUE
      DO 120 I=1,KU2
      DO 120 J=1,KU2
       WRITE(1) (EMT1(I,J,L) , L=1,LU)
       WRITE(1) (EMT2(I,J,L) , L=1,LU)
120   CONTINUE
******
*     Compute the time-stepping matrices for f.
*     Compute the radial collocation points, the KU zeros of T(KU) on (-1,1).
      PI=DACOS(-1.D0)
      DO 125 K=1,KU
       XU(K)=-DCOS(DFLOAT(2*K - 1)/DFLOAT(2*KU) * PI)
125   CONTINUE
      DO 135 I=1,KU4
       DO 130 J=1,KU4
        FFFF(I,J)=0.D0
130    CONTINUE
       FFFF(I,I)=1.D0
135   CONTINUE
      DO 140 I=1,KU
       DO 140 J=1,KU
        FFFF(I,J)=DT * TTT(0,J-1,XU(I))
140   CONTINUE
*
      DO 170 L=1,LU
       DO 145 I=1,KU
        R=(R2 + R1)/2.D0  +  (R2 - R1)/2.D0*XU(I)
        FAC=DFLOAT(L*(L+1))/R**2
        DO 145 J=1,KU4
         T0=TTT(0,J-1,XU(I))
         T1=TTT(1,J-1,XU(I)) * DXDR
         T2=TTT(2,J-1,XU(I)) * DXDR**2
         T4=TTT(4,J-1,XU(I)) * DXDR**4
         FMAT(I,J)=-FAC * ( (T2 - FAC*T0) - EKM*0.8D0*DT
     &                    *(T4 - 2.D0*FAC*T2 + 4.D0*FAC/R*T1
     &                           + FAC**2*T0 - 6.D0*FAC/R**2*T0) )
         FFMT(I,J)=-FAC * ( (T2 - FAC*T0) + EKM*0.2D0*DT
     &                    *(T4 - 2.D0*FAC*T2 + 4.D0*FAC/R*T1
     &                           + FAC**2*T0 - 6.D0*FAC/R**2*T0) )
145    CONTINUE
*      Boundary Conditions on f
       DO 150 J=1,KU4
*       no-slip
        FMAT(KU+1,J)=(-1.D0)**(J-1)
        FMAT(KU+2,J)=  1.D0
        FMAT(KU+3,J)=TTT(1,J-1,-1.D0)*DXDR
        FMAT(KU+4,J)=TTT(1,J-1, 1.D0)*DXDR
*       stress-free
*        FMAT(KU+1,J)=(-1.D0)**(J-1)
*        FMAT(KU+2,J)=  1.D0
*        FMAT(KU+3,J)=TTT(2,J-1,-1.D0)*DXDR**2/R1**2
*     &              -TTT(1,J-1,-1.D0)*DXDR*2.D0/R1**3
*        FMAT(KU+4,J)=TTT(2,J-1, 1.D0)*DXDR**2/R2**2
*     &              -TTT(1,J-1, 1.D0)*DXDR*2.D0/R2**3
        FFMT(KU+1,J)=  0.D0
        FFMT(KU+2,J)=  0.D0
        FFMT(KU+3,J)=  0.D0
        FFMT(KU+4,J)=  0.D0
150    CONTINUE
       DO 160 I=1,KU4
        DO 155 J=1,KU4
         FMT0(I,J)=0.D0
155     CONTINUE
        FMT0(I,I)=1.D0
160    CONTINUE
       CALL F04AEF(FMAT,KU4,FMT0,KU4,KU4,KU4,FMT3,
     &               KU4,WK,FMT4,KU4,FMT5,KU4,IFL)
       CALL MATMULT(FMT3,FFMT,FMT4,KU4)
       CALL MATMULT(FMT3,FFFF,FMT5,KU4)
       DO 165 I=1,KU4
       DO 165 J=1,KU4
        FMT1(I,J,L)=FMT4(J,I)
        FMT2(I,J,L)=FMT5(J,I)
165    CONTINUE
170   CONTINUE
      DO 175 I=1,KU4
      DO 175 J=1,KU4
       WRITE(1) (FMT1(I,J,L) , L=1,LU)
       WRITE(1) (FMT2(I,J,L) , L=1,LU)
175   CONTINUE
******
*     Compute the time-stepping matrices for theta.
*     Compute the radial collocation points, the KT zeros of T(KT) on (-1,1).
      PI=DACOS(-1.D0)
      DO 180 K=1,KT
       XT(K)=-DCOS(DFLOAT(2*K - 1)/DFLOAT(2*KT) * PI)
180   CONTINUE
      DO 190 I=1,KT2
       DO 185 J=1,KT2
        TTTT(I,J)=0.D0
185    CONTINUE
       TTTT(I,I)=1.D0
190   CONTINUE
      DO 195 I=1,KT
       DO 195 J=1,KT
        TTTT(I,J)=DT * TTT(0,J-1,XT(I))
195   CONTINUE
*
      DO 225 L=0,LT
       DO 200 I=1,KT
        R=(R2 + R1)/2.D0  +  (R2 - R1)/2.D0*XT(I)
        FAC=DFLOAT(L*(L+1))/R**2
        DO 200 J=1,KT2
         T0=TTT(0,J-1,XT(I))
         T1=TTT(1,J-1,XT(I)) * DXDR
         T2=TTT(2,J-1,XT(I)) * DXDR**2
         TMAT(I,J)=T0 - EKM/PR*0.5D0*DT*(T2 + 2.D0/R*T1 - FAC*T0)
         TTMT(I,J)=T0 + EKM/PR*0.5D0*DT*(T2 + 2.D0/R*T1 - FAC*T0)
200    CONTINUE
*      Boundary Conditions on theta
       DO 205 J=1,KT2
        TMAT(KT+1,J)=(-1.D0)**(J-1)
        TMAT(KT+2,J)=  1.D0
        TTMT(KT+1,J)=  0.D0
        TTMT(KT+2,J)=  0.D0
205    CONTINUE
       DO 215 I=1,KT2
        DO 210 J=1,KT2
         TMT0(I,J)=0.D0
210     CONTINUE
        TMT0(I,I)=1.D0
215    CONTINUE
       CALL F04AEF(TMAT,KT2,TMT0,KT2,KT2,KT2,TMT3,
     &               KT2,WK,TMT4,KT2,TMT5,KT2,IFL)
       CALL MATMULT(TMT3,TTMT,TMT4,KT2)
       CALL MATMULT(TMT3,TTTT,TMT5,KT2)
       DO 220 I=1,KT2
       DO 220 J=1,KT2
        TMT1(I,J,L)=TMT4(J,I)
        TMT2(I,J,L)=TMT5(J,I)
220    CONTINUE
225   CONTINUE
      DO 230 I=1,KT2
      DO 230 J=1,KT2
       WRITE(1) (TMT1(I,J,L) , L=0,LT)
       WRITE(1) (TMT2(I,J,L) , L=0,LT)
230   CONTINUE
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
