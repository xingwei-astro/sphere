      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
*
      COMMON /FIELD/ GC(KB2, LB, 0:MB) , GS(KB2, LB, 0:MB)
     &             , HC(KB2, LB, 0:MB) , HS(KB2, LB, 0:MB)
*
      DIMENSION TORB(LB,0:MB) , POLB(LB,0:MB)
      DIMENSION EMAGM(0:MB), EMAGL(LB) 
*
      READ(4) GC
      READ(4) GS
      READ(4) HC
      READ(4) HS
      CLOSE(4)
*
      do m=0, mb
       do l=1, lb
        torb(l,m)=0.d0
        polb(l,m)=0.d0
       enddo
      enddo      
******
      CALL INIT
      CALL ENERGY(TORB,POLB)
******
      open(1,file='emag.dat',form='formatted')
      do m=0,mb
       do l=1,lb
        enr=torb(l,m)+polb(l,m)
        write(1,*) l, m, enr
       enddo
      enddo
      close(1)      
******
      OPEN(1,FILE='emag-m.dat',FORM='formatted')
      enr0=0.d0
      DO M=0,MB
       EMAGM(M)=0.D0
       LS=MAX(M,1)
       DO L=LS,LB
        EMAGM(M)=EMAGM(M)+TORB(L,M)+POLB(L,M)
        enr=torb(l,m)+polb(l,m)
        if(enr.ge.enr0) then
         enr0=enr
         l0=l
         m0=m
        endif
       ENDDO
       WRITE(1,*) M, EMAGM(M)
      ENDDO
      CLOSE(1)
      write(6,*) l0, m0

      OPEN(2,FILE='emag-l.dat',FORM='formatted')
      DO L=1,LB
       EMAGL(L)=0.D0
       DO M=0,MB
        EMAGL(L)=EMAGL(L)+TORB(L,M)+POLB(L,M)
       ENDDO
       WRITE(2,*) L, EMAGL(L)
      ENDDO
      CLOSE(2)

      DO L=1,LB
       ENR=0.D0
       ENR=ENR+TORB(L,0)
      ENDDO
      WRITE(6,*) 'ZONAL MAGNETIC ENERGY=', ENR
******
      STOP
      END
************************************************************************
      SUBROUTINE ENERGY(TORB,POLB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
*
      COMMON /FIELD/ GC(KB2, LB, 0:MB) , GS(KB2, LB, 0:MB)
     &             , HC(KB2, LB, 0:MB) , HS(KB2, LB, 0:MB)
*
      COMMON /WORK/ XINT1(LB, 0:MB) , XINT2(LB, 0:MB) , WK(0:NKB)
     &           , TT1(KB2,0:NKB) , TT2(KB2,0:NKB) , TT3(KB2,0:NKB)
     &           , TT4(KB2,0:NKB) , TT5(KB2,0:NKB)
*
      DIMENSION TORB(LB,0:MB) , POLB(LB,0:MB)
******
      DO M=0,MB
       LS=MAX(M,1)
       DO L=LS,LB
        XN1=0.D0
        XN2=0.D0
        XN3=0.D0
        DO IK=0,NKB
         WHC1=0.D0
         WHS1=0.D0
         WHC2=0.D0
         WHS2=0.D0
         DO K=1,KB2
          WHC1=WHC1  +  HC(K,L,M) * TT1(K,IK)
          WHS1=WHS1  +  HS(K,L,M) * TT1(K,IK)
          WHC2=WHC2  +  HC(K,L,M) * TT2(K,IK)
          WHS2=WHS2  +  HS(K,L,M) * TT2(K,IK)
         ENDDO
         WGC3=0.D0
         WGS3=0.D0
         DO K=1,KB2
          WGC3=WGC3  +  GC(K,L,M) * TT3(K,IK)
          WGS3=WGS3  +  GS(K,L,M) * TT3(K,IK)
         ENDDO
         XN1=XN1  +  (WHC1**2 + WHS1**2) * WK(IK)
         XN2=XN2  +  (WHC2**2 + WHS2**2) * WK(IK)
         XN3=XN3  +  (WGC3**2 + WGS3**2) * WK(IK)
        ENDDO
        POLB(L,M)=XN1*XINT1(L,M) + XN2*XINT2(L,M)
        TORB(L,M)=XN3*XINT2(L,M)
       ENDDO
      ENDDO
*
      RETURN
      END
************************************************************************
      SUBROUTINE INIT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
*
      COMMON /WORK/ XINT1(LB, 0:MB) , XINT2(LB, 0:MB) , WK(0:NKB)
     &           , TT1(KB2,0:NKB) , TT2(KB2,0:NKB) , TT3(KB2,0:NKB)
     &           , TT4(KB2,0:NKB) , TT5(KB2,0:NKB)
******
      PI=DACOS(-1.D0)
      DO M=0,MB
       LS=MAX(M,1)
       DO L=LS,LB
        FAC=2.D0/DFLOAT(2*L+1)
        DO N=L-M+1,L+M
         FAC=FAC * DFLOAT(N)
        ENDDO
        XINT1(L,M)=FAC * DFLOAT(L*(L+1))**2
        XINT2(L,M)=FAC * DFLOAT(L*(L+1))
        XINT1(L,M)=XINT1(L,M) * PI/2.D0
        XINT2(L,M)=XINT2(L,M) * PI/2.D0
        IF ( M .EQ. 0 ) XINT1(L,M)=XINT1(L,M) * 2.D0
        IF ( M .EQ. 0 ) XINT2(L,M)=XINT2(L,M) * 2.D0
       ENDDO
      ENDDO
******
      DR=(R2 - R1)/DFLOAT(NKB)
      DO IK=0,NKB
       IF ( MOD(IK,2)  .EQ. 1 ) WK(IK)=4.D0/3.D0 * DR
       IF ( MOD(IK,2)  .EQ. 0 ) WK(IK)=2.D0/3.D0 * DR
       IF ( IK*(IK-NKB) .EQ. 0 ) WK(IK)=1.D0/3.D0 * DR
       XX=-1.D0 + DFLOAT(2*IK)/DFLOAT(NKB)
       R=(R2 + R1)/2.D0  +  (R2 - R1)/2.D0*XX
       DO K=1,KB2
        TT1(K,IK)=TTT(0,K-1,XX) / R
        TT2(K,IK)=TTT(1,K-1,XX) * DXDR
        TT3(K,IK)=TTT(0,K-1,XX)
       ENDDO
      ENDDO
******
      RETURN
      END
************************************************************************
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
