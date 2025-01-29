      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
*
      COMMON /FLOW/ EC(KU2, LU, 0:MU) , ES(KU2, LU, 0:MU)
     &            , FC(KU4, LU, 0:MU) , FS(KU4, LU, 0:MU)
*
      DIMENSION TORU(LU,0:MU) , POLU(LU,0:MU)
      DIMENSION EKINM(0:MU), EKINL(LU) 
*
      READ(3) EC
      READ(3) ES
      READ(3) FC
      READ(3) FS
      CLOSE(3)
*
      do m=0, mu
       do l=1, lu
        toru(l,m)=0.d0
        polu(l,m)=0.d0
       enddo
      enddo      
******
      CALL INIT
      CALL ENERGY(TORU,POLU)      
******
      open(1,file='ekin.dat',form='formatted')
      do m=0,mu
       do l=1,lu
        enr=toru(l,m)+polu(l,m)
        write(1,*) l, m, enr
       enddo
      enddo
      close(1)
******
      OPEN(1,FILE='ekin-m.dat',FORM='formatted')
      enr0=0.d0
      DO M=0,MU
       EKINM(M)=0.D0
       LS=MAX(M,1)
       DO L=LS,LU
        EKINM(M)=EKINM(M)+TORU(L,M)+POLU(L,M)
        enr=toru(l,m)+polu(l,m)
        if(enr.ge.enr0) then
         enr0=enr
         l0=l
         m0=m
        endif
       ENDDO
       WRITE(1,*) M, EKINM(M)
      ENDDO
      CLOSE(1)
      write(6,*) l0, m0

      OPEN(2,FILE='ekin-l.dat',FORM='formatted')
      DO L=1,LU
       EKINL(L)=0.D0
       DO M=0,MU
        EKINL(L)=EKINL(L)+TORU(L,M)+POLU(L,M)
       ENDDO
       WRITE(2,*) L, EKINL(L)
      ENDDO
      CLOSE(2)

      DO L=1,LU
       ENR=0.D0
       ENR=ENR+TORU(L,0)
      ENDDO
      WRITE(6,*) 'ZONAL KINETIC ENERGY=', ENR
******
      STOP
      END
************************************************************************
      SUBROUTINE ENERGY(TORU,POLU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
*
      COMMON /FLOW/ EC(KU2, LU, 0:MU) , ES(KU2, LU, 0:MU)
     &            , FC(KU4, LU, 0:MU) , FS(KU4, LU, 0:MU)
*
      COMMON /WORK/ XINT1(LU, 0:MU) , XINT2(LU, 0:MU) , WK(0:NKU)
     &           , TT1(KU4,0:NKU) , TT2(KU4,0:NKU) , TT3(KU4,0:NKU)
     &           , TT4(KU4,0:NKU) , TT5(KU4,0:NKU)
*
      DIMENSION TORU(LU,0:MU) , POLU(LU,0:MU)
******
      DO M=0,MU
       LS=MAX(M,1)
       DO L=LS,LU
        XN1=0.D0
        XN2=0.D0
        XN3=0.D0
        DO IK=0,NKU
         WFC1=0.D0
         WFS1=0.D0
         WFC2=0.D0
         WFS2=0.D0
         DO K=1,KU4
          WFC1=WFC1  +  FC(K,L,M) * TT1(K,IK)
          WFS1=WFS1  +  FS(K,L,M) * TT1(K,IK)
          WFC2=WFC2  +  FC(K,L,M) * TT2(K,IK)
          WFS2=WFS2  +  FS(K,L,M) * TT2(K,IK)
         ENDDO
         WEC3=0.D0
         WES3=0.D0
         DO K=1,KU2
          WEC3=WEC3  +  EC(K,L,M) * TT3(K,IK)
          WES3=WES3  +  ES(K,L,M) * TT3(K,IK)
         ENDDO
         XN1=XN1  +  (WFC1**2 + WFS1**2) * WK(IK)
         XN2=XN2  +  (WFC2**2 + WFS2**2) * WK(IK)
         XN3=XN3  +  (WEC3**2 + WES3**2) * WK(IK)
        ENDDO
        POLU(L,M)=XN1*XINT1(L,M) + XN2*XINT2(L,M)
        TORU(L,M)=XN3*XINT2(L,M)
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
      COMMON /WORK/ XINT1(LU, 0:MU) , XINT2(LU, 0:MU) , WK(0:NKU)
     &           , TT1(KU4,0:NKU) , TT2(KU4,0:NKU) , TT3(KU4,0:NKU)
     &           , TT4(KU4,0:NKU) , TT5(KU4,0:NKU)
******
      PI=DACOS(-1.D0)
      DO M=0,MU
       LS=MAX(M,1)
       DO L=LS,LU
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
      DR=(R2 - R1)/DFLOAT(NKU)
      DO IK=0,NKU
       IF ( MOD(IK,2)  .EQ. 1 ) WK(IK)=4.D0/3.D0 * DR
       IF ( MOD(IK,2)  .EQ. 0 ) WK(IK)=2.D0/3.D0 * DR
       IF ( IK*(IK-NKU) .EQ. 0 ) WK(IK)=1.D0/3.D0 * DR
       XX=-1.D0 + DFLOAT(2*IK)/DFLOAT(NKU)
       R=(R2 + R1)/2.D0  +  (R2 - R1)/2.D0*XX
       DO K=1,KU4
        TT1(K,IK)=TTT(0,K-1,XX) / R
        TT2(K,IK)=TTT(1,K-1,XX) * DXDR
        TT3(K,IK)=TTT(0,K-1,XX)
        TT4(K,IK)=TTT(0,K-1,XX) / R**2
        TT5(K,IK)=TTT(2,K-1,XX) * DXDR**2
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
