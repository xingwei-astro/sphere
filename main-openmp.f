      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
      include 'omp_lib.h'
*
      COMMON /FIELD/  GC(KB2, LB, 0:MB) ,  GS(KB2, LB, 0:MB)
     &             , DGC(KB2, LB, 0:MB) , DGS(KB2, LB, 0:MB)
     &             ,  HC(KB2, LB, 0:MB) ,  HS(KB2, LB, 0:MB)
     &             , DHC(KB2, LB, 0:MB) , DHS(KB2, LB, 0:MB)
*
      COMMON /FLOW/   EC(KU2, LU, 0:MU) ,  ES(KU2, LU, 0:MU)
     &             , DEC(KU2, LU, 0:MU) , DES(KU2, LU, 0:MU)
     &             ,  FC(KU4, LU, 0:MU) ,  FS(KU4, LU, 0:MU)
     &             , DFC(KU4, LU, 0:MU) , DFS(KU4, LU, 0:MU)
*
      COMMON /TEMP/   TC(KT2,0:LT,0:MT) ,  TS(KT2,0:LT,0:MT)
     &             , DTC(KT2,0:LT,0:MT) , DTS(KT2,0:LT,0:MT)
*
      DIMENSION TORB(LB,0:MB), POLB(LB,0:MB)
     &        , TORU(LU,0:MU), POLU(LU,0:MU)
      DIMENSION EM(0:MB), EMPHI(0:MB), EK(0:MU), EKPHI(0:MU)
      DIMENSION TORCB(LB,0:MB), POLCB(LB,0:MB)
     &        , TORCU(LU,0:MU), POLCU(LU,0:MU)
      DIMENSION DM(0:MB), DMPHI(0:MB), DK(0:MU), DKPHI(0:MU)
******
      nthreads=64
      call omp_set_num_threads(nthreads)
      s_time=omp_get_wtime()

      WRITE(6,*) 'EKM=',EKM,'AA=',AA,'CC=',CC,'nthreads=',nthreads
      CALL READIN
      CALL INIT
      DO 10 M=0,MB
       LS=MAX(M,1)
       DO 10 L=LS,LB
        DO 10 K=1,KB2
         GC(K,L,M)=0.D0
         GS(K,L,M)=0.D0
         HC(K,L,M)=0.D0
         HS(K,L,M)=0.D0
*
         DGC(K,L,M)=0.D0
         DGS(K,L,M)=0.D0
         DHC(K,L,M)=0.D0
         DHS(K,L,M)=0.D0
10    CONTINUE
      DO 20 M=0,MU
       LS=MAX(M,1)
       DO 20 L=LS,LU
        DO 20 K=1,KU2
         EC(K,L,M)=0.D0
         ES(K,L,M)=0.D0
*
         DEC(K,L,M)=0.D0
         DES(K,L,M)=0.D0
20    CONTINUE
      DO 30 M=0,MU
       LS=MAX(M,1)
       DO 30 L=LS,LU
        DO 30 K=1,KU4
         FC(K,L,M)=0.D0
         FS(K,L,M)=0.D0
*
         DFC(K,L,M)=0.D0
         DFS(K,L,M)=0.D0
30    CONTINUE
      DO 40 M=0,MT
       DO 40 L=M,LT
        DO 40 K=1,KT2
         TC(K,L,M)=0.D0
         TS(K,L,M)=0.D0
*
         DTC(K,L,M)=0.D0
         DTS(K,L,M)=0.D0
40    CONTINUE
******
      do m=0,mu
       ek(m)=0.d0
       ekphi(m)=0.d0
       dk(m)=0.d0
       dkphi(m)=0.d0
       do l=1,lu
        toru(l,m)=0.d0
        polu(l,m)=0.d0
        torcu(l,m)=0.d0
        polcu(l,m)=0.d0
       enddo
      enddo
      do m=0,mb
       em(m)=0.d0
       emphi(m)=0.d0
       dm(m)=0.d0
       dmphi(m)=0.d0
       do l=1,lb
        torb(l,m)=0.d0
        polb(l,m)=0.d0
        torcb(l,m)=0.d0
        polcb(l,m)=0.d0
       enddo
      enddo
******
*     Initialize appropriate slots of Dxx arrays to
*     implement inhomogeneous boundary conditions.
******
      IC=1
      IF(IC.EQ.0)  THEN
       TIME0=0.D0
       EC(1,1,0)=1.D-6
       ES(1,1,1)=1.D-6
       FC(1,1,0)=1.D-6
       FS(1,1,1)=1.D-6
       TC(1,1,0)=1.D-6
       TS(1,1,1)=1.D-6
       GC(1,1,0)=1.D-6
       GS(1,1,1)=1.D-6
       HC(1,1,0)=1.D-6
       HS(1,1,1)=1.D-6
      ENDIF
      IF(IC.NE.0) THEN
       OPEN(1,FILE='time',FORM='formatted')
       READ(1,*) TIME0
       CLOSE(1)
       OPEN(3,FORM='unformatted')
       READ(3) EC
       READ(3) ES
       READ(3) FC
       READ(3) FS
       CLOSE(3)
       OPEN(4,FORM='unformatted')
       READ(4) GC
       READ(4) GS
       READ(4) HC
       READ(4) HS
       CLOSE(4)
       OPEN(7,FORM='unformatted')
       READ(7) TC
       READ(7) TS
       CLOSE(7)
      ENDIF
*      
      OPEN(10,FILE='energy',FORM='FORMATTED')
      OPEN(11,FILE='gauss',FORM='FORMATTED')
      DO 150 NN=0,NT
       TIME=TIME0+DT*DFLOAT(NN)
*** collision EC for u_phi, ES for u_theta
       if(time.GE.time0) then
        EC(KU2,1,0)=CC*R2*exp(-(time-time0)/tau)
*        ES(KU2,1,1)=-CC*R2*exp(-(time-time0)/tau)
       endif
*** tidally induced bounary radial flow
*       FC(KU2,2,2)=xi*omega*R2**2/3.0*cos(2.0*omega*time)       
*       FS(KU2,2,2)=xi*omega*R2**2/3.0*sin(2.0*omega*time)
***       
       CALL STEP(TIME)
       IF ( DABS(EC(1,1,0)) .GT. 1.D50 ) STOP
       IF ( DABS(GC(1,1,0)) .GT. 1.D50 ) STOP
       CALL ENER(TORB,POLB,TORU,POLU,EM,EMPHI,EK,EKPHI)
       CALL DISS(TORCB,POLCB,TORCU,POLCU,DM,DMPHI,DK,DKPHI)
       ENERB =0.D0
       ENERBT=0.D0
       DISSB =0.D0
       DISSBT=0.D0
       ENERU =0.D0
       ENERUT=0.D0
       DISSU =0.D0
       DISSUT=0.D0
       DO M=0,MB
        ENERB =ENERB  + EM(M)
        ENERBT=ENERBT + EMPHI(M)
        DISSB =DISSB  + DM(M)
        DISSBT=DISSBT + DMPHI(M)
       ENDDO
       DO M=0,MU
        ENERU =ENERU  + EK(M)
        ENERUT=ENERUT + EKPHI(M)
        DISSU =DISSU  + DK(M)
        DISSUT=DISSUT + DKPHI(M)
       ENDDO
* magnetic moments, dipole moment = Psi^c_10(ro) = -ro*h^c_10(ro) = -HH1
       HH1=0.D0
       HH2=0.D0
       HH3=0.D0
       HH4=0.D0
       HH5=0.D0
       HH6=0.D0
       HH7=0.D0
       HH8=0.D0
       DO KK=1,KB2
        HH1=HH1+HC(KK,1,0)
        HH2=HH2+HC(KK,1,1)
        HH3=HH3+HS(KK,1,1)
        HH4=HH4+HC(KK,2,0)
        HH5=HH5+HC(KK,2,1)
        HH6=HH6+HS(KK,2,1)
        HH7=HH7+HC(KK,2,2)
        HH8=HH8+HS(KK,2,2)
       ENDDO
       WRITE(10,'(9E15.6)') TIME, ENERU, ENERUT, ENERB, ENERBT,
     &                            DISSU, DISSUT, DISSB, DISSBT
       WRITE(11,'(9E15.6)') TIME, HH1, HH2, HH3, HH4, HH5, HH6,
     &                            HH7, HH8
150   CONTINUE
      CLOSE(10)
      CLOSE(11)
*
      OPEN(1,FILE='time',FORM='formatted')
      WRITE(1,*) TIME
      CLOSE(1)
      OPEN(3,FORM='unformatted')
      WRITE(3) EC
      WRITE(3) ES
      WRITE(3) FC
      WRITE(3) FS
      CLOSE(3)
      OPEN(4,FORM='unformatted')
      WRITE(4) GC
      WRITE(4) GS
      WRITE(4) HC
      WRITE(4) HS
      CLOSE(4)
      OPEN(7,FORM='unformatted')
      WRITE(7) TC
      WRITE(7) TS
      CLOSE(7)
*
      e_time=omp_get_wtime()
*      write(6,*) 'cpu time=', e_time-s_time
******
      STOP
      END
*******************************************************************************
      SUBROUTINE READIN
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
*
      COMMON /STEP1/ GMT1(KB2, KB2, LB), HMT1(KB2, KB2, LB)
     &             , GMT2(KB2, KB2, LB), HMT2(KB2, KB2, LB)
*
      COMMON /STEP2/ EMT1(KU2, KU2, LU), FMT1(KU4, KU4, LU)
     &             , EMT2(KU2, KU2, LU), FMT2(KU4, KU4, LU)
*
      COMMON /STEP3/ TMT1(KT2, KT2, 0:LT)
     &             , TMT2(KT2, KT2, 0:LT)
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
******
      DO 10 K2=1,KB2
       DO 10 K1=1,KB2
        READ(1) (GMT1(K2,K1,L) , L=1,LB)
        READ(1) (GMT2(K2,K1,L) , L=1,LB)
        READ(1) (HMT1(K2,K1,L) , L=1,LB)
        READ(1) (HMT2(K2,K1,L) , L=1,LB)
10    CONTINUE
      DO 20 K2=1,KU2
       DO 20 K1=1,KU2
        READ(1) (EMT1(K2,K1,L) , L=1,LU)
        READ(1) (EMT2(K2,K1,L) , L=1,LU)
20    CONTINUE
      DO 30 K2=1,KU4
       DO 30 K1=1,KU4
        READ(1) (FMT1(K2,K1,L) , L=1,LU)
        READ(1) (FMT2(K2,K1,L) , L=1,LU)
30    CONTINUE
      DO 40 K2=1,KT2
       DO 40 K1=1,KT2
        READ(1) (TMT1(K2,K1,L) , L=0,LT)
        READ(1) (TMT2(K2,K1,L) , L=0,LT)
40    CONTINUE
*
      READ(2) TT1
      READ(2) TT2
      READ(2) TT3
      READ(2) TT4
      READ(2) TT5
      READ(2) RNV
      DO 50 M=0,MM
       DO 50 L=M,LL
        READ(2) (PP1(L,M,IL) , IL=1,LN)
        READ(2) (PP2(L,M,IL) , IL=1,LN)
        READ(2) (PP3(L,M,IL) , IL=1,LN)
50    CONTINUE
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
      DO 60 IL=1,LN
       DO 60 L=1,LN
        READ(2) (PPP1(IL,L,M) , M=0,MM)
        READ(2) (PPP2(IL,L,M) , M=0,MM)
        READ(2) (PPP4(IL,L,M) , M=0,MM)
60    CONTINUE
*
      READ(2) TTT
      DO 70 IL=1,LN
       DO 70 L=0,LN
        READ(2) (PPP(IL,L,M) , M=0,MM)
70    CONTINUE
*
      READ(2) CALPHA
      READ(2) SALPHA
*
      RETURN
      END
*******************************************************************************
      SUBROUTINE STEP(TIME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
      include 'omp_lib.h'
*
      COMMON /FIELD/  GC(KB2, LB, 0:MB) ,  GS(KB2, LB, 0:MB)
     &             , DGC(KB2, LB, 0:MB) , DGS(KB2, LB, 0:MB)
     &             ,  HC(KB2, LB, 0:MB) ,  HS(KB2, LB, 0:MB)
     &             , DHC(KB2, LB, 0:MB) , DHS(KB2, LB, 0:MB)
*
      COMMON /FLOW/   EC(KU2, LU, 0:MU) ,  ES(KU2, LU, 0:MU)
     &             , DEC(KU2, LU, 0:MU) , DES(KU2, LU, 0:MU)
     &             ,  FC(KU4, LU, 0:MU) ,  FS(KU4, LU, 0:MU)
     &             , DFC(KU4, LU, 0:MU) , DFS(KU4, LU, 0:MU)
*
      COMMON /TEMP/   TC(KT2,0:LT,0:MT) ,  TS(KT2,0:LT,0:MT)
     &             , DTC(KT2,0:LT,0:MT) , DTS(KT2,0:LT,0:MT)
*
      COMMON /STEP1/ GMT1(KB2, KB2, LB), HMT1(KB2, KB2, LB)
     &             , GMT2(KB2, KB2, LB), HMT2(KB2, KB2, LB)
*
      COMMON /STEP2/ EMT1(KU2, KU2, LU), FMT1(KU4, KU4, LU)
     &             , EMT2(KU2, KU2, LU), FMT2(KU4, KU4, LU)
*
      COMMON /STEP3/ TMT1(KT2, KT2, 0:LT)
     &             , TMT2(KT2, KT2, 0:LT)
*
      DIMENSION WGC1(KB2, LB, 0:MB) , WGC2(KB2, LB, 0:MB)
      DIMENSION WGS1(KB2, LB, 0:MB) , WGS2(KB2, LB, 0:MB)
      DIMENSION WHC1(KB2, LB, 0:MB) , WHC2(KB2, LB, 0:MB)
      DIMENSION WHS1(KB2, LB, 0:MB) , WHS2(KB2, LB, 0:MB)
*
      DIMENSION WEC1(KU2, LU, 0:MU) , WEC2(KU2, LU, 0:MU)
      DIMENSION WES1(KU2, LU, 0:MU) , WES2(KU2, LU, 0:MU)
      DIMENSION WFC1(KU4, LU, 0:MU) , WFC2(KU4, LU, 0:MU)
      DIMENSION WFS1(KU4, LU, 0:MU) , WFS2(KU4, LU, 0:MU)
*
      DIMENSION WTC1(KT2,0:LT,0:MT) , WTC2(KT2,0:LT,0:MT)
      DIMENSION WTS1(KT2,0:LT,0:MT) , WTS2(KT2,0:LT,0:MT)
******
!$omp parallel do private(M,LS,L,K1,K2) schedule(static)
      DO 10 M=0,MB
       LS=MAX(M,1)
       DO 10 L=LS,LB
        DO 10 K1=1,KB2
         WGC1(K1,L,M)=0.D0
         WGS1(K1,L,M)=0.D0
         WHC1(K1,L,M)=0.D0
         WHS1(K1,L,M)=0.D0
         DO 10 K2=1,KB2
          WGC1(K1,L,M)=WGC1(K1,L,M) + GMT1(K2,K1,L)*GC(K2,L,M)
          WGS1(K1,L,M)=WGS1(K1,L,M) + GMT1(K2,K1,L)*GS(K2,L,M)
          WHC1(K1,L,M)=WHC1(K1,L,M) + HMT1(K2,K1,L)*HC(K2,L,M)
          WHS1(K1,L,M)=WHS1(K1,L,M) + HMT1(K2,K1,L)*HS(K2,L,M)
10    CONTINUE
!$omp end parallel do
!$omp parallel do private(M,LS,L,K1,K2) schedule(static)
      DO 15 M=0,MU
       LS=MAX(M,1)
       DO 15 L=LS,LU
        DO 15 K1=1,KU2
         WEC1(K1,L,M)=0.D0
         WES1(K1,L,M)=0.D0
         DO 15 K2=1,KU2
          WEC1(K1,L,M)=WEC1(K1,L,M) + EMT1(K2,K1,L)*EC(K2,L,M)
          WES1(K1,L,M)=WES1(K1,L,M) + EMT1(K2,K1,L)*ES(K2,L,M)
15    CONTINUE
!$omp end parallel do
!$omp parallel do private(M,LS,L,K1,K2) schedule(static)
      DO 20 M=0,MU
       LS=MAX(M,1)
       DO 20 L=LS,LU
        DO 20 K1=1,KU4
         WFC1(K1,L,M)=0.D0
         WFS1(K1,L,M)=0.D0
         DO 20 K2=1,KU4
          WFC1(K1,L,M)=WFC1(K1,L,M) + FMT1(K2,K1,L)*FC(K2,L,M)
          WFS1(K1,L,M)=WFS1(K1,L,M) + FMT1(K2,K1,L)*FS(K2,L,M)
20    CONTINUE
!$omp end parallel do
!$omp parallel do private(M,L,K1,K2) schedule(static)
      DO 25 M=0,MT
       DO 25 L=M,LT
        DO 25 K1=1,KT2
         WTC1(K1,L,M)=0.D0
         WTS1(K1,L,M)=0.D0
         DO 25 K2=1,KT2
          WTC1(K1,L,M)=WTC1(K1,L,M) + TMT1(K2,K1,L)*TC(K2,L,M)
          WTS1(K1,L,M)=WTS1(K1,L,M) + TMT1(K2,K1,L)*TS(K2,L,M)
25    CONTINUE
!$omp end parallel do
*
      CALL NONLIN(TIME)

!$omp parallel do private(M,LS,L,K1,K2,K) schedule(static)
      DO 35 M=0,MB
       LS=MAX(M,1)
       DO 35 L=LS,LB
        DO 30 K1=1,KB2
         WGC2(K1,L,M)=0.D0
         WGS2(K1,L,M)=0.D0
         WHC2(K1,L,M)=0.D0
         WHS2(K1,L,M)=0.D0
         DO 30 K2=1,KB2
          WGC2(K1,L,M)=WGC2(K1,L,M) + GMT2(K2,K1,L)*DGC(K2,L,M)
          WGS2(K1,L,M)=WGS2(K1,L,M) + GMT2(K2,K1,L)*DGS(K2,L,M)
          WHC2(K1,L,M)=WHC2(K1,L,M) + HMT2(K2,K1,L)*DHC(K2,L,M)
          WHS2(K1,L,M)=WHS2(K1,L,M) + HMT2(K2,K1,L)*DHS(K2,L,M)
30      CONTINUE
        DO 35 K=1,KB2
         GC(K,L,M)=WGC1(K,L,M) + WGC2(K,L,M)
         GS(K,L,M)=WGS1(K,L,M) + WGS2(K,L,M)
         HC(K,L,M)=WHC1(K,L,M) + WHC2(K,L,M)
         HS(K,L,M)=WHS1(K,L,M) + WHS2(K,L,M)
35    CONTINUE
!$omp end parallel do
!$omp parallel do private(M,LS,L,K1,K2,K) schedule(static)
      DO 45 M=0,MU
       LS=MAX(M,1)
       DO 45 L=LS,LU
        DO 40 K1=1,KU2
         WEC2(K1,L,M)=0.D0
         WES2(K1,L,M)=0.D0
         DO 40 K2=1,KU2
          WEC2(K1,L,M)=WEC2(K1,L,M) + EMT2(K2,K1,L)*DEC(K2,L,M)
          WES2(K1,L,M)=WES2(K1,L,M) + EMT2(K2,K1,L)*DES(K2,L,M)
40      CONTINUE
        DO 45 K=1,KU2
         EC(K,L,M)=WEC1(K,L,M) + WEC2(K,L,M)
         ES(K,L,M)=WES1(K,L,M) + WES2(K,L,M)
45    CONTINUE
!$omp end parallel do
!$omp parallel do private(M,LS,L,K1,K2,K) schedule(static)
      DO 55 M=0,MU
       LS=MAX(M,1)
       DO 55 L=LS,LU
        DO 50 K1=1,KU4
         WFC2(K1,L,M)=0.D0
         WFS2(K1,L,M)=0.D0
         DO 50 K2=1,KU4
          WFC2(K1,L,M)=WFC2(K1,L,M) + FMT2(K2,K1,L)*DFC(K2,L,M)
          WFS2(K1,L,M)=WFS2(K1,L,M) + FMT2(K2,K1,L)*DFS(K2,L,M)
50      CONTINUE
        DO 55 K=1,KU4
         FC(K,L,M)=WFC1(K,L,M) + WFC2(K,L,M)
         FS(K,L,M)=WFS1(K,L,M) + WFS2(K,L,M)
55    CONTINUE
!$omp end parallel do
!$omp parallel do private(M,L,K1,K2,K) schedule(static)
      DO 65 M=0,MT
       DO 65 L=M,LT
        DO 60 K1=1,KT2
         WTC2(K1,L,M)=0.D0
         WTS2(K1,L,M)=0.D0
         DO 60 K2=1,KT2
          WTC2(K1,L,M)=WTC2(K1,L,M) + TMT2(K2,K1,L)*DTC(K2,L,M)
          WTS2(K1,L,M)=WTS2(K1,L,M) + TMT2(K2,K1,L)*DTS(K2,L,M)
60      CONTINUE
        DO 65 K=1,KT2
         TC(K,L,M)=WTC1(K,L,M) + WTC2(K,L,M)
         TS(K,L,M)=WTS1(K,L,M) + WTS2(K,L,M)
65    CONTINUE
!$omp end parallel do
*
      CALL NONLIN(TIME+DT)

!$omp parallel do private(M,LS,L,K1,K2,K) schedule(static)
      DO 80 M=0,MB
       LS=MAX(M,1)
       DO 80 L=LS,LB
        DO 75 K1=1,KB2
         DO 70 K2=1,KB2
          WGC2(K1,L,M)=WGC2(K1,L,M) + GMT2(K2,K1,L)*DGC(K2,L,M)
          WGS2(K1,L,M)=WGS2(K1,L,M) + GMT2(K2,K1,L)*DGS(K2,L,M)
          WHC2(K1,L,M)=WHC2(K1,L,M) + HMT2(K2,K1,L)*DHC(K2,L,M)
          WHS2(K1,L,M)=WHS2(K1,L,M) + HMT2(K2,K1,L)*DHS(K2,L,M)
70       CONTINUE
         WGC2(K1,L,M)=0.5D0*WGC2(K1,L,M)
         WGS2(K1,L,M)=0.5D0*WGS2(K1,L,M)
         WHC2(K1,L,M)=0.5D0*WHC2(K1,L,M)
         WHS2(K1,L,M)=0.5D0*WHS2(K1,L,M)
75      CONTINUE
        DO 80 K=1,KB2
         GC(K,L,M)=WGC1(K,L,M) + WGC2(K,L,M)
         GS(K,L,M)=WGS1(K,L,M) + WGS2(K,L,M)
         HC(K,L,M)=WHC1(K,L,M) + WHC2(K,L,M)
         HS(K,L,M)=WHS1(K,L,M) + WHS2(K,L,M)
80    CONTINUE
!$omp end parallel do
!$omp parallel do private(M,LS,L,K1,K2,K) schedule(static)
      DO 95 M=0,MU
       LS=MAX(M,1)
       DO 95 L=LS,LU
        DO 90 K1=1,KU2
         DO 85 K2=1,KU2
          WEC2(K1,L,M)=WEC2(K1,L,M) + EMT2(K2,K1,L)*DEC(K2,L,M)
          WES2(K1,L,M)=WES2(K1,L,M) + EMT2(K2,K1,L)*DES(K2,L,M)
85       CONTINUE
         WEC2(K1,L,M)=0.5D0*WEC2(K1,L,M)
         WES2(K1,L,M)=0.5D0*WES2(K1,L,M)
90      CONTINUE
        DO 95 K=1,KU2
         EC(K,L,M)=WEC1(K,L,M) + WEC2(K,L,M)
         ES(K,L,M)=WES1(K,L,M) + WES2(K,L,M)
95    CONTINUE
!$omp end parallel do
!$omp parallel do private(M,LS,L,K1,K2,K) schedule(static)
      DO 110 M=0,MU
       LS=MAX(M,1)
       DO 110 L=LS,LU
        DO 105 K1=1,KU4
         DO 100 K2=1,KU4
          WFC2(K1,L,M)=WFC2(K1,L,M) + FMT2(K2,K1,L)*DFC(K2,L,M)
          WFS2(K1,L,M)=WFS2(K1,L,M) + FMT2(K2,K1,L)*DFS(K2,L,M)
100      CONTINUE
         WFC2(K1,L,M)=0.5D0*WFC2(K1,L,M)
         WFS2(K1,L,M)=0.5D0*WFS2(K1,L,M)
105     CONTINUE
        DO 110 K=1,KU4
         FC(K,L,M)=WFC1(K,L,M) + WFC2(K,L,M)
         FS(K,L,M)=WFS1(K,L,M) + WFS2(K,L,M)
110   CONTINUE
!$omp end parallel do
!$omp parallel do private(M,L,K1,K2,K) schedule(static)
      DO 125 M=0,MT
       DO 125 L=M,LT
        DO 120 K1=1,KT2
         DO 115 K2=1,KT2
          WTC2(K1,L,M)=WTC2(K1,L,M) + TMT2(K2,K1,L)*DTC(K2,L,M)
          WTS2(K1,L,M)=WTS2(K1,L,M) + TMT2(K2,K1,L)*DTS(K2,L,M)
115      CONTINUE
         WTC2(K1,L,M)=0.5D0*WTC2(K1,L,M)
         WTS2(K1,L,M)=0.5D0*WTS2(K1,L,M)
120     CONTINUE
        DO 125 K=1,KT2
         TC(K,L,M)=WTC1(K,L,M) + WTC2(K,L,M)
         TS(K,L,M)=WTS1(K,L,M) + WTS2(K,L,M)
125   CONTINUE
!$omp end parallel do
*
      RETURN
      END
*******************************************************************************
      SUBROUTINE NONLIN(TIME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
      include 'omp_lib.h'
*
      COMMON /FIELD/  GC(KB2, LB, 0:MB) ,  GS(KB2, LB, 0:MB)
     &             , DGC(KB2, LB, 0:MB) , DGS(KB2, LB, 0:MB)
     &             ,  HC(KB2, LB, 0:MB) ,  HS(KB2, LB, 0:MB)
     &             , DHC(KB2, LB, 0:MB) , DHS(KB2, LB, 0:MB)
*
      COMMON /FLOW/   EC(KU2, LU, 0:MU) ,  ES(KU2, LU, 0:MU)
     &             , DEC(KU2, LU, 0:MU) , DES(KU2, LU, 0:MU)
     &             ,  FC(KU4, LU, 0:MU) ,  FS(KU4, LU, 0:MU)
     &             , DFC(KU4, LU, 0:MU) , DFS(KU4, LU, 0:MU)
*
      COMMON /TEMP/   TC(KT2,0:LT,0:MT) ,  TS(KT2,0:LT,0:MT)
     &             , DTC(KT2,0:LT,0:MT) , DTS(KT2,0:LT,0:MT)
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
      COMMON /DTCS/ TTT(KN, KN) , PPP(LN, 0:LN, 0:MM)
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
*
      DIMENSION WTC2(0:LT, 0:MT) , WTC3(0:LT, 0:MT) , WTC(LN)
      DIMENSION WTS2(0:LT, 0:MT) , WTS3(0:LT, 0:MT) , WTS(LN)
*
      DIMENSION WTC23(0:MT), WTS23(0:MT), WTC32(0:MT), WTS32(0:MT)
      DIMENSION WTC33(0:MT), WTS33(0:MT)
*
      DIMENSION B1(2*MN, LN, KN) , B2(2*MN, LN, KN) , B3(2*MN, LN, KN)
      DIMENSION U1(2*MN, LN, KN) , U2(2*MN, LN, KN) , U3(2*MN, LN, KN)
      DIMENSION F1(2*MN, LN, KN) , F2(2*MN, LN, KN) , F3(2*MN, LN, KN)
*
      DIMENSION UXB1(2*MN) , UXB2(2*MN) , UXB3(2*MN), UT(2*MN, LN, KN)
*
      COMMON /CURL1/ F1C(KN,LN,0:MM), F2C(KN,LN,0:MM), F3C(KN,LN,0:MM)
     &             , F1S(KN,LN,0:MM), F2S(KN,LN,0:MM), F3S(KN,LN,0:MM)
******
*     Evaluate b
!$omp parallel do private(IK,M,LS,L,K,IL,IM, 
!$omp& WHC1,WHS1,WHC2,WHS2,WGC3,WGS3,WHC11,WHS11,WHC22,WHS22,
!$omp& WHC23,WHS23,WGC32,WGS32,WGC33,WGS33) schedule(static)
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
!$omp end parallel do
*     Evaluate curl b
*     Note that the work arrays have been reused from the calcu-
*     lation of b, so the names no longer necessarily make sense.
!$omp parallel do private(IK,M,LS,L,K,IL,IM, 
!$omp& WHC1,WHS1,WHC2,WHS2,WHC4,WHS4,WHC5,WHS5,WGC3,WGS3,WHC11,WHS11,
!$omp& WHC22,WHS22,WHC23,WHS23,WGC32,WGS32,WGC33,WGS33,
!$omp& CB1,CB2,CB3) schedule(static)
      DO 45 IK=1,KN
       DO 30 M=0,MB
        LS=MAX(M,1)
        DO 30 L=LS,LB
         WHC1(L,M)=0.D0
         WHS1(L,M)=0.D0
         WHC2(L,M)=0.D0
         WHS2(L,M)=0.D0
         WHC4     =0.D0
         WHS4     =0.D0
         WHC5     =0.D0
         WHS5     =0.D0
         DO 25 K=1,KB2
          WHC1(L,M)=WHC1(L,M)  +  GC(K,L,M) * TT1(K,IK)
          WHS1(L,M)=WHS1(L,M)  +  GS(K,L,M) * TT1(K,IK)
          WHC2(L,M)=WHC2(L,M)  +  GC(K,L,M) * TT2(K,IK)
          WHS2(L,M)=WHS2(L,M)  +  GS(K,L,M) * TT2(K,IK)
          WHC4     =WHC4       +  HC(K,L,M) * TT4(K,IK)
          WHS4     =WHS4       +  HS(K,L,M) * TT4(K,IK)
          WHC5     =WHC5       +  HC(K,L,M) * TT5(K,IK)
          WHS5     =WHS5       +  HS(K,L,M) * TT5(K,IK)
25       CONTINUE
         WGC3(L,M)=XLL(L)*WHC4 - WHC5
         WGS3(L,M)=XLL(L)*WHS4 - WHS5
30     CONTINUE
       DO 45 IL=1,LN
        DO 35 M=0,MB
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
         DO 35 L=LS,LB
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
35      CONTINUE
        DO 45 IM=1,2*MN
         CB1=0.D0
         CB2=0.D0
         CB3=0.D0
         DO 40 M=0,MB
          CB1=CB1 + WHC11(M)* CS(M,IM) + WHS11(M)* SN(M,IM)
          CB2=CB2 + WHC22(M)* CS(M,IM) + WHS22(M)* SN(M,IM)
     &            - WGC33(M)*SNM(M,IM) + WGS33(M)*CSM(M,IM)
          CB3=CB3 - WHC23(M)*SNM(M,IM) + WHS23(M)*CSM(M,IM)
     &            - WGC32(M)* CS(M,IM) - WGS32(M)* SN(M,IM)
40       CONTINUE
* lorentz force
         F1(IM,IL,IK)=CB2*B3(IM,IL,IK) - CB3*B2(IM,IL,IK)
         F2(IM,IL,IK)=CB3*B1(IM,IL,IK) - CB1*B3(IM,IL,IK)
         F3(IM,IL,IK)=CB1*B2(IM,IL,IK) - CB2*B1(IM,IL,IK)
45    CONTINUE
!$omp end parallel do
*     Evaluate u
!$omp parallel do private(IK,M,LS,L,K,IL,IM, 
!$omp& WFC1,WFS1,WFC2,WFS2,WEC3,WES3,WFC11,WFS11,WFC22,WFS22,
!$omp& WFC23,WFS23,WEC32,WES32,WEC33,WES33) schedule(static)
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
!$omp end parallel do
* coriolis force
!$omp parallel do private(IK,IL,IM) schedule(static)
      DO IK=1,KN
       DO IL=1,LN
        DO IM=1,2*MN
         F1(IM,IL,IK)=F1(IM,IL,IK)+ 2.D0*STH(IL)*U3(IM,IL,IK)
         F2(IM,IL,IK)=F2(IM,IL,IK)+ 2.D0*CTH(IL)*U3(IM,IL,IK)
         F3(IM,IL,IK)=F3(IM,IL,IK)- 2.D0*STH(IL)*U1(IM,IL,IK)
     &                            - 2.D0*CTH(IL)*U2(IM,IL,IK)
        ENDDO
       ENDDO
      ENDDO
!$omp end parallel do
*     Evaluate curl u
*     Note that the work arrays have been reused from the calcu-
*     lation of u, so the names no longer necessarily make sense.
!$omp parallel do private(IK,M,LS,L,K,IL,IM, 
!$omp& WFC1,WFS1,WFC2,WFS2,WFC4,WFS4,WFC5,WFS5,WEC3,WES3,WFC11,WFS11,
!$omp& WFC22,WFS22,WFC23,WFS23,WEC32,WES32,WEC33,WES33,
!$omp& CU1,CU2,CU3) schedule(static)
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
         CU1=0.D0
         CU2=0.D0
         CU3=0.D0
         DO 75 M=0,MU
          CU1=CU1 + WFC11(M)* CS(M,IM) + WFS11(M)* SN(M,IM)
          CU2=CU2 + WFC22(M)* CS(M,IM) + WFS22(M)* SN(M,IM)
     &            - WEC33(M)*SNM(M,IM) + WES33(M)*CSM(M,IM)
          CU3=CU3 - WFC23(M)*SNM(M,IM) + WFS23(M)*CSM(M,IM)
     &            - WEC32(M)* CS(M,IM) - WES32(M)* SN(M,IM)
75       CONTINUE
* inertial force          
         F1(IM,IL,IK)=F1(IM,IL,IK)
     &              - CU2*U3(IM,IL,IK) + CU3*U2(IM,IL,IK)
         F2(IM,IL,IK)=F2(IM,IL,IK)
     &              - CU3*U1(IM,IL,IK) + CU1*U3(IM,IL,IK)
         F3(IM,IL,IK)=F3(IM,IL,IK)
     &              - CU1*U2(IM,IL,IK) + CU2*U1(IM,IL,IK)
80    CONTINUE
!$omp end parallel do
* precession: Poincare force
*      PRECX= SALPHA*DCOS(TIME)
*      PRECY=-SALPHA*DSIN(TIME)
*      PRECZ= CALPHA
*!$omp parallel do private(IK,IL,IM,R,
*!$omp& PREC1,PREC2,PREC3,POIN1,POIN2,POIN3) schedule(static)
*      DO IK=1,KN
*       R=1.D0/RNV(IK)
*       DO IL=1,LN
*        DO IM=1,2*MN
*         PREC1=PRECX*STH(IL)*CPHI(IM)+PRECY*STH(IL)*SPHI(IM)
*     &        +PRECZ*CTH(IL)
*         PREC2=PRECX*CTH(IL)*CPHI(IM)+PRECY*CTH(IL)*SPHI(IM)
*     &        -PRECZ*STH(IL)
*         PREC3=-PRECX*SPHI(IM)+PRECY*CPHI(IM)
*         POIN1=2.D0*(U2(IM,IL,IK)*PREC3-U3(IM,IL,IK)*PREC2)
*         POIN2=2.D0*(U3(IM,IL,IK)*PREC1-U1(IM,IL,IK)*PREC3)
*         POIN3=2.D0*(U1(IM,IL,IK)*PREC2-U2(IM,IL,IK)*PREC1)
*         POIN2=POIN2+R*(CTH(IL)*PREC2+STH(IL)*PREC1)
*         POIN3=POIN3+R*CTH(IL)*PREC3
*         F1(IM,IL,IK)=F1(IM,IL,IK)+ PO*POIN1
*         F2(IM,IL,IK)=F2(IM,IL,IK)+ PO*POIN2
*         F3(IM,IL,IK)=F3(IM,IL,IK)- PO*POIN3
*        ENDDO
*       ENDDO
*      ENDDO
*!$omp end parallel do
*     Evaluate @ and (u . grad) @
!$omp parallel do private(IK,IL,IM,K,L,M,R,
!$omp& WTC3,WTS3,WTC2,WTS2,WTC33,WTS33,WTC23,WTS23,WTC32,WTS32,
!$omp& TH,RT,TT,PT) schedule(static)
      DO 110 IK=1,KN
       R=1.D0/RNV(IK)
       DO 90 M=0,MT
        DO 90 L=M,LT
         WTC3(L,M)=0.D0
         WTS3(L,M)=0.D0
         WTC2(L,M)=0.D0
         WTS2(L,M)=0.D0
         DO 85 K=1,KT2
          WTC3(L,M)=WTC3(L,M)  +  TC(K,L,M) * TT3(K,IK)
          WTS3(L,M)=WTS3(L,M)  +  TS(K,L,M) * TT3(K,IK)
          WTC2(L,M)=WTC2(L,M)  +  TC(K,L,M) * TT2(K,IK)
          WTS2(L,M)=WTS2(L,M)  +  TS(K,L,M) * TT2(K,IK)
85       CONTINUE
         WTC3(L,M)=WTC3(L,M) * R
         WTS3(L,M)=WTS3(L,M) * R
         WTC2(L,M)=WTC2(L,M) * R
         WTS2(L,M)=WTS2(L,M) * R
90     CONTINUE
       DO 110 IL=1,LN
        DO 100 M=0,MT
         WTC33(M)=0.D0
         WTS33(M)=0.D0
         WTC23(M)=0.D0
         WTS23(M)=0.D0
         WTC32(M)=0.D0
         WTS32(M)=0.D0
         DO 95 L=M,LT
          WTC33(M)=WTC33(M)  +  WTC3(L,M) * PP3(L,M,IL)
          WTS33(M)=WTS33(M)  +  WTS3(L,M) * PP3(L,M,IL)
          WTC23(M)=WTC23(M)  +  WTC2(L,M) * PP3(L,M,IL)
          WTS23(M)=WTS23(M)  +  WTS2(L,M) * PP3(L,M,IL)
          WTC32(M)=WTC32(M)  +  WTC3(L,M) * PP2(L,M,IL)
          WTS32(M)=WTS32(M)  +  WTS3(L,M) * PP2(L,M,IL)
95       CONTINUE
         WTC33(M)=WTC33(M) * STH(IL)
         WTS33(M)=WTS33(M) * STH(IL)
         WTC23(M)=WTC23(M) * STH(IL)
         WTS23(M)=WTS23(M) * STH(IL)
100     CONTINUE
        DO 110 IM=1,2*MN
         TH=0.D0
         RT=0.D0
         TT=0.D0
         PT=0.D0
         DO 105 M=0,MT
          TH=TH + WTC33(M)* CS(M,IM) + WTS33(M)* SN(M,IM)
          RT=RT + WTC23(M)* CS(M,IM) + WTS23(M)* SN(M,IM)
          TT=TT + WTC32(M)* CS(M,IM) + WTS32(M)* SN(M,IM)
          PT=PT - WTC33(M)*SNM(M,IM) + WTS33(M)*CSM(M,IM)
105      CONTINUE
* buoyancy force
         F1(IM,IL,IK)=F1(IM,IL,IK)  +  AA*TH*R
         UT(IM,IL,IK)=- U1(IM,IL,IK) * RT
     &                - U2(IM,IL,IK) * TT * RNV(IK)
     &                - U3(IM,IL,IK) * PT * RNV(IK)*CSC(IL)
* stable stratification
*     &                - U1(IM,IL,IK)
* unstable stratification
     &                + aspect*U1(IM,IL,IK)
110   CONTINUE
!$omp end parallel do
*
*     Evaluate the spectral coefficients DHC, DHS, DGC, DGS.
*     Given u and b at the collocation points, compute u x b
*     and separate out its different m modes again.  Note that
*     this is a slow Fourier transform.
!$omp parallel do private(IK,IL,IM,M,UXB1,UXB2,UXB3) schedule(static)
      DO 120 IK=1,KN
       DO 120 IL=1,LN
        DO 115 IM=1,2*MN
         UXB1(IM)=U2(IM,IL,IK)*B3(IM,IL,IK) - U3(IM,IL,IK)*B2(IM,IL,IK)
         UXB2(IM)=U3(IM,IL,IK)*B1(IM,IL,IK) - U1(IM,IL,IK)*B3(IM,IL,IK)
         UXB3(IM)=U1(IM,IL,IK)*B2(IM,IL,IK) - U2(IM,IL,IK)*B1(IM,IL,IK)
115     CONTINUE
        DO 120 M=0,MB
         F1C(IK,IL,M)=0.D0
         F1S(IK,IL,M)=0.D0
         F2C(IK,IL,M)=0.D0
         F2S(IK,IL,M)=0.D0
         F3C(IK,IL,M)=0.D0
         F3S(IK,IL,M)=0.D0
         DO 120 IM=1,2*MN
          F1C(IK,IL,M)=F1C(IK,IL,M)  +  CT(IM,M) * UXB1(IM)
          F1S(IK,IL,M)=F1S(IK,IL,M)  +  ST(IM,M) * UXB1(IM)
          F2C(IK,IL,M)=F2C(IK,IL,M)  +  CT(IM,M) * UXB2(IM)
          F2S(IK,IL,M)=F2S(IK,IL,M)  +  ST(IM,M) * UXB2(IM)
          F3C(IK,IL,M)=F3C(IK,IL,M)  +  CT(IM,M) * UXB3(IM)
          F3S(IK,IL,M)=F3S(IK,IL,M)  +  ST(IM,M) * UXB3(IM)
120   CONTINUE
!$omp end parallel do
      CALL CURL(DHC, DHS, KB2, DGC, DGS, KB2, KB, LB, MB)
*
*     Evaluate the spectral coefficients DEC, DES, DFC, DFS.
!$omp parallel do private(IK,IL,IM,M) schedule(static)
      DO 125 IK=1,KN
       DO 125 IL=1,LN
        DO 125 M=0,MU
         F1C(IK,IL,M)=0.D0
         F1S(IK,IL,M)=0.D0
         F2C(IK,IL,M)=0.D0
         F2S(IK,IL,M)=0.D0
         F3C(IK,IL,M)=0.D0
         F3S(IK,IL,M)=0.D0
         DO 125 IM=1,2*MN
          F1C(IK,IL,M)=F1C(IK,IL,M)  +  CT(IM,M) * F1(IM,IL,IK)
          F1S(IK,IL,M)=F1S(IK,IL,M)  +  ST(IM,M) * F1(IM,IL,IK)
          F2C(IK,IL,M)=F2C(IK,IL,M)  +  CT(IM,M) * F2(IM,IL,IK)
          F2S(IK,IL,M)=F2S(IK,IL,M)  +  ST(IM,M) * F2(IM,IL,IK)
          F3C(IK,IL,M)=F3C(IK,IL,M)  +  CT(IM,M) * F3(IM,IL,IK)
          F3S(IK,IL,M)=F3S(IK,IL,M)  +  ST(IM,M) * F3(IM,IL,IK)
125   CONTINUE
!$omp end parallel do
      CALL CURL(DEC, DES, KU2, DFC, DFS, KU4, KU, LU, MU)
*
*     Evaluate the spectral coefficients DTC, DTS.
!$omp parallel do private(IK,IL,IM,M) schedule(static)
      DO 130 IK=1,KN
       DO 130 IL=1,LN
        DO 130 M=0,MT
         F1C(IK,IL,M)=0.D0
         F1S(IK,IL,M)=0.D0
         DO 130 IM=1,2*MN
          F1C(IK,IL,M)=F1C(IK,IL,M)  +  CT(IM,M) * UT(IM,IL,IK)
          F1S(IK,IL,M)=F1S(IK,IL,M)  +  ST(IM,M) * UT(IM,IL,IK)
130   CONTINUE
!$omp end parallel do
!$omp parallel do private(M,K,IL,IK,L,WTC,WTS) schedule(static)
      DO 140 M=0,MT
       DO 140 K=1,KT
        DO 135 IL=1,LN
         WTC(IL)=0.D0
         WTS(IL)=0.D0
         DO 135 IK=1,KN
          WTC(IL)=WTC(IL)  +  TTT(IK,K) * F1C(IK,IL,M)
          WTS(IL)=WTS(IL)  +  TTT(IK,K) * F1S(IK,IL,M)
135     CONTINUE
        DO 140 L=M,LT
         DTC(K,L,M)=0.D0
         DTS(K,L,M)=0.D0
         DO 140 IL=1,LN
          DTC(K,L,M)=DTC(K,L,M)  +  WTC(IL) * PPP(IL,L,M)
          DTS(K,L,M)=DTS(K,L,M)  +  WTS(IL) * PPP(IL,L,M)
140   CONTINUE
!$omp end parallel do
*
      RETURN
      END
*******************************************************************************
      SUBROUTINE CURL(DUC, DUS, KDU, DVC, DVS, KDV, KX, LX, MX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
*
      COMMON /CURL0/ TTT1(KN, KN) , PPP1(LN, LN, 0:MM)
     &             , TTT3(KN, KN) , PPP2(LN, LN, 0:MM)
     &             , TTT4(KN, KN) , PPP4(LN, LN, 0:MM)
*
      COMMON /CURL1/ F1C(KN,LN,0:MM), F2C(KN,LN,0:MM), F3C(KN,LN,0:MM)
     &             , F1S(KN,LN,0:MM), F2S(KN,LN,0:MM), F3S(KN,LN,0:MM)
*
      DIMENSION DUC(KDU, LX, 0:MX) , DUS(KDU, LX, 0:MX)
      DIMENSION DVC(KDV, LX, 0:MX) , DVS(KDV, LX, 0:MX)
*
      DIMENSION W13C(LN) , W12C(LN) , W42C(LN) , W31C(LN) , W43C(LN)
      DIMENSION W13S(LN) , W12S(LN) , W42S(LN) , W31S(LN) , W43S(LN)
******
*     On input F1C, F1S, F2C, F2S, F3C, F3S contain the three components
*     of some vector, separated out into the azimuthal modes M, and at the
*     collocation points IK, IL.  On output, DUC, DUS, DVC, DVS contain
*     the spectral coefficients of the r components of the curl and the
*     curl of the curl of the original vector.
******
!$omp parallel do private(M,LS,K,IL,IK,L,
!$omp& W13C,W12C,W42C,W31C,W43C,W13S,W12S,W42S,W31S,W43S) 
!$omp& schedule(static)
      DO 20 M=0,MX
       LS=MAX(M,1)
       DO 20 K=1,KX
        DO 10 IL=1,LN
         W13C(IL)=0.D0
         W12C(IL)=0.D0
         W42C(IL)=0.D0
         W31C(IL)=0.D0
         W43C(IL)=0.D0
         W13S(IL)=0.D0
         W12S(IL)=0.D0
         W42S(IL)=0.D0
         W31S(IL)=0.D0
         W43S(IL)=0.D0
         DO 10 IK=1,KN
          W13C(IL)=W13C(IL)  +  TTT1(IK,K) * F3C(IK,IL,M)
          W12C(IL)=W12C(IL)  +  TTT1(IK,K) * F2C(IK,IL,M)
          W42C(IL)=W42C(IL)  +  TTT4(IK,K) * F2C(IK,IL,M)
          W31C(IL)=W31C(IL)  +  TTT3(IK,K) * F1C(IK,IL,M)
          W43C(IL)=W43C(IL)  +  TTT4(IK,K) * F3C(IK,IL,M)
          W13S(IL)=W13S(IL)  +  TTT1(IK,K) * F3S(IK,IL,M)
          W12S(IL)=W12S(IL)  +  TTT1(IK,K) * F2S(IK,IL,M)
          W42S(IL)=W42S(IL)  +  TTT4(IK,K) * F2S(IK,IL,M)
          W31S(IL)=W31S(IL)  +  TTT3(IK,K) * F1S(IK,IL,M)
          W43S(IL)=W43S(IL)  +  TTT4(IK,K) * F3S(IK,IL,M)
10      CONTINUE
        DO 20 L=LS,LX
         DUC(K,L,M)=0.D0
         DUS(K,L,M)=0.D0
         DVC(K,L,M)=0.D0
         DVS(K,L,M)=0.D0
         DO 20 IL=1,LN
          DUC(K,L,M)=DUC(K,L,M)  +  W13C(IL) * PPP1(IL,L,M)
     &                           -  W12S(IL) * PPP2(IL,L,M)
          DUS(K,L,M)=DUS(K,L,M)  +  W13S(IL) * PPP1(IL,L,M)
     &                           +  W12C(IL) * PPP2(IL,L,M)
          DVC(K,L,M)=DVC(K,L,M)  +  W42C(IL) * PPP1(IL,L,M)
     &                           +  W31C(IL) * PPP4(IL,L,M)
     &                           +  W43S(IL) * PPP2(IL,L,M)
          DVS(K,L,M)=DVS(K,L,M)  +  W42S(IL) * PPP1(IL,L,M)
     &                           +  W31S(IL) * PPP4(IL,L,M)
     &                           -  W43C(IL) * PPP2(IL,L,M)
20    CONTINUE
!$omp end parallel do
      RETURN
      END
*****************************************************************************
      SUBROUTINE ENER(TORB,POLB,TORU,POLU,EM,EMPHI,EK,EKPHI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'      
*
      COMMON /FIELD/  GC(KB2, LB, 0:MB) ,  GS(KB2, LB, 0:MB)
     &             , DGC(KB2, LB, 0:MB) , DGS(KB2, LB, 0:MB)
     &             ,  HC(KB2, LB, 0:MB) ,  HS(KB2, LB, 0:MB)
     &             , DHC(KB2, LB, 0:MB) , DHS(KB2, LB, 0:MB)
*
      COMMON /FLOW/   EC(KU2, LU, 0:MU) ,  ES(KU2, LU, 0:MU)
     &             , DEC(KU2, LU, 0:MU) , DES(KU2, LU, 0:MU)
     &             ,  FC(KU4, LU, 0:MU) ,  FS(KU4, LU, 0:MU)
     &             , DFC(KU4, LU, 0:MU) , DFS(KU4, LU, 0:MU)
*
      COMMON /NRGB/ BINT1(LB, 0:MB) , BINT2(LB, 0:MB) , BWK(0:NKB)
     &           , BTT1(KB2,0:NKB) , BTT2(KB2,0:NKB) , BTT3(KB2,0:NKB)
     &           , BTT4(KB2,0:NKB) , BTT5(KB2,0:NKB)
*
      COMMON /NRGU/ UINT1(LU, 0:MU) , UINT2(LU, 0:MU) , UWK(0:NKU)
     &           , UTT1(KU4,0:NKU) , UTT2(KU4,0:NKU) , UTT3(KU4,0:NKU)
     &           , UTT4(KU4,0:NKU) , UTT5(KU4,0:NKU)
*
      DIMENSION TORB(LB,0:MB), POLB(LB,0:MB)
     &        , TORU(LU,0:MU), POLU(LU,0:MU)
      DIMENSION EM(0:MB), EMPHI(0:MB), EK(0:MU), EKPHI(0:MU)
******
!$omp parallel do private(M,LS,L,IK,K,XNRG1,XNRG2,XNRG3,
!$omp& WHC1,WHS1,WHC2,WHS2,WGC3,WGS3) schedule(dynamic)
      DO 40 M=0,MB
       EM(M)=0.D0
       EMPHI(M)=0.D0
       LS=MAX(M,1)
       DO 40 L=LS,LB
        XNRG1=0.D0
        XNRG2=0.D0
        XNRG3=0.D0
        DO 30 IK=0,NKB
         WHC1=0.D0
         WHS1=0.D0
         WHC2=0.D0
         WHS2=0.D0
         DO 10 K=1,KB2
          WHC1=WHC1  +  HC(K,L,M) * BTT1(K,IK)
          WHS1=WHS1  +  HS(K,L,M) * BTT1(K,IK)
          WHC2=WHC2  +  HC(K,L,M) * BTT2(K,IK)
          WHS2=WHS2  +  HS(K,L,M) * BTT2(K,IK)
10    CONTINUE
         WGC3=0.D0
         WGS3=0.D0
         DO 20 K=1,KB2
          WGC3=WGC3  +  GC(K,L,M) * BTT3(K,IK)
          WGS3=WGS3  +  GS(K,L,M) * BTT3(K,IK)
20    CONTINUE
        XNRG1=XNRG1  +  (WHC1**2 + WHS1**2) * BWK(IK)
        XNRG2=XNRG2  +  (WHC2**2 + WHS2**2) * BWK(IK)
        XNRG3=XNRG3  +  (WGC3**2 + WGS3**2) * BWK(IK)
30    CONTINUE
       XNRG1=XNRG1 * BINT1(L,M)
       XNRG2=XNRG2 * BINT2(L,M)
       XNRG3=XNRG3 * BINT2(L,M)
       POLB(L,M)=XNRG1 + XNRG2
       TORB(L,M)=XNRG3
       EM(M)=EM(M) + XNRG1 + XNRG2 + XNRG3
       EMPHI(M)=EMPHI(M) + XNRG3
40    CONTINUE
!$omp end parallel do
******
!$omp parallel do private(M,LS,L,IK,K,XNRG1,XNRG2,XNRG3,
!$omp& WFC1,WFS1,WFC2,WFS2,WEC3,WES3) schedule(dynamic)
      DO 80 M=0,MU
       EK(M)=0.D0
       EKPHI(M)=0.D0
       LS=MAX(M,1)
       DO 80 L=LS,LU
        XNRG1=0.D0
        XNRG2=0.D0
        XNRG3=0.D0
        DO 70 IK=0,NKU
         WFC1=0.D0
         WFS1=0.D0
         WFC2=0.D0
         WFS2=0.D0
         DO 50 K=1,KU4
          WFC1=WFC1  +  FC(K,L,M) * UTT1(K,IK)
          WFS1=WFS1  +  FS(K,L,M) * UTT1(K,IK)
          WFC2=WFC2  +  FC(K,L,M) * UTT2(K,IK)
          WFS2=WFS2  +  FS(K,L,M) * UTT2(K,IK)
50    CONTINUE
         WEC3=0.D0
         WES3=0.D0
         DO 60 K=1,KU2
          WEC3=WEC3  +  EC(K,L,M) * UTT3(K,IK)
          WES3=WES3  +  ES(K,L,M) * UTT3(K,IK)
60    CONTINUE
        XNRG1=XNRG1  +  (WFC1**2 + WFS1**2) * UWK(IK)
        XNRG2=XNRG2  +  (WFC2**2 + WFS2**2) * UWK(IK)
        XNRG3=XNRG3  +  (WEC3**2 + WES3**2) * UWK(IK)
70    CONTINUE
       XNRG1=XNRG1 * UINT1(L,M)
       XNRG2=XNRG2 * UINT2(L,M)
       XNRG3=XNRG3 * UINT2(L,M)
       POLU(L,M)=XNRG1 + XNRG2
       TORU(L,M)=XNRG3
       EK(M)=EK(M) + XNRG1 + XNRG2 + XNRG3
       EKPHI(M)=EKPHI(M) + XNRG3
80    CONTINUE
!$omp end parallel do
      RETURN
      END
************************************************************************
      SUBROUTINE DISS(TORCB,POLCB,TORCU,POLCU,DM,DMPHI,DK,DKPHI)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
*
      COMMON /FIELD/  GC(KB2, LB, 0:MB) ,  GS(KB2, LB, 0:MB)
     &             , DGC(KB2, LB, 0:MB) , DGS(KB2, LB, 0:MB)
     &             ,  HC(KB2, LB, 0:MB) ,  HS(KB2, LB, 0:MB)
     &             , DHC(KB2, LB, 0:MB) , DHS(KB2, LB, 0:MB)
*
      COMMON /FLOW/   EC(KU2, LU, 0:MU) ,  ES(KU2, LU, 0:MU)
     &             , DEC(KU2, LU, 0:MU) , DES(KU2, LU, 0:MU)
     &             ,  FC(KU4, LU, 0:MU) ,  FS(KU4, LU, 0:MU)
     &             , DFC(KU4, LU, 0:MU) , DFS(KU4, LU, 0:MU)
*
      COMMON /NRGB/ BINT1(LB, 0:MB) , BINT2(LB, 0:MB) , BWK(0:NKB)
     &           , BTT1(KB2,0:NKB) , BTT2(KB2,0:NKB) , BTT3(KB2,0:NKB)
     &           , BTT4(KB2,0:NKB) , BTT5(KB2,0:NKB)
*
      COMMON /NRGU/ UINT1(LU, 0:MU) , UINT2(LU, 0:MU) , UWK(0:NKU)
     &           , UTT1(KU4,0:NKU) , UTT2(KU4,0:NKU) , UTT3(KU4,0:NKU)
     &           , UTT4(KU4,0:NKU) , UTT5(KU4,0:NKU)
*
      DIMENSION TORCB(LB,0:MB), POLCB(LB,0:MB)
     &        , TORCU(LU,0:MU), POLCU(LU,0:MU)
      DIMENSION DM(0:MB), DMPHI(0:MB), DK(0:MU), DKPHI(0:MU)
******
!$omp parallel do private(M,LS,L,IK,K,XN1,XN2,XN3,
!$omp& WGC1,WGS1,WGC2,WGS2,WHC4,WHS4,WHC5,WHS5,WHC3,WHS3)
!$omp& schedule(dynamic)
      DO M=0,MB
       DM(M)=0.D0
       DMPHI(M)=0.D0
       LS=MAX(M,1)
       DO L=LS,LB
        XN1=0.D0
        XN2=0.D0
        XN3=0.D0
        DO IK=0,NKB
         WGC1=0.D0
         WGS1=0.D0
         WGC2=0.D0
         WGS2=0.D0
         WHC4=0.D0
         WHS4=0.D0
         WHC5=0.D0
         WHS5=0.D0
         DO K=1,KB2
          WGC1=WGC1  +  GC(K,L,M) * BTT1(K,IK)
          WGS1=WGS1  +  GS(K,L,M) * BTT1(K,IK)
          WGC2=WGC2  +  GC(K,L,M) * BTT2(K,IK)
          WGS2=WGS2  +  GS(K,L,M) * BTT2(K,IK)
          WHC4=WHC4  +  HC(K,L,M) * BTT4(K,IK)
          WHS4=WHS4  +  HS(K,L,M) * BTT4(K,IK)
          WHC5=WHC5  +  HC(K,L,M) * BTT5(K,IK)
          WHS5=WHS5  +  HS(K,L,M) * BTT5(K,IK)
         ENDDO
         WHC3=DFLOAT(L*(L+1))*WHC4 - WHC5
         WHS3=DFLOAT(L*(L+1))*WHS4 - WHS5
         XN1=XN1  +  (WGC1**2 + WGS1**2) * BWK(IK)
         XN2=XN2  +  (WGC2**2 + WGS2**2) * BWK(IK)
         XN3=XN3  +  (WHC3**2 + WHS3**2) * BWK(IK)
        ENDDO
       XN1=XN1 * BINT1(L,M)
       XN2=XN2 * BINT2(L,M)
       XN3=XN3 * BINT2(L,M)
       TORCB(L,M)=XN1 + XN2
       POLCB(L,M)=XN3
       DM(M)=DM(M) + XN1 + XN2 + XN3
       DMPHI(M)=DMPHI(M) + XN3
       ENDDO
      ENDDO
!$omp end parallel do
***
!$omp parallel do private(M,LS,L,IK,K,XN1,XN2,XN3,
!$omp& WEC1,WES1,WEC2,WES2,WFC4,WFS4,WFC5,WFS5,WFC3,WFS3)
!$omp& schedule(dynamic)
      DO M=0,MU
       DK(M)=0.D0
       DKPHI(M)=0.D0
       LS=MAX(M,1)
       DO L=LS,LU
        XN1=0.D0
        XN2=0.D0
        XN3=0.D0
        DO IK=0,NKU
         WEC1=0.D0
         WES1=0.D0
         WEC2=0.D0
         WES2=0.D0
         DO K=1,KU2
          WEC1=WEC1  +  EC(K,L,M) * UTT1(K,IK)
          WES1=WES1  +  ES(K,L,M) * UTT1(K,IK)
          WEC2=WEC2  +  EC(K,L,M) * UTT2(K,IK)
          WES2=WES2  +  ES(K,L,M) * UTT2(K,IK)
         ENDDO
         WFC4=0.D0
         WFS4=0.D0
         WFC5=0.D0
         WFS5=0.D0
         DO K=1,KU4
          WFC4=WFC4  +  FC(K,L,M) * UTT4(K,IK)
          WFS4=WFS4  +  FS(K,L,M) * UTT4(K,IK)
          WFC5=WFC5  +  FC(K,L,M) * UTT5(K,IK)
          WFS5=WFS5  +  FS(K,L,M) * UTT5(K,IK)
         ENDDO
         WFC3=DFLOAT(L*(L+1))*WFC4 - WFC5
         WFS3=DFLOAT(L*(L+1))*WFS4 - WFS5
         XN1=XN1  +  (WEC1**2 + WES1**2) * UWK(IK)
         XN2=XN2  +  (WEC2**2 + WES2**2) * UWK(IK)
         XN3=XN3  +  (WFC3**2 + WFS3**2) * UWK(IK)
        ENDDO
       XN1=XN1 * UINT1(L,M)
       XN2=XN2 * UINT2(L,M)
       XN3=XN3 * UINT2(L,M)
       TORCU(L,M)=XN1 + XN2
       POLCU(L,M)=XN3
       DK(M)=DK(M) + XN1 + XN2 + XN3
       DKPHI(M)=DKPHI(M) + XN3
       ENDDO
      ENDDO
!$omp end parallel do
      RETURN
      END
************************************************************************
      SUBROUTINE INIT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'      
*
      COMMON /NRGB/ BINT1(LB, 0:MB) , BINT2(LB, 0:MB) , BWK(0:NKB)
     &           , BTT1(KB2,0:NKB) , BTT2(KB2,0:NKB) , BTT3(KB2,0:NKB)
     &           , BTT4(KB2,0:NKB) , BTT5(KB2,0:NKB)
*
      COMMON /NRGU/ UINT1(LU, 0:MU) , UINT2(LU, 0:MU) , UWK(0:NKU)
     &           , UTT1(KU4,0:NKU) , UTT2(KU4,0:NKU) , UTT3(KU4,0:NKU)
     &           , UTT4(KU4,0:NKU) , UTT5(KU4,0:NKU)
******
      PI=DACOS(-1.D0)
*
      DO 20 M=0,MB
       LS=MAX(M,1)
       DO 20 L=LS,LB
        FAC=2.D0/DFLOAT(2*L+1)
        DO 10 N=L-M+1,L+M
         FAC=FAC * DFLOAT(N)
10    CONTINUE
        BINT1(L,M)=FAC * DFLOAT(L*(L+1))**2
        BINT2(L,M)=FAC * DFLOAT(L*(L+1))
        BINT1(L,M)=BINT1(L,M) * PI/2.D0
        BINT2(L,M)=BINT2(L,M) * PI/2.D0
        IF ( M .EQ. 0 ) BINT1(L,M)=BINT1(L,M) * 2.D0
        IF ( M .EQ. 0 ) BINT2(L,M)=BINT2(L,M) * 2.D0
20    CONTINUE
******
      DR=(R2 - R1)/DFLOAT(NKB)
      DO 30 IK=0,NKB
       IF ( MOD(IK,2)  .EQ. 1 ) BWK(IK)=4.D0/3.D0 * DR
       IF ( MOD(IK,2)  .EQ. 0 ) BWK(IK)=2.D0/3.D0 * DR
       IF ( IK*(IK-NKB) .EQ. 0 ) BWK(IK)=1.D0/3.D0 * DR
       XX=-1.D0 + DFLOAT(2*IK)/DFLOAT(NKB)
       R=(R2 + R1)/2.D0  +  (R2 - R1)/2.D0*XX
       DO 30 K=1,KB2
        BTT1(K,IK)=TTT(0,K-1,XX) / R
        BTT2(K,IK)=TTT(1,K-1,XX) * DXDR
        BTT3(K,IK)=TTT(0,K-1,XX)
        BTT4(K,IK)=TTT(0,K-1,XX) / R**2
        BTT5(K,IK)=TTT(2,K-1,XX) * DXDR**2
30    CONTINUE
******
      DO 50 M=0,MU
       LS=MAX(M,1)
       DO 50 L=LS,LU
        FAC=2.D0/DFLOAT(2*L+1)
        DO 40 N=L-M+1,L+M
         FAC=FAC * DFLOAT(N)
40    CONTINUE
        UINT1(L,M)=FAC * DFLOAT(L*(L+1))**2
        UINT2(L,M)=FAC * DFLOAT(L*(L+1))
        UINT1(L,M)=UINT1(L,M) * PI/2.D0
        UINT2(L,M)=UINT2(L,M) * PI/2.D0
        IF ( M .EQ. 0 ) UINT1(L,M)=UINT1(L,M) * 2.D0
        IF ( M .EQ. 0 ) UINT2(L,M)=UINT2(L,M) * 2.D0
50    CONTINUE
******
      DR=(R2 - R1)/DFLOAT(NKU)
      DO 60 IK=0,NKU
       IF ( MOD(IK,2)  .EQ. 1 ) UWK(IK)=4.D0/3.D0 * DR
       IF ( MOD(IK,2)  .EQ. 0 ) UWK(IK)=2.D0/3.D0 * DR
       IF ( IK*(IK-NKU) .EQ. 0 ) UWK(IK)=1.D0/3.D0 * DR
       XX=-1.D0 + DFLOAT(2*IK)/DFLOAT(NKU)
       R=(R2 + R1)/2.D0  +  (R2 - R1)/2.D0*XX
       DO 60 K=1,KU4
        UTT1(K,IK)=TTT(0,K-1,XX) / R
        UTT2(K,IK)=TTT(1,K-1,XX) * DXDR
        UTT3(K,IK)=TTT(0,K-1,XX)
        UTT4(K,IK)=TTT(0,K-1,XX) / R**2
        UTT5(K,IK)=TTT(2,K-1,XX) * DXDR**2
60    CONTINUE
      RETURN
      END
************************************************************************
*     TTT(K,M,X) = the K-th derivative of Tm(X),
*     the M-th Chebyshev polynomial evaLBated at X.
*
      DOUBLE PRECISION FUNCTION TTT(K,M,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(0:1000), B(0:1000)
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
      DOUBLE PRECISION A(0:1000), C(0:1000)
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
