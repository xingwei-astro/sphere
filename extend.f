      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (KB=50 , LB=50 , MB=20)
      PARAMETER (KB2=KB+2)
*
      PARAMETER (KU=50 , LU=50 , MU=20)
      PARAMETER (KU2=KU+2 , KU4=KU+4)
*
      PARAMETER (KT=50 , LT=50 , MT=20)
      PARAMETER (KT2=KT+2)
*
      PARAMETER (KBNEW=70 , LBNEW=70 , MBNEW=30)
      PARAMETER (KB2NEW=KBNEW+2)
*
      PARAMETER (KUNEW=70 , LUNEW=70 , MUNEW=30)
      PARAMETER (KU2NEW=KUNEW+2 , KU4NEW=KUNEW+4)
*
      PARAMETER (KTNEW=70 , LTNEW=70 , MTNEW=30)
      PARAMETER (KT2NEW=KTNEW+2)
*
      DIMENSION GC(KB2, LB, 0:MB) ,  GS(KB2, LB, 0:MB)
     &        , HC(KB2, LB, 0:MB) ,  HS(KB2, LB, 0:MB)
*
      DIMENSION EC(KU2, LU, 0:MU) ,  ES(KU2, LU, 0:MU)
     &        , FC(KU4, LU, 0:MU) ,  FS(KU4, LU, 0:MU)
*
      DIMENSION TC(KT2,0:LT,0:MT) ,  TS(KT2,0:LT,0:MT)
*
      DIMENSION GCNEW(KB2NEW, LBNEW, 0:MBNEW)
     &        , GSNEW(KB2NEW, LBNEW, 0:MBNEW)
     &        , HCNEW(KB2NEW, LBNEW, 0:MBNEW) 
     &        , HSNEW(KB2NEW, LBNEW, 0:MBNEW)
*
      DIMENSION ECNEW(KU2NEW, LUNEW, 0:MUNEW)
     &        , ESNEW(KU2NEW, LUNEW, 0:MUNEW)
     &        , FCNEW(KU4NEW, LUNEW, 0:MUNEW)
     &        , FSNEW(KU4NEW, LUNEW, 0:MUNEW)
*
      DIMENSION TCNEW(KT2NEW,0:LTNEW,0:MTNEW) 
     &        , TSNEW(KT2NEW,0:LTNEW,0:MTNEW)
******
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
      READ(7) TC
      READ(7) TS
      CLOSE(7)
*  
      DO M=0,MBNEW
       DO L=1,LBNEW
        DO K=1,KB2NEW
         GCNEW(K,L,M)=0.D0
         GSNEW(K,L,M)=0.D0
         HCNEW(K,L,M)=0.D0
         HSNEW(K,L,M)=0.D0
        ENDDO
       ENDDO
      ENDDO
*
      DO M=0,MIN(MB,MBNEW)
       DO L=1,MIN(LB,LBNEW)
        DO K=1,MIN(KB2,KB2NEW)
         GCNEW(K,L,M)=GC(K,L,M)
         GSNEW(K,L,M)=GS(K,L,M)
         HCNEW(K,L,M)=HC(K,L,M)
         HSNEW(K,L,M)=HS(K,L,M)
        ENDDO
       ENDDO
      ENDDO
*
      DO M=0,MUNEW
       DO L=1,LUNEW
        DO K=1,KU2NEW
         ECNEW(K,L,M)=0.D0
         ESNEW(K,L,M)=0.D0
        ENDDO
        DO K=1,KU4NEW
         FCNEW(K,L,M)=0.D0
         FSNEW(K,L,M)=0.D0
        ENDDO
       ENDDO
      ENDDO
*
      DO M=0,MIN(MU,MUNEW)
       DO L=1,MIN(LU,LUNEW)
        DO K=1,MIN(KU2,KU2NEW)
         ECNEW(K,L,M)=EC(K,L,M)
         ESNEW(K,L,M)=ES(K,L,M)
        ENDDO
        DO K=1,MIN(KU4,KU4NEW)
         FCNEW(K,L,M)=FC(K,L,M)
         FSNEW(K,L,M)=FS(K,L,M)
        ENDDO
       ENDDO
      ENDDO
*
      DO M=0,MTNEW
       DO L=0,LTNEW
        DO K=1,KT2NEW
         TCNEW(K,L,M)=0.D0
         TSNEW(K,L,M)=0.D0
        ENDDO
       ENDDO
      ENDDO
*
      DO M=0,MIN(MT,MTNEW)
       DO L=0,MIN(LT,LTNEW)
        DO K=1,MIN(KT2,KT2NEW)
         TCNEW(K,L,M)=TC(K,L,M)
         TSNEW(K,L,M)=TS(K,L,M)
        ENDDO
       ENDDO
      ENDDO
*
      WRITE(3) ECNEW
      WRITE(3) ESNEW
      WRITE(3) FCNEW
      WRITE(3) FSNEW
      CLOSE(3)
      WRITE(4) GCNEW
      WRITE(4) GSNEW
      WRITE(4) HCNEW
      WRITE(4) HSNEW
      CLOSE(4)
      WRITE(7) TCNEW
      WRITE(7) TSNEW
      CLOSE(7)
******
      STOP
      END
