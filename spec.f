      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      INCLUDE 'params.f'
*
      DIMENSION GC(KB2, LB, 0:MB) , GS(KB2, LB, 0:MB)
      DIMENSION HC(KB2, LB, 0:MB) , HS(KB2, LB, 0:MB)
      DIMENSION EC(KU2, LU, 0:MU) , ES(KU2, LU, 0:MU)
      DIMENSION FC(KU4, LU, 0:MU) , FS(KU4, LU, 0:MU)
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
      DO K=1,KU2
       E=0.D0
       DO L=1,LU
        DO M=0,MU
         E=E+EC(K,L,M)**2+ES(K,L,M)**2
        ENDDO
       ENDDO
       WRITE(31,10) K, DLOG10(E)
      ENDDO
*
      DO K=1,KU4
       F=0.D0
       DO L=1,LU
        DO M=0,MU
         F=F+FC(K,L,M)**2+FS(K,L,M)**2
        ENDDO
       ENDDO
       WRITE(32,10) K, DLOG10(F)
      ENDDO      
*
      DO L=1,LU
       E=0.D0
       F=0.D0
       DO M=0,MU
        DO K=1,KU2
         E=E+EC(K,L,M)**2+ES(K,L,M)**2
        ENDDO
        DO K=1,KU4
         F=F+FC(K,L,M)**2+FS(K,L,M)**2
        ENDDO
       ENDDO
       WRITE(33,20) L, DLOG10(E), DLOG10(F)
      ENDDO
*
      DO M=0,MU
       E=0.D0
       F=0.D0
       DO L=1,LU
        DO K=1,KU2
         E=E+EC(K,L,M)**2+ES(K,L,M)**2
        ENDDO
        DO K=1,KU4
         F=F+FC(K,L,M)**2+FS(K,L,M)**2
        ENDDO
       ENDDO
       WRITE(34,20) M, DLOG10(E), DLOG10(F)
      ENDDO
*
      DO K=1,KB2
       G=0.D0
       H=0.D0
       DO L=1,LB
        DO M=0,MB
         G=G+GC(K,L,M)**2+GS(K,L,M)**2
         H=H+HC(K,L,M)**2+HS(K,L,M)**2
        ENDDO
       ENDDO
       WRITE(35,20) K, DLOG10(G), DLOG10(H)
      ENDDO
*
      DO L=1,LB
       G=0.D0
       H=0.D0
       DO K=1,KB2
        DO M=0,MB
         G=G+GC(K,L,M)**2+GS(K,L,M)**2
         H=H+HC(K,L,M)**2+HS(K,L,M)**2
        ENDDO
       ENDDO
       WRITE(36,20) L, DLOG10(G), DLOG10(H)
      ENDDO
*
      DO M=0,MB
       G=0.D0
       H=0.D0
       DO K=1,KB2
        DO L=1,LB
         G=G+GC(K,L,M)**2+GS(K,L,M)**2
         H=H+HC(K,L,M)**2+HS(K,L,M)**2
        ENDDO
       ENDDO
       WRITE(37,20) M, DLOG10(G), DLOG10(H)
      ENDDO
*
      CLOSE(31)
      CLOSE(32)
      CLOSE(33)
      CLOSE(34)
      CLOSE(35)
      CLOSE(36)
      CLOSE(37)
*
10    FORMAT(I5, E16.8)
20    FORMAT(I5, 2E16.8)
      STOP
      END
