      program main
      implicit none
      integer KB, KB2, LB, MB, KU, KU2, KU4, LU, MU, KT, KT2, LT, MT
      integer KBn, KB2n, LBn, MBn, KUn, KU2n, KU4n, LUn, MUn, KTn, KT2n, LTn, MTn
      integer K, L, M

      PARAMETER (KB=50 , LB=50 , MB=50)
      PARAMETER (KB2=KB+2)
      PARAMETER (KU=50 , LU=50 , MU=50)
      PARAMETER (KU2=KU+2 , KU4=KU+4)
      PARAMETER (KT=50 , LT=50 , MT=50)
      PARAMETER (KT2=KT+2)
      
      PARAMETER (KBn=100 , LBn=100 , MBn=100)
      PARAMETER (KB2n=KBn+2)
      PARAMETER (KUn=100 , LUn=100 , MUn=100)
      PARAMETER (KU2n=KUn+2 , KU4n=KUn+4)
      PARAMETER (KTn=100 , LTn=100 , MTn=100)
      PARAMETER (KT2n=KTn+2)
      
      double precision  EC(KU2, LU, 0:MU) ,  ES(KU2, LU, 0:MU) &
     &               ,  FC(KU4, LU, 0:MU) ,  FS(KU4, LU, 0:MU) 
      double precision  GC(KB2, LB, 0:MB) ,  GS(KB2, LB, 0:MB) &
     &               ,  HC(KB2, LB, 0:MB) ,  HS(KB2, LB, 0:MB) 
      double precision  TC(KT2,0:LT,0:MT) ,  TS(KT2,0:LT,0:MT)
      
      double precision  ECn(KU2n, LUn, 0:MUn) ,  ESn(KU2n, LUn, 0:MUn) &
     &               ,  FCn(KU4n, LUn, 0:MUn) ,  FSn(KU4n, LUn, 0:MUn)
      double precision  GCn(KB2n, LBn, 0:MBn) ,  GSn(KB2n, LBn, 0:MBn) &
     &               ,  HCn(KB2n, LBn, 0:MBn) ,  HSn(KB2n, LBn, 0:MBn)     
      double precision  TCn(KT2n,0:LTn,0:MTn) ,  TSn(KT2n,0:LTn,0:MTn)
     
      do M=0, MUn
       do L=1, LUn
        do K=1, KU2n
         ECn(K,L,M)=0.d0
         ESn(K,L,M)=0.d0
        enddo
        do K=1, KU4n
         FCn(K,L,M)=0.d0
         FSn(K,L,M)=0.d0
        enddo
       enddo
      enddo      
      do M=0, MBn
       do L=1, LBn
        do K=1, KB2n
         GCn(K,L,M)=0.d0
         GSn(K,L,M)=0.d0
         HCn(K,L,M)=0.d0
         HSn(K,L,M)=0.d0
        enddo
       enddo
      enddo       
      do M=0, MTn
       do L=1, LTn
        do K=1, KT2n
         TCn(K,L,M)=0.d0
         TSn(K,L,M)=0.d0
        enddo
       enddo
      enddo 
      
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
      
      do M=0, MU
       do L=1, LU
        do K=1, KU2
         ECn(K,L,M)=EC(K,L,M)
         ESn(K,L,M)=ES(K,L,M)
        enddo
        do K=1, KU4
         FCn(K,L,M)=FC(K,L,M)
         FSn(K,L,M)=FS(K,L,M)
        enddo
       enddo
      enddo      
      do M=0, MB
       do L=1, LB
        do K=1, KB2
         GCn(K,L,M)=GC(K,L,M)
         GSn(K,L,M)=GS(K,L,M)
         HCn(K,L,M)=HC(K,L,M)
         HSn(K,L,M)=HS(K,L,M)
        enddo
       enddo
      enddo       
      do M=0, MT
       do L=1, LT
        do K=1, KT2
         TCn(K,L,M)=TC(K,L,M)
         TSn(K,L,M)=TS(K,L,M)
        enddo
       enddo
      enddo 
      
      OPEN(3,FORM='unformatted')
      write(3) ECn
      write(3) ESn
      write(3) FCn
      write(3) FSn
      CLOSE(3)
      OPEN(4,FORM='unformatted')
      write(4) GCn
      write(4) GSn
      write(4) HCn
      write(4) HSn
      CLOSE(4)
      OPEN(7,FORM='unformatted')
      write(7) TCn
      write(7) TSn
      CLOSE(7)
        
      end program main
