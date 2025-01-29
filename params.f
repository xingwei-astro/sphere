*********************************************************************
* 2KN >= MAX(3KU+1, KU+2KB+1)
* 2LN >= MAX(3LU+1, LU+2LB+1)
* 2MN >= MAX(3MU+1, MU+2MB+1)
* LL   = MAX(LB, LU)
* MM   = MAX(MB, MU)
*********************************************************************
      PARAMETER (R1=0.1D0 , R2=1.0D0)
      PARAMETER (DXDR=2.D0/(R2-R1))
      PARAMETER (aspect=R2/(R2-R1))
*
      PARAMETER (KB=50 , LB=50 , MB=50)
      PARAMETER (KB2=KB+2)
*
      PARAMETER (KU=50 , LU=50 , MU=50)
      PARAMETER (KU2=KU+2 , KU4=KU+4)
*
      PARAMETER (KT=50 , LT=50 , MT=50)
      PARAMETER (KT2=KT+2)
*
      PARAMETER (KN=77 , LN=77 , MN=77)
      PARAMETER (LL=50 , MM=50)
*
      PARAMETER (NKU=8*KU, NKB=8*KB)
*
      PARAMETER (NT=100000, DT=1.D-1)
*
      PARAMETER (EKM=1.D-4, RA=1.0D7, PR=1.D0, PM=1.D0)
* AA measures convection       
      PARAMETER (AA=RA*EKM**2/Pr*aspect**3)
* collision      
      PARAMETER (CC=1.D-2, tau=10.d0) 
* precession      
      PARAMETER (PO=-1.0D-1, ALPHA=23.5D0)     
* tide      
      PARAMETER (xi=1.d-6, omega=0.1d0)
