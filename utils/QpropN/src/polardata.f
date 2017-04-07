      MODULE POLARDATAMODULE
      
      TYPE POLARDATA
      REAL :: RE;
      REAL :: ASTAR;
      REAL :: ALPHA(101);
      REAL :: CL(101);
      REAL :: CD(101);
C
      INTEGER :: FIT_CL_N;
      REAL :: FIT_CL_T(105);
      REAL :: FIT_CL_C(105);
C
      INTEGER :: FIT_CD_N;
      REAL :: FIT_CD_T(105);
      REAL :: FIT_CD_C(105);

      END TYPE POLARDATA

      END MODULE
