      PROGRAM test

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

C     POLARDATA MODULE
      REAL :: RE
      REAL :: ALPHA(5)
      REAL :: CL(5)
      REAL :: CD(5)
      
C     Variaveis Locais
      REAL :: FA(5)
      INTEGER :: pos

C     INPUT  > POLAR,A
      RE = 10000
      ALPHA = [1,2,3,4,5]
      CL = [5,6,7,8,9]
      CD = [8,9,10,11,12]

      A = 2
C     OUTPUT > CD,CD_CL

C     Returns CD-related data from polar given AoA

C     CD_CL = dCD/dCL = dCD/dA * dA/DCL = dCD/dA * (dCL/dA)^-1

C     Finds the closest value of AoA from the polar daata
      pos = 1
      val = 1E3
      
      DO i = 1, 5
        FA(i) = abs(ALPHA(i) - A)
        IF(FA(i).LT.val) THEN
           val = FA(i)
           pos = i
         END IF
      END DO

C     If point is exact, the CL is given, otherwise, interpolated #IMPLEMENT
      CDX = CD(pos)

C     Finds the derivatives using the two surrounding points #IMPROVE
      DCL = (CL(pos+1)-CL(pos-1)) / (ALPHA(pos+i)-ALPHA(pos-i))
      DCD = (CD(pos+1)-CD(pos-1)) / (ALPHA(pos+i)-ALPHA(pos-i))
      
      CD_CL = DCD/DCL

      PRINT *,CDX, CD_CL
C      PAUSE
      
      END PROGRAM test
