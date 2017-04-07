      PROGRAM TST_INITPOL
      USE POLARDATAMODULE
      
      TYPE (POLARDATA) POLAR
      INTEGER IRR
      CHARACTER*80 POLARFILE

      POLARFILE = '../runs/be50_pol.dat'
      
      CALL INITPOLAR(POLARFILE, POLAR, IERR)
      IF(IERR.EQ.+1) GO TO 900
      IF(IERR.EQ.-1) GO TO 950
      
      WRITE(*,*) 'Read polar file: ', POLARFILE
      WRITE(*,*) 'Alpha - Cl - Cd'
      DO NA = 1, 101
         WRITE(*,*) POLAR%ALPHA(NA), POLAR%CL(NA), POLAR%CD(NA)
      ENDDO

      STOP

C     
 900  CONTINUE
      WRITE(*,9000) POLARFILE
 9000 FORMAT(/' Read error'
     &       /'   in file:  ', A)
      STOP
C
 950  CONTINUE
      WRITE(*,9500) POLARFILE
 9500 FORMAT(/' Unexpected end-of-file reached'
     &       /'   in file:  ', A)
      STOP

      END !TST_INITPOL
