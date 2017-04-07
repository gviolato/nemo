      PROGRAM TST_READPN

      IMPLICIT REAL (A-H,M,O-Z)
      
      CHARACTER*128 LINE
      CHARACTER*80 FILNAM, PNAME
      CHARACTER*80 ARGP1
C      
C---- Defined polars (from propfile)
      CHARACTER*80 POLFILES(10)

C---- get Unix command-line arguments, if any
      CALL GETARG0(1,ARGP1)
      
C---- read prop data file
      FILNAM = ARGP1
      IF(FILNAM.EQ.' ') GO TO 18
C
      LU = 1
      OPEN(LU,FILE=FILNAM,STATUS='OLD',ERR=18)
C
      ILINE = 0
C
C---- prop name
      CALL FREAD(LU,LINE,ILINE,IERR,PNAME)
      IF(IERR.EQ.+1) GO TO 900
      IF(IERR.EQ.-1) GO TO 950
C
C---- extract parameters on data lines
      NVAL = 2
      CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
      IF(IERR.EQ.+1) GO TO 900
      IF(IERR.EQ.-1) GO TO 950
      IF(NVAL.LT. 1) GO TO 980
C
      NPOLS = 0
      CALL READPOLNAMES(LU,LINE,ILINE,IERR,NPOLS,POLFILES)
      IF(IERR.EQ.+1) GO TO 900
      IF(IERR.EQ.-1) GO TO 950
      IF(NPOLS.LT.1) GO TO 980
      DO IP = 1, NPOLS
         WRITE(*,*) 'Polar file ', IP, ': ', POLFILES(IP)
      ENDDO

      CLOSE(LU)
      STOP
C
C
 18   CONTINUE
      WRITE(*,*)
      WRITE(*,*) 'Prop file not found:  ', FILNAM(1:48)
      STOP
C     
 900  CONTINUE
      WRITE(*,9000) FILNAM(1:64), ILINE, LINE
 9000 FORMAT(/' Read error'
     &       /'   in file:  ', A
     &       /'   line',I3,':  ', A)
      STOP
C
 950  CONTINUE
      WRITE(*,9500) FILNAM(1:64), ILINE
 9500 FORMAT(/' Unexpected end-of-file reached'
     &       /'   in file:  ', A
     &       /'   line',I3)
      STOP
C
C
 980  CONTINUE
      WRITE(*,9500) FILNAM(1:64), ILINE
 9800 FORMAT(/' Fewer parameters than required'
     &       /'   in file:  ', A
     &       /'   line',I3)
      STOP
C
      END ! TST_READPN      
