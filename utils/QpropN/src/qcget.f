
      SUBROUTINE QCGET(RHO,RMU,VSO)
      CHARACTER*80 FNAME,LINE
      LOGICAL ERROR
C
      LU = 11
      FNAME = 'qcon.def'
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=90)
      ILINE = 0
C
C---- extract parameters on data lines
      NVAL = 1
      CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
      IF(IERR.EQ.+1) GO TO 900
      IF(IERR.EQ.-1) GO TO 950
      IF(NVAL.LT. 1) GO TO 980
      RHO = RVAL
C
C---- extract parameters on data lines
      NVAL = 1
      CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
      IF(IERR.EQ.+1) GO TO 900
      IF(IERR.EQ.-1) GO TO 950
      IF(NVAL.LT. 1) GO TO 980
      RMU = RVAL
C
C---- extract parameters on data lines
      NVAL = 1
      CALL RREAD(LU,LINE,ILINE,IERR,NVAL,RVAL)
      IF(IERR.EQ.+1) GO TO 900
      IF(IERR.EQ.-1) GO TO 950
      IF(NVAL.LT. 1) GO TO 980
      VSO = RVAL
C
      CLOSE(LU)
      RETURN
C
 90   CONTINUE
      RETURN
C
C--------------------------------------------
 900  CONTINUE
      WRITE(*,9000) FNAME(1:64), ILINE, LINE
 9000 FORMAT(/' Read error'
     &       /'   in file:  ', A
     &       /'   line',I3,':  ', A)
      STOP
C
 950  CONTINUE
      WRITE(*,9500) FNAME(1:64), ILINE
 9500 FORMAT(/' Unexpected end-of-file reached'
     &       /'   in file:  ', A
     &       /'   line',I3)
      STOP
C
C
 980  CONTINUE
      WRITE(*,9500) FNAME(1:64), ILINE
 9800 FORMAT(/' Fewer parameters than required'
     &       /'   in file:  ', A
     &       /'   line',I3)
      STOP
C
      END ! QCGET
