c$$$  type mytype
c$$$  integer:: i
c$$$  real*8 :: a(3)
c$$$  end type mytype
c$$$  To create a variable of type mytype, use
c$$$  type (mytype) var
c$$$  An array of mytype can also be created.
c$$$  type (mytype) stuff(3)
c$$$  Elements of derived types are accessed with the "%" operator. For instance,
c$$$  var%i = 3
c$$$  var%a(1) = 4.0d0
c$$$  stuff(1)%a(2) = 8.0d0
C
C     
      SUBROUTINE INITPOLAR(POLARFILE, POLAR, IERR)
      USE POLARDATAMODULE
C----------------------------------------------------------------
C     Reads polar data in POLARFILE to a POLAR structure
C     defined type
C----------------------------------------------------------------
C
      TYPE (POLARDATA) POLAR
      
      CHARACTER*128 LINE
      CHARACTER*80 POLARFILE

      REAL RVAL(3)

      REAL w(101), wrk(2500)
      REAL fp, s, xb, xe
      INTEGER i,ier,iopt,k,lwrk,m,n,nest
      
      LU=999
      OPEN(LU, FILE=POLARFILE, STATUS='OLD', ERR=90)

      ILINE = 0
C
C---- Reynolds Number
      NVAL = 1
      CALL RREAD(LU, LINE, ILINE, IERR, NVAL, RVAL)
      IF(IERR.NE.0) THEN
         RETURN
      ELSE
         POLAR%RE = RVAL(1)
      ENDIF
C
C---- Stall AoA
      NVAL = 1
      CALL RREAD(LU, LINE, ILINE, IERR, NVAL, RVAL)
      IF(IERR.NE.0) THEN
         RETURN
      ELSE
         POLAR%ASTAR = RVAL(1)
      ENDIF
C
C---- Aerodynamic coefficients
      DO NA = 1, 101
         NVAL = 3
         CALL RREAD(LU, LINE, ILINE, IERR, NVAL, RVAL)
         IF(IERR.NE.0) THEN
            RETURN
         ELSE
            POLAR%ALPHA(NA) = RVAL(1)
            POLAR%CL(NA) = RVAL(2)
            POLAR%CD(NA) = RVAL(3)
         ENDIF
      ENDDO
      
      CLOSE(LU)
C
C     Setting up interpolation routines
C      
      iopt = 1
      s = 0.
      m = 101
      k = 3
      
      xb=POLAR%ALPHA(1)
      xe=POLAR%ALPHA(101)

      nest = 105
      lwrk = 2500

      DO i=1,m
         w(i) = 1.0
      ENDDO
C
C     Getting spline interpolation coefficients
C      
      CALL curfit(iopt,m,POLAR%ALPHA,POLAR%CL,w,xb,xe,k,s,nest,
     *     POLAR%FIT_CL_N,POLAR%FIT_CL_T,POLAR%FIT_CL_C,fp,
     *     wrk,lwrk,iwrk,ier)

      CALL curfit(iopt,m,POLAR%ALPHA,POLAR%CD,w,xb,xe,k,s,nest,
     *     POLAR%FIT_CD_N,POLAR%FIT_CD_T,POLAR%FIT_CD_C,fp,
     *     wrk,lwrk,iwrk,ier)
      
      RETURN

 90   CONTINUE
      WRITE(*,*)
      WRITE(*,*) 'Problems opening polar file: ', POLARFILE
      IERR = 1
      RETURN
      
      END
      
      SUBROUTINE GETPOLARCLINFO(POLAR,A,CL,DCLDA,CL0,STALL)
      USE POLARDATAMODULE
C----------------------------------------------------------------
C     Returns CL-related data from polar given AoA angle
C     
C----------------------------------------------------------------
C     
      REAL wrk(2500)
      INTEGER k, ier, nu
      REAL ADEG, CL, DCLDA, CL0, STALL

      REAL RAD2DEG
      TYPE (POLARDATA) POLAR

      RAD2DEG = 180./(2*3.14159265)
      k = 3
      nu = 1
      
      ADEG = A*RAD2DEG

C      WRITE(*,*) 'Fit N: ', POLAR%FIT_CL_N
      
      CALL splev(POLAR%FIT_CL_T,POLAR%FIT_CL_N,POLAR%FIT_CL_C,
     *     k,ADEG,CL,1,ier)
C      WRITE(*,*) 'CL IER: ', ier

C      WRITE(*,200) ADEG,CL
      
      CALL splder(POLAR%FIT_CL_T,POLAR%FIT_CL_N,POLAR%FIT_CL_C,
     *     k,nu,ADEG,DCLDA,1,wrk,ier)
      
C      WRITE(*,*) 'DCLDA IER: ', ier
C      WRITE(*,210) ADEG,DCLDA

      CALL splev(POLAR%FIT_CL_T,POLAR%FIT_CL_N,POLAR%FIT_CL_C,
     *     k,0.,CL0,1,ier)

C      WRITE(*,200) 0.,CL0

      IF(ADEG .GT. POLAR%ASTAR) THEN
         STALL=1.
      ELSE
         STALL=0.
      ENDIF
      
      RETURN

 200  FORMAT('CL for alpha ',f10.5,' is ',f7.3)
 210  FORMAT('DCLDA for alpha ',f10.5,' is ',f7.3)
 220  FORMAT('RAD2DEG ', f10.5)
      END

      SUBROUTINE GETPOLARCDINFO(POLAR,A,CD,CD_CL)
      USE POLARDATAMODULE
C----------------------------------------------------------------
C     Returns CD-related data from polar given AoA angle
C     
C----------------------------------------------------------------
C
      REAL wrk(2500)
      REAL dcdda, dclda
      REAL ADEG, RAD2DEG
      INTEGER k, ier, nu
      TYPE (POLARDATA) POLAR
      
      RAD2DEG = 180./(2*3.14159265)
      k = 3
      nu = 1

      ADEG = A*RAD2DEG
      
      CALL splev(POLAR%FIT_CD_T,POLAR%FIT_CD_N,POLAR%FIT_CD_C,
     *     k,ADEG,CD,1,ier)
      
      CALL splder(POLAR%FIT_CD_T,POLAR%FIT_CD_N,POLAR%FIT_CD_C,
     *     k,nu,ADEG,dcdda,1,wrk,ier)

      CALL splder(POLAR%FIT_CL_T,POLAR%FIT_CL_N,POLAR%FIT_CL_C,
     *     k,order,ADEG,dclda,1,wrk,ier)

      CD_CL = dcdda/dclda
      
      RETURN
      END
