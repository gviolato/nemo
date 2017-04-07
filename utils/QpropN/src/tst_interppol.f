      PROGRAM TST_INTERPPOL
      USE POLARDATAMODULE
      
      TYPE (POLARDATA) POLAR
      INTEGER IRR
      CHARACTER*80 POLARFILE
      REAL w(101), tl(105), cl(105), wrk(2500)
      REAL td(105), cd(105)
      REAL as(181), cli(181), cdi(181)
      REAL fp, s, xb, xe
      INTEGER i,ier,iopt,k,lwrk,m,nl,nd,nest,nk1

      POLARFILE = '../runs/be50_pol.dat'
      
      CALL INITPOLAR(POLARFILE, POLAR, IERR)
      IF(IERR.EQ.+1) GO TO 900
      IF(IERR.EQ.-1) GO TO 950
      
      WRITE(*,*) 'Read polar file: ', POLARFILE

      iopt = 1
      s = 0.
      m = 101
      k = 3
      
      xb=POLAR%ALPHA(1)
      xe=POLAR%ALPHA(101)

      nest = 105
      lwrk = 2500

      DO i=1,101
         w(i) = 1.0
      ENDDO
      
      CALL curfit(iopt,m,POLAR%ALPHA,POLAR%CL,w,xb,xe,k,s,nest,nl,tl,cl,
     *     fp,wrk,lwrk,iwrk,ier)

      CALL curfit(iopt,m,POLAR%ALPHA,POLAR%CD,w,xb,xe,k,s,nest,nd,td,cd,
     *     fp,wrk,lwrk,iwrk,ier)
      
      DO i=1,181
         as(i) = -45. + (i-1)*0.5
      ENDDO      
      
      CALL splev(tl,nl,cl,k,as,cli,181,ier)
      CALL splev(td,nd,cd,k,as,cdi,181,ier)
C      CALL splder(t,n,c,k,1,as,dcli,181,wrk,ier)

      WRITE(*,*) 'Interpolated values:'
      WRITE(*,*) 'alpha, cl, cd'
      DO i=1,181
         WRITE(*,*) as(i), cli(i), cdi(i)
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

      END !TST_INTERPPOL
