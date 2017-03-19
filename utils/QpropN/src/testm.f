
C---- SUBROUTINE MOTORQ derivative tester 

      REAL PARMOT(10)

      EPS = 1.0E-5
C
C---- set typical input parameters 
      OMEGA = 1000.0
      VOLT  = 20.0
C
      PARMOT(1) = 0.5
      PARMOT(2) = 1.5
      PARMOT(3) = 700.0
      NMPAR = 3


      CALL MOTORQ(OMEGA,VOLT, PARMOT,NMPAR, 
     &            QM,QM_OMEGA,QM_VOLT)

      write(*,*) 'Omega ', OMEGA
      write(*,*) 'Volt  ', VOLT
      write(*,*) 'Qmotor', QM

 2100 FORMAT(1X,8F16.10)


      DOMEGA = MAX( OMEGA*EPS , EPS )
      DVOLT  = MAX( VOLT *EPS , EPS )


      OMEGA$ = OMEGA + DOMEGA
      VOLT$  = VOLT
      CALL MOTORQ(OMEGA$,VOLT$, PARMOT,NMPAR, 
     &            QM$,QM_OMEGA$,QM_VOLT$)
      WRITE(*,*) 
      WRITE(*,*) 'OMEGA'
      WRITE(*,2100) (QM$-QM)/DOMEGA
      WRITE(*,2100) (QM_OMEGA$+QM_OMEGA)*0.5
      WRITE(*,*) 


      OMEGA$ = OMEGA
      VOLT$  = VOLT + DVOLT
      CALL MOTORQ(OMEGA$,VOLT$, PARMOT,NMPAR, 
     &            QM$,QM_OMEGA$,QM_VOLT$)
      WRITE(*,*) 
      WRITE(*,*) 'VOLT'
      WRITE(*,2100) (QM$-QM)/DVOLT
      WRITE(*,2100) (QM_VOLT$+QM_VOLT)*0.5
      WRITE(*,*) 

C
      STOP
      END


