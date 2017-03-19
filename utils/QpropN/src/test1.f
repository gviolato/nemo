
C---- SUBROUTINE GVCALC derivative tester

      LOGICAL STALL,LCONV
      REAL MCRIT

      EPS = 1.0E-5
C
C
C---- set typical input parameters 
      BLDS = 3.0

      C = 0.1
      B = 0.65
cc      B = 0.8   !  stalled
      R = 0.4
      RAD = 0.6
      VEL = 12.0
      OMG = 40.0
      VSO = 24.0
C
      CL0 = 0.5
      DCLDA = 5.8
      CLMIN = -0.4
      CLMAX = 1.3
      MCRIT = 0.6
C
      CALL GVCALC(C,B,R,
     &            BLDS,RAD,VEL,OMG,VSO,
     &            CL0,DCLDA,CLMIN,CLMAX,MCRIT,
     &            GAM,GAM_VEL,GAM_OMG,GAM_CH,GAM_BE,
     &             VA, VA_VEL, VA_OMG, VA_CH, VA_BE,
     &             VT, VT_VEL, VT_OMG, VT_CH, VT_BE,
     &             CL, CL_VEL, CL_OMG, CL_CH, CL_BE, STALL,LCONV)

      write(*,*) 'VA', VA
      write(*,*) 'VT', VT
      write(*,*) 'CL', CL
      write(*,*) 'GAM', GAM
      write(*,*) 'Mach', SQRT((VEL+VA)**2 + (OMG*R-VT)**2) / VSO

 2100 FORMAT(1X,8F16.10)

      VEL$ = VEL + EPS
      OMG$ = OMG
      C$ = C
      B$ = B
      CALL GVCALC(C$,B$,R,
     &            BLDS,RAD,VEL$,OMG$,VSO,
     &            CL0,DCLDA,CLMIN,CLMAX,MCRIT,
     &            GAM$,GAM_VEL$,GAM_OMG$,GAM_CH$,GAM_BE$,
     &             VA$, VA_VEL$, VA_OMG$, VA_CH$, VA_BE$,
     &             VT$, VT_VEL$, VT_OMG$, VT_CH$, VT_BE$,
     &             CL$, CL_VEL$, CL_OMG$, CL_CH$, CL_BE$, STALL,LCONV)
      WRITE(*,*) 
      WRITE(*,*) 'VEL'
      WRITE(*,2100) (VA$ - VA)/EPS,
     &           (VT$ - VT)/EPS,
     &           (CL$ - CL)/EPS,
     &           (GAM$ - GAM)/EPS
      WRITE(*,2100) (VA_VEL$+VA_VEL)*0.5,
     &           (VT_VEL$+VT_VEL)*0.5,
     &           (CL_VEL$+CL_VEL)*0.5,
     &           (GAM_VEL$+GAM_VEL)*0.5
      WRITE(*,*) 
C
C
      VEL$ = VEL
      OMG$ = OMG + EPS
      C$ = C
      B$ = B
      CALL GVCALC(C$,B$,R,
     &            BLDS,RAD,VEL$,OMG$,VSO,
     &            CL0,DCLDA,CLMIN,CLMAX,MCRIT,
     &            GAM$,GAM_VEL$,GAM_OMG$,GAM_CH$,GAM_BE$,
     &             VA$, VA_VEL$, VA_OMG$, VA_CH$, VA_BE$,
     &             VT$, VT_VEL$, VT_OMG$, VT_CH$, VT_BE$,
     &             CL$, CL_VEL$, CL_OMG$, CL_CH$, CL_BE$, STALL,LCONV)
      WRITE(*,*) 
      WRITE(*,*) 'OMG'
      WRITE(*,2100) (VA$ - VA)/EPS,
     &           (VT$ - VT)/EPS,
     &           (CL$ - CL)/EPS,
     &           (GAM$ - GAM)/EPS
      WRITE(*,2100) (VA_OMG$+VA_OMG)*0.5,
     &           (VT_OMG$+VT_OMG)*0.5,
     &           (CL_OMG$+CL_OMG)*0.5,
     &           (GAM_OMG$+GAM_OMG)*0.5
      WRITE(*,*) 

      VEL$ = VEL
      OMG$ = OMG
      C$ = C + EPS
      B$ = B
      CALL GVCALC(C$,B$,R,
     &            BLDS,RAD,VEL$,OMG$,VSO,
     &            CL0,DCLDA,CLMIN,CLMAX,MCRIT,
     &            GAM$,GAM_VEL$,GAM_OMG$,GAM_CH$,GAM_BE$,
     &             VA$, VA_VEL$, VA_OMG$, VA_CH$, VA_BE$,
     &             VT$, VT_VEL$, VT_OMG$, VT_CH$, VT_BE$,
     &             CL$, CL_VEL$, CL_OMG$, CL_CH$, CL_BE$, STALL,LCONV)
      WRITE(*,*) 
      WRITE(*,*) 'CH'
      WRITE(*,2100) (VA$ - VA)/EPS,
     &           (VT$ - VT)/EPS,
     &           (CL$ - CL)/EPS,
     &           (GAM$ - GAM)/EPS
      WRITE(*,2100) (VA_CH$+VA_CH)*0.5,
     &           (VT_CH$+VT_CH)*0.5,
     &           (CL_CH$+CL_CH)*0.5,
     &           (GAM_CH$+GAM_CH)*0.5
      WRITE(*,*) 

      VEL$ = VEL
      OMG$ = OMG
      C$ = C
      B$ = B + EPS
      CALL GVCALC(C$,B$,R,
     &            BLDS,RAD,VEL$,OMG$,VSO,
     &            CL0,DCLDA,CLMIN,CLMAX,MCRIT,
     &            GAM$,GAM_VEL$,GAM_OMG$,GAM_CH$,GAM_BE$,
     &             VA$, VA_VEL$, VA_OMG$, VA_CH$, VA_BE$,
     &             VT$, VT_VEL$, VT_OMG$, VT_CH$, VT_BE$,
     &             CL$, CL_VEL$, CL_OMG$, CL_CH$, CL_BE$, STALL,LCONV)
      WRITE(*,*) 
      WRITE(*,*) 'BE'
      WRITE(*,2100) (VA$ - VA)/EPS,
     &           (VT$ - VT)/EPS,
     &           (CL$ - CL)/EPS,
     &           (GAM$ - GAM)/EPS
      WRITE(*,2100) (VA_BE$+VA_BE)*0.5,
     &           (VT_BE$+VT_BE)*0.5,
     &           (CL_BE$+CL_BE)*0.5,
     &           (GAM_BE$+GAM_BE)*0.5
      WRITE(*,*) 

C
      STOP
      END


