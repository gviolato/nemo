
C---- SUBROUTINE CDFUN and DCDFUN tester
      REAL MCRIT, MA, MA$

      EPS = 1.0E-5

      DCLDA = 6.2
      CLCD0 = 0.5
      CL0 = 0.3
      CD0 = 0.09
      CD2 = 0.12
      REREF = 1.0E3
      REEXP = -0.7
      MCRIT = 0.75
C
C
      CL = 1.7
      RE = 3.0E3
      MA = 0.85
C
      BE = 0.7
      WT = 2.0
      WA = WT*TAN(BE-0.15)
C
      WRITE(*,*)
C
 2100 FORMAT(1X,2F16.10)
C
      CALL CDFUN(   CL,   RE,   MA, CLCD0,CD0,CD2,REREF,REEXP,MCRIT,
     &        CD,CD_CL,CD_RE,CD_MA )
C
      write(*,*) 'CD', CD
C-----------------------
      CL$ = CL + EPS
      RE$ = RE
      MA$ = MA
      CALL CDFUN(  CL$,   RE$,    MA$, CLCD0,CD0,CD2,REREF,REEXP,MCRIT,
     &      CD$,CD_CL$,CD_RE$, CD_MA$ )
C
      WRITE(*,*) 
      WRITE(*,*) 'CL'
      WRITE(*,2100) (CD$ - CD)/EPS
      WRITE(*,2100) (CD_CL$+CD_CL)*0.5
      WRITE(*,*) 
C-----------------------
      DRE = RE*EPS
      CL$ = CL 
      RE$ = RE + DRE
      MA$ = MA
      CALL CDFUN(  CL$,   RE$,    MA$, CLCD0,CD0,CD2,REREF,REEXP,MCRIT,
     &      CD$,CD_CL$,CD_RE$, CD_MA$ )
C
      WRITE(*,*) 
      WRITE(*,*) 'RE'
      WRITE(*,2100) (CD$ - CD)/DRE
      WRITE(*,2100) (CD_RE$+CD_RE)*0.5
      WRITE(*,*) 
C
C-----------------------
      DMA = EPS
      CL$ = CL 
      RE$ = RE
      MA$ = MA + DMA
      CALL CDFUN(  CL$,   RE$,    MA$, CLCD0,CD0,CD2,REREF,REEXP,MCRIT,
     &      CD$,CD_CL$,CD_RE$, CD_MA$ )
C
      WRITE(*,*) 
      WRITE(*,*) 'MA'
      WRITE(*,2100) (CD$ - CD)/DMA
      WRITE(*,2100) (CD_MA$+CD_MA)*0.5
      WRITE(*,*) 
C
C
C--------------------------------------------------------------------
      CALL DCDFUN(    BE,    WA,    WT, CLCD0,CL0,DCLDA,
     &        DCD,DCD_BE,DCD_WA,DCD_WT )
C
      write(*,*) 'DCD', DCD

C-----------------------
      BE$ = BE + EPS
      WA$ = WA
      WT$ = WT
      CALL DCDFUN(    BE$,    WA$,    WT$, CLCD0,CL0,DCLDA,
     &       DCD$,DCD_BE$,DCD_WA$,DCD_WT$ )
C
      WRITE(*,*) 
      WRITE(*,*) 'BE'
      WRITE(*,2100) (DCD$ - DCD)/EPS
      WRITE(*,2100) (DCD_BE$+DCD_BE)*0.5
      WRITE(*,*) 
C-----------------------
      BE$ = BE
      WA$ = WA + EPS
      WT$ = WT
      CALL DCDFUN(    BE$,    WA$,    WT$, CLCD0,CL0,DCLDA,
     &       DCD$,DCD_BE$,DCD_WA$,DCD_WT$ )
C
      WRITE(*,*) 
      WRITE(*,*) 'WA'
      WRITE(*,2100) (DCD$ - DCD)/EPS
      WRITE(*,2100) (DCD_WA$+DCD_WA)*0.5
      WRITE(*,*) 
C-----------------------
      BE$ = BE
      WA$ = WA
      WT$ = WT + EPS
      CALL DCDFUN(    BE$,    WA$,    WT$, CLCD0,CL0,DCLDA,
     &       DCD$,DCD_BE$,DCD_WA$,DCD_WT$ )
C
      WRITE(*,*) 
      WRITE(*,*) 'WT'
      WRITE(*,2100) (DCD$ - DCD)/EPS
      WRITE(*,2100) (DCD_WT$+DCD_WT)*0.5
      WRITE(*,*) 
C
      STOP
      END


