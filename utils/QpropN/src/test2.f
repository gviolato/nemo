
C---- SUBROUTINE TQCALC derivative tester

      PARAMETER (NDIM=10)
      REAL C(NDIM), B(NDIM), R(NDIM), DR(NDIM)
      REAL CL0(NDIM), DCLDA(NDIM), CLMIN(NDIM), CLMAX(NDIM), MCRIT(NDIM)
      REAL CD0(NDIM), CD2U(NDIM), CD2L(NDIM), CLCD0(NDIM)
      REAL REREF(NDIM), REEXP(NDIM)
C
      REAL C$(NDIM), B$(NDIM)
C
      REAL VA(NDIM), VT(NDIM), CL(NDIM), CD(NDIM)
      REAL TP_B(NDIM), TP_C(NDIM),
     &     QP_B(NDIM), QP_C(NDIM)

      REAL TP_B$(NDIM), TP_C$(NDIM),
     &     QP_B$(NDIM), QP_C$(NDIM)

      LOGICAL STALL(NDIM)
C
      EPS = 1.0E-6
C
      N = 2

C---- set typical input parameters 
      C(1) = 0.1
      C(2) = 0.07
      B(1) = 0.40
      B(2) = 0.08
      R(1) = 0.4
      R(2) = 0.7
      DR(1) = 0.3
      DR(2) = 0.2
C
      BLDS = 3.0
      RAD = 0.9
      VEL = 2.0
      VSO = 34.0
      OMG = 40.0
      DBE = 0.02
C
      DO I = 1, N
        CL0(I) = 0.5
        DCLDA(I) = 5.8
        CLMIN(I) = -0.4
        CLMAX(I) = 1.0
C
        CD0(I) = 0.02
        CD2U(I) = 0.08
        CD2L(I) = 0.12
        CLCD0(I) = 0.5
        REREF(I) = 1.0E5
        REEXP(I) = -0.5
C
        MCRIT(I) = 0.7
      ENDDO
C
      RHO = 1.2
      RMU = 1.8E-5
C
      DO I = 1, N
        C$(I) = C(I)
        B$(I) = B(I)
      ENDDO
C
C
C
      CALL TQCALC(N,C,B,R,DR,
     &              VA,VT,CL,CD,STALL,
     &              BLDS,RAD,VEL,OMG,DBE,
     &              RHO,RMU,VSO,
     &              CL0,DCLDA,CLMIN,CLMAX,MCRIT,
     &              CD0,CD2U,CD2L,CLCD0,REREF,REEXP,
     &              TP, TP_VEL, TP_OMG, TP_DBE, TP_C, TP_B,
     &              QP, QP_VEL, QP_OMG, QP_DBE, QP_C, QP_B )

      AM = SQRT((VEL+VA(I))**2 + (OMG*R(I)-VT(I))**2)
      write(*,*) 'VA', (VA(I), I=1, N)
      write(*,*) 'VT', (VT(I), I=1, N)
      write(*,*) 'Ma', (SQRT((VEL+VA(I))**2+(OMG*R(I)-VT(I))**2)/VSO,
     &                          I=1,N)
      write(*,*) 'CL', (CL(I), I=1, N)
      write(*,*) 'st', (STALL(I), I=1, N)
      write(*,*) 'TP', TP
      write(*,*) 'QP', QP

 2100 FORMAT(1X,8F16.10)

      I = N
C
C---------------------------------------
      VEL$ = VEL + EPS
      OMG$ = OMG
      DBE$ = DBE
      B$(I) = B(I)
      C$(I) = C(I)
      CALL TQCALC(N,C$,B$,R,DR,
     &              VA,VT,CL,CD,STALL,
     &              BLDS,RAD,VEL$,OMG$,DBE$,
     &              RHO,RMU,VSO,
     &              CL0,DCLDA,CLMIN,CLMAX,MCRIT,
     &              CD0,CD2U,CD2L,CLCD0,REREF,REEXP,
     &              TP$, TP_VEL$, TP_OMG$, TP_DBE$, TP_C$, TP_B$,
     &              QP$, QP_VEL$, QP_OMG$, QP_DBE$, QP_C$, QP_B$ )

      WRITE(*,*) 
      WRITE(*,*) 'VEL'
      WRITE(*,2100)
     &           (TP$ - TP)/EPS,
     &           (QP$ - QP)/EPS
      WRITE(*,2100)
     &           (TP_VEL$+TP_VEL)*0.5,
     &           (QP_VEL$+QP_VEL)*0.5
      WRITE(*,*) 
C
C
C---------------------------------------
      VEL$ = VEL
      OMG$ = OMG + EPS
      DBE$ = DBE
      B$(I) = B(I)
      C$(I) = C(I)
      CALL TQCALC(N,C$,B$,R,DR,
     &              VA,VT,CL,CD,STALL,
     &              BLDS,RAD,VEL$,OMG$,DBE$,
     &              RHO,RMU,VSO,
     &              CL0,DCLDA,CLMIN,CLMAX,MCRIT,
     &              CD0,CD2U,CD2L,CLCD0,REREF,REEXP,
     &              TP$, TP_VEL$, TP_OMG$, TP_DBE$, TP_C$, TP_B$,
     &              QP$, QP_VEL$, QP_OMG$, QP_DBE$, QP_C$, QP_B$ )

      WRITE(*,*) 
      WRITE(*,*) 'OMG'
      WRITE(*,2100)
     &           (TP$ - TP)/EPS,
     &           (QP$ - QP)/EPS
      WRITE(*,2100)
     &           (TP_OMG$+TP_OMG)*0.5,
     &           (QP_OMG$+QP_OMG)*0.5
      WRITE(*,*) 
C
C---------------------------------------
      VEL$ = VEL
      OMG$ = OMG
      DBE$ = DBE + EPS
      B$(I) = B(I)
      C$(I) = C(I)
      CALL TQCALC(N,C$,B$,R,DR,
     &              VA,VT,CL,CD,STALL,
     &              BLDS,RAD,VEL$,OMG$,DBE$,
     &              RHO,RMU,VSO,
     &              CL0,DCLDA,CLMIN,CLMAX,MCRIT,
     &              CD0,CD2U,CD2L,CLCD0,REREF,REEXP,
     &              TP$, TP_VEL$, TP_OMG$, TP_DBE$, TP_C$, TP_B$,
     &              QP$, QP_VEL$, QP_OMG$, QP_DBE$, QP_C$, QP_B$ )


      WRITE(*,*) 
      WRITE(*,*) 'DBE'
      WRITE(*,2100)
     &           (TP$ - TP)/EPS,
     &           (QP$ - QP)/EPS
      WRITE(*,2100)
     &           (TP_DBE$+TP_DBE)*0.5,
     &           (QP_DBE$+QP_DBE)*0.5
      WRITE(*,*) 
C
C---------------------------------------
      VEL$ = VEL
      OMG$ = OMG
      DBE$ = DBE
      B$(I) = B(I) + EPS
      C$(I) = C(I)
      CALL TQCALC(N,C$,B$,R,DR,
     &              VA,VT,CL,CD,STALL,
     &              BLDS,RAD,VEL$,OMG$,DBE$,
     &              RHO,RMU,VSO,
     &              CL0,DCLDA,CLMIN,CLMAX,MCRIT,
     &              CD0,CD2U,CD2L,CLCD0,REREF,REEXP,
     &              TP$, TP_VEL$, TP_OMG$, TP_DBE$, TP_C$, TP_B$,
     &              QP$, QP_VEL$, QP_OMG$, QP_DBE$, QP_C$, QP_B$ )


      WRITE(*,*) 
      WRITE(*,*) 'B'
      WRITE(*,2100)
     &           (TP$ - TP)/EPS,
     &           (QP$ - QP)/EPS
      WRITE(*,2100)
     &           (TP_B$(I)+TP_B(I))*0.5,
     &           (QP_B$(I)+QP_B(I))*0.5
      WRITE(*,*) 
C
C---------------------------------------
      VEL$ = VEL
      OMG$ = OMG
      DBE$ = DBE
      B$(I) = B(I)
      C$(I) = C(I) + EPS
      CALL TQCALC(N,C$,B$,R,DR,
     &              VA,VT,CL,CD,STALL,
     &              BLDS,RAD,VEL$,OMG$,DBE$,
     &              RHO,RMU,VSO,
     &              CL0,DCLDA,CLMIN,CLMAX,MCRIT,
     &              CD0,CD2U,CD2L,CLCD0,REREF,REEXP,
     &              TP$, TP_VEL$, TP_OMG$, TP_DBE$, TP_C$, TP_B$,
     &              QP$, QP_VEL$, QP_OMG$, QP_DBE$, QP_C$, QP_B$ )


      WRITE(*,*) 
      WRITE(*,*) 'C'
      WRITE(*,2100)
     &           (TP$ - TP)/EPS,
     &           (QP$ - QP)/EPS
      WRITE(*,2100)
     &           (TP_C$(I)+TP_C(I))*0.5,
     &           (QP_C$(I)+QP_C(I))*0.5
      WRITE(*,*) 
C
C
      STOP
      END

