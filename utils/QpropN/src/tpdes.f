
      SUBROUTINE TPDES(TDES1,PDES1, RHO,VEL,OMG, RAD,BLDS,
     &                 N, R,DR, CL,CD, CL0,DCLDA,
     &                 GAM, CH, BETA, ADW)
C---------------------------------------------------------------
C     Computes approximate propeller geometry 
C     using a modified Larrabee algorithm
C
C  Input:
C      TDES1    prescribed thrust (0 if power prescribed)
C      PDES1    prescribed power  (0 if thrust prescribed)
C                (if both TDES1=PDES1=0.0, impose maximum windmill power)
C      RHO      density
C      VEL      freestream velocity
C      OMG      rotation rate,  in radians/time-unit
C      RAD      tip radius
C      BLDS     number of blades
C
C      N        number of radial stations
C      R(.)     radii where geometry is to be defined
C      DR(.)    integration radius intervals,   Int [ ] dr  =  Sum [ ]_i DR(i)
C      CL(.)    prescribed CL
C      CD(.)    prescribed CD
C      CL0(.)   zero-lift angle
C      DCLDA(.) 2D lift-curve slope  dcl/da
C    
C  Output:
C     GAM(.)    circulation
C     CH(.)     chord
C     BETA(.)   pitch angle (radians)
C     ADW       wake advance ratio
C
C---------------------------------------------------------------
      REAL R(N), DR(N), CL(N), CD(N), CL0(N), DCLDA(N)
      REAL GAM(N), CH(N), BETA(N)
C
      DATA PI / 3.14159265 /
C
      TDES = TDES1
      PDES = PDES1
C
C---- set VD = velocity through disk, using actuator-disk theory
C
C-----------------------------------------------------------------
      IF(TDES1 .NE. 0.0) THEN

C----- Betz limit estimates
       DVDMIN = -0.50*VEL
       TMIN = RHO*PI*RAD**2
     &      * 2.0*( DVDMIN**2 + VEL*DVDMIN )
C
C----- specified thrust -- used thrust relation for actuator disk
       IF(TDES .LT. TMIN) THEN
        WRITE(*,*) 'Limiting windmill thrust to Betz limit:', TMIN
        TDES = 0.99 * TMIN
        DVD  = 0.90 * DVDMIN
       ELSE
        VT2 = TDES / (0.5*RHO * PI*RAD**2)
        DVD = 0.5*( SQRT(VEL**2 + VT2) - VEL )
       ENDIF
C
       VD = VEL + DVD
       TC = TDES / (0.5*RHO*VD**2 * PI*RAD**2)
C
C-----------------------------------------------------------------
      ELSE
C----- r^2 weighted cd/cl ratio
       EPSR3 = 0.
       DO I = 1, N
         EPSR3 = EPSR3 + CD(I)/CL(I) * R(I)**2 * DR(I)
       ENDDO
C
C----- additional effective velocity through disk accounting for viscous power loss
       VVIS = 2.0*OMG*EPSR3/RAD**2
C
C----- Betz limit estimates
ccc    DVDMIN = -0.33*VEL
       DVDMIN =  0.33*(SQRT(VEL**2+VEL*VVIS+VVIS**2) - 2.0*VEL - VVIS)
       PMIN = RHO*PI*RAD**2
     &      * 2.0*(                 DVDMIN**3
     &             + (2.0*VEL+VVIS)*DVDMIN**2
     &             + VEL*(VEL+VVIS)*DVDMIN   )
C
C----- specified power -- use power relation for actuator disk
       IF(PDES1 .EQ. 0.0 .OR. PDES .LT. PMIN) THEN
C------ guess power as nearly the actuator-disk Betz limit
        WRITE(*,*) 'Limiting windmill power to Betz limit:', PMIN
        PDES = 0.99 * PMIN
        DVD  = 0.90 * DVDMIN
C
       ELSE
C------ estimate DVD from power disk loading
        VP3 = PDES / (0.5*RHO * PI*RAD**2)
        DVDA = ABS(VP3)**(1.0/3.0)
        DVD = SIGN( DVDA , VP3 )
C
        DVD = MAX( DVD , DVDMIN+0.1*DVD0 )
       ENDIF
C
       DVD0 = ABS(DVD)
C
C----- Newton iteration for actual DVD
       DO ITER = 1, 20
         RES     =                DVD**3
     &           + (2.0*VEL+VVIS)*DVD**2
     &           + VEL*(VEL+VVIS)*DVD
     &           - 0.25 * PDES/(0.5*RHO * PI*RAD**2)
         RES_DVD =            3.0*DVD**2
     &           + (2.0*VEL+VVIS)*DVD * 2.0
     &           + VEL*(VEL+VVIS)
C
         DEL = -RES/RES_DVD

cc         write(*,'(1x,i3,4f12.4)') iter, dvd, res, res_dvd, del

         IF(DVD.EQ.DVDMIN .AND. DEL.LT.0.0) THEN
          WRITE(*,*) '? Windmill power exceeds Betz limit:', PMIN
          GO TO 8
         ENDIF
         DVD = DVD + DEL
         IF(ABS(DEL/DVD0) .LT. 1.0E-5) GO TO 8
         DVD = MAX( DVD , DVDMIN )
       ENDDO
       WRITE(*,*) '***  Vd convergence failed.  dVd =', DEL
 8     CONTINUE
C
       VD = VEL + DVD
       PC = PDES/(0.5*RHO*VD**3 * PI*RAD**2)
      ENDIF
C-----------------------------------------------------------------
C
      ADV = VEL/(OMG*RAD)
      ADW = VD /(OMG*RAD)
C
C
      WRITE(*,1050) INT(BLDS), ADV, ADW
 1050 FORMAT(/' -----------------------------------------------------'
     &       /' Larrabee-method initialization...'
     &      //'   blades  B =', I2,
     &       /'   true advance ratio  V /wR =', F8.4
     &       /'   wake advance ratio  Vd/wR =', F8.4 )
C
C---- go over radial stations, summing I1, I2, J1, J2 integrals
      WRITE(*,*)
      WRITE(*,*)
     & '  r/R    wr/Vd     F       G    dI1/dx  dI2/dx  dJ1/dx  dJ2/dx'
      RI1 = 0.
      RI2 = 0.
      RJ1 = 0.
      RJ2 = 0.
      DO 10 I = 1, N
        XI = R(I)/RAD
        DXI = DR(I)/RAD
C
        XX = XI/ADW
C
        ARG = MIN( 20.0 , 0.5*BLDS*(1.0-XI)/ADW )
        EK = EXP(-ARG)
        F = ATAN2( SQRT(1.0 - EK*EK) , EK )*2.0/PI
C
        FMOD = F * SQRT(1.0 + (4.0*ADW/(PI*BLDS*XI))**2)
C
        G = FMOD * XX*XX / (1.0 + XX*XX)
        DOL = CD(I)/CL(I)
        RIX1 = 4.0*XI*G*(1.0 - DOL/XX)
        RIX2 = 2.0*XI*G*(1.0 - DOL/XX) / (1.0 + XX*XX)
        RJX1 = 4.0*XI*G*(1.0 + DOL*XX)
        RJX2 = 2.0*XI*G*(1.0 + DOL*XX) * XX*XX / (1.0 + XX*XX)
        RI1 = RI1 + RIX1*DXI
        RI2 = RI2 + RIX2*DXI
        RJ1 = RJ1 + RJX1*DXI
        RJ2 = RJ2 + RJX2*DXI
        WRITE(*,1100) XI, XX, F, G, RIX1,RIX2,RJX1,RJX2
 1100   FORMAT(2X,F5.3, 7F8.4)
C
C------ save G for re-use later
        GAM(I) = G
   10 CONTINUE
C
      WRITE(*,1120) RI1,RI2,RJ1,RJ2
 1120 FORMAT(// ' I1 = ', F8.4, '    I2 = ', F8.4,
     &    '       J1 = ', F8.4, '    J2 = ', F8.4 )
C
C---- Calculate displacement velocity and total efficiency
      IF(TDES1 .NE. 0.0) THEN
       DISCR = 1.0 - 4.0*TC*RI2/RI1**2
       IF(DISCR.LT.0.0) THEN
        WRITE(*,*) 'Setting thrust to minimum'
        DISCR = 0.0
        TC = 0.25*RI1**2 / RI2
       ENDIF
       ZETA = 0.5*RI1/RI2 * (1.0 - SQRT(DISCR))
       PC = RJ1*ZETA + RJ2*ZETA**2
       ETA = TC/PC * VEL/VD
       T = 0.5*RHO*VD**2 * PI*RAD**2 * TC
       P = 0.5*RHO*VD**3 * PI*RAD**2 * PC
       WRITE(*,1150) TC, T, ZETA, PC, P, ETA
 1150  FORMAT(/ ' Tc =', F9.4,'  T =',G13.5,' =>   zeta =', F9.4
     &        / ' Pc =', F9.4,'  P =',G13.5,'      eff. =', F9.4)
C
      ELSE
       DISCR = 1.0 + 4.0*PC*RJ2/RJ1**2
       IF(DISCR .LT. 0.0 .OR. PDES1 .EQ. 0.0) THEN
        WRITE(*,*) 'Setting power to minimum'
        DISCR = 0.0
        PC = -0.25*RJ1**2 / RJ2
       ENDIF
       ZETA = 0.5*RJ1/RJ2 * (SQRT(DISCR) - 1.0)
       TC = RI1*ZETA - RI2*ZETA**2
       ETA = TC/PC * VEL/VD
       T = 0.5*RHO*VD**2 * PI*RAD**2 * TC
       P = 0.5*RHO*VD**3 * PI*RAD**2 * PC
       WRITE(*,1160) PC, P, ZETA, TC, T, ETA
 1160  FORMAT(/ ' Pc =', F9.4,'  P =',G13.5,' =>   zeta =', F9.4
     &        / ' Tc =', F9.4,'  T =',G13.5,'      eff. =', F9.4)
      ENDIF
C
C---------------------------------------------------------
C---- Calculate geometry
      DO I = 1, N
        G = GAM(I)
C
        PHI = ATAN2( VEL + 0.5*ZETA*VD , OMG*R(I) )
        WSQ = VEL**2 + (OMG*R(I))**2 - (VD*0.5*ZETA*COS(PHI))**2
        W = SQRT(WSQ)

        BETA(I) = PHI + (CL(I)-CL0(I))/DCLDA(I)
C
        GAM(I) = G * (2.0*PI/BLDS) * ZETA * VD**2/OMG
        CH(I) = 2.0*GAM(I)/(W*CL(I))
      ENDDO
C
      RETURN
      END ! TPDES

