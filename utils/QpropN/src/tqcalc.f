C***********************************************************************
C    Module:  tqcalc.f
C 
C    Copyright (C) 2003 Mark Drela 
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************

      SUBROUTINE TQCALC(N,C,B,R,DR,
     &     VA,VT,CL,CD,STALL,
     &     BLDS,RAD,VEL,OMG,DBE,
     &     RHO,RMU,VSO,
     &     POLAR,MCRIT,REREF,REEXP,
     &     TP, TP_VEL, TP_OMG, TP_DBE, TP_C, TP_B,
     &     QP, QP_VEL, QP_OMG, QP_DBE, QP_C, QP_B )
C
      USE POLARDATAMODULE
C      
      IMPLICIT REAL(A-H,M,O-Z)
      REAL C(N), B(N), R(N), DR(N)
      TYPE (POLARDATA) POLAR(N)
      REAL MCRIT(N)
      REAL REREF(N), REEXP(N)
C
      REAL VA(N), VT(N), CL(N), CD(N)
      REAL TP_C(N), TP_B(N),
     &     QP_C(N), QP_B(N), RES_A(N)
      LOGICAL STALL(N)
C
C-----------------------------------------------------
C     Computes propeller thrust and torque
C     by integrating forces along radius.
C
C
C  Input   N      number of radial stations  i = 1..N
C          C(i)   chords
C          B(i)   angles (radians)
C          R(i)   radii
C          DR(i)  radius interval for i station
C
C          BLDS    number of blades
C          RAD     tip radius  (for Prandtl's factor)
C          VEL     forward flight speed
C          OMG     rotational speed  (radians/time)
C          DBE     delta(B) to be added to all B(i) angles
C
C          RHO     fluid density
C          RMU     fluid viscosity
C
C          CL0(i)     constants for CL(alpha,M) function
C          DCLDA(i)
C          CLMIN(i)
C          CLMAX(i)
C          MCRIT(i)
C
C          CD0(i)     constant  coefficient for CD(CL) function
C          CD2(i)     quadratic coefficient for CD(CL) function
C          CLCD0(i)   CL at minimum drag point (where CD=CD0)
C          REREF(i)   reference Reynolds number where CD(CL) applies
C          REEXP(i)   Re-scaling exponent:   CD ~ (Re/REREF)^REEXP
C
C
C  Output  VA(i)    axial induced velocity
C          VT(i)    tangential induced velocity
C          CL(i)    section lift coefficient
C          CD(i)    section drag coefficient
C          STALL(i) T if alpha is outside stall limits
C
C          TP       prop thrust
C          QP       prop torque
C          ()_VEL   derivatives  d()/dVEL
C          ()_OMG   derivatives  d()/dOMG
C          ()_DBE   derivatives  d()/dDBE
C          ()_C(i)  derivatives  d()/dC(i)
C          ()_B(i)  derivatives  d()/dB(i)
C
C-----------------------------------------------------
      LOGICAL LCONV
C
      DATA MSQMAX    / 0.9    /
C
C---- clear for accumulation
      TP     = 0.
      TP_VEL = 0.
      TP_OMG = 0.
      TP_DBE = 0.
C
      QP     = 0.
      QP_VEL = 0.
      QP_OMG = 0.
      QP_DBE = 0.
C
C---- go over radial stations
      DO I = 1, N
        BTOT = B(I) + DBE
        CALL GVCALC(C(I),BTOT,R(I),POLAR(I),
     &              BLDS,RAD,VEL,OMG,VSO,
     &              MCRIT(I),
     &        GAM  ,GAM_VEL,GAM_OMG,GAM_C,GAM_B, 
     &        VA(I), VA_VEL, VA_OMG, VA_C, VA_B,
     &        VT(I), VT_VEL, VT_OMG, VT_C, VT_B,
     &       CL(I), CL_VEL, CL_OMG, CL_C, CL_B, STALL(I), LCONV,
     &       RES_A(I))
C
C------ perturbation imposed axial and tangential velocities
        U0A = 0.
        U0T = 0.
        GAM_U0A = 0.
        GAM_U0T = 0.
        VA_U0A = 0.
        VA_U0T = 0.
        VT_U0A = 0.
        VT_U0T = 0.
        CL_U0A = 0.
        CL_U0T = 0.
C
C------ total imposed axial and tangential velocities
        UA     = VEL   + U0A
        UA_VEL = 1.0
        UA_U0A =         1.0
C
        UT     = OMG*R(I) - U0T
        UT_OMG =     R(I)
        UT_U0T =          - 1.0
C
C------ geometric velocity
        WZ = SQRT(UA**2 + UT**2)
        WZ_VEL = (UA/WZ)*UA_VEL
        WZ_OMG = (UT/WZ)*UT_OMG
C
C------ total axial and tangential velocities
        WA     = UA     + VA(I)
        WA_VEL = UA_VEL + VA_VEL
        WA_OMG =          VA_OMG
        WA_B   =          VA_B  
        WA_C   =          VA_C
        WA_U0A = UA_U0A + VA_U0A
        WA_U0T =          VA_U0T
C
        WT     = UT     - VT(I)
        WT_VEL =        - VT_VEL
        WT_OMG = UT_OMG - VT_OMG
        WT_B   =        - VT_B  
        WT_C   =        - VT_C 
        WT_U0A =        - VT_U0A
        WT_U0T = UT_U0T - VT_U0T
C
C------ total velocity^2
        WSQ = WA**2 + WT**2
        WSQ_VEL = 2.0*WA*WA_VEL + 2.0*WT*WT_VEL
        WSQ_OMG = 2.0*WA*WA_OMG + 2.0*WT*WT_OMG
        WSQ_B   = 2.0*WA*WA_B   + 2.0*WT*WT_B  
        WSQ_C   = 2.0*WA*WA_C   + 2.0*WT*WT_C  
        WSQ_U0A = 2.0*WA*WA_U0A + 2.0*WT*WT_U0A
        WSQ_U0T = 2.0*WA*WA_U0T + 2.0*WT*WT_U0T
C
C------ total velocity
        W = SQRT(WSQ)
        W_VEL = 0.5*WSQ_VEL/W
        W_OMG = 0.5*WSQ_OMG/W
        W_B   = 0.5*WSQ_B  /W
        W_C   = 0.5*WSQ_C  /W
        W_U0A = 0.5*WSQ_U0A/W
        W_U0T = 0.5*WSQ_U0T/W
C
C------ chord Reynolds number
        RE     = RHO*C(I)*W    /RMU
        RE_VEL = RHO*C(I)*W_VEL/RMU
        RE_OMG = RHO*C(I)*W_OMG/RMU
        RE_B   = RHO*C(I)*W_B  /RMU
        RE_C   = RHO*C(I)*W_C  /RMU
     &         + RHO     *W    /RMU
        RE_U0A = RHO*C(I)*W_U0A/RMU
        RE_U0T = RHO*C(I)*W_U0T/RMU
C
C------ local Mach and PG factor
        MSQ     = WSQ     / VSO**2
        MSQ_VEL = WSQ_VEL / VSO**2
        MSQ_OMG = WSQ_OMG / VSO**2
        MSQ_B   = WSQ_B   / VSO**2
        MSQ_C   = WSQ_C   / VSO**2
        MSQ_U0A = WSQ_U0A / VSO**2
        MSQ_U0T = WSQ_U0T / VSO**2
        IF(MSQ .GT. MSQMAX) THEN
         MSQ = MSQMAX
         MSQ_VEL = 0.
         MSQ_OMG = 0.
         MSQ_B   = 0.
         MSQ_C   = 0.
         MSQ_U0A = 0.
         MSQ_U0T = 0.
        ENDIF
C
        PG = 1.0 / SQRT(1.0 - MSQ)
        PG_MSQ = 0.5*PG / (1.0 - MSQ)
C
        PG_VEL = PG_MSQ*MSQ_VEL
        PG_OMG = PG_MSQ*MSQ_OMG
        PG_B   = PG_MSQ*MSQ_B
        PG_C   = PG_MSQ*MSQ_C
        PG_U0A = PG_MSQ*MSQ_U0A
        PG_U0T = PG_MSQ*MSQ_U0T
C
        MA = SQRT(MSQ)
        MA = MAX( MA , 1.0E-8 )
        MA_VEL = (0.5/MA)*MSQ_VEL
        MA_OMG = (0.5/MA)*MSQ_OMG
        MA_B   = (0.5/MA)*MSQ_B  
        MA_C   = (0.5/MA)*MSQ_C  
        MA_U0A = (0.5/MA)*MSQ_U0A
        MA_U0T = (0.5/MA)*MSQ_U0T
C
C------set unstalled CD
        CALL CDFUN(RES_A(I),RE,MA, 
     &       POLAR(I), REREF, REEXP(I), MCRIT,
     &       CD(I), CD_CL, CD_RE, CD_MA)
        CD_VEL = CD_CL*CL_VEL + CD_RE*RE_VEL + CD_MA*MA_VEL
        CD_OMG = CD_CL*CL_OMG + CD_RE*RE_OMG + CD_MA*MA_OMG
        CD_B   = CD_CL*CL_B   + CD_RE*RE_B   + CD_MA*MA_B   
        CD_C   = CD_CL*CL_C   + CD_RE*RE_C   + CD_MA*MA_C   
        CD_U0A = CD_CL*CL_U0A + CD_RE*RE_U0A + CD_MA*MA_U0A
        CD_U0T = CD_CL*CL_U0T + CD_RE*RE_U0T + CD_MA*MA_U0T
C
c$$$        IF(STALL(I)) THEN
c$$$C------- additional CD with normal-force stall model
c$$$         CALL DCDFUN(BTOT,WA,WT, CLCD0(I),CL0(I),DCLDA(I),
c$$$     &               DCD,DCD_B,DCD_WA,DCD_WT)
c$$$C
c$$$         CD(I)  = CD(I)  + DCD
c$$$         CD_VEL = CD_VEL + DCD_WA*WA_VEL + DCD_WT*WT_VEL
c$$$         CD_OMG = CD_OMG + DCD_WA*WA_OMG + DCD_WT*WT_OMG
c$$$         CD_B   = CD_B   + DCD_WA*WA_B   + DCD_WT*WT_B   + DCD_B
c$$$         CD_C   = CD_C   + DCD_WA*WA_C   + DCD_WT*WT_C  
c$$$         CD_U0A = CD_U0A + DCD_WA*WA_U0A + DCD_WT*WT_U0A
c$$$         CD_U0T = CD_U0T + DCD_WA*WA_U0T + DCD_WT*WT_U0T
c$$$C
c$$$        ENDIF
C
C------ Axial and tangential load/span
        HRC = 0.5*RHO*C(I)
        FA     = HRC*W    *(  CL(I) *WT     - CD(I) *WA    )
        FA_VEL = HRC*W_VEL*(  CL(I) *WT     - CD(I) *WA    )
     &         + HRC*W    *(  CL_VEL*WT     - CD_VEL*WA
     &                      + CL(I) *WT_VEL - CD(I) *WA_VEL)
        FA_OMG = HRC*W_OMG*(  CL(I) *WT     - CD(I) *WA    )
     &         + HRC*W    *(  CL_OMG*WT     - CD_OMG*WA
     &                      + CL(I) *WT_OMG - CD(I) *WA_OMG)
        FA_B   = HRC*W_B  *(  CL(I) *WT     - CD(I) *WA    )
     &         + HRC*W    *(  CL_B  *WT     - CD_B  *WA
     &                      + CL(I) *WT_B   - CD(I) *WA_B  )
        FA_C   = HRC*W_C  *(  CL(I) *WT     - CD(I) *WA    )
     &         + HRC*W    *(  CL_C  *WT     - CD_C  *WA
     &                      + CL(I) *WT_C   - CD(I) *WA_C  )
     &     + 0.5*RHO*W    *(  CL(I) *WT     - CD(I) *WA    )
C
        FT     = HRC*W    *(  CL(I) *WA     + CD(I) *WT    )
        FT_VEL = HRC*W_VEL*(  CL(I) *WA     + CD(I) *WT    )
     &         + HRC*W    *(  CL_VEL*WA     + CD_VEL*WT
     &                      + CL(I) *WA_VEL + CD(I) *WT_VEL)
        FT_OMG = HRC*W_OMG*(  CL(I) *WA     + CD(I) *WT    )
     &         + HRC*W    *(  CL_OMG*WA     + CD_OMG*WT
     &                      + CL(I) *WA_OMG + CD(I) *WT_OMG)
        FT_B   = HRC*W_B  *(  CL(I) *WA     + CD(I) *WT    )
     &         + HRC*W    *(  CL_B  *WA     + CD_B  *WT
     &                      + CL(I) *WA_B   + CD(I) *WT_B  )
        FT_C   = HRC*W_C  *(  CL(I) *WA     + CD(I) *WT    )
     &         + HRC*W    *(  CL_C  *WA     + CD_C  *WT
     &                      + CL(I) *WA_C   + CD(I) *WT_C  )
     &     + 0.5*RHO*W    *(  CL(I) *WA     + CD(I) *WT    )
C
C------ sum thrust and torque
        TP     = TP     + BLDS*FA    *DR(I)
        TP_VEL = TP_VEL + BLDS*FA_VEL*DR(I)
        TP_OMG = TP_OMG + BLDS*FA_OMG*DR(I)
        TP_DBE = TP_DBE + BLDS*FA_B  *DR(I)
        TP_B(I) =         BLDS*FA_B  *DR(I)
        TP_C(I) =         BLDS*FA_C  *DR(I)
C
        QP     = QP     + BLDS*FT    *DR(I)*R(I)
        QP_VEL = QP_VEL + BLDS*FT_VEL*DR(I)*R(I)
        QP_OMG = QP_OMG + BLDS*FT_OMG*DR(I)*R(I)
        QP_DBE = QP_DBE + BLDS*FT_B  *DR(I)*R(I)
        QP_B(I) =         BLDS*FT_B  *DR(I)*R(I)
        QP_C(I) =         BLDS*FT_C  *DR(I)*R(I)

c       write(*,*) i
c       write(*,*) tp, qp
c       write(*,*) tp_vel, qp_vel
c       write(*,*) tp_omg, qp_omg
c       write(*,*) tp_b(i),qp_b(i)
c       write(*,*) tp_c(i),qp_c(i)
c       pause
c
C
      ENDDO
C
      RETURN
      END
