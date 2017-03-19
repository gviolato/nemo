C***********************************************************************
C    Module:  cdfun.f
C 
C    Copyright (C) 2005 Mark Drela 
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


      SUBROUTINE CDFUN(A, RE, MA, POLAR, REREF, REEXP, MCRIT,
     &     CD,CD_CL,CD_RE,CD_MA )
C
      USE POLARDATAMODULE
C      

      REAL MA, MCRIT
C
      DATA CDMF / 10.0 /
      DATA IEXP / 3 /
C
C---- CD-scaling factor
      FAC = (RE/REREF)**REEXP
      FAC_RE = REEXP*FAC/RE
C
C---- CD(CL;Re) function
c$$$  CLB = CL - CLCD0
c$$$  CD    = (CD0 + CD2*CLB**2 )*FAC
c$$$  CD_CL = (      CD2*CLB*2.0)*FAC
c$$$  CD_RE = (CD0 + CD2*CLB**2 )*FAC_RE
      CALL GETPOLARCDINFO(POLAR,A,CD,CD_CL)
      CD_CL = CD_CL*FAC
      CD_RE = CD*FAC_RE
      CD = CD*FAC
      CD_MA = 0.
C
      IF(MA .GT. MCRIT) THEN
       CD    = CD    + CDMF*(MA-MCRIT)** IEXP
       CD_MA = CD_MA + CDMF*(MA-MCRIT)**(IEXP-1) * FLOAT(IEXP)
      ENDIF
C
      RETURN
      END ! CDFUN



c$$$      SUBROUTINE DCDFUN(    BE,    WA,    WT, CLCD0,CL0,DCLDA,
c$$$     &              DCD,DCD_BE,DCD_WA,DCD_WT )
c$$$C
c$$$      A    = BE - ATAN2(WA,WT)
c$$$      A_BE = 1.0
c$$$      A_WA = -WT/(WA**2 + WT**2)
c$$$      A_WT =  WA/(WA**2 + WT**2)
c$$$C
c$$$      ACD0 = (CLCD0-CL0)/DCLDA
c$$$      DCD   = 2.0*SIN(A-ACD0)**2
c$$$      DCD_A = 4.0*SIN(A-ACD0)*COS(A-ACD0)
c$$$C
c$$$      DCD_BE = DCD_A*A_BE
c$$$      DCD_WA = DCD_A*A_WA
c$$$      DCD_WT = DCD_A*A_WT
c$$$C
c$$$      RETURN
c$$$      END ! DCDFUN
