C***********************************************************************
C    Module:  motor.f
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

      SUBROUTINE MOTORQ(  OMEGA,   VOLT, IMOTYPE, PARMOT,NMPAR, 
     &               Q, Q_OMEGA, Q_VOLT,
     &               I, I_OMEGA, I_VOLT )
C-----------------------------------------------------------
C     Motor+gearbox   torque(rpm,Voltage)  function.
C
C Input:  OMEGA      output shaft rotation rate  radians/s
C         VOLT       terminal voltage (or throttle for IC motors)
C         IMOTYPE    specifies type of motor model to be used
C         PARMOT(.)  motor parameters  (lines 3,4... in motor file)
C         NMPAR      number of motor parameters in motor file
C
C Output: Q          output shaft torque   (N-m)
C         Q_OMEGA    dQ/dOMEGA  function derivative
C         Q_VOLT     dQ/dVOLT   function derivative
C         I          motor current (or fuel flow for IC motors)
C         I_OMEGA    dI/dOMEGA  function derivative
C         I_VOLT     dI/dVOLT   function derivative
C
C-----------------------------------------------------------
      REAL PARMOT(*)
      REAL I, I_OMEGA, I_VOLT
C
      REAL KVRPM, KQRPM, KVRAD, KQRAD
      DATA PI /3.1415926/
      DATA EPS /1.0E-6/
C
C----------------------------------------------------
      IF(IMOTYPE.EQ.1) THEN
C----- Brushed DC motor - 1st-order model
       IF(NMPAR.LT.3) THEN
        WRITE(*,*) 'MOTORQ: Motor model 1 needs  3  parameters.',
     &             '  Number passed in:', IMOTYPE
        STOP
       ENDIF
C
       RMOTOR = PARMOT(1)    ! R   (Ohms)      motor resistance
       ZLOADI = PARMOT(2)    ! Io  (Amps)      zero-load current
       KVRPM  = PARMOT(3)    ! Kv  (rpm/Volt)  motor constant
C
       KVRAD = KVRPM * PI/30.0
       KQRAD = KVRAD
C     
       VM       = OMEGA/KVRAD
       VM_OMEGA = 1.0  /KVRAD
C
       I       = (VOLT - VM      )/RMOTOR
       I_OMEGA =       - VM_OMEGA /RMOTOR
       I_VOLT  =  1.0             /RMOTOR
C
       Q       = (I       - ZLOADI)/KQRAD
       Q_OMEGA =  I_OMEGA          /KQRAD
       Q_VOLT  =  I_VOLT           /KQRAD
C
      ELSEIF(IMOTYPE.EQ.2) THEN
C----- Brushed DC motor - 2nd-order model
       IF(NMPAR.LT.3) THEN
        WRITE(*,*) 'MOTORQ: Motor model 2 needs at least 3 parameters.',
     &             '  Number passed in:', IMOTYPE
        STOP
       ENDIF
C
       RMOTOR0 = PARMOT(1)    ! R0  (Ohms)      motor resistance
       ZLOADI0 = PARMOT(2)    ! Io0 (Amps)      zero-load current
       KVRPM   = PARMOT(3)    ! Kv  (rpm/Volt)  motor constant
       KQRPM   = PARMOT(4)    ! Kq  (rpm/Volt)  motor constant
       TAU     = PARMOT(5)    ! tau
       ZLOADI1 = PARMOT(6)    ! Io1
       ZLOADI2 = PARMOT(7)    ! Io2
       RMOTOR2 = PARMOT(8)    ! R2  (Ohms/Amp^2)
C
C----- default case for Kq is Kq=Kv
       IF(KQRPM .EQ. 0.0) KQRPM = KVRPM
C
       KVRAD = KVRPM * PI/30.0
       KQRAD = KQRPM * PI/30.0
C
       VM       = (1.0 + TAU*OMEGA    )*OMEGA/KVRAD
       VM_OMEGA = (1.0 + TAU*OMEGA*2.0)      /KVRAD
C
       I = (VOLT - VM)/RMOTOR0
       DO ITER = 1, 10
         RES = I*(RMOTOR0 + RMOTOR2*I**2) + VM - VOLT
         RES_I = RMOTOR0 + 3.0*RMOTOR2*I**2
         I = I - RES/RES_I
         IF(ABS(RES) .LT. EPS*MAX(1.0,ABS(VOLT))) GO TO 11
       ENDDO
       WRITE(*,*) 'MOTOR: Current convergence failed'
 11    CONTINUE
       RES_OMEGA = VM_OMEGA
       RES_VOLT  = -1.0
       I_OMEGA = -RES_OMEGA/RES_I
       I_VOLT  = -RES_VOLT /RES_I
C
       ZLOADI       = ZLOADI0 + ZLOADI1*OMEGA + ZLOADI2*OMEGA**2
       ZLOADI_OMEGA =           ZLOADI1       + ZLOADI2*OMEGA * 2.0
C
       Q       = (I       - ZLOADI      ) / KQRAD
       Q_OMEGA = (I_OMEGA - ZLOADI_OMEGA) / KQRAD
       Q_VOLT  =  I_VOLT                  / KQRAD
C
      ELSE
C----- Other motor models would go here
       WRITE(*,*) 'MOTORQ: Undefined motor type index:', IMOTYPE
       STOP
      ENDIF
C
      RETURN
      END ! MOTORQ



      SUBROUTINE VOLTM(  OMEGA,      Q, IMOTYPE, PARMOT,NMPAR, 
     &        VOLT, VOLT_OMEGA, VOLT_Q, 
     &        AMPS, AMPS_OMEGA, AMPS_Q )
C-----------------------------------------------------------
C     Motor+gearbox   Voltage(rpm,torque)  function.
C     Inverts MOTORQ's torque(rpm,Voltage)  function
C       via Newton iteration.
C
C Input:  OMEGA      output shaft rotation rate  (radians/s)
C         Q          output shaft torque   (N-m)
C         IMOTYPE    specifies type of motor model to be used
C         PARMOT(.)  motor parameters  (lines 3,4... in motor file)
C         NMPAR      number of motor parameters in motor file
C
C Output: VOLT       terminal voltage
C         VOLT_OMEGA dVOLT/dOMEGA  function derivative
C         VOLT_Q     dVOLT/dQ      function derivative
C         AMPS       current
C         AMPS_OMEGA dAMPS/dOMEGA  function derivative
C         AMPS_Q     dAMPS/dQ      function derivative
C
C-----------------------------------------------------------
      REAL PARMOT(*)
      REAL KVRPM, KQRPM, KVRAD, KQRAD

      DATA PI /3.1415926/
      DATA EPS /1.0E-6/
C
C---- Initial guess for Newton iteration, for each IMOTYPE
C-    May need to be different for different IMOTYPE.
      IF(IMOTYPE.EQ.1) THEN
       RMOTOR = PARMOT(1)    ! R   (Ohms)      motor resistance
       ZLOADI = PARMOT(2)    ! Io  (Amps)      zero-load current
       KVRPM  = PARMOT(3)    ! Kv  (rpm/Volt)  motor constant
       KVRAD = KVRPM*PI/30.0
       KQRAD = KVRAD
       AMPS = Q*KQRAD + ZLOADI
       VOLT = AMPS*RMOTOR + OMEGA/KVRAD
      ELSEIF(IMOTYPE.EQ.2) THEN
       RMOTOR = PARMOT(1)    ! R   (Ohms)      motor resistance
       ZLOADI = PARMOT(2)    ! Io  (Amps)      zero-load current
       KVRPM  = PARMOT(3)    ! Kv  (rpm/Volt)  motor constant
       KVRAD = KVRPM*PI/30.0
       KQRAD = KVRAD
       AMPS = Q*KQRAD + ZLOADI
       VOLT = AMPS*RMOTOR + OMEGA/KVRAD
      ELSE
C----- default case
       VOLT = 1.0
      ENDIF
C
      DO ITER = 1, 20
        CALL MOTORQ(OMEGA,VOLT, IMOTYPE, PARMOT,NMPAR, 
     &              QM,QM_OMEGA,QM_VOLT,
     &              AM,AM_OMEGA,AM_VOLT )
C
        RES      = QM - Q
        RES_VOLT = QM_VOLT
C
        DVOLT = -RES/RES_VOLT
ccc        write(*,*) iter, dvolt
        IF(ABS(DVOLT) .LT. EPS*MAX(1.0,ABS(VOLT))) GO TO 10
C
        VOLT = VOLT + DVOLT
      ENDDO
      WRITE(*,*) '** VOLTM: Voltage convergence failed. Res =', RES
C
 10   CONTINUE
      RES_OMEGA = QM_OMEGA
      RES_Q     = -1.0
C
      VOLT_OMEGA = -RES_OMEGA / RES_VOLT
      VOLT_Q     = -RES_Q     / RES_VOLT
C
      AMPS_OMEGA = AM_VOLT*VOLT_OMEGA + AM_OMEGA
      AMPS_Q     = AM_VOLT*VOLT_Q
C
      RETURN
      END ! VOLTM

