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
      SUBROUTINE INITPOLAR(POLARFILE, POLAR)
      USE POLARDATAMODULE
C----------------------------------------------------------------
C     Reads polar data in POLARFILE to a POLAR structure
C     defined type
C----------------------------------------------------------------
C     
      CONTINUE
      RETURN
      END
      
      SUBROUTINE GETPOLARCLINFO(POLAR,A,CL,DCLDA,CL0,STALL)
      USE POLARDATAMODULE
C----------------------------------------------------------------
C     Returns CL-related data from polar given AoA angle
C     
C----------------------------------------------------------------
C     
      CONTINUE
      RETURN
      END

      SUBROUTINE GETPOLARCDINFO(POLAR,A,CD,CD_CL)
      USE POLARDATAMODULE
C----------------------------------------------------------------
C     Returns CD-related data from polar given AoA angle
C     
C----------------------------------------------------------------
C     
      CONTINUE
      RETURN
      END
      

