C***********************************************************************
C    Module:  bnsolv.f
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

      SUBROUTINE BNSOLV(A,B,C,R,NB,N,NRHS,NRMAX)
      DIMENSION A(NB,NB,N), B(NB,NB,N), C(NB,NB,N), R(NB,NRMAX,N)
C     **********************************************************************
C      This routine solves an N-long block-tridiagonal system with NBxNB
C      blocks, and with NRHS righthand sides by a standard block elimination
C      scheme.  The solutions are returned in the Rj vectors.
C
C      |A C      ||d|   |R..   |
C      |B A C    ||d|   |R..   |
C      |  B . .  ||.| = |R.....|
C      |    . . C||.|   |R..   |
C      |      B A||d|   |R..   |
C                                                  Mark Drela   10 June 89
C     **********************************************************************
C
CCC** Forward sweep: Elimination of lower block diagonal (B's).
      DO 1 I=1, N
C
        IM = I-1
C
C------ don't eliminate first B block because it doesn't exist
        IF(I.EQ.1) GO TO 12
C
C------ eliminate Bi block, thus modifying Ai and Ci blocks
        DO 11 K=1, NB
          DO 111 J=1, NB
            BTMP = B(K,J,I)
            DO 1111 L=1, NB
              A(K,L,I) = A(K,L,I) - BTMP*C(J,L,IM)
 1111       CONTINUE
            DO 1112 L=1, NRHS
              R(K,L,I) = R(K,L,I) - BTMP*R(J,L,IM)
 1112       CONTINUE
  111     CONTINUE
   11   CONTINUE
C
C                                                              -1
CCC---- multiply Ci block and righthand side Ri vectors by (Ai)
C       using Gaussian elimination.
C
   12   DO 13 KPIV=1, NB-1
          KP1 = KPIV+1
C
C-------- find max pivot index KX
          KX = KPIV
          DO 131 K=KP1, NB
            IF(ABS(A(K,KPIV,I))-ABS(A(KX,KPIV,I))) 131,131,1311
 1311        KX = K
  131     CONTINUE
C
          IF(A(KX,KPIV,I).EQ.0.0) THEN
           WRITE(*,*) 'Singular A block, i = ',I
           STOP
          ENDIF
C
          PIVOT = 1.0/A(KX,KPIV,I)
C
C-------- switch pivots
          A(KX,KPIV,I) = A(KPIV,KPIV,I)
C
C-------- switch rows & normalize pivot row
          DO 132 L=KP1, NB
            TEMP = A(KX,L,I)*PIVOT
            A(KX,L,I) = A(KPIV,L,I)
            A(KPIV,L,I) = TEMP
  132     CONTINUE
C
          DO 133 L=1, NB
            TEMP = C(KX,L,I)*PIVOT
            C(KX,L,I) = C(KPIV,L,I)
            C(KPIV,L,I) = TEMP
  133     CONTINUE
C
          DO 134 L=1, NRHS
            TEMP = R(KX,L,I)*PIVOT
            R(KX,L,I) = R(KPIV,L,I)
            R(KPIV,L,I) = TEMP
  134     CONTINUE
C
C-------- forward eliminate everything
          DO 135 K=KP1, NB
            ATMP = A(K,KPIV,I)
            DO 1351 L=KP1, NB
              A(K,L,I) = A(K,L,I) - ATMP*A(KPIV,L,I)
 1351       CONTINUE
            DO 1352 L=1, NB
              C(K,L,I) = C(K,L,I) - ATMP*C(KPIV,L,I)
 1352       CONTINUE
            DO 1353 L=1, NRHS
              R(K,L,I) = R(K,L,I) - ATMP*R(KPIV,L,I)
 1353       CONTINUE
  135     CONTINUE
C
   13   CONTINUE
C
C------ solve for last row
        IF(A(NB,NB,I).EQ.0.0) THEN
         WRITE(*,*) 'Singular A block, i = ',I
         STOP
        ENDIF
        PIVOT = 1.0/A(NB,NB,I)
        DO 141 L=1, NB
          C(NB,L,I) = C(NB,L,I)*PIVOT
  141   CONTINUE
        DO 142 L=1, NRHS
          R(NB,L,I) = R(NB,L,I)*PIVOT
  142   CONTINUE
C
C------ back substitute everything
        DO 15 KPIV=NB-1, 1, -1
          KP1 = KPIV+1
          DO 151 K=KP1, NB
            ATMP = A(KPIV,K,I)
            DO 1511 L=1, NB
              C(KPIV,L,I) = C(KPIV,L,I) - ATMP*C(K,L,I)
 1511       CONTINUE
            DO 1512 L=1, NRHS
              R(KPIV,L,I) = R(KPIV,L,I) - ATMP*R(K,L,I)
 1512       CONTINUE
  151     CONTINUE
   15   CONTINUE
    1 CONTINUE
C
CCC** Backward sweep: Back substitution using upper block diagonal (Ci's).
      DO 2 I=N-1, 1, -1
        IP = I+1
        DO 21 L=1, NRHS
          DO 211 K=1, NB
            DO 2111 J=1, NB
              R(K,L,I) = R(K,L,I) - R(J,L,IP)*C(K,J,I)
 2111       CONTINUE
  211     CONTINUE
   21   CONTINUE
    2 CONTINUE
C
      RETURN
      END ! BNSOLV
