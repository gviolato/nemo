C***********************************************************************
C    Module:  io.f
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

      SUBROUTINE FREAD(LU,LINE,ILINE,IERR, CVAL)
      CHARACTER*(*) LINE, CVAL
C
 10   CONTINUE
      ILINE = ILINE + 1
C
      READ(LU,1000,END=90) LINE
 1000 FORMAT(A)
C
      IF(INDEX('#!',LINE(1:1)) .NE. 0) GO TO 10
      IF(LINE.EQ.' ') GO TO 10
C
      KB = INDEX(LINE,'!') - 1
      IF(KB.LE.0) KB = LEN(LINE)
C
      CVAL = LINE(1:KB)
      IERR = 0
      RETURN
C
 80   IERR = 1
      RETURN
C
 90   IERR = -1
      RETURN
      END


      SUBROUTINE RREAD(LU,LINE,ILINE,IERR, NVAL,VAL)
      CHARACTER*(*) LINE
      REAL VAL(*)
C
      LOGICAL ERROR
C
 10   CONTINUE
      ILINE = ILINE + 1
C
      READ(LU,1000,END=90) LINE
 1000 FORMAT(A)
C
      IF(INDEX('#!',LINE(1:1)) .NE. 0) GO TO 10
      IF(LINE.EQ.' ') GO TO 10
C
      KB = INDEX(LINE,'!') - 1
      IF(KB.LE.0) KB = LEN(LINE)
C
      CALL GETFLT(LINE(1:KB),VAL,NVAL,ERROR)
      IF(ERROR) THEN
       IERR = 1
      ELSE
       IERR = 0
      ENDIF
      RETURN
C
 90   IERR = -1
      RETURN
      END



      SUBROUTINE IREAD(LU,LINE,ILINE,IERR, NVAL,IVAL)
      CHARACTER*(*) LINE
      INTEGER IVAL(*)
C
      LOGICAL ERROR
C
 10   CONTINUE
      ILINE = ILINE + 1
C
      READ(LU,1000,END=90) LINE
 1000 FORMAT(A)
C
      IF(INDEX('#!',LINE(1:1)) .NE. 0) GO TO 10
      IF(LINE.EQ.' ') GO TO 10
C
      KB = INDEX(LINE,'!') - 1
      IF(KB.LE.0) KB = LEN(LINE)
C
      CALL GETINT(LINE(1:KB),IVAL,NVAL,ERROR)
      IF(ERROR) THEN
       IERR = 1
      ELSE
       IERR = 0
      ENDIF
      RETURN
C
 90   IERR = -1
      RETURN
      END



      SUBROUTINE READPOLNAMES(LU,LINE,ILINE,IERR,POLFILES,NPOLS)
      CHARACTER*(*) LINE
C----------------------------------------------------------------
C     Parses character string INPUT into an array
C     of strings POLFILES(1 NPOLS)
C
C     Each filename must be separated by string
C
C     NPOLS returns how many numbers were actually extracted.
C----------------------------------------------------------------
C     
 10   CONTINUE
      ILINE = ILINE + 1
C     
      READ(LU,1000,END=90) LINE
 1000 FORMAT(A)
C     
      IF(INDEX('#!',LINE(1:1)) .NE. 0) GO TO 10
      IF(LINE.EQ.' ') GO TO 10
C
C      
C     Magic happens here
C
      
      IERR = 0
      RETURN
C     
 80   IERR = 1
      RETURN
C     
 90   IERR = -1
      RETURN
      END


      SUBROUTINE GETFLT(INPUT,RNUM,NR,ERROR)
      CHARACTER*(*) INPUT
      REAL RNUM(*)
      LOGICAL ERROR
C----------------------------------------------------------------
C     Parses character string INPUT into an array
C     of real numbers returned in RNUM(1..NR).
C
C     Will attempt to extract no more than NR numbers,
C     unless NR = 0, in which case all numbers present 
C     in INPUT will be extracted.
C
C     NR returns how many numbers were actually extracted.
C----------------------------------------------------------------
      CHARACTER*1 TAB
C
      TAB = CHAR(9)
C
C---- number of characters to be examined
      ILEN = LEN(INPUT)
C
C---- ignore everything after a "!" character
      K = INDEX(INPUT,'!')
      IF(K.GT.0) ILEN = K-1
C
C---- set limit on numbers to be read
      NINP = NR
      IF(NINP.EQ.0) NINP = ILEN/2 + 1
C
      NR = 0
C
      IF(ILEN.EQ.0) RETURN
C
C---- extract numbers
      N = 0
      K = 1
      DO 10 IPASS=1, ILEN
C------ find next blank or tab (pretend there's one after the end of the string)
        KBLK = INDEX(INPUT(K:ILEN),' ') + K - 1
        KTAB = INDEX(INPUT(K:ILEN),TAB) + K - 1
C
        IF(KBLK.EQ.K-1) KBLK = ILEN + 1
        IF(KTAB.EQ.K-1) KTAB = ILEN + 1
C
        KSPACE = MIN( KBLK , KTAB )
C
        IF(KSPACE.EQ.K) THEN
C------- just skip this space
         K = K+1
         GO TO 9
        ENDIF
C
C------ also find next comma
        KCOMMA = INDEX(INPUT(K:ILEN),',') + K - 1
        IF(KCOMMA.EQ.K-1) KCOMMA = ILEN + 1
C
C------ space is farther down, so we ran into something...
        N = N+1
C
C------ bug out early if no more numbers are to be read
        IF(N.GT.NINP) GO TO 11
C
C------ set ending delimiter position for this number
        KDELIM = MIN(KSPACE,KCOMMA)
C
        IF(K.EQ.KDELIM) THEN
C------- nothing but a comma... just keep looking
         K = K+1
         GO TO 9
        ENDIF
C
C------ whatever we have, it is in substring K:KEND
        KEND = KDELIM - 1
        READ(INPUT(K:KEND),*,ERR=20) RNUM(N)
        NR = N
C
C------ keep looking after delimiter
        K = KDELIM + 1
C
  9     IF(K.GE.ILEN) GO TO 11
 10   CONTINUE
C
C---- normal return
 11   CONTINUE
      ERROR = .FALSE.
      RETURN
C
C---- bzzzt !!!
 20   CONTINUE
ccc   WRITE(*,*) 'GETFLT: List-directed read error.'
      ERROR = .TRUE.
      RETURN
      END ! GETFLT



      SUBROUTINE GETINT(INPUT,A,N,ERROR)
      CHARACTER*(*) INPUT
      INTEGER A(*)
      LOGICAL ERROR
C----------------------------------------------------------
C     Parses character string INPUT into an array
C     of integer numbers returned in A(1...N)
C
C     Will attempt to extract no more than N numbers, 
C     unless N = 0, in which case all numbers present
C     in INPUT will be extracted.
C
C     N returns how many numbers were actually extracted.
C----------------------------------------------------------
      CHARACTER*130 REC
C
C---- only first 128 characters in INPUT will be parsed
      ILEN = MIN( LEN(INPUT) , 128 )
      ILENP = ILEN + 2
C
C---- put input into local work string (which will be munched)
      REC(1:ILENP) = INPUT(1:ILEN) // ' ,'
C
C---- ignore everything after a "!" character
      K = INDEX(REC,'!')
      IF(K.GT.0) REC(1:ILEN) = REC(1:K-1)
C
      NINP = N
C
C---- count up how many numbers are to be extracted
      N = 0
      K = 1
      DO 10 IPASS=1, ILEN
C------ search for next space or comma starting with current index K
        KSPACE = INDEX(REC(K:ILENP),' ') + K - 1
        KCOMMA = INDEX(REC(K:ILENP),',') + K - 1
C
        IF(K.EQ.KSPACE) THEN
C------- just skip this space
         K = K+1
         GO TO 9
        ENDIF
C
        IF(K.EQ.KCOMMA) THEN
C------- comma found.. increment number count and keep looking
         N = N+1
         K = K+1
         GO TO 9
        ENDIF
C
C------ neither space nor comma found, so we ran into a number...
C-    ...increment number counter and keep looking after next space or comma
        N = N+1
        K = MIN(KSPACE,KCOMMA) + 1
C
  9     IF(K.GE.ILEN) GO TO 11
 10   CONTINUE
C
C---- decide on how many numbers to read, and go ahead and read them
 11   IF(NINP.GT.0) N = MIN( N, NINP )
      READ(REC(1:ILEN),*,ERR=20) (A(I),I=1,N)
      ERROR = .FALSE.
      RETURN
C
C---- bzzzt !!!
 20   CONTINUE
ccc   WRITE(*,*) 'GETINT: String-to-integer conversion error.'
      N = 0
      ERROR = .TRUE.
      RETURN
      END





      SUBROUTINE PPARSE(LINE,P1,P2,NP,IERR)
      CHARACTER*(*) LINE
C----------------------------------------------------
C     Extracts sequence parameters from string LINE,
C     which must have one of the following formats:
C
C       P1
C       P1,P2
C       P1,P2,DP
C       P1,P2/NP
C----------------------------------------------------
C
      N = LEN(LINE)
C
      KCOMMA1 = INDEX(LINE,',')
      KSLASH  = INDEX(LINE,'/')
      IF(KCOMMA1.GT.0 .AND. KSLASH.EQ.0) THEN
       K = KCOMMA1 + 1
       KCOMMA2 = INDEX(LINE(K:N),',') + K - 1
      ELSE
       KCOMMA2 = 0
      ENDIF
C
C
      IERR = 0
C
      IF(KCOMMA1.EQ.0) THEN
C----- only a single parameter is present
       READ(LINE,*,ERR=80,END=90) P1
       P2 = P1
       NP = 1
       RETURN
C
      ELSE
C----- two or three sequence parameters
       K = KCOMMA1 - 1
       READ(LINE(1:K),*,ERR=80,END=90) P1
C
       IF(KCOMMA2.EQ.0 .AND. KSLASH.EQ.0) THEN
C------ no step size parameter given
        K = KCOMMA1 + 1
        READ(LINE(K:N),*,ERR=80,END=90) P2
        IF(P1.EQ.P2) THEN
C------- same two endpoints are given -- zero interval with no steps
         NP = 1
        ELSE
C------- only two endpoints are given -- set default number of steps
         NP = 5
        ENDIF
        RETURN
C
       ELSEIF(KCOMMA2.NE.0) THEN
C------ interval step size is given
        KA = KCOMMA1 + 1
        KB = KCOMMA2 - 1
        READ(LINE(KA:KB),*,ERR=80,END=90) P2
        K = KCOMMA2 + 1
        READ(LINE(K:N),*,ERR=80,END=90) PDEL
        RNUM = (P2-P1)/PDEL
        NP = INT( ABS(RNUM) + 1.49)
        RETURN
C         
       ELSEIF(KSLASH.NE.0) THEN
        KA = KCOMMA1 + 1
        KB = KSLASH  - 1

        READ(LINE(KA:KB),*,ERR=80,END=90) P2
        K = KSLASH + 1
        READ(LINE(K:N),*,ERR=80,END=90) RNUM
        NP = INT( ABS(RNUM) + 0.01 )
        RETURN
       ENDIF
C
      ENDIF
C
 80   CONTINUE
      IERR = +1
      RETURN
C
 90   CONTINUE
      IERR = -1
      RETURN
C
      END ! PPARSE


      SUBROUTINE GETARG0(IARG,ARG)
C------------------------------------------------
C     Same as GETARG, but...
C
C     ...in the case of Intel Fortran, this one
C     doesn't barf if there's no Unix argument 
C      (just returns blank string instead)
C------------------------------------------------
      CHARACTER*(*) ARG
C
      NARG = IARGC()
      IF(NARG.GE.IARG) THEN
       CALL GETARG(IARG,ARG)
      ELSE
       ARG = ' '
      ENDIF
C
      RETURN
      END ! GETARG0



      SUBROUTINE STRIP(STRING,NS)
      CHARACTER*(*) STRING
C---------------------------------------------------
C     Strips all blanks off string STRING
C     and returns non-blank part in STRING.
C     Returns length of the non-blank part in NS.
C---------------------------------------------------
      N = LEN(STRING)
C
C---- find last non-blank character
      DO NS = N, 1, -1
        IF(STRING(NS:NS).NE.' ') GO TO 11
      ENDDO
      NS = 0
      RETURN
C
C---- start blank-stripping loop
   11 CONTINUE
C
C---- find first blank character before NS
      DO K0 = 1, NS
        IF(STRING(K0:K0).EQ.' ') GO TO 21
      ENDDO
      RETURN
   21 CONTINUE
C
C---- find first non-blank character following first blank
      DO K1 = K0+1, NS
        IF(STRING(K1:K1).NE.' ') GO TO 31
      ENDDO
   31 CONTINUE
C
C---- shift STRING to remove blank segment
      NSDEL = K1 - K0
      STRING(K0:NS-NSDEL) = STRING(K1:NS)
      NS = NS - NSDEL
C
      GO TO 11
C
      END ! STRIP
