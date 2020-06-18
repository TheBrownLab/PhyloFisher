* COPYRIGHT (c) 1988 AEA Technology
*######DATE 4 Feb 1993
C       Toolpack tool decs employed.
C       SAVE statements added.
C 16th October 2002: STOP and WRITE statements removed.
C
      INTEGER FUNCTION ID05A(INUM)
C----------------------------------------------------------------
C  Integer constants for: IBM Mainframe
C  single precision (4-byte arithmetic)
C
C  Obtained from H.S.L. subroutine ZE02AM.
C  Nick Gould and Sid Marlow, Harwell Laboratory, July 1988.
C----------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER INUM
C     ..
C     .. Local Arrays ..
      INTEGER IC(10)
C     ..
C     .. Save statement ..
      SAVE IC
C     ..
C     .. Data statements ..
C  IC(1) THE BASE (RADIX) OF THE FLOATING-POINT ARITHMETIC.
C  IC(2) THE NUMBER OF BASE IC(1) DIGITS IN THE SIGNIFICAND.
C  IC(3) THE NUMBER OF BITS USED FOR THE EXPONENT
C  IC(4) = 0 FLOATING-POINT ADDITION CHOPS, = 1 IT ROUNDS.
C  IC(5) = 0 IEEE-STANDARD ROUNDING IS NOT USED, = 1 IT IS.
C  IC(6) LARGEST -VE INTEGER:1.0 + REAL(IC(1))**IC(6) > 1.0.
C  IC(7) LARGEST -VE INTEGER:1.0 - REAL(IC(1))**IC(7) < 1.0.
C  IC(8) LARGEST -VE INTEGER: REAL(IC(1))**IC(8) > 0.0.
C  IC(9) LARGEST -VE INTEGER: REAL(IC(1))**IC(9) IS NORMAL.
C  IC(10) LARGEST +VE INTEGER: REAL(IC(1))**IC(10) FINITE.
C
      DATA IC(1)/16/
      DATA IC(2)/6/
      DATA IC(3)/7/
      DATA IC(4)/0/
      DATA IC(5)/0/
      DATA IC(6)/-5/
      DATA IC(7)/-6/
      DATA IC(8)/-65/
      DATA IC(9)/-65/
      DATA IC(10)/62/
C     ..
C     .. Executable Statements ..

      IF ( INUM .LE. 0 ) THEN
        ID05A = IC( 1 )
      ELSE IF ( INUM .GE. 11 ) THEN
        ID05A = IC( 10 )
      ELSE
        ID05A = IC( INUM )
      ENDIF
      RETURN
      END
