************************************************************************
      DOUBLE PRECISION FUNCTION URAND2(IX)
*      FUNCTION URAND2(IX)
************************************************************************
* UNIFORM RANDOM NUMBER GENERATOR (MULTIPLICATIVE CONGRUENTIAL METHOD) *
*    FAST BUT NOT COMPLETELY PORTABLE.  INTEGERS SHOULD BE 32 BIT LONG.*
*    OVERFLOW IN INTERGER ARITHMETIC SHOULD NOT BE SENSED.             *
*    BITWISE AND FUNCTION 'IAND(IA,IB)' SHOULD BE SUPPORTED.           *
* PARAMETERS                                                           *
*   (1) N      (I) THE NUMBER OF RANDOM NUMBERS TO BE GENERATED        *
*                  (INPUT)                                             *
*   (2) X      (D) UNIFORM RANDOM NUMBERS (OUTPUT)                     *
*   (3) IR     (I) THE INITIAL SEED  (INPUT)                           *
*                  THE SEED FOR THE NEXT CALL (OUTPUT)                 *
* COPYRIGHT: Y. OYANAGI, JUNE 30, 1989  V.1                            *
************************************************************************
       Implicit Double Precision (A-H,O-Z)
*       COMMON/RANDOM1/IX
*       common/rancheck/NRancheck,NRanspc(1000000)
       DOUBLE PRECISION XX, INVM
       PARAMETER (LAMBDA = 48828125, MASK=2**30+(2**30-1) )
       PARAMETER (INVM = 0.5D0 ** 31)
*PAREMETER CHECK
C      IF( N .LE. 0) THEN
C       WRITE(6,*) '(SUBR.URAND2) PARAMETER ERROR. N = ', N
C       WRITE(6,*) 'RETURN WITH NO FURTHER CALCULATION.'
C       RETURN
C      END IF
      IF( MOD(IX, 2) .NE. 1) THEN
       WRITE(6,*) '(SUBR.URAND2) PARAMETER ERROR. IX = ', IX,
     >            ' IS EVEN.'
       WRITE(6,*) 'RETURN WITH NO FURTHER CALCULATION.'
C       RETURN
      END IF
*MAIN LOOP
C      DO 10 I = 1, N
       IX = IAND(LAMBDA * IX, MASK)
       XX = IX * INVM
C 10     CONTINUE
      URAND2 = XX
*      NRancheck=NRancheck+1
      RETURN
      END

