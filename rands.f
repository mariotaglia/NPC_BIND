   
C****************************************************************
C **********************************************************************
        double precision FUNCTION RANDS (SEED)
C **********************************************************************

C-----  this is a special function for random number generation
C        on 32-bit machines that do not support long integer
C        multiplication and truncation.  the technique used is to do
C        the multiplication and addition in parts, by splitting all
C       integers in a 'high' and a 'low' part.  the algorithm is
C-----        exact, and should give machine-independent results.

C-----        the algorithm implemented is (following d.e. knuth):
C        seed = seed*1592653589 + 453816691
C        if (seed.lt.0) seed = seed + 1 + 2147483647
C-----        note that 1592653589 = 48603*2**15 + 30485

C 32768 = 2^15, 65536 = 2^16, 4.65661287308E-10 = 2^(-31)

        INTEGER SEED, I1, I2

        I2 = MOD (SEED, 32768) * 30485
        I1 = MOD (SEED / 32768 * 30485, 65536) + MOD (MOD (SEED, 32768)
     X    * 48603, 65536) + I2 / 32768 + MOD (SEED / 32768, 2) * 32768 +
     X     13849
        I2 = MOD (I2, 32768) + 12659
        I1 = I1 + I2 / 32768
        I2 = MOD (I2, 32768)
        I1 = MOD (I1, 65536)
        SEED = I1 * 32768 + I2
        RANDS = REAL(I1 * 32768 + I2) * 4.65661287308E-10

        RETURN
        END


