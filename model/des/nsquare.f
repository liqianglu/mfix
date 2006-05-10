!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: NSQUARE(PARTS)                                         C
!  Purpose: DES - N-Square neighbor search                             C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE NSQUARE(PARTS)

      USE param1
      USE discretelement
      IMPLICIT NONE

      INTEGER L, I, K, KK, II, TC1, TC2, TCR, TCM, PARTS, J
      DOUBLE PRECISION CT
      DOUBLE PRECISION DIST, R_LM

      IF (DO_NSQUARE) THEN
 !        CALL SYSTEM_CLOCK(TC1,TCR,TCM)
         DO L = 1, PARTICLES 
            DO I = 1, PARTICLES
               DIST = ZERO           
               IF (I.GT.L) THEN
                  DO II = 1, DIMN
                     DIST = DIST + (DES_POS_NEW(I,II)-DES_POS_NEW(L,II))**2 
                  END DO
                  DIST = SQRT(DIST)
                  R_LM = DES_RADIUS(L) + DES_RADIUS(I)
                  IF (DIST.LE.R_LM) THEN
                     NEIGHBOURS(L,1) = NEIGHBOURS(L,1)+1
                     NEIGHBOURS(I,1) = NEIGHBOURS(I,1)+1
                     K = NEIGHBOURS(L,1) 
                     KK = NEIGHBOURS(I,1)
                     IF (K.LE.MN) THEN
                        NEIGHBOURS(L,K+1) = I
                     ELSE 
                        PRINT *,'NSQUARE - NEIGHBORS GT MN'
                        PRINT *, L,':',(NEIGHBOURS(L,II), II=1,MAXNEIGHBORS) 
                        STOP
                     END IF
                     IF (KK.LE.MN) THEN
                        NEIGHBOURS(I,KK+1) = L
                     ELSE 
                        PRINT *,'NSQUARE - NEIGHBORS GT MN'
                        PRINT *, I,':',(NEIGHBOURS(I,II), II=1,MAXNEIGHBORS) 
                        STOP
                     END IF
                  END IF
               END IF
            END DO
         END DO

!         CALL SYSTEM_CLOCK(TC2,TCR,TCM)
!         CT = TC2-TC1
!         IF(CT.LE.0) THEN
!            CT = TC2 + TCM - TC1
!         END IF
!         CT = CT/TCR
!         N2CT = CT
!     PRINT *,'N2:- CPU TIME TAKEN:',N2CT
!      ELSE
!         N2CT = ZERO
      END IF

      RETURN
      END SUBROUTINE NSQUARE


