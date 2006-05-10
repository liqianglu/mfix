!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_VDW_WALL_INTERACTION(I,J)                        C
!  Purpose: This module calculates the attractive force between        C
!           a wall and a particle using the Hamaker van der Waals modelC
!                                                                      C
!   Author: Mike Weber                              Date: 9/3ZERO4      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


      SUBROUTINE CHECK_VDW_WALL_INTERACTION(I,J)

!-----MODULES USED
      USE param1
      USE discretelement
      IMPLICIT NONE

!-----LOCAL DUMMY VARIABLES
      INTEGER I,J,N,K, II
      INTEGER LINK
      INTEGER CHECK_LINK
      DOUBLE PRECISION RADIUS
      DOUBLE PRECISION FORCE, RELPOS(2)
      DOUBLE PRECISION DIST
      LOGICAL ALREADY_EXISTS

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START CHECK VAN DER WAALS WALL INTERACTION'
      END IF

      IF(J.eq.PARTICLES+1)THEN !West wall
         DES_POS_NEW(J,1)=WX1
         DES_POS_NEW(J,2)=DES_POS_NEW(I,2)
         IF(DIMN.EQ.3) DES_POS_NEW(J,3)=DES_POS_NEW(I,3)
         DES_VEL_NEW(J,1)=ZERO
         DES_VEL_NEW(J,2)=ZERO
         IF(DIMN.EQ.3) DES_VEL_NEW(J,3)=ZERO
      END IF         

      IF(J.eq.PARTICLES+2)THEN !Bottom wall
         DES_POS_NEW(J,1)=DES_POS_NEW(I,1)
         DES_POS_NEW(J,2)=BY1
         IF(DIMN.EQ.3) DES_POS_NEW(J,3)=DES_POS_NEW(I,3)
         DES_VEL_NEW(J,1)=ZERO
         DES_VEL_NEW(J,2)=ZERO
         IF(DIMN.EQ.3) DES_VEL_NEW(J,3)=ZERO
      END IF  

      IF(J.eq.PARTICLES+3)THEN !East wall
         DES_POS_NEW(J,1)=EX2
         DES_POS_NEW(J,2)=DES_POS_NEW(I,2)
         IF(DIMN.EQ.3) DES_POS_NEW(J,3)=DES_POS_NEW(I,3)
         DES_VEL_NEW(J,1)=ZERO
         DES_VEL_NEW(J,2)=ZERO
         IF(DIMN.EQ.3) DES_VEL_NEW(J,3)=ZERO
      END IF  

      IF(J.eq.PARTICLES+4)THEN !Top wall
         DES_POS_NEW(J,1)=DES_POS_NEW(I,1)
         DES_POS_NEW(J,2)=TY2
         IF(DIMN.EQ.3) DES_POS_NEW(J,3)=DES_POS_NEW(I,3)
         DES_VEL_NEW(J,1)=ZERO
         DES_VEL_NEW(J,2)=ZERO
         IF(DIMN.EQ.3) DES_VEL_NEW(J,3)=ZERO
      END IF  

      IF(J.eq.PARTICLES+5)THEN !North wall
         DES_POS_NEW(J,1)=DES_POS_NEW(I,1)
         DES_POS_NEW(J,2)=DES_POS_NEW(I,2)
         IF(DIMN.EQ.3) DES_POS_NEW(J,3)=NZ2
         DES_VEL_NEW(J,1)=ZERO
         DES_VEL_NEW(J,2)=ZERO
         IF(DIMN.EQ.3) DES_VEL_NEW(J,3)=ZERO
      END IF 

      IF(J.eq.PARTICLES+6)THEN !South wall
         DES_POS_NEW(J,1)=DES_POS_NEW(I,1)
         DES_POS_NEW(J,2)=DES_POS_NEW(I,2)
         IF(DIMN.EQ.3) DES_POS_NEW(J,3)=SZ1
         DES_VEL_NEW(J,1)=ZERO
         DES_VEL_NEW(J,2)=ZERO
         IF(DIMN.EQ.3) DES_VEL_NEW(J,3)=ZERO
      END IF 

      IF(DIMN.EQ.3) THEN
         RADIUS=SQRT((DES_POS_NEW(J,1)-DES_POS_NEW(I,1))**2+&
                (DES_POS_NEW(J,2)-DES_POS_NEW(I,2))**2+&
                (DES_POS_NEW(J,3)-DES_POS_NEW(I,3))**2)
      ELSE
         RADIUS=SQRT((DES_POS_NEW(J,1)-DES_POS_NEW(I,1))**2+&
                (DES_POS_NEW(J,2)-DES_POS_NEW(I,2))**2)
      END IF

      DIST=RADIUS-DES_RADIUS(I)

      IF(DIST.lt.WALL_VDW_OUTER_CUTOFF)THEN
         DO II=1,2
            RELPOS(II)=DES_POS_NEW(J,II)-DES_POS_NEW(I,II)
         END DO
         IF(DIST.gt.WALL_VDW_INNER_CUTOFF)THEN
            FORCE=WALL_HAMAKER_CONSTANT*DES_RADIUS(I)/(6*DIST*DIST)
         ELSE
            FORCE=4*3.14*WALL_SURFACE_ENERGY*DES_RADIUS(I)
         END IF !Long range or surface?

         DO K=1,2
           FC(I,K)=FC(I,K)+RELPOS(K)/RADIUS*FORCE
           FC(J,K)=FC(J,K)-RELPOS(K)/RADIUS*FORCE
         END DO 
                    
      END IF !Is particle within cutoff?


      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**END CHECK VAN DER WAALS WALL INTERACTION'
      END IF

      END SUBROUTINE CHECK_VDW_WALL_INTERACTION 
