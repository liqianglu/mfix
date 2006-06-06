!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFRELVEL(L, II, VRl)                                   C
!  Purpose: DES - Calculate relative velocity between a particle pair  C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Comments: Relative (translational) velocity required                C
!  for Eqn 6  from the following paper                                 C
!  Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical         C
!  simulation of plug glow of cohesionless particles in a              C
!  horizontal pipe", Powder technology, 71, 239-250, 1992              C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFRELVEL(L, II, VRl)
      
      USE discretelement
      USE param1
      IMPLICIT NONE

      INTEGER L, KK, II
      DOUBLE PRECISION VRl(DIMN)

!-----------------------------------------------------------------------


      VRl(:) = (DES_VEL_NEW(L,:) - DES_VEL_NEW(II,:))

      RETURN
      END SUBROUTINE CFRELVEL


