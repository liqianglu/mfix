!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFUPDATEOLD
!  Purpose: Update old arrays to store current time step information
!           before the new arrays are updated
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer:                                          Date:
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE CFUPDATEOLD

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      USE run
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters (no. particles)
      INTEGER LL
!-----------------------------------------------

!!$omp parallel do if(max_pip .ge. 10000) default(shared)         &
!!$omp private(ll)                    &
!!$omp schedule (guided,50)
      DO LL = 1, MAX_PIP
         IF(.NOT.PEA(LL,1) .or. pea(ll,4)) CYCLE

         DES_POS_OLD(:,LL)  = DES_POS_NEW(:,LL)
         DES_VEL_OLD(:,LL)  = DES_VEL_NEW(:,LL)
         OMEGA_OLD(:,LL)    = OMEGA_NEW(:,LL)

      ENDDO
!!$omp end parallel do

      RETURN
      END SUBROUTINE CFUPDATEOLD
