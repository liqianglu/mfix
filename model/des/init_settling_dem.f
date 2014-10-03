!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INIT_SETTLING_DEM

!      USE param1
!      USE funits
      USE run
!      USE compar
      USE discretelement
!      USE cutcell
!      use desmpi
!      use mpi_utility
!      USE geometry
!      USE des_rxns
!      USE des_thermo
!      USE des_stl_functions

      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: FACTOR

!-----------------------------------------------
! Include statement functions
!-----------------------------------------------


! Skip this routine if there are no particles.
      IF(PARTICLES == 0) RETURN
! Skip this routine if not a new run.
      IF(RUN_TYPE /= 'NEW') RETURN
! Skip if not coupled.
      IF(.NOT.DES_CONTINUUM_COUPLED) RETURN
! Skip if using cohesion. (Why?)
      IF(USE_COHESION) RETURN

      WRITE(ERR_MSG, 1100) trim(iVal(NFACTOR))
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
 1100 FORMAT('Begining DEM settling period: ',A,' steps.')


! Disable the coupling flag.
      DES_CONTINUUM_COUPLED = .FALSE.


      DO FACTOR = 1, NFACTOR
! calculate forces
         CALL CALC_FORCE_DEM
! update particle position/velocity
         CALL CFNEWVALUES
! set the flag do_nsearch before calling particle in cell (for mpi)
         DO_NSEARCH = (MOD(FACTOR,NEIGHBOR_SEARCH_N)==0)
! find particles on grid
         CALL PARTICLES_IN_CELL
! perform neighbor search
         IF(DO_NSEARCH) CALL NEIGHBOUR
      ENDDO

! Reset the comoupling flag.
      DES_CONTINUUM_COUPLED = .TRUE.

      WRITE(ERR_MSG, 1200)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
 1200 FORMAT('DEM settling period complete.')

! Calculate the average solids temperature in each fluid cell
      CALL SET_INIT_avgTs

! this write_des_data is needed to properly show the initial state of
! the simulation (granular or coupled). In the coupled case, the
! particles may have 'settled' according to above.  In the granular
! case, the initial state won't be written until after the particles
! have moved without this call.
      IF(PRINT_DES_DATA) CALL WRITE_DES_DATA

      RETURN
      END SUBROUTINE INIT_SETTLING_DEM

