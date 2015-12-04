!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CALC_THERMO_DES                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_THERMO_DES

      USE compar
      USE des_rxns
      USE des_thermo
      USE derived_types, only: pic
      USE discretelement
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE interpolation
      USE param1
      USE run

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! Index of neighbor particle of particle I such that I < J
      INTEGER IJK
! Loop index for particles.
      INTEGER NP, lNP
! Phase index for particle NP
      INTEGER M
! Identifies that the indicated particle is of interest for debugging
      LOGICAL FOCUS
! Variables needed for calculating new interpolation quantities for
! species and energy equations
      INTEGER INTERP_IJK(2**3)
      DOUBLE PRECISION INTERP_WEIGHTS(2**3)

! Functions
!---------------------------------------------------------------------//

! This is a quick work-around to keep the thermo routines from causes
! issues while the "check_data" routines are rewritten. Moving forward
! this routine should be split apart to avoid the particle loops for
! cold-flow, non-reacting cases.
      IF(.NOT.ENERGY_EQ .AND. .NOT.ANY_SPECIES_EQ) RETURN

! Loop over fluid cells.
!---------------------------------------------------------------------//
      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE IJK_LP
         IF(PINC(IJK) == 0) CYCLE IJK_LP

! Interpolation: Removed J.Musser 11/8/2012
!---------------------------------------------------------------------//
!     IF(DES_INTERP_ON .AND. (ANY_SPECIES_EQ .OR. DES_CONV_EQ)) THEN
!         INTERP_IJK(:) = -1
!         INTERP_WEIGHTS(:) = ZERO
!         CALL INTERPOLATE_CC(NP, INTERP_IJK, INTERP_WEIGHTS, FOCUS)
!      ENDIF

! Preform user-defined calculations from fluid grid.
         IF(CALL_USR) CALL USR4_DES(IJK)

! Loop over all particles in cell IJK.
!---------------------------------------------------------------------//
         lNP_LP: DO lNP = 1, PINC(IJK)
            NP = PIC(IJK)%p(lNP)

! Skip indices that do not represent particles
            IF(.NOT.IS_NORMAL(NP)) CYCLE lNP_LP

! Reset the debug flag
            FOCUS = .FALSE.

! Calculate time dependent physical properties
            CALL DES_PHYSICAL_PROP(NP, FOCUS)

! Identify the solid phases of each particle
            M = PIJK(NP,5)

! calculate heat transfer via convection
            IF(CALC_CONV_DES) CALL DES_CONVECTION(NP, M, IJK, &
               INTERP_IJK, INTERP_WEIGHTS, FOCUS)

! calculate heat transfer via radiation
            IF(CALC_RADT_DES(M)) CALL DES_RADIATION(NP, M, IJK, FOCUS)

! Calculate reaction rates and interphase mass transfer
            IF(ANY_SPECIES_EQ) CALL DES_RRATES0(NP, M, IJK, &
               INTERP_IJK, INTERP_WEIGHTS, FOCUS)

         ENDDO lNP_LP ! End loop over all particles
      ENDDO IJK_LP ! End loop over fluid cells

      END SUBROUTINE CALC_THERMO_DES
