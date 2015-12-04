!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_REACTION_MODEL                                     !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_REACTION_MODEL

      USE compar
      Use constant
      Use des_rxns
      Use des_thermo
      Use derived_types, only: pic
      Use discretelement
      USE geometry
      USE indices
      Use param, only: dimension_n_s
      use run, only: ANY_SPECIES_EQ, SPECIES_EQ
      use physprop, only: SMAX, NMAX
      use run, only: SOLVE_ROs
      use functions

      IMPLICIT NONE

! Passed variables
!-----------------------------------------------
! None

! Local variables
!-----------------------------------------------
! Index of neighbor particle of particle I such that I < J
      INTEGER IJK
! Index value of particle
      INTEGER lNP, NP
! index of solids phase and species
      INTEGER M, NN
! total rate of consumption/production of species (g/sec)
      DOUBLE PRECISION SUM_DES_R_sc, SUM_DES_Rs
! masses of species comprising the particle (g)
      DOUBLE PRECISION S_MASS( DIMENSION_N_s )

      LOGICAL ALL_GONE( DIMENSION_N_s )

      DOUBLE PRECISION, PARAMETER :: P43 = 4.0d0/3.0d0

      DOUBLE PRECISION dMdt, dXdt, dXMdt, dXMdt_DTs

! Logical for Adams-Bashfort integration.
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.

!---------------------------------------------------------------------//

      IF(.NOT.ANY_SPECIES_EQ) RETURN

! Loop over fluid cells.
!---------------------------------------------------------------------//
      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE IJK_LP
         IF(PINC(IJK) == 0) CYCLE IJK_LP

! Loop over all particles in cell IJK.
!---------------------------------------------------------------------//
         lNP_LP: DO lNP = 1, PINC(IJK)
            NP = PIC(IJK)%p(lNP)
! Skip indices that do not represent particles
            IF(IS_NONEXISTENT(NP)) CYCLE lNP_LP
! Skip indices that represent ghost particles
            IF(IS_GHOST(NP) .OR. IS_ENTERING_GHOST(NP) .OR. IS_EXITING_GHOST(NP)) CYCLE lNP_LP

! Set the particle phase index
            M = PIJK(NP,5) + SMAX

! Skip particles when not solving species equations
            IF(.NOT.SPECIES_EQ(M)) CYCLE lNP_LP

! Check to see if any of the reaction rate will consume more than the
! available species in the particle.
!---------------------------------------------------------------------//
! Initialize local variables
            SUM_DES_Rs = ZERO
            SUM_DES_R_sc = ZERO
            ALL_GONE(:) = .FALSE.

! Reset the flag
            DO NN=1,NMAX(M)
! Calculate the current mass of species NN in the particle
               S_MASS(NN) = DES_X_s(NP,NN) * PMASS(NP)
! Calculte the amount of species N mass that is produced/consumed.
!    dXMdt_DTs > 0 :: Net production in mass units
!    dXMdt_DTs < 0 :: Net consumption in mass units
               dXMdt_DTs = (DES_R_sp(NP,NN)-DES_R_sc(NP,NN))*DTSOLID
! Check to see if the amount of solids phase species N consumed will
! be larger than the amount of available solids.
               IF((S_MASS(NN) + dXMdt_DTs) < ZERO) THEN
! Indicate that all of species is consumed.
                  ALL_GONE(NN) = .TRUE.
! Limit the consumption rate so that only the amount of species mass
! that the particle contains is consumed.
                  dXMdt = S_MASS(NN)/DTSOLID
! Calculate the total rate of change in particle mass
                  SUM_DES_Rs = SUM_DES_Rs + dXMdt
! Calculate the total rate of consumption
                  SUM_DES_R_sc = SUM_DES_R_sc + dXMdt
               ELSE
! Calculate the total particle mass rate of change (mass/time)
                  SUM_DES_Rs = SUM_DES_Rs + (DES_R_sp(NP,NN)-DES_R_sc(NP,NN))
! Calculate the total rate of consumption (mass/time)
                  SUM_DES_R_sc = SUM_DES_R_sc + DES_R_sc(NP,NN)
               ENDIF
            ENDDO

! Update the particle's mass.
! The mass of the particle is updated first so that it can be used in
! updating the species mass percent of the particle.
!---------------------------------------------------------------------//
            dMdt = SUM_DES_Rs
! First-order method: Euler
            IF(INTG_EULER) THEN
               PMASS(NP) = PMASS(NP) + DTSOLID*dMdt
! Second-order Adams-Bashforth scheme
            ELSEIF(INTG_ADAMS_BASHFORTH) THEN
               IF(FIRST_PASS) dMdt_OLD(NP) = dMdt
               PMASS(NP) = PMASS(NP) + DTSOLID * &
                  (1.50d0*dMdt - 0.5d0*dMdt_OLD(NP))
! Store the current value as the old value.
               dMdt_OLD(NP) = dMdt
            ENDIF

! Update the species mass percent
!---------------------------------------------------------------------//
            DO NN=1,NMAX(M)
               IF(ALL_GONE(NN))THEN
                  DES_X_s(NP,NN) = ZERO
               ELSE
                  dXdt = ((DES_R_sp(NP,NN) - DES_R_sc(NP,NN)) - &
                     DES_X_s(NP,NN) * SUM_DES_Rs)/PMASS(NP)
! First-order method: Euler
                  IF(INTG_EULER) THEN
                     DES_X_s(NP,NN) = DES_X_s(NP,nn) + DTSOLID*dXdt
! Second-order Adams-Bashforth scheme
                  ELSEIF(INTG_ADAMS_BASHFORTH) THEN
                     IF(FIRST_PASS) dXdt_OLD(NP,nn) = dXdt
                     DES_X_s(NP,nn) = DES_X_s(NP,nn) + DTSOLID * &
                        (1.50d0*dXdt - 0.5d0*dXdt_OLD(NP,nn))
! Store the current value as the old value.
                     dXdt_OLD(NP,nn) = dXdt
                  ENDIF
               ENDIF
            ENDDO

! Variable density
!---------------------------------------------------------------------//
            IF(SOLVE_ROs(M)) THEN
! update variables
               RO_Sol(NP) = PMASS(NP)/PVOL(NP)

! Shrinking particle
!---------------------------------------------------------------------//
            ELSE
               DES_RADIUS(NP) = &
                  (PMASS(NP)/(P43*Pi*Ro_Sol(NP)))**(1.0d0/3.0d0)
! update variables
               PVOL(NP) = PMASS(NP)/Ro_Sol(NP)
            ENDIF

! Update one over the particle's moment of inertia
            OMOI(NP) = 5.0d0 / (2.0d0 * PMASS(NP) * DES_RADIUS(NP)**2)

         ENDDO lNP_LP ! End loop over all particles
      ENDDO IJK_LP ! End loop over fluid cells

! Clear the necessary variables.
      DES_R_sp(:,:) = ZERO
      DES_R_sc(:,:) = ZERO

! Flag that the first pass is over
      FIRST_PASS = .FALSE.

      RETURN
      END SUBROUTINE DES_REACTION_MODEL
