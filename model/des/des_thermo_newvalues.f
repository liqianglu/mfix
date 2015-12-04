!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_THERMO_NEWVALUES                                   !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_THERMO_NEWVALUES

      USE compar
      Use des_thermo
      Use des_rxns
      Use derived_types, only: pic
      Use discretelement
      USE geometry
      USE indices
      Use param1
      Use physprop
      use run, only: ENERGY_EQ
      use functions
      use funits, only: dmp_log

      IMPLICIT NONE

! Passed variables
!-----------------------------------------------
! NONE

! Local variables
!---------------------------------------------------------------------//
! Index of neighbor particle of particle I such that I < J
      INTEGER IJK
! Loop index for particles.
      INTEGER NP, lNP
! Logical for Adams-Bashfort integration.
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
! Sum of particle temperatures in fluid cell.
      DOUBLE PRECISION SUM_T_s
!---------------------------------------------------------------------//

      IF(.NOT.ENERGY_EQ) RETURN

! Second-order Adams-Bashforth scheme defaults to Euler on first pass.
      IF(FIRST_PASS) THEN
         IF(INTG_ADAMS_BASHFORTH) &
            Q_Source0(:) = Q_Source(:)/ (PMASS(:) * DES_C_ps(:))
         FIRST_PASS = .FALSE.
      ENDIF

! Clear the average solids temperature for all fluid cells.
      avgDES_T_s(:) = ZERO

! Loop over fluid cells.
!---------------------------------------------------------------------//
      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE IJK_LP
         IF(PINC(IJK) == 0) CYCLE IJK_LP

! Initialize local solids temperature.
         SUM_T_s = ZERO

! Loop over all particles in cell IJK.
!---------------------------------------------------------------------//
         lNP_LP: DO lNP = 1, PINC(IJK)
            NP = PIC(IJK)%p(lNP)
! Skip indices that do not represent particles
            IF(IS_NONEXISTENT(NP)) CYCLE lNP_LP
! Skip indices that represent ghost particles
            IF(IS_GHOST(NP) .OR. IS_ENTERING_GHOST(NP) .OR. IS_EXITING_GHOST(NP)) CYCLE lNP_LP
! Advance particle position, velocity
            IF (INTG_EULER) THEN

! First-order method
               DES_T_s_NEW(NP) = DES_T_s_NEW(NP) + &
                  DTSOLID*(Q_Source(NP) / (PMASS(NP) * DES_C_ps(NP)))
            ELSE
! Second-order Adams-Bashforth scheme
               DES_T_s_NEW(NP) = DES_T_s_OLD(NP) + &
                  (1.5d0 * (Q_Source(NP)/(PMASS(NP)*DES_C_ps(NP))) - &
                   0.5d0 * Q_Source0(NP)) * DTSOLID
               Q_Source0(NP) = Q_Source(NP) / (PMASS(NP) * DES_C_ps(NP))
            ENDIF

! Write out the debugging information.
            IF(DEBUG_DES) THEN
               IF(DMP_LOG) THEN
                  IF(NP == FOCUS_PARTICLE) THEN
                     WRITE(*,"(//5X,A)")'From: DES_THERMO_NEWVALUES -'
                     WRITE(*,"(8X,A,D13.6)")'Tp:  ',DES_T_s_NEW(NP)
                     WRITE(*,"(8X,A,D13.6)")'Tp0: ',DES_T_s_OLD(NP)
                     WRITE(*,"(8X,A,D13.6)")'Qsrc:',Q_Source(NP)
                     WRITE(*,"(5X,25('-')/)")
                  ENDIF
               ENDIF
            ENDIF

! Update the old temperature value
            DES_T_s_OLD(NP) = DES_T_s_NEW(NP)
! Update the sum of particle temperatures in fluid cell IJK.
            SUM_T_s = SUM_T_s + DES_T_s_NEW(NP)
         ENDDO lNP_LP ! End loop over all particles
! Average solids temperature in fluid cell IJK. The average method
! (over particles) will need changed for Hybrid model (area? volume?).
         avgDES_T_s(IJK) = SUM_T_s/PINC(IJK)

      ENDDO IJK_LP ! End loop over fluid cells

      Q_Source(:) = ZERO


      RETURN

      END SUBROUTINE DES_THERMO_NEWVALUES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: SET_INIT_avgTs                                         !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-NOV-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_INIT_avgTs

      USE compar
      Use des_thermo
      Use des_rxns
      Use derived_types, only: pic
      Use discretelement
      USE geometry
      USE indices
      Use param1
      Use physprop
      use run, only: ENERGY_EQ
      use functions
      IMPLICIT NONE

! Passed variables
!-----------------------------------------------
! NONE

! Local variables
!---------------------------------------------------------------------//
! Index of neighbor particle of particle I such that I < J
      INTEGER IJK
! Loop index for particles.
      INTEGER NP, lNP
! Sum of particle temperatures in fluid cell.
      DOUBLE PRECISION SUM_T_s
!---------------------------------------------------------------------//

      IF(.NOT.ENERGY_EQ) RETURN

! Clear the average solids temperature for all fluid cells.
      avgDES_T_s(:) = ZERO

! Loop over fluid cells.
!---------------------------------------------------------------------//
      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE IJK_LP
         IF(PINC(IJK) == 0) CYCLE IJK_LP

! Initialize local solids temperature.
         SUM_T_s = ZERO

! Loop over all particles in cell IJK.
!---------------------------------------------------------------------//
         lNP_LP: DO lNP = 1, PINC(IJK)
            NP = PIC(IJK)%p(lNP)
! Skip indices that do not represent particles
            IF(IS_NONEXISTENT(NP)) CYCLE lNP_LP
! Skip indices that represent ghost particles
            IF(IS_GHOST(NP) .OR. IS_ENTERING_GHOST(NP) .OR. IS_EXITING_GHOST(NP)) CYCLE lNP_LP
! Update the sum of particle temperatures in fluid cell IJK.
            SUM_T_s = SUM_T_s + DES_T_s_NEW(NP)
         ENDDO lNP_LP ! End loop over all particles
! Average solids temperature in fluid cell IJK. The average method
! (over particles) will need changed for Hybrid model (area? volume?).
         avgDES_T_s(IJK) = SUM_T_s/PINC(IJK)

      ENDDO IJK_LP ! End loop over fluid cells

      RETURN

      END SUBROUTINE SET_INIT_avgTs
