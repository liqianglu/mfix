!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                         !
!  Subrourtine: DES_INIT_ARRAYS                                           !
!  Author: Jay Boyalakuntla                              Date: 12-Jun-04  !
!                                                                         !
!  Purpose: Initialize arrays at the start of the simulation. Note that   !
!  arrays based on the number of particles (MAX_PIP) should be added to   !
!  the DES_INIT_PARTICE_ARRAYS as they need to be reinitialized after     !
!  particle arrays are grown (see PARTICLE_GROW).                         !
!                                                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_INIT_ARRAYS

      USE param
      USE param1
      USE discretelement
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE des_bc
      USE run
      use desgrid
      use desmpi
      USE des_thermo
      USE des_rxns

      IMPLICIT NONE

      INTEGER :: II
!-----------------------------------------------


      PINC(:) = ZERO

      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO
      DES_ROP_S(:,:) = ZERO
      DES_ROP_SO(:,:) = ZERO

      P_FORCE(:,:) = ZERO

      IF(allocated(DRAG_AM)) DRAG_AM = ZERO
      IF(allocated(DRAG_BM)) DRAG_BM = ZERO

      F_GDS = ZERO
      VXF_GDS = ZERO

      IF (DES_CONTINUUM_HYBRID) THEN
         F_SDS = ZERO
         VXF_SDS = ZERO
         SDRAG_AM = ZERO
         SDRAG_BM = ZERO
      ENDIF

      GRAV(:) = ZERO

      DO II = 1, SIZE(particle_wall_collisions)
         nullify(particle_wall_collisions(II)%pp)
      ENDDO


      IF(ENERGY_EQ)THEN
         avgDES_T_s(:) = ZERO
         DES_ENERGY_SOURCE(:) = ZERO
      ENDIF

      CALL DES_INIT_PARTICLE_ARRAYS(1,MAX_PIP)

      RETURN
      END SUBROUTINE DES_INIT_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                         !
!  Subrourtine: DES_INIT_PARTICLE_ARRAYS                                  !
!  Author: Jay Boyalakuntla                              Date: 12-Jun-04  !
!                                                                         !
!  Purpose: Initialize particle arrays. The upper and lower bounds are    !
!  passed so that after resizing particle arrays (see GROW_PARTICLE) the  !
!  new portions of the arrays can be initialized.                         !
!                                                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_INIT_PARTICLE_ARRAYS(LB,UB)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use discretelement
      use desgrid
      use desmpi
      use des_thermo
      use des_rxns

      use run, only: ENERGY_EQ
      use run, only: ANY_SPECIES_EQ

      use particle_filter, only: FILTER_SIZE
      use particle_filter, only: FILTER_CELL, FILTER_WEIGHT
      use mfix_pic, only: MPPIC, DES_STAT_WT, PS_GRAD
      use mfix_pic, only: AVGSOLVEL_P, EPG_P

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LB, UB

      IGLOBAL_ID(LB:UB) = 0

! Physical properties:
      DES_RADIUS(LB:UB) = ZERO
      RO_Sol(LB:UB) = ZERO
      PVOL(LB:UB) = ZERO
      PMASS(LB:UB) = ZERO
      OMOI(LB:UB) = ZERO

! Particle position, velocity, etc
      DES_POS_NEW(:,LB:UB) = ZERO
      DES_VEL_NEW(:,LB:UB) = ZERO
      OMEGA_NEW(:,LB:UB) = ZERO

! Particle state flag
      PEA(LB:UB,:) = .FALSE.

! DES grid bin information
      DG_PIJK(LB:UB) = -1
      DG_PIJKPRV(LB:UB) = -1
      IGHOST_UPDATED(LB:UB) = -1

! Fluid cell bin information
      PIJK(LB:UB,:) = 0

! Translation and rotational forces
      FC(:,LB:UB) = ZERO
      TOW(:,LB:UB) = ZERO

! Collision data
      WALL_COLLISION_FACET_ID(:,LB:UB) = -1
      WALL_COLLISION_PFT(:,:,LB:UB) = ZERO

! Initializing user defined array
      DES_USR_VAR(:,LB:UB) = ZERO

! Paritcle center drag coefficient and explit drag force
      F_GP(LB:UB) = ZERO
      DRAG_FC(:,LB:UB) = ZERO


! Interpolation variables.
      IF(FILTER_SIZE > 0)THEN
         FILTER_CELL(:,LB:UB) = -1
         FILTER_WEIGHT(:,LB:UB) = ZERO
      ENDIF

! MPPIC variables
      IF(MPPIC) THEN
         DES_STAT_WT(LB:UB) = ZERO
         PS_GRAD(LB:UB,:) = ZERO
         AVGSOLVEL_P(:,LB:UB) = ZERO
         EPG_P(LB:UB) = ZERO
      ENDIF

! Higher order time integration variables.
      IF (DO_OLD) THEN
         DES_POS_OLD(:,LB:UB) = ZERO
         DES_VEL_OLD(:,LB:UB) = ZERO
         DES_ACC_OLD(:,LB:UB) = ZERO
         OMEGA_OLD(:,LB:UB) = ZERO
         ROT_ACC_OLD(:,LB:UB) = ZERO
      ENDIF

! Energy equation variables.
      IF(ENERGY_EQ)THEN
         DES_T_s_OLD(LB:UB) = ZERO
         DES_T_s_NEW(LB:UB) = ZERO
         DES_C_PS(LB:UB) = ZERO
         DES_X_s(LB:UB,:) = ZERO
         Q_Source(LB:UB) = ZERO
         IF (INTG_ADAMS_BASHFORTH) &
            Q_Source0(LB:UB) = ZERO
      ENDIF

! Chemical reaction variables.
      IF(ANY_SPECIES_EQ)THEN
         DES_R_sp(:,LB:UB) = ZERO
         DES_R_sc(:,LB:UB) = ZERO
         IF (INTG_ADAMS_BASHFORTH) THEN
            dMdt_OLD(LB:UB) = ZERO
            dXdt_OLD(:,LB:UB) = ZERO
         ENDIF
         Qint(LB:UB) = ZERO
      ENDIF


! Cohesion VDW forces
      IF(USE_COHESION) THEN
         PostCohesive (LB:UB) = ZERO
      ENDIF


      RETURN
      END SUBROUTINE DES_INIT_PARTICLE_ARRAYS



