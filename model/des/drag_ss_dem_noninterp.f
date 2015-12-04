!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DRAG_SS_DEM_NONINTERP                                   !
!                                                                      !
!  Purpose: This routine is called from the DISCRETE side to calculate !
!  the solids-solids drag force acting on each particle using cell-    !
!  center continuum solids velocity. The total contact force is also   !
!  updated to include P* for over-packing.                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_SS_DEM_NONINTERP

! Global Variables:
!---------------------------------------------------------------------//
! The count and a list of particles in IJK
      use discretelement, only: PINC
      use derived_types, only: PIC
! Particle velocity and density
      use discretelement, only: DES_VEL_NEW, RO_SOL
! Particle radius and volume.
      use discretelement, only: DES_RADIUS, PVOL
! Total forces acting on particle
      use discretelement, only: FC
      use functions, only: IS_NONEXISTENT, IS_GHOST, IS_ENTERING_GHOST, IS_EXITING_GHOST
! Number of continuum solids phases
      use physprop, only: MMAX
! Diameter of continuum solids phases
      use fldvar, only: D_P
! Material (ROs) and bulk density (ROPs) of continuum solids
      use fldvar, only: ROP_s, RO_s
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Function to deterine if a cell contains fluid.
      use functions, only: FLUID_AT
! Gas phase volume fraction
      use fldvar, only: EP_g
! Solids pressure
      use fldvar, only: P_STAR
! Flag that a continuum solids can pack
      use physprop, only: CLOSE_PACKED
! Fudge factor for SS drag (c.f. Gera, 2004)
      use constant, only: SEGREGATION_SLOPE_COEFFICIENT

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO, ONE
! Array sizes for solids
      use param, only: DIMENSION_M

      IMPLICIT NONE

! Local variables:
!---------------------------------------------------------------------//
! Loop indices:
      INTEGER :: IJK, M, NINDX, NP
! average solids velocity at scalar cell center in array form
      DOUBLE PRECISION :: VELCS(3, DIMENSION_M)
! relative velocity between solids phase m and l
      DOUBLE PRECISION :: VSLP(3), VREL
! radial distribution function between phase M and L
      DOUBLE PRECISION :: G0_ML
! Sum over all phases of ratio volume fraction over particle diameter
      DOUBLE PRECISION :: EPSoDP
! DEM particle diameter calculated from radius.
      DOUBLE PRECISION :: lDP
! solid-solid drag coefficient
      DOUBLE PRECISION :: lDss
! Intermediate calculations for volume fraction.
      DOUBLE PRECISION :: OoEPg, EPg_2
! Drag force acting on each phase.
      DOUBLE PRECISION :: D_FORCE(3)

!......................................................................!


! Calculate the solid-solid drag for each particle.
!---------------------------------------------------------------------//
!!$omp parallel do schedule(guided, 50) default(none)              &
!!$omp shared(IJKSTART3, IJKEND3, PINC, PIC, DES_VEL_NEW,          &
!!$omp   MMAX, D_P, RO_s, ROP_s, EP_G, DES_RADIUS, RO_SOL, FC,     &
!!$omp   PVOL, SEGREGATION_SLOPE_COEFFICIENT, CLOSE_PACKED, P_STAR)&
!!$omp private(IJK, OoEPg, EPg_2, EPSoDP, NP, lDP, D_FORCE, G0_ML, &
!!$omp   VELCS, VSLP, VREL, lDss)
      DO IJK = IJKSTART3, IJKEND3

         IF(.NOT.FLUID_AT(IJK)) CYCLE
         IF(PINC(IJK) == 0) CYCLE

         OoEPg = ONE/EP_g(IJK)
         EPg_2 = EP_g(IJK)*EP_g(IJK)

! Calculate solids volume fraction over particle diameter.
         CALL CALC_EPSoDP(IJK, EPSoDP)

! Calculate the continuum solids velocity at scalar cell center.
         CALL CALC_NONINTERP_VELFP_SOLIDS(IJK, VelCS)

! Calculate the solids drag for each particle in the current cell.
         DO NINDX = 1,PINC(IJK)
            NP = PIC(IJK)%P(NINDX)
! skipping indices that do not represent particles and ghost particles
            IF(IS_NONEXISTENT(NP)) CYCLE
            IF(IS_GHOST(NP) .OR. IS_ENTERING_GHOST(NP) .OR. IS_EXITING_GHOST(NP)) CYCLE

! Diameter of particle (not a phase diameter).
            lDP = 2.0d0*DES_RADIUS(NP)

            D_FORCE = ZERO
            DO M = 1, MMAX

! evaluating g0 - taken from G_0.f subroutine (lebowitz form)
! this section is needed to account for all solids phases until g0 for
! multiple solids types (i.e. discrete & continuum) can be addressed
! more effectively.
               G0_ML = OoEPg + 3.0d0*EPSoDP*D_P(IJK,M)*lDP /           &
                  (EPg_2 *(D_P(IJK,M) + lDP))


! Relative (slip) velocity.
               VSLP = DES_VEL_NEW(NP,:) - VelCS(:,M)
! Relative velocity magnitude.
               VREL = sqrt(dot_product(VSLP,VSLP))

               CALL DRAG_SS_SYAM(lDss, D_P(IJK,M), lDP, RO_S(IJK,M),   &
                  RO_SOL(NP), G0_ML, VREL)

               lDss = lDss*ROP_S(IJK,M)*RO_Sol(NP)

! accounting for particle-particle drag due to enduring contact in a
! close-packed system
               IF(CLOSE_PACKED(M)) lDss = lDss +                       &
                  SEGREGATION_SLOPE_COEFFICIENT*P_star(IJK)

! Calculating the accumulated solids-solids drag force.
               D_FORCE(:) = D_FORCE(:) - lDss*VSLP(:)

            ENDDO ! end do loop (M=1,MMAX)

            FC(NP,:) = FC(NP,:) + D_FORCE(:)*PVOL(NP)

         ENDDO ! END DO LOOP (NP=1,MAX_PIP)
      ENDDO ! END DO LOOP (IJK=IJKSTART3, IJKEND3)
!!$omp end parallel do



      RETURN
      END SUBROUTINE DRAG_SS_DEM_NONINTERP



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_DRAG_SS                                             !
!                                                                      !
!  Purpose: This subroutine is called from the CONTINUUM side. It      !
!  calculate the solids-solids drag force coefficient between          !
!  continuum and discrete solids.                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_SS_TFM_NONINTERP

! Global Variables:
!---------------------------------------------------------------------//
! The count and a list of particles in IJK
      use derived_types, only: PIC
      use discretelement, only: PINC
! Particle velocity and density
      use discretelement, only: DES_VEL_NEW, RO_SOL
! Particle radius and volume.
      use discretelement, only: DES_RADIUS
! Total forces acting on particle
      use discretelement, only: SDRAG_AM, SDRAG_BM, F_SDS
! Number of continuum solids phases
      use physprop, only: MMAX
! Diameter of continuum solids phases
      use fldvar, only: D_P
! Material (ROs) and bulk density (ROPs) of continuum solids
      use fldvar, only: ROP_s, RO_s
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Function to deterine if a cell contains fluid.
      use functions, only: FLUID_AT
! Gas phase volume fraction
      use fldvar, only: EP_g
! Solids pressure
      use fldvar, only: P_STAR
! Flag that a continuum solids can pack
      use physprop, only: CLOSE_PACKED
! Fudge factor for SS drag (c.f. Gera, 2004)
      use constant, only: SEGREGATION_SLOPE_COEFFICIENT

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO, ONE
! Array sizes for solids
      use param, only: DIMENSION_M

      use functions, only: IS_NONEXISTENT, IS_GHOST, IS_ENTERING_GHOST, IS_EXITING_GHOST

      IMPLICIT NONE

! Local variables:
!---------------------------------------------------------------------//
! Loop indices:
      INTEGER :: IJK, M, NINDX, NP
! average solids velocity at scalar cell center in array form
      DOUBLE PRECISION :: VELCS(3, DIMENSION_M)
! relative velocity between solids phase m and l
      DOUBLE PRECISION :: VSLP(3), VREL
! radial distribution function between phase M and L
      DOUBLE PRECISION :: G0_ML
! Sum over all phases of ratio volume fraction over particle diameter
      DOUBLE PRECISION :: EPSoDP
! DEM particle diameter calculated from radius.
      DOUBLE PRECISION :: lDP
! solid-solid drag coefficient
      DOUBLE PRECISION :: lDss
! Intermediate calculations for volume fraction.
      DOUBLE PRECISION :: OoEPg, EPg_2
! Drag force acting on each phase.
      DOUBLE PRECISION :: lFORCE
! One divided by fluid cell volume
      DOUBLE PRECISION :: OoVOL

!......................................................................!

      DO IJK = IJKSTART3, IJKEND3

         F_SDS(IJK,:) = ZERO
         SDRAG_AM(IJK,:) = ZERO
         SDRAG_BM(IJK,:,:) = ZERO

         IF(.NOT.FLUID_AT(IJK)) CYCLE
         IF(PINC(IJK) == 0) CYCLE

         OoEPg = ONE/EP_g(IJK)
         EPg_2 = EP_g(IJK)*EP_g(IJK)

! Calculate solids volume fraction over particle diameter.
         CALL CALC_EPSoDP(IJK, EPSoDP)

! Calculate the continuum solids velocity at scalar cell center.
         CALL CALC_NONINTERP_VELFP_SOLIDS(IJK, VelCS)

! Calculate the solids drag for each particle in the current cell.
         DO NINDX = 1,PINC(IJK)
            NP = PIC(IJK)%P(NINDX)
! skipping indices that do not represent particles and ghost particles
            IF(IS_NONEXISTENT(NP)) CYCLE
            IF(IS_GHOST(NP) .OR. IS_ENTERING_GHOST(NP) .OR. IS_EXITING_GHOST(NP)) CYCLE

! Diameter of particle (not a phase diameter).
            lDP = 2.0d0*DES_RADIUS(NP)

            DO M = 1, MMAX

! Relative (slip) velocity.
               VSLP = DES_VEL_NEW(NP,:) - VelCS(:,M)
! Relative velocity magnitude.
               VREL = sqrt(dot_product(VSLP,VSLP))

! evaluating g0 - taken from G_0.f subroutine (lebowitz form)
! this section is needed to account for all solids phases until g0 for
! multiple solids types (i.e. discrete & continuum) can be addressed
! more effectively.
               G0_ML = OoEPg + 3.0d0*EPSoDP*D_P(IJK,M)*lDP /           &
                  (EPg_2 *(D_P(IJK,M) + lDP))

               CALL DRAG_SS_SYAM(lDss, D_P(IJK,M), lDP, RO_S(IJK,M),   &
                  RO_SOL(NP), G0_ML, VREL)

               lDss = lDss*ROP_S(IJK,M)*RO_Sol(NP)

! accounting for particle-particle drag due to enduring contact in a
! close-packed system
               IF(CLOSE_PACKED(M)) lDss = lDss +                       &
                  SEGREGATION_SLOPE_COEFFICIENT*P_star(IJK)

! Calculating the accumulated solids-solids drag force.
               lFORCE = OoVOL*lDss

              SDRAG_AM(IJK,M) = SDRAG_AM(IJK,M) + lFORCE
              SDRAG_BM(IJK,:,M) = SDRAG_BM(IJK,:,M) +                  &
                 lFORCE*DES_VEL_NEW(NP,:)

            ENDDO ! end do loop (M=1,MMAX)

         ENDDO ! END DO LOOP (NP=1,MAX_PIP)

        F_SDS(IJK,:) = SDRAG_AM(IJK,:)

      ENDDO ! END DO LOOP (IJK=IJKSTART3, IJKEND3)



      RETURN
      END SUBROUTINE DRAG_SS_TFM_NONINTERP




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_NONINTER_VELFP                                     !
!  Author: J.Musser                                   Date: 07-NOV-14  !
!                                                                      !
!  Purpose: Calculate the scalar cell center gas velocity. This code   !
!  is common to the DEM and GAS calls for non-interpolated drag        !
!  routines.                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_NONINTERP_VELFP_SOLIDS(IJK, lVELFP)

! Global Variables:
!---------------------------------------------------------------------//
! Number of solid phases.
      use physprop, only: MMAX
! Solids velocities.
      use fldvar, only: U_S, V_S, W_S
! Functions to average momentum to scalar cell center.
      use fun_avg, only: AVG_X_E, AVG_Y_N, AVG_Z_T
! Flags and correction factors for cut momentum cells.
      use cutcell, only: CUT_U_TREATMENT_AT, THETA_UE, THETA_UE_BAR
      use cutcell, only: CUT_V_TREATMENT_AT, THETA_VN, THETA_VN_BAR
      use cutcell, only: CUT_W_TREATMENT_AT, THETA_WT, THETA_WT_BAR
! Functions to lookup adjacent cells by index.
      use functions, only: IM_OF, JM_OF, KM_OF
      use indices, only: I_OF
! Flag for 3D simulatoins.
      use geometry, only: DO_K

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision parameters.
      use param1, only: ZERO
! Array sizes for solids
      use param, only: DIMENSION_M

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Fluid cell and solids phase indices
      INTEGER, INTENT(IN) :: IJK
! Fluid velocity vector at IJK cell center.
      DOUBLE PRECISION, INTENT(OUT) :: lVELFP(3, DIMENSION_M)

! Local variables:
!---------------------------------------------------------------------//
! Phase loop counter
      INTEGER :: M
! Indices of adjacent cells
      INTEGER :: IMJK, IJMK, IJKM

! Calculate the average fluid velocity at scalar cell center.
      DO M=1,MMAX

         IMJK = IM_OF(IJK)
         IF(CUT_U_TREATMENT_AT(IMJK)) THEN
            lVELFP(1,M) = (THETA_UE_BAR(IMJK)*U_S(IMJK,M) +            &
               THETA_UE(IMJK)*U_S(IJK,M))
         ELSE
            lVELFP(1,M) = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I_OF(IJK))
         ENDIF

         IJMK = JM_OF(IJK)
         IF(CUT_V_TREATMENT_AT(IJMK)) THEN
            lVELFP(2,M) = (THETA_VN_BAR(IJMK)*V_S(IJMK,M) +            &
               THETA_VN(IJMK)*V_S(IJK,M))
         ELSE
            lVELFP(2,M) = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M))
         ENDIF

         IF(DO_K) THEN
            IJKM = KM_OF(IJK)
            IF(CUT_W_TREATMENT_AT(IJKM)) THEN
               lVELFP(3,M) = (THETA_WT_BAR(IJKM)*W_S(IJKM,M) +         &
                  THETA_WT(IJKM)* W_S(IJK,M))
            ELSE
               lVELFP(3,M) = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M))
            ENDIF
         ELSE
            lVELFP(3,M) = ZERO
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CALC_NONINTERP_VELFP_SOLIDS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_NONINTER_VELFP                                     !
!  Author: J.Musser                                   Date: 07-NOV-14  !
!                                                                      !
!  Purpose: Calculate the scalar cell center gas velocity. This code   !
!  is common to the DEM and GAS calls for non-interpolated drag        !
!  routines.                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_EPSoDP(IJK, lEPSoDP)

! Global Variables:
!---------------------------------------------------------------------//
! The count and a list of particles in IJK
      use derived_types, only: PIC
      use discretelement, only: PINC
! Particle volume and radius
      use discretelement, only: PVOL, DES_RADIUS
! Volume of scalar cell
      use geometry, only: VOL
! Continuum solids diamter
      use fldvar, only: D_P
! Function to calculate continuum solids volume fraction
      use fldvar, only: EP_s
! Number of continuum solids phases
      use physprop, only: MMAX

      use functions

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision parameters.
      use param1, only: ZERO

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Fluid cell and solids phase indices
      INTEGER, INTENT(IN) :: IJK
! Fluid velocity vector at IJK cell center.
      DOUBLE PRECISION, INTENT(OUT) :: lEPSoDP

! Local variables:
!---------------------------------------------------------------------//
! Loop counters
      INTEGER :: NP, NINDX, M

! Calculate the sum of the particle volumes divided by radius.
      lEPSoDP = ZERO
      DO NINDX = 1,PINC(IJK)
         NP = PIC(IJK)%P(NINDX)
         IF(IS_NONEXISTENT(NP)) CYCLE
         IF(IS_GHOST(NP) .OR. IS_ENTERING_GHOST(NP) .OR. IS_EXITING_GHOST(NP)) CYCLE
         lEPSoDP = lEPSoDP + PVOL(NP)/DES_RADIUS(NP)
      ENDDO
! Convert radius to diameter and divide by cell volume.
      lEPSoDP = lEPSoDP/(2.0d0*VOL(IJK))

! Add contributions from continuum solids.
      DO M = 1, MMAX
         lEPSoDP = lEPSoDP + EP_s(IJK,M)/D_p(IJK,M)
      ENDDO

      RETURN
      END SUBROUTINE CALC_EPSoDP
