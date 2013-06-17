!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_SPECIES_EQ                                        C
!  Purpose: Solve species mass balance equations in matrix equation    C
!     form Ax=b. The center coefficient (ap) and source vector (b)     C
!     are negative.  The off-diagonal coefficients are positive.       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 11-FEB-98  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLVE_SPECIES_EQ(IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      USE run
      USE physprop
      USE geometry
      USE fldvar
      USE output
      USE indices
      USE drag
      USE residual
      USE ur_facs 
      USE pgcor
      USE pscor
      USE leqsol 
      USE bc
      USE energy
      USE rxns
      Use ambm
      USE matrix 
      USE ChiScheme
      Use tmp_array, S_p => Array1, S_c => Array2, EPs => Array3, VxGama => Array4
      USE compar       
      USE mpi_utility      
      USE sendrecv 
      USE mflux    
      use ps
 
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index 
      INTEGER, INTENT(INOUT) :: IER 
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! phase index 
      INTEGER :: M
! species index 
      INTEGER :: LN
! previous time step term
      DOUBLE PRECISION :: apo
! Indices 
      INTEGER :: IJK, IMJK, IJMK, IJKM 
! linear equation solver method and iterations 
      INTEGER :: LEQM, LEQI 
! tmp array to pass to set_chi
      DOUBLE PRECISION :: X_s_temp(DIMENSION_3, DIMENSION_N_s) 

! Local/Global error flags.
      LOGICAL :: lDiverged, gDiverged

      DOUBLE PRECISION EP_SS_L_TOT
! temporary use of global arrays:
! array1 (locally s_p) 
! source lhs: coefficient of dependent variable
! becomes part of a_m matrix; must be positive
!      DOUBLE PRECISION :: S_P(DIMENSION_3)
! array2 (locally s_c) 
! source rhs vector: constant part becomes part of b_m vector
!      DOUBLE PRECISION :: S_C(DIMENSION_3)
! array3 (locally eps)
! alias for solids volume fraction      
!      DOUBLE PRECISION :: eps(DIMENSION_3)
! array4 (locally vxgama)
!      DOUBLE PRECISION :: vxgama(DIMENSION_3)
! Septadiagonal matrix A_m, vector b_m
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)

!-----------------------------------------------
! External functions
!-----------------------------------------------
      LOGICAL , EXTERNAL :: IS_SMALL 
      DOUBLE PRECISION , EXTERNAL :: Check_conservation 
!-----------------------------------------------
! Local statement functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

! Initialize error flag.
      lDiverged = .FALSE.
     
      call lock_ambm       ! locks arrys a_m and b_m
      call lock_tmp_array  ! locks array1,array2,array3,array4
                           ! (locally s_p, s_c, eps, vxgama) 

! Fluid phase species mass balance equations
! ---------------------------------------------------------------->>>
      IF (SPECIES_EQ(0)) THEN 
         IF(chi_scheme) call set_chi(DISCRETIZE(7), X_g, NMAX(0), &
                                     U_g, V_g, W_g, IER)

! looping over species         
         DO LN = 1, NMAX(0) 
            CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0, IER) 
!!$omp    parallel do private(IJK, APO)
            DO IJK = ijkstart3, ijkend3 
               IF (FLUID_AT(IJK)) THEN
! calculate the source terms to be used in the a matrix and b vector
                   APO = ROP_GO(IJK)*VOL(IJK)*ODT 
                   S_P(IJK) = APO+&
                      (ZMAX(SUM_R_G(IJK))+ROX_GC(IJK,LN))*VOL(IJK) 
                   S_C(IJK) = APO*X_GO(IJK,LN) + &
                      X_G(IJK,LN)*ZMAX((-SUM_R_G(IJK)))*VOL(IJK) +&
                      R_GP(IJK,LN)*VOL(IJK)
               ELSE 
                  S_P(IJK) = ZERO 
                  S_C(IJK) = ZERO 
               ENDIF 
            ENDDO

! calculate the convection-diffusion terms 
            IF(.NOT.ADDED_MASS) THEN
               CALL CONV_DIF_PHI (X_G(1,LN), DIF_G(1,LN), &
                  DISCRETIZE(7), U_G, V_G, W_G, &
                  Flux_gE, Flux_gN, Flux_gT, 0, A_M, B_M, IER)
            ELSE
               CALL CONV_DIF_PHI (X_G(1,LN), DIF_G(1,LN),&
                  DISCRETIZE(7), U_G, V_G, W_G, &
                  Flux_gSE, Flux_gSN, Flux_gST, 0, A_M, B_M, IER)
            ENDIF

! calculate standard bc
            CALL BC_PHI (X_G(1,LN), BC_X_G(1,LN), BC_XW_G(1,LN), &
               BC_HW_X_G(1,LN), BC_C_X_G(1,LN), 0, A_M, B_M, IER) 

! set the source terms in a and b matrix form
            CALL SOURCE_PHI (S_P, S_C, EP_G, X_G(1,LN), 0, A_M, B_M, IER)

! Add point souce contriubtions.
            IF(POINT_SOURCE) CALL POINT_SOURCE_SPECIES_EQ (X_G(1,LN), &
               BC_X_G(:,LN), BC_MASSFLOW_G(:), 0, A_M, B_M, IER)

            CALL CALC_RESID_S (X_G(1,LN), A_M, B_M, 0, &
               NUM_RESID(RESID_X+(LN-1),0), &
               DEN_RESID(RESID_X+(LN-1),0), RESID(RESID_X+(LN-1),0), &
               MAX_RESID(RESID_X+(LN-1),0), IJK_RESID(RESID_X+(LN-1),0), &
               ZERO_X_GS, IER) 

            CALL UNDER_RELAX_S (X_G(1,LN), A_M, B_M, 0, UR_FAC(7), IER) 

            CALL ADJUST_LEQ (RESID(RESID_X+(LN-1),0), LEQ_IT(7), &
               LEQ_METHOD(7), LEQI, LEQM, IER)

! solve the phi equation 
            CALL SOLVE_LIN_EQ ('X_g', 7, X_G(1,LN), A_M, B_M, 0, &
               LEQI, LEQM, LEQ_SWEEP(7), LEQ_TOL(7), LEQ_PC(7), IER) 

! Check for linear solver divergence.
            IF(ier == -2) lDiverged = .TRUE.
            CALL BOUND_X (X_G(1,LN), IJKMAX2, IER) 

         ENDDO    ! end do loop (ln = 1, nmax(0)
         IF(chi_scheme) call unset_chi(IER)
      ENDIF 
! end fluid phase species equations      
! ----------------------------------------------------------------<<<

! Solids phase species balance equations
! ---------------------------------------------------------------->>>
      DO M = 1, SMAX 
         IF (SPECIES_EQ(M)) THEN 
            IF(chi_scheme) THEN
               DO LN = 1, NMAX(M)
                  DO IJK = ijkstart3, ijkend3
                     X_S_temp(IJK, LN) = X_S(IJK,M,LN)
                  ENDDO
                ENDDO 
              call set_chi(DISCRETIZE(7), X_S_temp, NMAX(M), &
                           U_S(1,M), V_S(1,M), W_S(1,M), IER)
            ENDIF ! for chi_scheme

            DO LN = 1, NMAX(M) 
               CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 

!!$omp    parallel do private(IJK, APO)
               DO IJK = ijkstart3, ijkend3 
                  IF (FLUID_AT(IJK)) THEN 
                    APO = ROP_SO(IJK,M)*VOL(IJK)*ODT 
                    S_P(IJK) = APO + &
                       (ZMAX(SUM_R_S(IJK,M))+ROX_SC(IJK,M,LN))*VOL(IJK) 
                    S_C(IJK) = APO*X_SO(IJK,M,LN) + &
                       X_S(IJK,M,LN)*ZMAX((-SUM_R_S(IJK,M)))*VOL(IJK) + &
                       R_SP(IJK,M,LN)*VOL(IJK)
                    EPS(IJK) = EP_S(IJK,M) 
                  ELSE 
                     S_P(IJK) = ZERO 
                     S_C(IJK) = ZERO 
                     EPS(IJK) = ZERO 
                  ENDIF 
               ENDDO 

               IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                  CALL CONV_DIF_PHI (X_S(1,M,LN), DIF_S(1,M,LN), &
                    DISCRETIZE(7), U_S(1,M), V_S(1,M), W_S(1,M), &
                    Flux_sE(1,M), Flux_sN(1,M), Flux_sT(1,M), M, A_M, B_M, IER)
               ELSE
                  CALL CONV_DIF_PHI (X_S(1,M,LN), DIF_S(1,M,LN), &
                    DISCRETIZE(7), U_S(1,M), V_S(1,M), W_S(1,M), &
                    Flux_sSE, Flux_sSN, Flux_sST, M, A_M, B_M, IER)
               ENDIF

               CALL BC_PHI (X_S(1,M,LN), BC_X_S(1,M,LN), &
                  BC_XW_S(1,M,LN), BC_HW_X_S(1,M,LN), &
                  BC_C_X_S(1,M,LN), M, A_M, B_M, IER) 

               CALL SOURCE_PHI (S_P, S_C, EPS, X_S(1,M,LN), M, A_M, B_M, IER)

               IF(POINT_SOURCE) CALL POINT_SOURCE_SPECIES_EQ(X_S(1,M,LN), &
                  BC_X_S(:,M,LN), BC_MASSFLOW_S(:,M), M, A_M, B_M, IER)

               CALL CALC_RESID_S (X_S(1,M,LN), A_M, B_M, M, &
                  NUM_RESID(RESID_X+(LN-1),M), &
                  DEN_RESID(RESID_X+(LN-1),M), RESID(RESID_X+(LN-1),&
                  M), MAX_RESID(RESID_X+(LN-1),M), IJK_RESID(RESID_X+(LN-1),M), &
                  ZERO_X_GS, IER) 

               CALL UNDER_RELAX_S (X_S(1,M,LN), A_M, B_M, M, UR_FAC(7), IER) 

!               call check_ab_m(a_m, b_m, m, .false., ier)
!               write(*,*) resid(resid_x+(LN-1), m), &
!                  max_resid(resid_x+(LN-1), m), &
!                  ijk_resid(resid_x+(LN-1), m)
!               call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!
!               call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1,-3,M),&
!                  1, DO_K, ier)

               CALL ADJUST_LEQ (RESID(RESID_X+(LN-1),M), LEQ_IT(7), &
                  LEQ_METHOD(7), LEQI, LEQM, IER) 

               CALL SOLVE_LIN_EQ ('X_s', 7, X_S(1,M,LN), A_M, B_M, M,&
                  LEQI, LEQM, LEQ_SWEEP(7), LEQ_TOL(7), LEQ_PC(7), IER) 

! Check for linear solver divergence.
               IF(ier == -2) lDiverged = .TRUE.

               CALL BOUND_X (X_S(1,M,LN), IJKMAX2, IER) 
!               call out_array(X_s(1,m,LN), 'X_s')

            END DO 

            IF (SOLID_RO_V) THEN
               DO IJK = ijkstart3, ijkend3
                  IF (FLUID_AT(IJK)) THEN 
                     EP_SS_L_TOT = ZERO
                     DO LN = 1,NMAX(M)
                        EP_SS_L_TOT = EP_SS_L_TOT + X_S(IJK,M,LN) / RO_SS(M,LN)
                     ENDDO
                     IF(EP_SS_L_TOT .GT. ZERO) THEN
!                       phase density
                        RO_SV(IJK,M) = ONE/EP_SS_L_TOT
                     ENDIF

                     if(ROP_s(IJK,M) .eq. 0.E0) then
                        if(RO_SV(IJK,M) .eq. 0.d0) then
                           write(*,*) 'WARNING: zero ro_s in solve_species_eq.f'
                        endif
                     endif

                  ENDIF  ! Fluid cell
               ENDDO  ! IJK Loop



            ENDIF
!end

            if(chi_scheme) call unset_chi(IER)
         ENDIF ! check for any species in phase m
      END DO ! for m = 1, mmax
! end solids phases species equations      
! ----------------------------------------------------------------<<<
      

! If the linear solver diverged, species mass fractions may take on
! unphysical values. To prevent them from propogating through the domain
! or causing failure in other routines, force an exit from iterate and
! reduce the time step.
      CALL GLOBAL_ALL_OR(lDiverged, gDiverged)
      if(gDiverged) IER = -100

      call unlock_ambm
      call unlock_tmp_array

      RETURN  
      END SUBROUTINE SOLVE_SPECIES_EQ 
      

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_SPECIES_EQ                                 C
!  Purpose: Adds point sources to the species equations.               C
!                                                                      C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_SPECIES_EQ(X_x, BC_X, BC_FLOW, M, &
         A_M, B_M, IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use compar
      use fldvar
      use funits
      use geometry
      use indices
      use physprop
      use ps
      use run
      use toleranc

! To be removed upon complete integration of point source routines.
      use bc
      use usr

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m 
      DOUBLE PRECISION, intent(in) :: X_x(DIMENSION_3)

! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(IN) :: BC_X(DIMENSION_BC)

! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(IN) :: BC_FLOW(DIMENSION_BC) 

      INTEGER, intent(in) :: M

! Septadiagonal matrix A_m 
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 

! Vector b_m 
      DOUBLE PRECISION, INTENT(INOUT) :: B_M(DIMENSION_3, 0:DIMENSION_M) 

! Error index 
      INTEGER, INTENT(INOUT) :: IER 

!----------------------------------------------- 
! Local Variables
!----------------------------------------------- 

! Indices 
      INTEGER :: IJK, I, J, K
      INTEGER :: BCV

! terms of bm expression
      DOUBLE PRECISION pSource

!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------
      BC_LP: do BCV = 50, DIMENSION_BC
         if(POINT_SOURCES(BCV) == 0) cycle BC_LP

         do k = BC_K_B(BCV), BC_K_T(BCV)
         do j = BC_J_S(BCV), BC_J_N(BCV)
         do i = BC_I_W(BCV), BC_I_E(BCV)

            if(.NOT.IS_ON_myPE_plus2layers(I,J,K)) cycle

            ijk = funijk(i,j,k)
            if(.NOT.fluid_at(ijk)) cycle

            if(A_M(IJK,0,M) == -ONE .AND. &
               B_M(IJK,M) == -X_x(IJK)) then
               B_M(IJK,M) = -BC_X(BCV)

            else
               pSource = BC_FLOW(BCV) * (VOL(IJK)/PS_VOLUME(BCV))
               A_M(IJK,0,M) = A_M(IJK,0,M) - pSource
               B_M(IJK,M) = B_M(IJK,M) - BC_X(BCV) * pSource 
            endif

         enddo
         enddo
         enddo

      enddo BC_LP


      RETURN
      END SUBROUTINE POINT_SOURCE_SPECIES_EQ
