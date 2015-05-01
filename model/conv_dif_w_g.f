!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_DIF_W_g                                            C
!  Purpose: Determine convection diffusion terms for w_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_w_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_DIF_W_G(A_M, B_M, IER)

! Modules
!---------------------------------------------------------------------//
      USE param, only: dimension_3, dimension_m
      USE run, only: momentum_z_eq
      USE run, only: discretize
      USE run, only: def_cor
      USE visc_g, only: mu_gt
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!---------------------------------------------------------------------//


      IF (.NOT.MOMENTUM_Z_EQ(0)) RETURN

      IF (DEF_COR) THEN
! USE DEFERRED CORRECTION TO SOLVE W_G
         CALL STORE_A_W_G0 (A_M(1,-3,0), IER)
         IF (DISCRETIZE(5) > 1) CALL STORE_A_W_GDC(B_M(1,0))

      ELSE
! DO NOT USE DEFERRED CORRECTOIN TO SOLVE FOR W_G
         IF (DISCRETIZE(5) == 0) THEN               ! 0 & 1 => FOUP
            CALL STORE_A_W_G0 (A_M(1,-3,0), IER)
         ELSE
            CALL STORE_A_W_G1 (A_M(1,-3,0))
         ENDIF
      ENDIF

      CALL DIF_W_IS(MU_GT, A_M, 0, IER)

      RETURN
      END SUBROUTINE CONV_DIF_W_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of velocity on the east, north,   C
!  and top face of a w-momentum cell                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_WCELL_GVTERMS(U, V, WW)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE cutcell, only: cut_w_treatment_at
      USE cutcell, only: theta_wt, theta_wt_bar
      USE cutcell, only: theta_w_be, theta_w_te
      USE cutcell, only: theta_w_bn, theta_w_tn
      USE cutcell, only: alpha_we_c, alpha_wn_c, alpha_wt_c

      USE fldvar, only: u_g, v_g, w_g

      USE fun_avg, only: avg_z_t, avg_z
      USE functions, only: kp_of
      USE indices, only: k_of

      USE param, only: dimension_3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(OUT) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: WW(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: IJK, K, IJKP
! for cartesian grid
      DOUBLE PRECISION :: AW, HW, VELW
!---------------------------------------------------------------------//


!!!$omp parallel do private(IJK,K,IJKT,IJKP)
      DO IJK = ijkstart3, ijkend3
         K = K_OF(IJK)
         IJKP = KP_OF(IJK)

         IF(CUT_W_TREATMENT_AT(IJK)) THEN

! East face (i+1/2, j, k+1/2)
            U(IJK) = (Theta_W_be(IJK) * U_G(IJK) + &
                      Theta_W_te(IJK) * U_G(IJKP))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM', &
               ALPHA_We_c(IJK), AW, HW, VELW)
            U(IJK) = U(IJK) * AW

! North face (i, j+1/2, k+1/2)
            V(IJK) = (Theta_W_bn(IJK) * V_G(IJK) + &
                      Theta_W_tn(IJK) * V_G(IJKP))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM', &
               ALPHA_Wn_c(IJK), AW, HW, VELW)
            V(IJK) = V(IJK) * AW

! Top face (i, j, k+1)
            WW(IJK) = (Theta_Wt_bar(IJK) * W_G(IJK) + &
                       Theta_Wt(IJK) * W_G(IJKP))
            CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM', &
               alpha_Wt_c(IJK), AW, HW, VELW)
            WW(IJK) = WW(IJK) * AW

         ELSE
            U(IJK) = AVG_Z(U_G(IJK),U_G(IJKP),K)
            V(IJK) = AVG_Z(V_G(IJK),V_G(IJKP),K)
            WW(IJK) = AVG_Z_T(W_G(IJK),W_G(IJKP))
         ENDIF   ! end if/else cut_w_treatment_at
      ENDDO   ! end do ijk

      RETURN
      END SUBROUTINE GET_WCELL_GVTERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the convective fluxes through the faces of a     C
!  w-momentum cell. Note the fluxes are calculated at all faces of     C
!  regardless of flow_at_t of condition of the west, south, or         C
!  bottom face.                                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_WCELL_GCFLUX_TERMS(FLUX_E, FLUX_W, FLUX_N, &
         FLUX_S, FLUX_T, FLUX_B, IJK)

! Modules
!---------------------------------------------------------------------//
      USE cutcell, only: cut_w_treatment_at
      USE cutcell, only: theta_wt, theta_wt_bar
      USE cutcell, only: theta_w_be, theta_w_te
      USE cutcell, only: theta_w_bn, theta_w_tn
      USE cutcell, only: alpha_we_c, alpha_wn_c, alpha_wt_c

      USE functions, only: kp_of, im_of, jm_of, km_of

      USE mflux, only: flux_ge, flux_gn, flux_gt

      USE param1, only: half
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! fluxes through faces of given ijk u-momentum cell
      DOUBLE PRECISION, INTENT(OUT) :: flux_e, flux_w
      DOUBLE PRECISION, INTENT(OUT) :: flux_n, flux_s
      DOUBLE PRECISION, INTENT(OUT) :: flux_t, flux_b
! ijk index
      INTEGER, INTENT(IN) :: ijk

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: imjk, ijmk, ijkm
      INTEGER :: ijkp, imjkp, ijmkp

! for cartesian grid
      DOUBLE PRECISION :: AW, HW, VELW
!---------------------------------------------------------------------//

      IJKP = KP_OF(IJK)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)
      IMJKP = KP_OF(IMJK)
      IJMKP = KP_OF(IJMK)

      IF(CUT_W_TREATMENT_AT(IJK)) THEN
! East face (i+1/2, j, k+1/2)
         Flux_e = (Theta_W_be(IJK) * Flux_gE(IJK) + &
                 Theta_W_te(IJK) * Flux_gE(IJKP))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',&
            ALPHA_We_c(IJK), AW, HW, VELW)
         Flux_e = Flux_e * AW
! West face (i-1/2, j, k+1/2)
         Flux_w = (Theta_W_be(IMJK) * Flux_gE(IMJK) + &
                   Theta_W_te(IMJK) * Flux_gE(IMJKP))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',&
            ALPHA_We_c(IMJK), AW, HW, VELW)
         Flux_w = Flux_w * AW


! North face (i, j+1/2, k+1/2)
         Flux_n = (Theta_W_bn(IJK) * Flux_gN(IJK) + &
                 Theta_W_tn(IJK) * Flux_gN(IJKP))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',&
            ALPHA_Wn_c(IJK), AW, HW, VELW)
         Flux_n = Flux_n * AW
! South face (i, j-1/2, k+1/2)
         Flux_s = (Theta_W_bn(IJMK) * Flux_gN(IJMK) + &
                   Theta_W_tn(IJMK) * Flux_gN(IJMKP))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',&
            ALPHA_Wn_c(IJMK), AW, HW, VELW)
         Flux_s = Flux_s * AW


! Top face (i, j, k+1)
         Flux_t = (Theta_Wt_bar(IJK) * Flux_gT(IJK) + &
                   Theta_Wt(IJK) * Flux_gT(IJKP))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',&
            alpha_Wt_c(IJK),AW,HW,VELW)
         Flux_t = Flux_t * AW
! Bottom face (i, j, k)
         Flux_b = (Theta_Wt_bar(IJKM) * Flux_gT(IJKM) + &
                   Theta_Wt(IJKM) * Flux_gT(IJK))
         CALL GET_INTERPOLATION_TERMS_G(IJK,'W_MOMENTUM',&
            alpha_Wt_c(IJKM), AW, HW, VELW)
         Flux_b = Flux_b * AW
      ELSE
         Flux_e = HALF * (Flux_gE(IJK) + Flux_gE(IJKP))
         Flux_w = HALF * (Flux_gE(IMJK) + Flux_gE(IMJKP))
         Flux_n = HALF * (Flux_gN(IJK) + Flux_gN(IJKP))
         Flux_s = HALF * (Flux_gN(IJMK) + Flux_gN(IJMKP))
         Flux_t = HALF * (Flux_gT(IJK) + Flux_gT(IJKP))
         Flux_b = HALF * (Flux_gT(IJKM) + Flux_gT(IJK))
      ENDIF   ! end if/else cut_w_treatment_at

      RETURN
      END SUBROUTINE GET_WCELL_GCFLUX_TERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of diffusive flux through the     C
!  faces of a w-momentum cell. Note the fluxes are calculated at       C
!  all faces regardless of flow_at_t condition of the west, south      C
!  or bottom face.                                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_WCELL_GDIFF_TERMS(D_FE, D_FW, D_FN, D_FS, &
         D_FT, D_FB, IJK)

! Modules
!---------------------------------------------------------------------//
      USE cutcell, only: cut_w_treatment_at
      USE cutcell, only: oneodx_e_w, oneody_n_w, oneodz_t_w

      USE functions, only: wall_at
      USE functions, only: east_of, north_of, top_of
      USE functions, only: west_of, south_of
      USE functions, only: im_of, jm_of, km_of

      USE fun_avg, only: avg_x_h, avg_y_h, avg_z_h

      USE geometry, only: odx_e, ody_n, odz
      USE geometry, only: ox
      USE geometry, only: ayz_w, axz_w, axy_w

      USE indices, only: i_of, j_of, k_of
      USE indices, only: kp1, im1, jm1

      USE visc_g, only: mu_gt
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! diffusion through faces of given ijk w-momentum cell
      DOUBLE PRECISION, INTENT(OUT) :: d_fe, d_fw
      DOUBLE PRECISION, INTENT(OUT) :: d_fn, d_fs
      DOUBLE PRECISION, INTENT(OUT) :: d_ft, d_fb
! ijk idnex
      INTEGER, INTENT(IN) :: ijk

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: imjk, ijmk, ijkm
      INTEGER :: i, j, k, kp, im, jm
      INTEGER :: ijkc, ijkt, ijke, ijkte, ijkw, ijkwt
      INTEGER :: ijkn, ijktn, ijks, ijkst

! length terms
      DOUBLE PRECISION :: C_AE, C_AW, C_AN, C_AS, C_AT, C_AB
!---------------------------------------------------------------------//

      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      KP = KP1(K)
      IM = IM1(I)
      JM = JM1(J)

      IJKT = TOP_OF(IJK)
      IF (WALL_AT(IJK)) THEN
         IJKC = IJKT
      ELSE
         IJKC = IJK
      ENDIF
      IJKE = EAST_OF(IJK)
      IJKTE = EAST_OF(IJKT)
      IJKN = NORTH_OF(IJK)
      IJKTN = NORTH_OF(IJKT)
      IJKW = WEST_OF(IJK)
      IJKWT = TOP_OF(IJKW)
      IJKS = SOUTH_OF(IJK)
      IJKST = TOP_OF(IJKS)

      IF(CUT_W_TREATMENT_AT(IJK)) THEN
         C_AE = ONEoDX_E_W(IJK)
         C_AW = ONEoDX_E_W(IMJK)
         C_AN = ONEoDY_N_W(IJK)
         C_AS = ONEoDY_N_W(IJMK)
         C_AT = ONEoDZ_T_W(IJK)
         C_AB = ONEoDZ_T_W(IJKM)
      ELSE
         C_AE = ODX_E(I)
         C_AW = ODX_E(IM)
         C_AN = ODY_N(J)
         C_AS = ODY_N(JM)
         C_AT = ODZ(KP)
         C_AB = ODZ(K)
      ENDIF

! East face (i+1/2, j, k+1/2)
      D_Fe = AVG_Z_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),&
                       AVG_X_H(MU_GT(IJKT),MU_GT(IJKTE),I),K)*&
               C_AE*AYZ_W(IJK)
! West face (i-1/2, j, k+1/2)
      D_Fw = AVG_Z_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),&
                      AVG_X_H(MU_GT(IJKWT),MU_GT(IJKT),IM),K)*&
                C_AW*AYZ_W(IMJK)

! North face (i, j+1/2, k+1/2)
      D_Fn = AVG_Z_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),&
                        AVG_Y_H(MU_GT(IJKT),MU_GT(IJKTN),J),K)*&
                C_AN*AXZ_W(IJK)
! South face (i, j-1/2, k+1/2)
      D_Fs = AVG_Z_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),&
                        AVG_Y_H(MU_GT(IJKST),MU_GT(IJKT),JM),K)*&
                C_AS*AXZ_W(IJMK)


! Top face (i, j, k+1)
      D_Ft = MU_GT(IJKT)*OX(I)*C_AT*AXY_W(IJK)
! Bottom face (i, j, k)
      D_Fb = MU_GT(IJK)*OX(I)*C_AB*AXY_W(IJKM)

      RETURN
      END SUBROUTINE GET_WCELL_GDIFF_TERMS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_W_g0                                            C
!  Purpose: Determine convection diffusion terms for W_g momentum eqs. C
!  The off-diagonal coefficients calculated here must be positive.     C
!  The center coefficient and the source vector are negative. See      C
!  source_w_g.                                                         C
!  Implement FOUP discretization                                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-APR-96  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STORE_A_W_G0(A_W_G, IER)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE functions, only: flow_at_t
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of

      USE param, only: dimension_3
      USE param1, only: zero
      USE matrix, only: e, w, n, s, t, b

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_U_g
      DOUBLE PRECISION, INTENT(INOUT) :: A_W_g(DIMENSION_3, -3:3)
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
      INTEGER :: IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
! Face mass flux
      DOUBLE PRECISION :: flux_e, flux_w, flux_n, flux_s
      DOUBLE PRECISION :: flux_t, flux_b
! Diffusion parameter
      DOUBLE PRECISION :: D_fe, d_fw, d_fn, d_fs, d_ft, d_fb

!---------------------------------------------------------------------//

!$omp     parallel do default(none)                                &
!$omp     private(IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,         &
!$omp             D_fe, d_fw, d_fn, d_fs, d_ft, d_fb,              &
!$omp             flux_e, flux_w, flux_n, flux_s, flux_t, flux_b)  &
!$omp     shared(ijkstart3, ijkend3, a_w_g)

      DO IJK = ijkstart3, ijkend3

         IF (FLOW_AT_T(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_WCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, ijk)

            CALL GET_WCELL_GDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, ijk)

            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

! East face (i+1/2, j, k+1/2)
            IF (Flux_e >= ZERO) THEN
               A_W_G(IJK,E) = D_Fe
               A_W_G(IPJK,W) = D_Fe + Flux_e
            ELSE
               A_W_G(IJK,E) = D_Fe - Flux_e
               A_W_G(IPJK,W) = D_Fe
            ENDIF
! West face (i-1/2, j, k+1/2)
            IF (.NOT.FLOW_AT_T(IMJK)) THEN
               IF (Flux_w >= ZERO) THEN
                  A_W_G(IJK,W) = D_Fw + Flux_w
               ELSE
                  A_W_G(IJK,W) = D_Fw
               ENDIF
            ENDIF


! North face (i, j+1/2, k+1/2)
            IF (Flux_n >= ZERO) THEN
               A_W_G(IJK,N) = D_Fn
               A_W_G(IJPK,S) = D_Fn + Flux_n
            ELSE
               A_W_G(IJK,N) = D_Fn - Flux_n
               A_W_G(IJPK,S) = D_Fn
            ENDIF
! South face (i, j-1/2, k+1/2)
            IF (.NOT.FLOW_AT_T(IJMK)) THEN
              IF (Flux_s >= ZERO) THEN
                  A_W_G(IJK,S) = D_Fs + Flux_s
               ELSE
                  A_W_G(IJK,S) = D_Fs
               ENDIF
            ENDIF


! Top face (i, j, k+1)
            IF (Flux_T >= ZERO) THEN
               A_W_G(IJK,T) = D_Ft
               A_W_G(IJKP,B) = D_Ft + Flux_t
            ELSE
               A_W_G(IJK,T) = D_Ft - Flux_t
               A_W_G(IJKP,B) = D_Ft
            ENDIF
! Bottom face (i, j, k)
            IF (.NOT.FLOW_AT_T(IJKM)) THEN
               IF (Flux_b >= ZERO) THEN
                  A_W_G(IJK,B) = D_Fb + Flux_b
               ELSE
                  A_W_G(IJK,B) = D_Fb
               ENDIF
            ENDIF
         ENDIF   ! end if (flow_at_t)
      ENDDO   ! end do (ijk)
!$omp end parallel do

      RETURN
      END SUBROUTINE STORE_A_W_G0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_W_GDC                                           C
!  Purpose: Use deferred correction method to solve the w-momentum     C
!  equation. This method combines first order upwind and a user        C
!  specified higher order method                                       C
!                                                                      C
!  Author: C. GUENTHER                                Date: 8-APR-99   C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STORE_A_W_GDC(B_M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE discretization, only: fpfoi_of

      USE fldvar, only: w_g

      USE function3, only: funijk3
      USE functions, only: flow_at_t
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of

      USE indices, only: i_of, j_of, k_of

      USE matrix, only: e, w, n, s, t, b

      USE param, only: dimension_3
      USE param1, only: zero, half

      USE run, only: discretize, fpfoi
      USE sendrecv3, only: send_recv3

      USE tmp_array, only: U => Array1, V => Array2, WW => Array3
      USE tmp_array, only: tmp4
      USE tmp_array, only: lock_tmp_array, unlock_tmp_array
      USE tmp_array, only: lock_tmp4_array, unlock_tmp4_array

      USE xsi, only: calc_xsi
      USE xsi_array, only: xsi_e, xsi_n, xsi_t
      USE xsi_array, only: lock_xsi_array, unlock_xsi_array
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IPJK, IMJK, IJMK, IJPK, IJKM, IJKP
      INTEGER :: IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER :: IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
! indication for shear
      INTEGER :: incr
! deferred corrction contribution form high order method
      DOUBLE PRECISION :: MOM_HO
! low order approximation
      DOUBLE PRECISION :: MOM_LO
! convection factor at each face
      DOUBLE PRECISION :: flux_e, flux_w, flux_n, flux_s
      DOUBLE PRECISION :: flux_t, flux_b
! deferred correction contributions from each face
      DOUBLE PRECISION :: EAST_DC
      DOUBLE PRECISION :: WEST_DC
      DOUBLE PRECISION :: NORTH_DC
      DOUBLE PRECISION :: SOUTH_DC
      DOUBLE PRECISION :: TOP_DC
      DOUBLE PRECISION :: BOTTOM_DC

! temporary use of global arrays:
! array1 (locally u)  - the x directional velocity
!      DOUBLE PRECISION :: U(DIMENSION_3)
! array2 (locally v)  - the y directional velocity
!      DOUBLE PRECISION :: V(DIMENSION_3)
! array3 (locally ww) - the z directional velocity
!      DOUBLE PRECISION :: WW(DIMENSION_3)
!---------------------------------------------------------------------//

      call lock_tmp4_array
      call lock_tmp_array
      call lock_xsi_array

      CALL GET_WCELL_GVTERMS(U, V, WW)

! Send recv the third ghost layer
      IF (FPFOI) THEN
         Do IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = W_G(IJK)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF


! shear indicator: 
      incr=0
      CALL CALC_XSI (DISCRETIZE(5), W_G, U, V, WW, XSI_E, XSI_N,&
         XSI_T, incr)


!!!$omp      parallel do                                             &
!!!$omp&     private( I, J, K, IJK, IPJK, IMJK, IJPK, IJMK,          &
!!!$omp&             flux_e, flux_w, flux_n, flux_s, flux_b, flux_t, &
!!!$omp&             MOM_HO, MOM_LO, EAST_DC, WEST_DC, NORTH_DC,     &
!!!$omp&             SOUTH_DC, TOP_DC, BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3

         IF (FLOW_AT_T(IJK)) THEN

! Calculate convection fluxes through each of the faces
            CALL GET_WCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, ijk)

            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)

! Third Ghost layer information
            IPPP  = IP_OF(IP_OF(IPJK))
            IPPP4 = funijk3(I_OF(IPPP), J_OF(IPPP), K_OF(IPPP))
            IMMM  = IM_OF(IM_OF(IMJK))
            IMMM4 = funijk3(I_OF(IMMM), J_OF(IMMM), K_OF(IMMM))
            JPPP  = JP_OF(JP_OF(IJPK))
            JPPP4 = funijk3(I_OF(JPPP), J_OF(JPPP), K_OF(JPPP))
            JMMM  = JM_OF(JM_OF(IJMK))
            JMMM4 = funijk3(I_OF(JMMM), J_OF(JMMM), K_OF(JMMM))
            KPPP  = KP_OF(KP_OF(IJKP))
            KPPP4 = funijk3(I_OF(KPPP), J_OF(KPPP), K_OF(KPPP))
            KMMM  = KM_OF(KM_OF(IJKM))
            KMMM4 = funijk3(I_OF(KMMM), J_OF(KMMM), K_OF(KMMM))


! East face (i+1/2, j, k+1/2)
            IF(U(IJK) >= ZERO)THEN
               MOM_LO = W_G(IJK)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_G(IPJK), W_G(IJK), &
                                   W_G(IMJK), W_G(IM_OF(IMJK)))
            ELSE
               MOM_LO = W_G(IPJK)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_G(IJK), W_G(IPJK), &
                                   W_G(IP_OF(IPJK)), TMP4(IPPP4))
            ENDIF
            IF (.NOT. FPFOI) MOM_HO = XSI_E(IJK)*W_G(IPJK)+ &
                                      (1.0-XSI_E(IJK))*W_G(IJK)
            EAST_DC = Flux_e*(MOM_LO-MOM_HO)


! North face (i, j+1/2, k+1/2)
            IF(V(IJK) >= ZERO)THEN
               MOM_LO = W_G(IJK)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_G(IJPK), W_G(IJK), &
                                   W_G(IJMK), W_G(JM_OF(IJMK)))
            ELSE
               MOM_LO = W_G(IJPK)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_G(IJK), W_G(IJPK), &
                                   W_G(JP_OF(IJPK)), TMP4(JPPP4))
            ENDIF
            IF (.NOT. FPFOI) MOM_HO = XSI_N(IJK)*W_G(IJPK)+ &
                                   (1.0-XSI_N(IJK))*W_G(IJK)
            NORTH_DC = Flux_n*(MOM_LO-MOM_HO)


! Top face (i, j, k+1)
            IF(WW(IJK) >= ZERO)THEN
               MOM_LO = W_G(IJK)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_G(IJKP), W_G(IJK), &
                                   W_G(IJKM), W_G(KM_OF(IJKM)))
            ELSE
               MOM_LO = W_G(IJKP)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_G(IJK), W_G(IJKP), &
                                   W_G(KP_OF(IJKP)), TMP4(KPPP4))
            ENDIF
            IF (.NOT. FPFOI) MOM_HO = XSI_T(IJK)*W_G(IJKP)+ &
                                      (1.0-XSI_T(IJK))*W_G(IJK)
            TOP_DC = Flux_t*(MOM_LO-MOM_HO)


! West face (i-1/2, j, k+1/2)
            IF(U(IMJK) >= ZERO)THEN
               MOM_LO = W_G(IMJK)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_G(IJK), W_G(IMJK), &
                                   W_G(IM_OF(IMJK)), TMP4(IMMM4))
            ELSE
               MOM_LO = W_G(IJK)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_G(IMJK), W_G(IJK), &
                                   W_G(IPJK), W_G(IP_OF(IPJK)))
            ENDIF
            IF (.NOT. FPFOI) MOM_HO = XSI_E(IMJK)*W_G(IJK)+ &
                                      (1.0-XSI_E(IMJK))*W_G(IMJK)
            WEST_DC = Flux_w*(MOM_LO-MOM_HO)


! South face (i, j-1/2, k+1/2)
            IF(V(IJMK) >= ZERO)THEN
               MOM_LO = W_G(IJMK)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_G(IJK), W_G(IJMK), &
                                   W_G(JM_OF(IJMK)), TMP4(JMMM4))
            ELSE
               MOM_LO = W_G(IJK)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_G(IJMK), W_G(IJK), &
                                   W_G(IJPK), W_G(JP_OF(IJPK)))
            ENDIF
            IF (.NOT. FPFOI) MOM_HO = XSI_N(IJMK)*W_G(IJK)+ &
                                      (1.0-XSI_N(IJMK))*W_G(IJMK)
            SOUTH_DC = Flux_s*(MOM_LO-MOM_HO)


! Bottom face (i, j, k)
            IF(WW(IJK) >= ZERO)THEN
               MOM_LO = W_G(IJKM)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_G(IJK), W_G(IJKM), &
                                   W_G(KM_OF(IJKM)), TMP4(KMMM4))
            ELSE
               MOM_LO = W_G(IJK)
               IF (FPFOI) MOM_HO = FPFOI_OF(W_G(IJKM), W_G(IJK), &
                                   W_G(IJKP), W_G(KP_OF(IJKP)))
            ENDIF
            IF (.NOT. FPFOI) MOM_HO = XSI_T(IJKM)*W_G(IJK)+ &
                                      (1.0-XSI_T(IJKM))*W_G(IJKM)
            BOTTOM_DC = Flux_b*(MOM_LO-MOM_HO)


! CONTRIBUTION DUE TO DEFERRED CORRECTION
            B_M(IJK) = B_M(IJK)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC&
                                +BOTTOM_DC-TOP_DC

         ENDIF   ! end if flow_at_t
      ENDDO   ! end do ijk

      call unlock_tmp4_array
      call unlock_tmp_array
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_W_GDC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_W_g1                                            C
!  Purpose: Determine convection diffusion terms for W_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive.     C
!  The center coefficient and the source vector are negative.          C
!  Implements higher order discretization.                             C
!  See source_w_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAR-97  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STORE_A_W_G1(A_W_G)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE fldvar, only: w_g

      USE functions, only: flow_at_t
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of

      USE param, only: dimension_3
      USE param1, only: one

      USE matrix, only: e, w, n, s, t, b

      USE run, only: discretize

      USE tmp_array, only: U => Array1, V => Array2, WW => Array3
      USE tmp_array, only: lock_tmp_array, unlock_tmp_array

      USE xsi, only: calc_xsi
      USE xsi_array, only: xsi_e, xsi_n, xsi_t
      USE xsi_array, only: lock_xsi_array, unlock_xsi_array

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_W_g
      DOUBLE PRECISION, INTENT(INOUT) :: A_W_g(DIMENSION_3, -3:3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, IPJK, IMJK, IJPK, IJMK, IJKP, IJKM
! indicator for shear
      INTEGER :: incr
! Diffusion parameter
      DOUBLE PRECISION :: d_fe, d_fw, d_fn, d_fs, d_ft, d_fb
! Face mass flux
      DOUBLE PRECISION :: Flux_e, flux_w, flux_n, flux_s
      DOUBLE PRECISION :: flux_t, flux_b

! temporary use of global arrays:
! array1 (locally u)  - the x directional velocity
!      DOUBLE PRECISION :: U(DIMENSION_3)
! array2 (locally v)  - the y directional velocity
!      DOUBLE PRECISION :: V(DIMENSION_3)
! array3 (locally ww) - the z directional velocity
!      DOUBLE PRECISION :: WW(DIMENSION_3)
!---------------------------------------------------------------------//

      call lock_tmp_array
      call lock_xsi_array

      CALL GET_WCELL_GVTERMS(U, V, WW)


! shear indicator:
      incr=0
      CALL CALC_XSI (DISCRETIZE(5), W_G, U, V, WW, XSI_E, XSI_N,&
         XSI_T, incr)

!!!$omp      parallel do                                                 &
!!!$omp&     private(IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,            &
!!!$omp&             d_fe, d_fw, d_fn, d_fs, d_ft, d_fb,                 &
!!!$omp&             flux_e, flux_w, flux_n, flux_s, flux_t, flux_b)
      DO IJK = ijkstart3, ijkend3

         IF (FLOW_AT_T(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_WCELL_GCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, ijk)
            CALL GET_WCELL_GDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, ijk)

            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

! East face (i+1/2, j, k+1/2)
            A_W_G(IJK,E) = D_Fe - XSI_E(IJK)*Flux_e
            A_W_G(IPJK,W) = D_Fe + (ONE - XSI_E(IJK))*Flux_e
! West face (i-1/2, j, k+1/2)
            IF (.NOT.FLOW_AT_T(IMJK)) THEN
               A_W_G(IJK,W) = D_Fw + (ONE - XSI_E(IMJK))*Flux_w
            ENDIF


! North face (i, j+1/2, k+1/2)
            A_W_G(IJK,N) = D_Fn - XSI_N(IJK)*Flux_n
            A_W_G(IJPK,S) = D_Fn + (ONE - XSI_N(IJK))*Flux_n
! South face (i, j-1/2, k+1/2)
            IF (.NOT.FLOW_AT_T(IJMK)) THEN
               A_W_G(IJK,S) = D_Fs + (ONE - XSI_N(IJMK))*Flux_s
            ENDIF


! Top face (i, j, k+1)
            A_W_G(IJK,T) = D_Ft - XSI_T(IJK)*Flux_t
            A_W_G(IJKP,B) = D_Ft + (ONE - XSI_T(IJK))*Flux_t
! Bottom face (i, j, k)
            IF (.NOT.FLOW_AT_T(IJKM)) THEN
              A_W_G(IJK,B) = D_Fb + (ONE - XSI_T(IJKM))*Flux_b
            ENDIF

         ENDIF   ! end if flow_at_t
      ENDDO   ! end do ijk

      call unlock_tmp_array
      call unlock_xsi_array

      RETURN
      END SUBROUTINE STORE_A_W_G1
