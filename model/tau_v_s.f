!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_Tau_V_s                                            C
!  Purpose: Cross terms in the gradient of stress in V_s momentum      c
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-DEC-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_TAU_V_S(TAU_V_S, IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
     USE param1, only: zero, half, one

! number of solids phases
      USE physprop, only: smax, mmax

! x,y,z-components of solids velocity
      USE fldvar, only: u_s, v_s, w_s
! solids phase particle bulk density and material density
! (needed for ep_s2.inc)
      USE fldvar, only: rop_s, ro_s

! viscous solids transport coefficients
      USE visc_s, only: mu_s, lambda_s
! trace of rate of strain tensor
      USE visc_s, only: trd_s

! dilute threshold      
      USE toleranc, only: dil_ep_s 
      
! kinetic theories
      USE run, only: kt_type_enum, ghd_2007
! runtime flag for treating the system as if shearing
      USE run, only: shear, V_sh

! primarily needed for function.inc
      USE compar
      USE geometry
      USE indices

! for sendrecv calls
      USE sendrecv

! for cutcell:
! wall velocity for slip conditions
      USE bc, only: bc_hw_s, bc_uw_s, bc_vw_s, bc_ww_s
      USE bc, only: bc_type
      USE quadric
      USE cutcell

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! TAU_V_s 
      DOUBLE PRECISION, INTENT(OUT) :: TAU_V_s(DIMENSION_3, DIMENSION_M)
! Error index 
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices 
      INTEGER :: I, J, K, IJK, IJKN, JP, IM,  KM, IJPK, IJMK,& 
                 IJKE, IJKNE, IJKW, IJKNW, IMJPK, IMJK, IJKT,& 
                 IJKTN, IJKB, IJKBN, IJKM, IJPKM 
! Phase index 
      INTEGER :: M, L 

! Average volume fraction 
      DOUBLE PRECISION :: EPSA, EPStmp
 
! Source terms (Surface) 
      DOUBLE PRECISION :: Sbv, Ssx, Ssy, Ssz 
! Shearing variables
      DOUBLE PRECISION :: Source_diff,Diffco_e,Diffco_w
! error message 
      CHARACTER*80     LINE

! for cartesian grid implementation:
      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      LOGICAL :: U_NODE_AT_NE,U_NODE_AT_NW,U_NODE_AT_SE,U_NODE_AT_SW
      LOGICAL :: W_NODE_AT_TN,W_NODE_AT_TS,W_NODE_AT_BN,W_NODE_AT_BS
      DOUBLE PRECISION :: U_SUM,W_SUM,X_SUM,Y_SUM,Z_SUM,Uc,Wc
      DOUBLE PRECISION :: Xuc,Yuc,Zuc,Xwc,Ywc,Zwc,Nxu,Nyu,Nzu,Nxw,Nyw,Nzw
      DOUBLE PRECISION :: dudy_at_E,dudy_at_W
      DOUBLE PRECISION :: dwdy_at_T,dwdy_at_B
      DOUBLE PRECISION :: Xi,Yi,Zi,Ui,Vi,Wi,Sx,Sy,Sz
      DOUBLE PRECISION :: x_circle,y_circle,angle
      DOUBLE PRECISION :: MU_S_CUT,SSX_CUT,SSZ_CUT
      DOUBLE PRECISION :: UW_s,VW_s,WW_s
      INTEGER :: N_SUM
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT  
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------
      DO M = 1, MMAX 
        IF(KT_TYPE_ENUM == GHD_2007 .AND. M /= MMAX) CYCLE


!!$omp  parallel do private( IJK, I, IJKE, EPSA, EPStmp,  J,  K, KM,  &
!!$omp& JP,IM,IJPK,IJKW,IJKNW,IMJPK,IJKTN,IJKBN,IJPKM, &
!!$omp&  IMJK,IJKN,IJKNE,IJMK,IJKT,  &
!!$omp&  IJKB,IJKM, &
!!$omp&  SBV,  SSX,SSY,   SSZ,&
!!$omp&  Source_diff, Diffco_e,Diffco_w) 
!!$omp&  schedule(static)

        DO IJK = IJKSTART3, IJKEND3

! Skip walls where some values are undefined.
          IF(WALL_AT(IJK)) cycle

          J = J_OF(IJK) 
          IJKN = NORTH_OF(IJK)

          IF (KT_TYPE_ENUM == GHD_2007) THEN
            EPStmp = ZERO
            DO L = 1, SMAX
              EPStmp = EPStmp + AVG_Y(EP_S(IJK,L),EP_S(IJKN,L),J) 
            ENDDO
            EPSA = EPStmp
          ELSE
            EPSA = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J) 
          ENDIF

          IF ( .NOT.SIP_AT_N(IJK) .AND. EPSA>DIL_EP_S) THEN 
            JP = JP1(J) 
            I = I_OF(IJK) 
            IM = IM1(I) 
            K = K_OF(IJK) 
            KM = KM1(K) 
            IJPK = JP_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
            IJKW = WEST_OF(IJK) 
            IJKNW = NORTH_OF(IJKW) 
            IMJPK = IM_OF(IJPK) 
            IMJK = IM_OF(IJK) 
            IJKT = TOP_OF(IJK) 
            IJKTN = NORTH_OF(IJKT) 
            IJKB = BOTTOM_OF(IJK) 
            IJKBN = NORTH_OF(IJKB) 
            IJKM = KM_OF(IJK) 
            IJPKM = JP_OF(IJKM) 

            IF((.NOT.CARTESIAN_GRID).OR.(CG_SAFE_MODE(4)==1)) THEN
! NON CARTESIAN GRID CASE
! ---------------------------------------------------------------->>>
! Surface forces

! bulk viscosity term
              SBV = (LAMBDA_S(IJKN,M)*TRD_S(IJKN,M)-&
                     LAMBDA_S(IJK,M)*TRD_S(IJK,M))*AXZ(IJK) 

! shear stress terms
              SSX = AVG_Y_H(AVG_X_H(MU_S(IJK,M),MU_S(IJKE,M),I),&
                            AVG_X_H(MU_S(IJKN,M),MU_S(IJKNE,M),I),J)*&
                    (U_S(IJPK,M)-U_S(IJK,M))*ODY_N(J)*AYZ_V(IJK) - &
                    AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJK,M),IM),&
                            AVG_X_H(MU_S(IJKNW,M),MU_S(IJKN,M),IM),J)*&
                    (U_S(IMJPK,M)-U_S(IMJK,M))*ODY_N(J)*AYZ_V(IMJK) 
              SSY = MU_S(IJKN,M)*(V_S(IJPK,M)-V_S(IJK,M))*ODY(JP)*AXZ_V(IJK) - &
                    MU_S(IJK,M)*(V_S(IJK,M)-V_S(IJMK,M))*ODY(J)*AXZ_V(IJMK)
!?              IF(DO_K) THEN
              SSZ = AVG_Y_H(AVG_Z_H(MU_S(IJK,M),MU_S(IJKT,M),K),&
                            AVG_Z_H(MU_S(IJKN,M),MU_S(IJKTN,M),K),J)*&
                    (W_S(IJPK,M)-W_S(IJK,M))*ODY_N(J)*AXY_V(IJK) - &
                    AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJK,M),KM),&
                            AVG_Z_H(MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)*&
                    (W_S(IJPKM,M)-W_S(IJKM,M))*ODY_N(J)*AXY_V(IJKM) 
!?              ELSE
!?                SSZ = ZERO
!?              ENDIF
! ----------------------------------------------------------------<<<

            ELSE

! CARTESIAN GRID CASE
! ---------------------------------------------------------------->>>
! Surface forces

! bulk viscosity term
              SBV = (LAMBDA_S(IJKN,M)*TRD_S(IJKN,M)) * AXZ_V(IJK) - & 
                    (LAMBDA_S(IJK,M) *TRD_S(IJK,M) ) * AXZ_V(IJMK) 

! shear stress terms
              IF(.NOT.CUT_V_CELL_AT(IJK)) THEN
                SSX = AVG_Y_H(AVG_X_H(MU_S(IJK,M),MU_S(IJKE,M),I),&
                              AVG_X_H(MU_S(IJKN,M),MU_S(IJKNE,M),I),J)*&
                      (U_S(IJPK,M)-U_S(IJK,M))*ONEoDY_N_U(IJK)*AYZ_V(IJK) - &
                      AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJK,M),IM),&
                              AVG_X_H(MU_S(IJKNW,M),MU_S(IJKN,M),IM),J)*&
                      (U_S(IMJPK,M)-U_S(IMJK,M))*ONEoDY_N_U(IMJK)*AYZ_V(IMJK)
                SSY = MU_S(IJKN,M)*(V_S(IJPK,M)-V_S(IJK,M))*&
                      ONEoDY_N_V(IJK)*AXZ_V(IJK) - &
                      MU_S(IJK,M)*(V_S(IJK,M)-V_S(IJMK,M))*&
                      ONEoDY_N_V(IJMK)*AXZ_V(IJMK) 
                IF(DO_K) THEN
                  SSZ = AVG_Y_H(AVG_Z_H(MU_S(IJK,M),MU_S(IJKT,M),K),&
                                AVG_Z_H(MU_S(IJKN,M),MU_S(IJKTN,M),K),J)*&
                        (W_S(IJPK,M)-W_S(IJK,M))*ONEoDY_N_W(IJK)*AXY_V(IJK) - &
                        AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJK,M),KM),&
                               AVG_Z_H(MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)*&
                        (W_S(IJPKM,M)-W_S(IJKM,M))*ONEoDY_N_W(IJKM)*AXY_V(IJKM) 
                ELSE
                  SSZ = ZERO
                ENDIF


              ELSE   ! CUT_CELL
                BCV = BC_V_ID(IJK)
                IF(BCV > 0 ) THEN
                  BCT = BC_TYPE(BCV)
                ELSE
                  BCT = 'NONE'
                ENDIF

                SELECT CASE (BCT) 
                  CASE ('CG_NSW')
                     CUT_TAU_VS = .TRUE.
                     NOC_VS     = .TRUE.
                     UW_s = ZERO
                     VW_s = ZERO
                     WW_s = ZERO
                  CASE ('CG_FSW')
                     CUT_TAU_VS = .FALSE.
                     NOC_VS     = .FALSE.
                     UW_s = ZERO
                     VW_s = ZERO
                     WW_s = ZERO
                  CASE('CG_PSW')
                     IF(BC_HW_S(BC_V_ID(IJK),M)==UNDEFINED) THEN   ! same as NSW
                        CUT_TAU_VS = .TRUE.
                        NOC_VS     = .TRUE.
                        UW_s = BC_UW_S(BCV,M)
                        VW_s = BC_VW_S(BCV,M)
                        WW_s = BC_WW_S(BCV,M)
                     ELSEIF(BC_HW_S(BC_V_ID(IJK),M)==ZERO) THEN   ! same as FSW
                        CUT_TAU_VS = .FALSE.
                        NOC_VS     = .FALSE.
                        UW_s = ZERO
                        VW_s = ZERO
                        WW_s = ZERO
                     ELSE                              ! partial slip
                        CUT_TAU_VS = .FALSE.
                        NOC_VS     = .FALSE.
                     ENDIF
                  CASE ('NONE')
                     TAU_V_S(IJK,M) = ZERO 
                     CYCLE    
                END SELECT 

                IF(CUT_TAU_VS) THEN
                  MU_S_CUT = (VOL(IJK)*MU_S(IJK,M) + &
                     VOL(IJPK)*MU_S(IJKN,M))/(VOL(IJK) + VOL(IJPK))
                ELSE
                  MU_S_CUT = ZERO
                ENDIF

! SSX:
                U_NODE_AT_NE = ((.NOT.BLOCKED_U_CELL_AT(IJPK)).AND.(.NOT.WALL_U_AT(IJPK)))
                U_NODE_AT_SE = ((.NOT.BLOCKED_U_CELL_AT(IJK)).AND.(.NOT.WALL_U_AT(IJK)))
                U_NODE_AT_NW = ((.NOT.BLOCKED_U_CELL_AT(IMJPK)).AND.(.NOT.WALL_U_AT(IMJPK)))
                U_NODE_AT_SW = ((.NOT.BLOCKED_U_CELL_AT(IMJK)).AND.(.NOT.WALL_U_AT(IMJK)))

                IF(U_NODE_AT_NE.AND.U_NODE_AT_SE) THEN
                   Ui = HALF * (U_S(IJPK,M) + U_S(IJK,M))
                   Xi = HALF * (X_U(IJPK) + X_U(IJK))
                   Yi = HALF * (Y_U(IJPK) + Y_U(IJK))
                   Zi = HALF * (Z_U(IJPK) + Z_U(IJK))
                   Sx = X_U(IJPK) - X_U(IJK)
                   Sy = Y_U(IJPK) - Y_U(IJK)
                   Sz = Z_U(IJPK) - Z_U(IJK)
                   CALL GET_DEL_H(IJK,'V_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
                   dudy_at_E =  (U_S(IJPK,M) - U_S(IJK,M)) * ONEoDY_N_U(IJK)
                   IF(NOC_VS) dudy_at_E = dudy_at_E - ((Ui-UW_s) * &
                                 ONEoDY_N_U(IJK)/DEL_H*(Sx*Nx+Sz*Nz))
                ELSE
                   dudy_at_E =  ZERO
                ENDIF

                IF(U_NODE_AT_NW.AND.U_NODE_AT_SW) THEN
                   Ui = HALF * (U_S(IMJPK,M) + U_S(IMJK,M))
                   Xi = HALF * (X_U(IMJPK) + X_U(IMJK))
                   Yi = HALF * (Y_U(IMJPK) + Y_U(IMJK))
                   Zi = HALF * (Z_U(IMJPK) + Z_U(IMJK))
                   Sx = X_U(IMJPK) - X_U(IMJK)
                   Sy = Y_U(IMJPK) - Y_U(IMJK)
                   Sz = Z_U(IMJPK) - Z_U(IMJK)
                   CALL GET_DEL_H(IJK,'V_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
                   dudy_at_W =  (U_S(IMJPK,M)-U_S(IMJK,M))*ONEoDY_N_U(IMJK)
                   IF(NOC_VS) dudy_at_W = dudy_at_W - ((Ui-UW_s) * &
                                 ONEoDY_N_U(IMJK)/DEL_H*(Sx*Nx+Sz*Nz))
                ELSE
                   dudy_at_W =  ZERO
                ENDIF

                IF(U_NODE_AT_SE) THEN
                   CALL GET_DEL_H(IJK,'V_MOMENTUM',X_U(IJK),Y_U(IJK),&
                                  Z_U(IJK),Del_H,Nx,Ny,Nz)
                   SSX_CUT = - MU_S_CUT * (U_S(IJK,M) - UW_s) / &
                      DEL_H * (Ny*Nx) * Area_V_CUT(IJK)        
                ELSE
                   SSX_CUT =  ZERO
                ENDIF

                SSX = AVG_Y_H(AVG_X_H(MU_S(IJK,M),MU_S(IJKE,M),I),&
                              AVG_X_H(MU_S(IJKN,M),MU_S(IJKNE,M),I),J)*&
                      dudy_at_E*AYZ_V(IJK) - &
                      AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJK,M),IM),&
                              AVG_X_H(MU_S(IJKNW,M),MU_S(IJKN,M),IM),J)*&
                      dudy_at_W*AYZ_V(IMJK) + SSX_CUT

! SSY:
                CALL GET_DEL_H(IJK,'V_MOMENTUM',X_V(IJK),Y_V(IJK),&
                               Z_V(IJK),Del_H,Nx,Ny,Nz)
                SSY = MU_S(IJKN,M)*(V_S(IJPK,M)-V_S(IJK,M))*ONEoDY_N_V(IJK)*AXZ_V(IJK) - &
                      MU_S(IJK,M)*(V_S(IJK,M)-V_S(IJMK,M))*ONEoDY_N_V(IJMK)*AXZ_V(IJMK) &
                      - MU_S_CUT * (V_S(IJK,M) - VW_s) / DEL_H * (Ny**2) * Area_V_CUT(IJK)

! SSZ:
                IF(DO_K) THEN  
                   W_NODE_AT_TN = ((.NOT.BLOCKED_W_CELL_AT(IJPK)).AND.(.NOT.WALL_W_AT(IJPK)))
                   W_NODE_AT_TS = ((.NOT.BLOCKED_W_CELL_AT(IJK)).AND.(.NOT.WALL_W_AT(IJK)))
                   W_NODE_AT_BN = ((.NOT.BLOCKED_W_CELL_AT(IJPKM)).AND.(.NOT.WALL_W_AT(IJPKM)))
                   W_NODE_AT_BS = ((.NOT.BLOCKED_W_CELL_AT(IJKM)).AND.(.NOT.WALL_W_AT(IJKM)))

                   IF(W_NODE_AT_TN.AND.W_NODE_AT_TS) THEN
                      Wi = HALF * (W_S(IJPK,M) + W_S(IJK,M))
                      Xi = HALF * (X_W(IJPK) + X_W(IJK))
                      Yi = HALF * (Y_W(IJPK) + Y_W(IJK))
                      Zi = HALF * (Z_W(IJPK) + Z_W(IJK))
                      Sx = X_W(IJPK) - X_W(IJK)
                      Sy = Y_W(IJPK) - Y_W(IJK)
                      Sz = Z_W(IJPK) - Z_W(IJK)
                      CALL GET_DEL_H(IJK,'V_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
                      dwdy_at_T =  (W_S(IJPK,M) - W_S(IJK,M)) * ONEoDY_N_W(IJK)
                      IF(NOC_VS) dwdy_at_T = dwdy_at_T - ((Wi-WW_s) * &
                                    ONEoDY_N_W(IJK)/DEL_H*(Sx*Nx+Sz*Nz))
                   ELSE
                      dwdy_at_T =  ZERO
                   ENDIF

                   IF(W_NODE_AT_BN.AND.W_NODE_AT_BS) THEN
                      Wi = HALF * (W_S(IJPKM,M) + W_S(IJKM,M))
                      Xi = HALF * (X_W(IJPKM) + X_W(IJKM))
                      Yi = HALF * (Y_W(IJPKM) + Y_W(IJKM))
                      Zi = HALF * (Z_W(IJPKM) + Z_W(IJKM))
                      Sx = X_W(IJPKM) - X_W(IJKM)
                      Sy = Y_W(IJPKM) - Y_W(IJKM)
                      Sz = Z_W(IJPKM) - Z_W(IJKM)
                      CALL GET_DEL_H(IJK,'V_MOMENTUM',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
                      dwdy_at_B =  (W_S(IJPKM,M) - W_S(IJKM,M)) * ONEoDY_N_W(IJKM)
                      IF(NOC_VS) dwdy_at_B = dwdy_at_B - ((Wi-WW_s) * &
                                    ONEoDY_N_W(IJKM)/DEL_H*(Sx*Nx+Sz*Nz))
                   ELSE
                      dwdy_at_B =  ZERO
                   ENDIF

                   IF(W_NODE_AT_TS) THEN
                      CALL GET_DEL_H(IJK,'V_MOMENTUM',X_W(IJK),Y_W(IJK),&
                                     Z_W(IJK),Del_H,Nx,Ny,Nz)
                      SSZ_CUT = - MU_S_CUT * (W_S(IJK,M) - WW_s) / &
                                DEL_H * (Ny*Nz) * Area_V_CUT(IJK)
                   ELSE
                      SSZ_CUT =  ZERO
                   ENDIF

                   SSZ = AVG_Y_H(AVG_Z_H(MU_S(IJK,M),MU_S(IJKT,M),K),&
                                 AVG_Z_H(MU_S(IJKN,M),MU_S(IJKTN,M),K),J)*&
                         dwdy_at_T*AXY_V(IJK) - &
                         AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJK,M),KM),&
                         AVG_Z_H(MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)*&
                         dwdy_at_B*AXY_V(IJKM) + SSZ_CUT  
                ELSE

                       SSZ = ZERO

                ENDIF  ! end if do_k

              ENDIF  ! end if/else cut_cell
! ----------------------------------------------------------------<<<
            ENDIF   ! end if/else cartesian grid

! Source terms from shear stress
            IF (SHEAR) THEN
               Diffco_e = AVG_Y_H((AVG_X_H(MU_S(IJK,m),MU_S(IJKE,m),I)),&
                                  (AVG_X_H(MU_S(IJKN,m),MU_S(IJKNE,m),I)),J)*AYZ_V(IJK)
               Diffco_w=AVG_Y_H((AVG_X_H(MU_S(IJK,m),MU_S(IJKW,m),I)),&
                                (AVG_X_H(MU_S(IJKN,m),MU_S(IJKNW,m),I)),J)*AYZ_V(IJKW)
               Source_diff=(2d0*V_sh/XLENGTH)*(Diffco_e-Diffco_w)
            ELSE
               Source_diff=0d0
            ENDIF

! Add the terms
            TAU_V_S(IJK,M) = SBV + SSX + SSY + SSZ + Source_diff

          ELSE 
            TAU_V_S(IJK,M) = ZERO 
          ENDIF   ! end if/else .not.sip_at_n .and. epsa>dil_ep_s
        ENDDO   ! end do ijk

      ENDDO 

      call send_recv(tau_v_s,2)
      RETURN  
      END SUBROUTINE CALC_TAU_V_S 


