calc_mflux.f                                                                                        0100644 0002444 0000146 00000014415 10247141615 011513  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_MFLUX(IER)                                        C
!  Purpose: Calculate the convection fluxes. Master routine.           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_MFLUX(IER) 
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE mflux
      USE physprop
      USE run
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      solids phase index 
      INTEGER          M 
! 
! 
!
!
      CALL CALC_MFLUX0 (U_g, V_g, W_g, ROP_gE, ROP_gN, ROP_gT, &
                        Flux_gE, Flux_gN, Flux_gT, IER) 
      DO M = 1, MMAX
        CALL CALC_MFLUX0 (U_s(1, M), V_s(1, M), W_s(1, M), &
	                  ROP_sE(1, M), ROP_sN(1, M), ROP_sT(1, M), &
                          Flux_sE(1, M), Flux_sN(1, M), Flux_sT(1, M), IER)  
      ENDDO
      
      RETURN  
      END SUBROUTINE CALC_MFLUX 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!  Module name: CALC_MFLUX0(U, V, W, ROP_E, ROP_N, ROP_T, Flux_E,      C
!                                         Flux_N, Flux_T, IER)         C
!  Purpose: Calculate the convection fluxes.                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_MFLUX0(U, V, W, ROP_E, ROP_N, ROP_T, Flux_E, Flux_N, Flux_T, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE geometry
      USE indices
      USE compar   
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
!
!                      Velocity components
      DOUBLE PRECISION U(DIMENSION_3), V(DIMENSION_3), W(DIMENSION_3) 
!
!                      Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION ROP_E(DIMENSION_3), ROP_N(DIMENSION_3), ROP_T(DIMENSION_3) 
!
!                      Convective mass fluxes
      DOUBLE PRECISION Flux_E(DIMENSION_3), Flux_N(DIMENSION_3), Flux_T(DIMENSION_3) 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, IMJK, IJMK, IJKM 
!-----------------------------------------------
      INCLUDE 'function.inc'
      

!
!  Interpolate the face value of density for calculating the convection fluxes 
!$omp  parallel do private( IJK, IMJK, IJMK, IJKM) &
!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLUID_AT(IJK)) THEN 
!
!         East face (i+1/2, j, k)
            Flux_E(IJK) = ROP_E(IJK)*AYZ(IJK)*U(IJK)
!
!         North face (i, j+1/2, k)
            Flux_N(IJK) = ROP_N(IJK)*AXZ(IJK)*V(IJK)
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               Flux_T(IJK) = ROP_T(IJK)*AXY(IJK)*W(IJK)
            ENDIF 
!
!         West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               Flux_E(IMJK) = ROP_E(IMJK)*AYZ(IMJK)*U(IMJK)
            ENDIF 
!
!         South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
              Flux_N(IJMK) = ROP_N(IJMK)*AXZ(IJMK)*V(IJMK)
            ENDIF 
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                 Flux_T(IJKM) = ROP_T(IJKM)*AXY(IJKM)*W(IJKM)
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE CALC_MFLUX0
                                                                                                                                                                                                                                                   conv_dif_phi.f                                                                                      0100644 0002444 0000146 00000105412 10247135154 012024  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:
!  CONV_DIF_Phi(Phi, Dif, Disc, Uf, Vf, Wf, Flux_E, Flux_N, Flux_T, M, A_m, B_m, IER)    C
!  Purpose: Determine convection diffusion terms for a sclar phi       C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C

!  The diffusion at the flow boundaries is prevented by setting the 
!  diffusion coefficients at boundary cells to zero and then using a 
!  harmonic average to calculate the boundary diffusivity.  The value
!  diffusivities at the boundaries are checked in check_data_30.  Ensure
!  that harmonic avergaing is used in this routine. 
!  See source_phi                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_DIF_PHI(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1
      USE run 
      USE geometry
      USE compar
      USE sendrecv
      Use xsi_array
      USE mpi_utility
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Scalar
      DOUBLE PRECISION Phi(DIMENSION_3)
!
!                      Gamma -- diffusion coefficient
      DOUBLE PRECISION Dif(DIMENSION_3)
!
!                      Discretizationindex
      INTEGER          Disc
!
!                      Velocity components
      DOUBLE PRECISION Uf(DIMENSION_3), Vf(DIMENSION_3), Wf(DIMENSION_3) 
!
!                      Mass flux components
      DOUBLE PRECISION Flux_E(DIMENSION_3), Flux_N(DIMENSION_3), Flux_T(DIMENSION_3) 
!
!                      Phase index
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Error index
      INTEGER          IER

!

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	IF DEFERRED CORRECTION IS USED WITH THE SCALAR TRANSPORT EQN.
!
	IF(DEF_COR)THEN
	  CALL CONV_DIF_PHI0(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER)
	  if (DISC > 1) CALL CONV_DIF_PHI_DC(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER)
	ELSE
!
!	NO DEFERRED CORRECTION IS USED WITH THE SCALAR TRANSPORT EQN.
!
	  IF (DISC == 0) THEN                        
            CALL CONV_DIF_PHI0(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER)
	  ELSE
            CALL CONV_DIF_PHI1(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER) 
          ENDIF
	ENDIF 
	
        CALL DIF_PHI_IS (DIF, A_M, B_M, M, IER)

        RETURN  
      END SUBROUTINE CONV_DIF_PHI 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:
!   CONV_DIF_Phi0(Phi, Dif, Disc, Uf, Vf, Wf, Flux_E,Flux_N,Flux_T, M, A_m, B_m, IER)
!  Purpose: Determine convection diffusion terms for Phi balance       C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_phi                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_DIF_PHI0(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE toleranc 
      USE run
      USE geometry
      USE compar
      USE sendrecv
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Scalar
      DOUBLE PRECISION Phi(DIMENSION_3)
!
!                      Gamma -- diffusion coefficient
      DOUBLE PRECISION Dif(DIMENSION_3)
!
!                      Discretizationindex
      INTEGER          Disc
!
!                      Velocity components
      DOUBLE PRECISION Uf(DIMENSION_3), Vf(DIMENSION_3), Wf(DIMENSION_3) 
!
!                      Mass flux components
      DOUBLE PRECISION Flux_E(DIMENSION_3), Flux_N(DIMENSION_3), Flux_T(DIMENSION_3) 
!
!                      Phase index
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          I,  J, K, IJK, IPJK, IJPK, IJKE, IJKN,&
                       IJKP, IJKT

      INTEGER          IMJK, IM, IJKW
      INTEGER          IJMK, JM, IJKS
      INTEGER          IJKM, KM, IJKB
!
!                      Face velocity
      DOUBLE PRECISION V_f
!
!                      Difusion parameter
      DOUBLE PRECISION D_f
!
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do                                              &
!$omp&     private(I,  J, K,  IJK,  IPJK, IJPK, IJKE, IJKN,         &
!$omp&             IJKP, IJKT,  V_f, D_f,                    &
!$omp&             IMJK, IM, IJKW,                                  &
!$omp&             IJMK, JM, IJKS,                                  &
!$omp&             IJKM, KM,  IJKB)                     
      DO IJK = ijkstart3, ijkend3
!
       I = I_OF(IJK)
       J = J_OF(IJK)
       K = K_OF(IJK)
!
         IF (FLUID_AT(IJK)) THEN 
!
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
!
!
!           East face (i+1/2, j, k)
            V_F = UF(IJK) 
            D_F = AVG_X_H(DIF(IJK),DIF(IJKE),I)*ODX_E(I)*AYZ(IJK) 
            IF (V_F >= ZERO) THEN 
               A_M(IJK,E,M) = D_F 
               A_M(IPJK,W,M) = D_F + FLUX_E(IJK) 
            ELSE 
               A_M(IJK,E,M) = D_F - FLUX_E(IJK) 
               A_M(IPJK,W,M) = D_F 
            ENDIF 
!
!
!           North face (i, j+1/2, k)
            V_F = VF(IJK) 
            D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*AXZ(IJK) 
            IF (V_F >= ZERO) THEN 
               A_M(IJK,N,M) = D_F 
               A_M(IJPK,S,M) = D_F + FLUX_N(IJK) 
            ELSE 
               A_M(IJK,N,M) = D_F - FLUX_N(IJK) 
               A_M(IJPK,S,M) = D_F 
            ENDIF 
!
!           Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               V_F = WF(IJK) 
               D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*OX(I)*ODZ_T(K)*AXY(IJK) 
               IF (V_F >= ZERO) THEN 
                  A_M(IJK,T,M) = D_F 
                  A_M(IJKP,B,M) = D_F + FLUX_T(IJK) 
               ELSE 
                  A_M(IJK,T,M) = D_F - FLUX_T(IJK) 
                  A_M(IJKP,B,M) = D_F 
               ENDIF 
            ENDIF 
!
!
!           West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
               V_F = UF(IMJK) 
               D_F = AVG_X_H(DIF(IJKW),DIF(IJK),IM)*ODX_E(IM)*AYZ(IMJK) 
               IF (V_F >= ZERO) THEN 
                  A_M(IJK,W,M) = D_F + FLUX_E(IMJK) 
               ELSE 
                  A_M(IJK,W,M) = D_F 
               ENDIF 
            ENDIF 
!
!           South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               JM = JM1(J) 
               IJKS = SOUTH_OF(IJK) 
               V_F = VF(IJMK) 
               D_F = AVG_Y_H(DIF(IJKS),DIF(IJK),JM)*ODY_N(JM)*AXZ(IJMK) 
               IF (V_F >= ZERO) THEN 
                  A_M(IJK,S,M) = D_F + FLUX_N(IJMK) 
               ELSE 
                  A_M(IJK,S,M) = D_F 
               ENDIF 
            ENDIF 
!
!           Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  KM = KM1(K) 
                  IJKB = BOTTOM_OF(IJK) 
                  V_F = WF(IJKM) 
                  D_F = AVG_Z_H(DIF(IJKB),DIF(IJK),KM)*OX_E(I)*ODZ_T(KM)*AXY(&
                     IJKM) 
                  IF (V_F >= ZERO) THEN 
                     A_M(IJK,B,M) = D_F + FLUX_T(IJKM) 
                  ELSE 
                     A_M(IJK,B,M) = D_F 
                  ENDIF 
               ENDIF 
            ENDIF 
!
         ENDIF
      END DO 
!
      RETURN  
      END SUBROUTINE CONV_DIF_PHI0 

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:
!    CONV_DIF_Phi_DC(Phi, Dif, Disc, Uf, Vf, Wf, Flux_E,Flux_N,Flux_T, M, A_m, B_m, IER)
!  Purpose: TO USE DEFERRED CORRECTION IN SOLVING THE SCALAR TRANSPORT C
!  EQN. THIS METHOD COMBINES FIRST ORDER UPWIND AND A USER SPECIFIED   C
!  HIGH ORDER METHOD TO SOLVE FOR THE SCALAR PHI.
!  See source_Phi                                                  C
!                                                                      C
!  Author: C. GUENTHER                                Date: 1-ARP-99   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_DIF_PHI_DC(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE toleranc 
      USE run
      USE geometry
      USE compar
      USE sendrecv
      USE sendrecv3
      USE indices
      Use xsi_array
      Use tmp_array
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Scalar
      DOUBLE PRECISION Phi(DIMENSION_3)
!
!                      Gamma -- diffusion coefficient
      DOUBLE PRECISION Dif(DIMENSION_3)
!
!                      Discretizationindex
      INTEGER          Disc
!
!                      Velocity components
      DOUBLE PRECISION Uf(DIMENSION_3), Vf(DIMENSION_3), Wf(DIMENSION_3) 
!
!                      Mass flux components
      DOUBLE PRECISION Flux_E(DIMENSION_3), Flux_N(DIMENSION_3), Flux_T(DIMENSION_3) 
!
!                      Phase index
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          I,  J, K, IJK, IPJK, IJPK, IJKE, IJKN,&
                       IJKP, IJKT

      INTEGER          IMJK, IJKW
      INTEGER          IJMK, IJKS
      INTEGER          IJKM, IJKB
      INTEGER          IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER          IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4

! loezos
      INTEGER  incr
! loezos

!
!                      Difusion parameter
      DOUBLE PRECISION D_f
!
!	FACE VELOCITY
	DOUBLE PRECISION V_F
!
!	DEFERRED CORRCTION CONTRIBUTION FORM HIGH ORDER METHOD
	DOUBLE PRECISION PHI_HO
!
!	LOW ORDER APPROXIMATION 
	DOUBLE PRECISION PHI_LO
!
!	DEFERRED CORRECTION CONTRIBUTIONS FROM EACH FACE
	DOUBLE PRECISION  EAST_DC
	DOUBLE PRECISION  WEST_DC
	DOUBLE PRECISION  NORTH_DC
	DOUBLE PRECISION  SOUTH_DC
        DOUBLE PRECISION  TOP_DC
        DOUBLE PRECISION  BOTTOM_DC
!
!
!---------------------------------------------------------------
!	EXTERNAL FUNCTIONS
!---------------------------------------------------------------
	DOUBLE PRECISION , EXTERNAL :: FPFOI_OF
!---------------------------------------------------------------
!
!
!
!---------------------------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'function3.inc'
      INCLUDE 'fun_avg2.inc'


      call lock_xsi_array
      call lock_tmp4_array
!
!  Calculate convection factors
!
!
! Send recv the third ghost layer
      IF ( FPFOI ) THEN
         Do IJK = ijkstart3, ijkend3
	    I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = PHI(IJK)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

! loezos
	incr=0		
! loezos

      CALL CALC_XSI (DISC, PHI, UF, VF, WF, XSI_E, XSI_N, XSI_T,incr) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do                                              &
!$omp&     private(I,  J, K,  IJK,  IPJK, IJPK, IJKE, IJKN,         &
!$omp&             IJKP, IJKT,  V_f, D_f,                    &
!$omp&             IMJK, IJKW,                                  &
!$omp&             IJMK, IJKS,                                  &
!$omp&             IJKM, IJKB, PHI_HO, PHI_LO, CONV_FAC,       &
!$omp&             EAST_DC, WEST_DC, NORTH_DC, SOUTH_DC, TOP_DC, BOTTOM_DC)                     
!
      DO IJK = ijkstart3, ijkend3
!
! Determine whether IJK falls within 1 ghost layer........
       I = I_OF(IJK)
       J = J_OF(IJK)
       K = K_OF(IJK)
!
         IF (FLUID_AT(IJK)) THEN 
!
!
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
!
!           Third Ghost layer information
            IPPP  = IP_OF(IP_OF(IPJK))
            IPPP4 = funijk3(I_OF(IPPP), J_OF(IPPP), K_OF(IPPP))
            IMMM  = IM_OF(IM_OF(IMJK))
            IMMM4 = funijk3(I_OF(IMMM), J_OF(IMMM), K_OF(IMMM))
            JPPP  = JP_OF(JP_OF(IJPK))
            JPPP4 = funijk3(I_OF(JPPP), J_OF(JPPP), K_OF(JPPP))
            JMMM  = JM_OF(JM_OF(IJMK))
            JMMM4 = funijk3(I_OF(JMMM), J_OF(JMMM), K_OF(JMMM))
            KPPP  = KP_OF(KP_OF(IJKP))
            KPPP4 = funijk3(K_OF(IPPP), J_OF(KPPP), K_OF(KPPP))
            KMMM  = KM_OF(KM_OF(IJKM))
            KMMM4 = funijk3(I_OF(KMMM), J_OF(KMMM), K_OF(KMMM))
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1/2, j, k)
!
		V_F = UF(IJK)
		IF(V_F >= ZERO)THEN
		   PHI_LO = PHI(IJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IPJK), PHI(IJK), & 
                            PHI(IMJK), PHI(IM_OF(IMJK)))
		ELSE
		   PHI_LO = PHI(IPJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJK), PHI(IPJK), & 
                            PHI(IP_OF(IPJK)), TMP4(IPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
                      PHI_HO = XSI_E(IJK)*PHI(IPJK)+(1.0-XSI_E(IJK))*PHI(IJK)
		EAST_DC = FLUX_E(IJK)*(PHI_LO - PHI_HO)
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1/2, k)
!
		V_F = VF(IJK)
		IF(V_F >= ZERO)THEN
		   PHI_LO = PHI(IJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJPK), PHI(IJK), & 
                            PHI(IJMK), PHI(JM_OF(IJMK)))
		ELSE
		   PHI_LO = PHI(IJPK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJK), PHI(IJPK), & 
                            PHI(JP_OF(IJPK)), TMP4(JPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		     PHI_HO = XSI_N(IJK)*PHI(IJPK)+(1.0-XSI_N(IJK))*PHI(IJK)
		NORTH_DC = FLUX_N(IJK)*(PHI_LO - PHI_HO)
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j, k+1/2)
!
              IF (DO_K) THEN
                IJKP = KP_OF(IJK) 
                IJKT = TOP_OF(IJK)
	        V_F = WF(IJK)
		IF(V_F >= ZERO)THEN
                   PHI_LO = PHI(IJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJKP),  PHI(IJK), &
                            PHI(IJKM), PHI(KM_OF(IJKM)))
	        ELSE
		   PHI_LO = PHI(IJKP)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJK), PHI(IJKP),  &
                            PHI(KP_OF(IJKP)), TMP4(KPPP4))
	        ENDIF
                IF (.NOT. FPFOI ) &
                     PHI_HO = XSI_T(IJK)*PHI(IJKP)+(1.0-XSI_T(IJK))*PHI(IJK)
                TOP_DC = FLUX_T(IJK)*(PHI_LO - PHI_HO)
	      ELSE
	          TOP_DC = ZERO	    
              ENDIF
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE West face (i-1/2, j, k)
!
	    	IMJK = IM_OF(IJK)
	    	IJKW = WEST_OF(IJK)
	    	V_F = UF(IMJK)
	    	IF(V_F >= ZERO)THEN
		   PHI_LO = PHI(IMJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJK), PHI(IMJK), &
                            PHI(IM_OF(IMJK)), TMP4(IMMM4))
		ELSE
		   PHI_LO = PHI(IJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IMJK), PHI(IJK), &
                            PHI(IPJK), PHI(IP_OF(IPJK)))
                ENDIF
                IF (.NOT. FPFOI ) &
		      PHI_HO = XSI_E(IMJK)*PHI(IJK)+(ONE-XSI_E(IMJK))*PHI(IMJK)
		WEST_DC = FLUX_E(IMJK)*(PHI_LO - PHI_HO)
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i, j-1/2, k)
!
            	IJMK = JM_OF(IJK) 
            	IJKS = SOUTH_OF(IJK)
		V_F = VF(IJMK)
		IF(V_F >= ZERO)THEN
		   PHI_LO = PHI(IJMK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJK), PHI(IJMK), & 
                            PHI(JM_OF(IJMK)), TMP4(JMMM4))
		ELSE
		   PHI_LO = PHI(IJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJMK), PHI(IJK), & 
                            PHI(IJPK), PHI(JP_OF(IJPK)))
                ENDIF
                IF (.NOT. FPFOI ) &
            	      PHI_HO = XSI_N(IJMK)*PHI(IJK)+(ONE-XSI_N(IJMK))*PHI(IJMK)
		SOUTH_DC = FLUX_N(IJMK)*(PHI_LO - PHI_HO)
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i, j, k-1/2)
              IF (DO_K) THEN 
                 IJKM = KM_OF(IJK) 
                 IJKB = BOTTOM_OF(IJK)
		 V_F = WF(IJKM)
		 IF(V_F >= ZERO)THEN
                   PHI_LO = PHI(IJKM)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJK), PHI(IJKM), &
                            PHI(KM_OF(IJKM)), TMP4(KMMM4))
		 ELSE
                   PHI_LO = PHI(IJK)
                   IF ( FPFOI ) &
                      PHI_HO = FPFOI_OF(PHI(IJKM), PHI(IJK), &
                            PHI(IJKP), PHI(KP_OF(IJKP)))
                 ENDIF
                 IF (.NOT. FPFOI ) &
                      PHI_HO = XSI_T(IJKM)*PHI(IJK)+(1.0-XSI_T(IJKM))*PHI(IJKM)
		 BOTTOM_DC = FLUX_T(IJKM)*(PHI_LO - PHI_HO)
	      ELSE
		   BOTTOM_DC = ZERO
	      ENDIF

!
!	    CONTRIBUTION DUE TO DEFERRED CORRECTION
!
		B_M(IJK,M) = B_M(IJK,M)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC&
				+BOTTOM_DC-TOP_DC
!
         ENDIF 
      END DO 
      call unlock_tmp4_array
      call unlock_xsi_array
!
!
      RETURN  
      END SUBROUTINE CONV_DIF_PHI_DC 
!
!

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:
!    CONV_DIF_Phi1(Phi, Dif, Disc, Uf, Vf, Wf, Flux_E,Flux_N,Flux_T, M, A_m, B_m, IER)
!  Purpose: Determine convection diffusion terms for gas energy eq Phi C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_Phi                                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_DIF_PHI1(PHI,DIF,DISC,UF,VF,WF,Flux_E,Flux_N,Flux_T,M,A_M,B_M,IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE toleranc 
      USE run
      USE geometry
      USE compar
      USE sendrecv
      USE indices
      USE vshear
      Use xsi_array

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Scalar
      DOUBLE PRECISION Phi(DIMENSION_3)
!
!                      Gamma -- diffusion coefficient
      DOUBLE PRECISION Dif(DIMENSION_3)
!
!                      Discretizationindex
      INTEGER          Disc
!
!                      Velocity components
      DOUBLE PRECISION Uf(DIMENSION_3), Vf(DIMENSION_3), Wf(DIMENSION_3) 
!
!                      Mass flux components
      DOUBLE PRECISION Flux_E(DIMENSION_3), Flux_N(DIMENSION_3), Flux_T(DIMENSION_3) 
!
!                      Phase index
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          I,  J, K, IJK, IPJK, IJPK, IJKE, IJKN,&
                       IJKP, IJKT

      INTEGER          IMJK, IM, IJKW
      INTEGER          IJMK, JM, IJKS
      INTEGER          IJKM, KM, IJKB
! start loezos
      INTEGER          I1, J1
      INTEGER incr
! end loezos

!
!                      Difusion parameter
      DOUBLE PRECISION D_f
!
!                      Convection weighting factors
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),&
!                       XSI_t(DIMENSION_3)
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      call lock_xsi_array
!
!  Calculate convection factors
!
!

! loezos
	incr=0		
! loezos	
 
      CALL CALC_XSI (DISC, PHI, UF, VF, WF, XSI_E, XSI_N, XSI_T,incr) 

! loezos
!update V to true velocity      

      IF (SHEAR) THEN
	 DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
	   VF(IJK)=VF(IJK)+VSH(IJK)	
          END IF
        END DO
      END IF
! loezos

!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!$omp      parallel do                                               &
!$omp&     private(I,  J, K,  IJK,  IPJK, IJPK, IJKE, IJKN,          &
!$omp&             IJKP, IJKT,    D_f,                        &
!$omp&             IMJK, IM, IJKW,                                   &
!$omp&             IJMK, JM,  IJKS,                                  &
!$omp&             IJKM, KM,  IJKB )                      
!
!
!
      DO IJK = ijkstart3, ijkend3
!
       I = I_OF(IJK) 
       J = J_OF(IJK) 
       K = K_OF(IJK) 
!
         IF (FLUID_AT(IJK)) THEN 
!
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
!
!
!           East face (i+1/2, j, k)
            D_F = AVG_X_H(DIF(IJK),DIF(IJKE),I)*ODX_E(I)*AYZ(IJK) 
!
            A_M(IJK,E,M) = D_F - XSI_E(IJK)*FLUX_E(IJK) 
!
            A_M(IPJK,W,M) = D_F + (ONE - XSI_E(IJK))*FLUX_E(IJK) 
!
!
!           North face (i, j+1/2, k)
            D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*AXZ(IJK) 
!
            A_M(IJK,N,M) = D_F - XSI_N(IJK)*FLUX_N(IJK) 
!
            A_M(IJPK,S,M) = D_F + (ONE - XSI_N(IJK))*FLUX_N(IJK) 
!
!
!           Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
!
               D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*OX(I)*ODZ_T(K)*AXY(IJK) 
!
               A_M(IJK,T,M) = D_F - XSI_T(IJK)*FLUX_T(IJK) 
!
               A_M(IJKP,B,M)=D_F+(ONE-XSI_T(IJK))*FLUX_T(IJK) 
            ENDIF 
!
!           West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
!
               D_F = AVG_X_H(DIF(IJKW),DIF(IJK),IM)*ODX_E(IM)*AYZ(IMJK) 
!
               A_M(IJK,W,M) = D_F + (ONE - XSI_E(IMJK))*FLUX_E(IMJK) 
            ENDIF 
!
!           South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               JM = JM1(J) 
               IJKS = SOUTH_OF(IJK) 
!
               D_F = AVG_Y_H(DIF(IJKS),DIF(IJK),JM)*ODY_N(JM)*AXZ(IJMK) 
!
               A_M(IJK,S,M) = D_F + (ONE - XSI_N(IJMK))*FLUX_N(IJMK) 
            ENDIF 
!
!           Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  KM = KM1(K) 
                  IJKB = BOTTOM_OF(IJK) 
!
                  D_F = AVG_Z_H(DIF(IJKB),DIF(IJK),KM)*OX_E(I)*ODZ_T(KM)*AXY(&
                     IJKM) 
!
                  A_M(IJK,B,M) = D_F + (ONE - XSI_T(IJKM))*FLUX_T(IJKM) 
               ENDIF 
            ENDIF 
!
         ENDIF 
      END DO 

! loezos 
       IF (SHEAR) THEN
	 DO IJK = ijkstart3, ijkend3
          IF (FLUID_AT(IJK)) THEN  	 
	   VF(IJK)=VF(IJK)-VSH(IJK)	
	  END IF
         END DO 	
        END IF
! loezos      
      call unlock_xsi_array

!
      RETURN  
      END SUBROUTINE CONV_DIF_PHI1 
!
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DIF_phi_IS(Dif, A_m, B_m, M, IER)                      C
!  Purpose: Remove diffusive fluxes across internal surfaces.          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE DIF_PHI_IS(DIF, A_M, B_M, M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE toleranc 
      USE run
      USE geometry
      USE compar
      USE sendrecv
      USE indices
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE output
      USE is
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Internal surface
      INTEGER          L
!
!                      Indices
      INTEGER          I,  J, K, I1, I2, J1, J2, K1, K2, IJK,&
                       IJKE, IJKN, IJKT, IPJK, IJPK, IJKP
!
!                      Solids phase
      INTEGER          M
!
!                      Gamma -- diffusion coefficient
      DOUBLE PRECISION Dif(DIMENSION_3)
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Difusion parameter
      DOUBLE PRECISION D_f
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
! Make user defined internal surfaces non-conducting
!
      DO L = 1, DIMENSION_IS 
         IF (IS_DEFINED(L)) THEN 
            I1 = IS_I_W(L) 
            I2 = IS_I_E(L) 
            J1 = IS_J_S(L) 
            J2 = IS_J_N(L) 
            K1 = IS_K_B(L) 
            K2 = IS_K_T(L) 

!       Limit I1, I2 and all to local processor first ghost layer

	    IF(I1.LE.IEND2)   I1 = MAX(I1, ISTART2)

            IF(J1.LE.JEND2)   J1 = MAX(J1, JSTART2)

            IF(K1.LE.KEND2)   K1 = MAX(K1, KSTART2)

            IF(I2.GE.ISTART2) I2 = MIN(I2, IEND2)

            IF(J2.GE.JSTART2) J2 = MIN(J2, JEND2)

            IF(K2.GE.KSTART2) K2 = MIN(K2, KEND2)

!     End of limiting to the first ghost cells of the processor....
            DO K = K1, K2 
               DO J = J1, J2 
                  DO I = I1, I2 
                     IJK = FUNIJK(I,J,K) 
!
                     SELECT CASE (TRIM(IS_PLANE(L)))  
                     CASE ('E')  
                        IJKE = EAST_OF(IJK) 
                        IPJK = IP_OF(IJK) 
!
                        D_F = AVG_X_H(DIF(IJK),DIF(IJKE),I)*ODX_E(I)*AYZ(IJK) 
!
                        A_M(IJK,E,M) = A_M(IJK,E,M) - D_F 
                        A_M(IPJK,W,M) = A_M(IPJK,W,M) - D_F 
!
                     CASE ('N')  
                        IJKN = NORTH_OF(IJK) 
                        IJPK = JP_OF(IJK) 
!
                        D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*AXZ(IJK) 
!
                        A_M(IJK,N,M) = A_M(IJK,N,M) - D_F 
                        A_M(IJPK,S,M) = A_M(IJPK,S,M) - D_F 
!
                     CASE ('T')  
                        IF (DO_K) THEN 
                           IJKT = TOP_OF(IJK) 
                           IJKP = KP_OF(IJK) 
!
                           D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*OX(I)*ODZ_T(K)*&
                              AXY(IJK) 
!
                           A_M(IJK,T,M) = A_M(IJK,T,M) - D_F 
                           A_M(IJKP,B,M) = A_M(IJKP,B,M) - D_F 
!
                        ENDIF 
                     CASE DEFAULT 
!
                     END SELECT 
                  END DO 
               END DO 
            END DO 
         ENDIF 
      END DO 
!
      RETURN  
      END SUBROUTINE DIF_PHI_IS 
                                                                                                                                                                                                                                                      conv_dif_u_g.f                                                                                      0100644 0002444 0000146 00000075521 10247670146 012032  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_U_g(A_m, B_m, IER)                            C
!  Purpose: Determine convection diffusion terms for U_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_u_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_DIF_U_G(A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE compar 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Error index
      INTEGER          IER
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
 
!-----------------------------------------------
!
!
!
      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN

            
!	IF DEFERRED CORRECTION IS USED TO SOLVE U_G
!
      IF(DEF_COR)THEN
	CALL STORE_A_U_G0(A_M(1,-3,0), IER)
        IF (DISCRETIZE(3) > 1) CALL STORE_A_U_GDC (A_M(1,-3,0),B_M(1,0), IER)
      ELSE
!
!	NO DEFERRED CORRECTION IS TO BE USED TO SOLVE FOR U_G   
!
        IF (DISCRETIZE(3) == 0) THEN               ! 0 & 1 => FOUP 
          CALL STORE_A_U_G0 (A_M(1,-3,0), IER) 
        ELSE 
          CALL STORE_A_U_G1 (A_M(1,-3,0), IER) 
        ENDIF 
!
      ENDIF
            
      CALL DIF_U_IS (MU_GT, A_M, B_M, 0, IER) 
            
      RETURN  
      END SUBROUTINE CONV_DIF_U_G 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_U_g0(A_U_g, IER)                               C
!  Purpose: Determine convection diffusion terms for U_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_u_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_U_G0(A_U_G, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE toleranc 
      USE physprop
      USE fldvar
      USE output
      USE compar
      USE mflux 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,&
                      IJKNE, IJKP, IJKT, IJKTE
 
      INTEGER          IMJK, IM, IJKW
      INTEGER          IJMK, JM, IPJMK, IJKS, IJKSE
      INTEGER          IJKM, KM, IPJKM, IJKB, IJKBE
!
!                      Solids phase
      INTEGER          M
!
!                      Face mass flux
      DOUBLE PRECISION Flux
!
!                      Diffusion parameter
      DOUBLE PRECISION D_f
!
!                      Septadiagonal matrix A_U_g
      DOUBLE PRECISION A_U_g(DIMENSION_3, -3:3)
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'

            

!
!  Calculate convection-diffusion fluxes through each of the faces
!
!     Fluid phase
      M = 0 
      
!$omp      parallel do                                                  &
!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,     &
!$omp&                    IJKNE, IJKP, IJKT, IJKTE, V_f, D_f,   &
!$omp&                    IMJK, IM, IJKW,                               &
!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,                 &
!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE)                 
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_E(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 	      
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKE 
            ELSE 
               IJKC = IJK 
            ENDIF 
            IP = IP1(I) 
            IJKN = NORTH_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!
!           East face (i+1, j, k)
	    Flux = HALF * (Flux_gE(IJK) + Flux_gE(IPJK))
            D_F = MU_GT(IJKE)*ODX(IP)*AYZ_U(IJK) 
            IF (Flux >= ZERO) THEN 
               A_U_G(IJK,E) = D_F 
               A_U_G(IPJK,W) = D_F + Flux
            ELSE 
               A_U_G(IJK,E) = D_F - Flux
               A_U_G(IPJK,W) = D_F 
            ENDIF 
!
!
!           North face (i+1/2, j+1/2, k)
	    Flux = HALF * (Flux_gN(IJK) + Flux_gN(IPJK))
            D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKE&
               ),MU_GT(IJKNE),J),I)*ODY_N(J)*AXZ_U(IJK) 
            IF (Flux >= ZERO) THEN 
               A_U_G(IJK,N) = D_F 
               A_U_G(IJPK,S) = D_F + Flux
            ELSE 
               A_U_G(IJK,N) = D_F - Flux
               A_U_G(IJPK,S) = D_F 
            ENDIF 
!
!           Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 
	       Flux = HALF * (Flux_gT(IJK) + Flux_gT(IPJK))
               D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                  IJKE),MU_GT(IJKTE),K),I)*OX_E(I)*ODZ_T(K)*AXY_U(IJK) 
               IF (Flux >= ZERO) THEN 
                  A_U_G(IJK,T) = D_F 
                  A_U_G(IJKP,B) = D_F + Flux
               ELSE 
                  A_U_G(IJK,T) = D_F - Flux
                  A_U_G(IJKP,B) = D_F 
               ENDIF 
            ENDIF 
!
!
!           West face (i, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
	       Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IJK))
               D_F = MU_GT(IJKC)*ODX(I)*AYZ_U(IMJK) 
               IF (Flux >= ZERO) THEN 
                  A_U_G(IJK,W) = D_F + Flux
               ELSE 
                  A_U_G(IJK,W) = D_F 
               ENDIF 
            ENDIF 
!
!           South face (i+1/2, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IJMK)) THEN 
               JM = JM1(J) 
               IPJMK = IP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKSE = EAST_OF(IJKS) 
	       Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IPJMK))
               D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                  IJKSE),MU_GT(IJKE),JM),I)*ODY_N(JM)*AXZ_U(IJMK) 
               IF (Flux >= ZERO) THEN 
                  A_U_G(IJK,S) = D_F + Flux
               ELSE 
                  A_U_G(IJK,S) = D_F 
               ENDIF 
            ENDIF 
!
!           Bottom face (i+1/2, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLOW_AT_E(IJKM)) THEN 
                  KM = KM1(K) 
                  IPJKM = IP_OF(IJKM) 
                  IJKB = BOTTOM_OF(IJK) 
                  IJKBE = EAST_OF(IJKB) 
	          Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IPJKM))
                  D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                     MU_GT(IJKBE),MU_GT(IJKE),KM),I)*OX_E(I)*ODZ_T(KM)*AXY_U(&
                     IJKM) 
                  IF (Flux >= ZERO) THEN 
                     A_U_G(IJK,B) = D_F + Flux
                  ELSE 
                     A_U_G(IJK,B) = D_F 
                  ENDIF 
               ENDIF 
            ENDIF 
!
         ENDIF 
      END DO 
      
    
      RETURN  
      END SUBROUTINE STORE_A_U_G0 

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_U_GDC(A_U_g, B_M, IER)                         C
!  Purpose: TO USE DEFERRED CORRECTION METHOD TO SOLVE THE U-MOMENTUM  C
!  EQUATION. THIS METHOD COMBINES FIRST ORDER UPWIND AND A USER        C
!  SPECIFIED HIGH ORDER METHOD                                         C
!                                                                      C
!  Author: C. GUENTHER                                Date: 8-APR-99   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_U_GDC(A_U_G, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE toleranc 
      USE physprop
      USE fldvar
      USE output
      USE vshear
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar   
      USE sendrecv
      USE sendrecv3
      USE mflux 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

!
!
!                      Error index
      INTEGER          IER
!
!                      Indices
!
!                      Indices
      INTEGER          I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,&
                       IJKNE, IJKP, IJKT, IJKTE
      INTEGER          IMJK, IM, IJKW
      INTEGER          IJMK, JM, IPJMK, IJKS, IJKSE
      INTEGER          IJKM, KM, IPJKM, IJKB, IJKBE
      INTEGER          IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER          IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
!
! loezos
	INTEGER incr
! loezos

!                      Diffusion parameter
      DOUBLE PRECISION D_f

!                      Septadiagonal matrix A_U_g
      DOUBLE PRECISION A_U_g(DIMENSION_3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3)
!
!	DEFERRED CORRCTION CONTRIBUTION FORM HIGH ORDER METHOD
	DOUBLE PRECISION MOM_HO
!
!	LOW ORDER APPROXIMATION 
	DOUBLE PRECISION MOM_LO
!
!	CONVECTION FACTOR AT THE FACE
	DOUBLE PRECISION Flux
!
!	DEFERRED CORRECTION CONTRIBUTIONS FROM EACH FACE
	DOUBLE PRECISION 	EAST_DC
	DOUBLE PRECISION 	WEST_DC
	DOUBLE PRECISION 	NORTH_DC
	DOUBLE PRECISION 	SOUTH_DC
      DOUBLE PRECISION  TOP_DC
      DOUBLE PRECISION  BOTTOM_DC
!
!
!
!-----------------------------------------------
!
!---------------------------------------------------------------
!	EXTERNAL FUNCTIONS
!---------------------------------------------------------------
	DOUBLE PRECISION , EXTERNAL :: FPFOI_OF
!---------------------------------------------------------------
!
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'function3.inc'


      call lock_tmp4_array
      call lock_tmp_array
      call lock_xsi_array


!
!  Calculate convection factors
!
!
! Send recv the third ghost layer
      IF ( FPFOI ) THEN
         Do IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = U_G(IJK)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

!$omp parallel do private(IJK,I,IP,IPJK,IJKE)
       DO IJK = ijkstart3, ijkend3
!
         I = I_OF(IJK) 
         IP = IP1(I) 
         IPJK = IP_OF(IJK) 
         IJKE = EAST_OF(IJK) 
!
!
!           East face (i+1, j, k)
         U(IJK) = AVG_X_E(U_G(IJK),U_G(IPJK),IP) 
!
!
!           North face (i+1/2, j+1/2, k)
         V(IJK) = AVG_X(V_G(IJK),V_G(IPJK),I) 
!
!
!           Top face (i+1/2, j, k+1/2)
         IF (DO_K) WW(IJK) = AVG_X(W_G(IJK),W_G(IPJK),I) 
      END DO 

! loezos
	incr=1		
! loezos

      CALL CALC_XSI (DISCRETIZE(3), U_G, U, V, WW, XSI_E, XSI_N, XSI_T,incr) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!$omp      parallel do                                                 &
!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,    &
!$omp&                    IJKNE, IJKP, IJKT, IJKTE,   D_f,  &
!$omp&                    IMJK, IM, IJKW,                              &
!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,                &
!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE, &
!$omp&              MOM_HO, MOM_LO, CONV_FAC,EAST_DC,WEST_DC,NORTH_DC,&
!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_E(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK) 
            IJKE = EAST_OF(IJK) 
!
!           Third Ghost layer information
            IPPP  = IP_OF(IP_OF(IPJK))
            IPPP4 = funijk3(I_OF(IPPP), J_OF(IPPP), K_OF(IPPP))
            IMMM  = IM_OF(IM_OF(IMJK))
            IMMM4 = funijk3(I_OF(IMMM), J_OF(IMMM), K_OF(IMMM))
            JPPP  = JP_OF(JP_OF(IJPK))
            JPPP4 = funijk3(I_OF(JPPP), J_OF(JPPP), K_OF(JPPP))
            JMMM  = JM_OF(JM_OF(IJMK))
            JMMM4 = funijk3(I_OF(JMMM), J_OF(JMMM), K_OF(JMMM))
            KPPP  = KP_OF(KP_OF(IJKP))
            KPPP4 = funijk3(K_OF(IPPP), J_OF(KPPP), K_OF(KPPP))
            KMMM  = KM_OF(KM_OF(IJKM))
            KMMM4 = funijk3(I_OF(KMMM), J_OF(KMMM), K_OF(KMMM))
!
!
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKE 
            ELSE 
               IJKC = IJK 
            ENDIF 
            IP = IP1(I) 
            IJKN = NORTH_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1, j, k)
!
		IF(U(IJK) >= ZERO)THEN
		    MOM_LO = U_G(IJK)
                     IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IPJK), U_G(IJK), & 
                            U_G(IMJK), U_G(IM_OF(IMJK)))
		ELSE
		    MOM_LO = U_G(IPJK)
                     IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJK), U_G(IPJK), & 
                            U_G(IP_OF(IPJK)), TMP4(IPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_E(IJK)*U_G(IPJK)+ &
                               (1.0-XSI_E(IJK))*U_G(IJK)
			      
	        Flux = HALF * (Flux_gE(IJK) + Flux_gE(IPJK))
		EAST_DC = Flux *(MOM_LO - MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i+1/2, j+1/2, k)
!            
		IF(V(IJK) >= ZERO)THEN
		    MOM_LO = U_G(IJK)
                     IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJPK), U_G(IJK), & 
                            U_G(IJMK), U_G(JM_OF(IJMK)))
		ELSE
		    MOM_LO = U_G(IJPK)
                     IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJK), U_G(IJPK), & 
                            U_G(JP_OF(IJPK)), TMP4(JPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_N(IJK)*U_G(IJPK)+ &
                                 (1.0-XSI_N(IJK))*U_G(IJK)
	        Flux = HALF * (Flux_gN(IJK) + Flux_gN(IPJK))
		NORTH_DC = Flux *(MOM_LO - MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i+1/2, j, k+1/2)
!
            IF (DO_K) THEN 
                IJKP = KP_OF(IJK) 
                IJKT = TOP_OF(IJK) 
                IJKTE = EAST_OF(IJKT) 
		IF(WW(IJK) >= ZERO)THEN
		   MOM_LO = U_G(IJK)
                   IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJKP), U_G(IJK), & 
                            U_G(IJKM), U_G(KM_OF(IJKM)))
		ELSE
		   MOM_LO = U_G(IJKP)
                   IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJK), U_G(IJKP), & 
                            U_G(KP_OF(IJKP)), TMP4(KPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_T(IJK)*U_G(IJKP)+ &
                               (1.0-XSI_T(IJK))*U_G(IJK)
	        Flux = HALF * (Flux_gT(IJK) + Flux_gT(IPJK))
		TOP_DC = Flux *(MOM_LO - MOM_HO)
		
	    ELSE
		TOP_DC = ZERO
		
            ENDIF
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE West face (i, j, k)
!
            IMJK = IM_OF(IJK) 
            IM = IM1(I) 
            IJKW = WEST_OF(IJK)
	    IF(U(IMJK) >= ZERO)THEN
	      MOM_LO = U_G(IMJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJK), U_G(IMJK), & 
                            U_G(IM_OF(IMJK)), TMP4(IMMM4))
 	    ELSE
	      MOM_LO = U_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IMJK), U_G(IJK), & 
                            U_G(IPJK), U_G(IP_OF(IPJK)))
	    ENDIF
            IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_E(IMJK)*U_G(IJK)+ &
                               (1.0-XSI_E(IMJK))*U_G(IMJK)
	    Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IJK))
	    WEST_DC = Flux * (MOM_LO - MOM_HO)

!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i+1/2, j-1/2, k)
!
            IJMK = JM_OF(IJK) 
            JM = JM1(J) 
            IPJMK = IP_OF(IJMK) 
            IJKS = SOUTH_OF(IJK) 
            IJKSE = EAST_OF(IJKS) 
 	    IF(V(IJMK) >= ZERO)THEN
	      MOM_LO = U_G(IJMK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJK), U_G(IJMK), & 
                            U_G(JM_OF(IJMK)), TMP4(JMMM4))
	    ELSE
	      MOM_LO = U_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJMK), U_G(IJK), & 
                            U_G(IJPK), U_G(JP_OF(IJPK)))
	    ENDIF
            IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_N(IJMK)*U_G(IJK)+ &
                              (1.0-XSI_N(IJMK))*U_G(IJMK)
	    Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IPJMK))
	    SOUTH_DC = Flux *(MOM_LO - MOM_HO)
	    
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i+1/2, j, k-1/2)
!
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               KM = KM1(K) 
               IPJKM = IP_OF(IJKM) 
               IJKB = BOTTOM_OF(IJK) 
               IJKBE = EAST_OF(IJKB)
	       IF(WW(IJK) >= ZERO)THEN
		 MOM_LO = U_G(IJKM)
                 IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJK), U_G(IJKM), & 
                            U_G(KM_OF(IJKM)), TMP4(KMMM4))
	       ELSE
		 MOM_LO = U_G(IJK)
                 IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_G(IJKM), U_G(IJK), & 
                            U_G(IJKP), U_G(KP_OF(IJKP)))
	       ENDIF
               IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_T(IJKM)*U_G(IJK)+ &
                               (1.0-XSI_T(IJKM))*U_G(IJKM)
	       Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IPJKM))
	       BOTTOM_DC = Flux * (MOM_LO - MOM_HO)
            ELSE
	       BOTTOM_DC = ZERO
	      
            ENDIF
!
!	    CONTRIBUTION DUE TO DEFERRED CORRECTION
!
	    B_M(IJK) = B_M(IJK)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC&
				+BOTTOM_DC-TOP_DC
!

         ENDIF 
      END DO 

      call unlock_tmp4_array
      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_U_GDC 


!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_U_g1(A_U_g, IER)                               C
!  Purpose: Determine convection diffusion terms for U_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; Higher order C
!  See source_u_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_U_G1(A_U_G, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE toleranc 
      USE physprop
      USE fldvar
      USE output
      USE vshear
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar  
      USE mflux 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

!
!
!                      Error index
      INTEGER          IER
!
!                      Indices
!
!                      Indices
      INTEGER          I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,&
                       IJKNE, IJKP, IJKT, IJKTE
      INTEGER          IMJK, IM, IJKW
      INTEGER          IJMK, JM, IPJMK, IJKS, IJKSE
      INTEGER          IJKM, KM, IPJKM, IJKB, IJKBE

! loezos                     SHEAR VELOCITY
      INTEGER incr    
!loezos

!
!                      Diffusion parameter
      DOUBLE PRECISION D_f
!
!                      Face mass flux
      DOUBLE PRECISION Flux
!
!                      Septadiagonal matrix A_U_g
      DOUBLE PRECISION A_U_g(DIMENSION_3, -3:3)
!
!                      Convection weighting factors
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),&
!                      XSI_t(DIMENSION_3)
!      DOUBLE PRECISION U(DIMENSION_3),&
!                      V(DIMENSION_3), WW(DIMENSION_3)
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'


      call lock_tmp_array
      call lock_xsi_array


!
!  Calculate convection factors
!
!$omp parallel do private(IJK,I,IP,IPJK,IJKE)
      DO IJK = ijkstart3, ijkend3
!
         I = I_OF(IJK) 
	 J=J_OF(IJK)
         IP = IP1(I) 
         IPJK = IP_OF(IJK) 
         IJKE = EAST_OF(IJK) 


!
!
!           East face (i+1, j, k)
         U(IJK) = AVG_X_E(U_G(IJK),U_G(IPJK),IP) 
!
!
!           North face (i+1/2, j+1/2, k)
         V(IJK) = AVG_X(V_G(IJK),V_G(IPJK),I) 
!
!
!           Top face (i+1/2, j, k+1/2)
         IF (DO_K) WW(IJK) = AVG_X(W_G(IJK),W_G(IPJK),I) 
      END DO 

! loezos
	incr=1		
! loezos

      CALL CALC_XSI (DISCRETIZE(3), U_G, U, V, WW, XSI_E, XSI_N, XSI_T,incr) 

! loezos      
!update V to true velocity
      IF (SHEAR) THEN
!$omp parallel do private(IJK)
	 DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
	   V(IJK)=V(IJK)+VSHE(IJK)	
          END IF
        END DO
      END IF
! loezos


!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!$omp      parallel do                                                 &
!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,    &
!$omp&                    IJKNE, IJKP, IJKT, IJKTE,   D_f,  &
!$omp&                    IMJK, IM, IJKW,                              &
!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,                &
!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_E(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKE 
            ELSE 
               IJKC = IJK 
            ENDIF 
            IP = IP1(I) 
            IJKN = NORTH_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           East face (i+1, j, k)
	    Flux = HALF * (Flux_gE(IJK) + Flux_gE(IPJK))
            D_F = MU_GT(IJKE)*ODX(IP)*AYZ_U(IJK) 
!
            A_U_G(IJK,E) = D_F - XSI_E(IJK) * Flux
!
            A_U_G(IPJK,W) = D_F + (ONE - XSI_E(IJK)) * Flux
!
!
!           North face (i+1/2, j+1/2, k)
	    Flux = HALF * (Flux_gN(IJK) + Flux_gN(IPJK))
            D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKE&
               ),MU_GT(IJKNE),J),I)*ODY_N(J)*AXZ_U(IJK) 
!
            A_U_G(IJK,N) = D_F - XSI_N(IJK) * Flux 
!
            A_U_G(IJPK,S) = D_F + (ONE - XSI_N(IJK)) * Flux 
!
!
!           Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 
!
	       Flux = HALF * (Flux_gT(IJK) + Flux_gT(IPJK))
               D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                  IJKE),MU_GT(IJKTE),K),I)*OX_E(I)*ODZ_T(K)*AXY_U(IJK) 
!
               A_U_G(IJK,T) = D_F - XSI_T(IJK) * Flux
!
               A_U_G(IJKP,B) = D_F + (ONE - XSI_T(IJK)) * Flux 
            ENDIF 
!
!           West face (i, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
!
	       Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IJK))
               D_F = MU_GT(IJKC)*ODX(I)*AYZ_U(IMJK) 
!
               A_U_G(IJK,W) = D_F + (ONE - XSI_E(IMJK)) * Flux
            ENDIF 
!
!           South face (i+1/2, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IJMK)) THEN 
               JM = JM1(J) 
               IPJMK = IP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKSE = EAST_OF(IJKS) 
!
	       Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IPJMK))
               D_F = AVG_X_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                  IJKSE),MU_GT(IJKE),JM),I)*ODY_N(JM)*AXZ_U(IJMK) 
!
               A_U_G(IJK,S) = D_F + (ONE - XSI_N(IJMK)) * Flux
            ENDIF 
!
!           Bottom face (i+1/2, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLOW_AT_E(IJKM)) THEN 
                  KM = KM1(K) 
                  IPJKM = IP_OF(IJKM) 
                  IJKB = BOTTOM_OF(IJK) 
                  IJKBE = EAST_OF(IJKB) 
!
	          Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IPJKM))
                  D_F = AVG_X_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                     MU_GT(IJKBE),MU_GT(IJKE),KM),I)*OX_E(I)*ODZ_T(KM)*AXY_U(&
                     IJKM) 
!
                  A_U_G(IJK,B) = D_F + (ONE - XSI_T(IJKM)) * Flux
               ENDIF 
            ENDIF 
!
         ENDIF 
      END DO 

      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_U_G1 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
                                                                                                                                                                               conv_dif_u_s.f                                                                                      0100644 0002444 0000146 00000076667 10247670663 012066  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_U_s(A_m, B_m, IER)                            C
!  Purpose: Determine convection diffusion terms for U_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_u_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_DIF_U_S(A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Solids phase index 
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
 
!-----------------------------------------------
!
!
      DO M = 1, MMAX 
        IF (MOMENTUM_X_EQ(M)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!	  IF DEFERRED CORRECTION IS USED TO SOLVE U_S
!
	  IF(DEF_COR)THEN
	    CALL STORE_A_U_S0 (A_M(1,-3,M), M, IER)
	    IF (DISCRETIZE(3) > 1) CALL STORE_A_U_SDC (A_M(1,-3,M), M, B_M, IER)
	  ELSE
!
!	  NO DEFERRED CORRECTION IS TO BE USED TO SOLVE FOR U_S
!  
            IF (DISCRETIZE(3) == 0) THEN         ! 0 & 1 => FOUP 
               CALL STORE_A_U_S0 (A_M(1,-3,M), M, IER) 
            ELSE 
               CALL STORE_A_U_S1 (A_M(1,-3,M), M, IER) 
            ENDIF 
!
          ENDIF
	  	  
	  CALL DIF_U_IS (MU_S(1,M), A_M, B_M, M, IER)
	ENDIF 
      END DO 
      
      RETURN  
      END SUBROUTINE CONV_DIF_U_S 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_U_s0(A_U_s, M, IER)                            C
!  Purpose: Determine convection diffusion terms for U_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; FOUP         C
!  See source_u_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_U_S0(A_U_S, M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc 
      USE fldvar
      USE output
      USE compar 
      USE mflux 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
! 
!                      Indices 
      INTEGER          I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,& 
                       IJKNE, IJKP, IJKT, IJKTE 
      INTEGER          IMJK, IM, IJKW 
      INTEGER          IJMK, JM, IPJMK, IJKS, IJKSE 
      INTEGER          IJKM, KM, IPJKM, IJKB, IJKBE 
! 
!                      Solids phase 
      INTEGER          M 
!
!                      Face mass flux
      DOUBLE PRECISION Flux
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_U_s 
      DOUBLE PRECISION A_U_s(DIMENSION_3, -3:3, M:M) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do 	&
!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,	&
!$omp&                    IJKNE, IJKP, IJKT, IJKTE,  V_f, D_f,	&
!$omp&                    IMJK, IM, IJKW,	&
!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,	&
!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_E(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKE 
            ELSE 
               IJKC = IJK 
            ENDIF 
            IP = IP1(I) 
            IJKN = NORTH_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           East face (i+1, j, k)
	    Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IPJK,M))
            D_F = MU_S(IJKE,M)*ODX(IP)*AYZ_U(IJK) 
            IF (Flux >= ZERO) THEN 
               A_U_S(IJK,E,M) = D_F 
               A_U_S(IPJK,W,M) = D_F + Flux
            ELSE 
               A_U_S(IJK,E,M) = D_F - Flux 
               A_U_S(IPJK,W,M) = D_F 
            ENDIF 
!
!           North face (i+1/2, j+1/2, k)
	    Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IPJK,M))
            D_F = AVG_X_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(&
               IJKE,M),MU_S(IJKNE,M),J),I)*ODY_N(J)*AXZ_U(IJK) 
            IF (Flux >= ZERO) THEN 
               A_U_S(IJK,N,M) = D_F 
               A_U_S(IJPK,S,M) = D_F + Flux
            ELSE 
               A_U_S(IJK,N,M) = D_F - Flux
               A_U_S(IJPK,S,M) = D_F 
            ENDIF 
!
!           Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 
	       Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IPJK,M))
               D_F = AVG_X_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),AVG_Z_H(MU_S(&
                  IJKE,M),MU_S(IJKTE,M),K),I)*OX_E(I)*ODZ_T(K)*AXY_U(IJK) 
               IF (Flux >= ZERO) THEN 
                  A_U_S(IJK,T,M) = D_F 
                  A_U_S(IJKP,B,M) = D_F + Flux
               ELSE 
                  A_U_S(IJK,T,M) = D_F - Flux
                  A_U_S(IJKP,B,M) = D_F 
               ENDIF 
            ENDIF 
!
!           West face (i, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
!
	       Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IJK,M))
               D_F = MU_S(IJKC,M)*ODX(I)*AYZ_U(IMJK) 
               IF (Flux >= ZERO) THEN 
                  A_U_S(IJK,W,M) = D_F + Flux
               ELSE 
                  A_U_S(IJK,W,M) = D_F 
               ENDIF 
            ENDIF 
!
!           South face (i+1/2, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IJMK)) THEN 
               JM = JM1(J) 
               IPJMK = IP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKSE = EAST_OF(IJKS) 
	       Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IPJMK,M))
               D_F = AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),AVG_Y_H(MU_S&
                  (IJKSE,M),MU_S(IJKE,M),JM),I)*ODY_N(JM)*AXZ_U(IJMK) 
               IF (Flux >= ZERO) THEN 
                  A_U_S(IJK,S,M) = D_F + Flux
               ELSE 
                  A_U_S(IJK,S,M) = D_F 
               ENDIF 
            ENDIF 
!
!           Bottom face (i+1/2, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLOW_AT_E(IJKM)) THEN 
                  KM = KM1(K) 
                  IPJKM = IP_OF(IJKM) 
                  IJKB = BOTTOM_OF(IJK) 
                  IJKBE = EAST_OF(IJKB) 
	          Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IPJKM,M))
                  D_F = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),AVG_Z_H(&
                     MU_S(IJKBE,M),MU_S(IJKE,M),KM),I)*OX_E(I)*ODZ_T(KM)*AXY_U(&
                     IJKM) 
                  IF (Flux >= ZERO) THEN 
                     A_U_S(IJK,B,M) = D_F + Flux
                  ELSE 
                     A_U_S(IJK,B,M) = D_F 
                  ENDIF 
               ENDIF 
            ENDIF 
!
         ENDIF 
      END DO 
      
      RETURN  
      END SUBROUTINE STORE_A_U_S0 

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_U_sdc(A_U_s, M, B_M, IER)                           C
!  Purpose:TO USE DEFERRED CORRECTION METHOD TO SOLVE THE U-MOMENTUM   C
!  EQUATION. THIS METHOD COMBINES FIRST ORDER UPWIND AND A USER        C
!  SPECIFIED HIGH ORDER METHOD                                         C
!                                                                      C
!  Author: C. GUENTHER                                 Date: 8-APR-99  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_U_SDC(A_U_S, M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc 
      USE fldvar
      USE output
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar  
      USE sendrecv
      USE sendrecv3
      USE mflux 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
! 
!                      Indices 
      INTEGER          I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,& 
                       IJKNE, IJKP, IJKT, IJKTE 
      INTEGER          IMJK, IM, IJKW 
      INTEGER          IJMK, JM, IPJMK, IJKS, IJKSE 
      INTEGER          IJKM, KM, IPJKM, IJKB, IJKBE 
      INTEGER          IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER          IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
! 
!                      Solids phase 
      INTEGER          M 
! 
! loezos
	INTEGER  incr
! loezos

!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
!
!                      Septadiagonal matrix A_U_s 
      DOUBLE PRECISION A_U_s(DIMENSION_3, -3:3, M:M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
!
!	DEFERRED CORRCTION CONTRIBUTION FORM HIGH ORDER METHOD
	DOUBLE PRECISION MOM_HO
!
!	LOW ORDER APPROXIMATION 
	DOUBLE PRECISION MOM_LO
!
!	CONVECTION FACTOR AT THE FACE
	DOUBLE PRECISION Flux
!
!	DEFERRED CORRECTION CONTRIBUTIONS FROM EACH FACE
	DOUBLE PRECISION  EAST_DC
	DOUBLE PRECISION  WEST_DC
	DOUBLE PRECISION  NORTH_DC
	DOUBLE PRECISION  SOUTH_DC
        DOUBLE PRECISION  TOP_DC
        DOUBLE PRECISION  BOTTOM_DC
!
!
!
!---------------------------------------------------------------
!	EXTERNAL FUNCTIONS
!---------------------------------------------------------------
	DOUBLE PRECISION , EXTERNAL :: FPFOI_OF
!---------------------------------------------------------------
! 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'function3.inc'


      call lock_tmp4_array

      
      call lock_tmp_array
      call lock_xsi_array      
!
!  Calculate convection factors
!
!
! Send recv the third ghost layer
      IF ( FPFOI ) THEN
         Do IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = U_S(IJK,M)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

!
!$omp parallel do private(IJK,I,IP,IPJK,IJKE)
      DO IJK = ijkstart3, ijkend3
!
         I = I_OF(IJK) 
         IP = IP1(I) 
         IPJK = IP_OF(IJK) 
         IJKE = EAST_OF(IJK) 
!
!
!           East face (i+1, j, k)
         U(IJK) = AVG_X_E(U_S(IJK,M),U_S(IPJK,M),IP) 
!
!
!           North face (i+1/2, j+1/2, k)
         V(IJK) = AVG_X(V_S(IJK,M),V_S(IPJK,M),I) 
!
!
!           Top face (i+1/2, j, k+1/2)
         IF (DO_K) WW(IJK) = AVG_X(W_S(IJK,M),W_S(IPJK,M),I) 
      END DO 

! loezos
	incr=1		
! loezos

      CALL CALC_XSI (DISCRETIZE(3), U_S(1,M), U, V, WW, XSI_E, XSI_N,&
	XSI_T,incr) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!$omp      parallel do 	&
!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,	&
!$omp&                    IJKNE, IJKP, IJKT, IJKTE,  D_f,	&
!$omp&                    IMJK, IM, IJKW,	&
!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,	&
!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE, &
!$omp&              MOM_HO, MOM_LO, CONV_FAC,EAST_DC,WEST_DC,NORTH_DC,&
!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_E(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKE = EAST_OF(IJK) 
!
!           Third Ghost layer information
            IPPP  = IP_OF(IP_OF(IPJK))
            IPPP4 = funijk3(I_OF(IPPP), J_OF(IPPP), K_OF(IPPP))
            IMMM  = IM_OF(IM_OF(IMJK))
            IMMM4 = funijk3(I_OF(IMMM), J_OF(IMMM), K_OF(IMMM))
            JPPP  = JP_OF(JP_OF(IJPK))
            JPPP4 = funijk3(I_OF(JPPP), J_OF(JPPP), K_OF(JPPP))
            JMMM  = JM_OF(JM_OF(IJMK))
            JMMM4 = funijk3(I_OF(JMMM), J_OF(JMMM), K_OF(JMMM))
            KPPP  = KP_OF(KP_OF(IJKP))
            KPPP4 = funijk3(K_OF(IPPP), J_OF(KPPP), K_OF(KPPP))
            KMMM  = KM_OF(KM_OF(IJKM))
            KMMM4 = funijk3(I_OF(KMMM), J_OF(KMMM), K_OF(KMMM))
!
!
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKE 
            ELSE 
               IJKC = IJK 
            ENDIF 
            IP = IP1(I) 
            IJKN = NORTH_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1, j, k)
!           
		IF(U(IJK) >= ZERO)THEN
		   MOM_LO = U_S(IJK,M)
                   IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_S(IPJK,M), U_S(IJK,M), & 
                            U_S(IMJK,M), U_S(IM_OF(IMJK),M))
		ELSE
		   MOM_LO = U_S(IPJK,M)
                   IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IPJK,M), & 
                            U_S(IP_OF(IPJK),M), TMP4(IPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_E(IJK)*U_S(IPJK,M)+ &
                                (1.0-XSI_E(IJK))*U_S(IJK,M)
	        Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IPJK,M))
		EAST_DC = Flux *(MOM_LO - MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i+1/2, j+1/2, k)
!            
		IF(V(IJK) >= ZERO)THEN
		    MOM_LO = U_S(IJK,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_S(IJPK,M), U_S(IJK,M), & 
                            U_S(IJMK,M), U_S(JM_OF(IJMK),M))
		ELSE
		    MOM_LO = U_S(IJPK,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJPK,M), & 
                            U_S(JP_OF(IJPK),M), TMP4(JPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_N(IJK)*U_S(IJPK,M)+ &
                              (1.0-XSI_N(IJK))*U_S(IJK,M)
	        Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IPJK,M))
		NORTH_DC = Flux * (MOM_LO - MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i+1/2, j, k+1/2)
!
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 
	       IF(WW(IJK) >= ZERO)THEN
		  MOM_LO = U_S(IJK,M)
                  IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_S(IJKP,M), U_S(IJK,M), & 
                            U_S(IJKM,M), U_S(KM_OF(IJKM),M))
	       ELSE
		 MOM_LO = U_S(IJKP,M)
                 IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJKP,M), & 
                            U_S(KP_OF(IJKP),M), TMP4(KPPP4))
	       ENDIF
               IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_T(IJK)*U_S(IJKP,M)+ &
                                 (1.0-XSI_T(IJK))*U_S(IJK,M)
	       Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IPJK,M))
	       TOP_DC = Flux *(MOM_LO - MOM_HO)
	    ELSE
	       TOP_DC = ZERO
	    
            ENDIF
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE West face (i, j, k)
!
            IMJK = IM_OF(IJK) 
            IM = IM1(I) 
            IJKW = WEST_OF(IJK) 
	    IF(U(IMJK) >= ZERO)THEN
	      MOM_LO = U_S(IMJK,M)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IMJK,M), & 
                            U_S(IM_OF(IMJK),M), TMP4(IMMM4))
	    ELSE
	      MOM_LO = U_S(IJK,M)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_S(IMJK,M), U_S(IJK,M), & 
                            U_S(IPJK,M), U_S(IP_OF(IPJK),M))
	    ENDIF
            IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_E(IMJK)*U_S(IJK,M)+ &
                               (1.0-XSI_E(IMJK))*U_S(IMJK,M)
	    Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IJK,M))
	    WEST_DC = Flux * (MOM_LO - MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i+1/2, j-1/2, k)
!
            IJMK = JM_OF(IJK) 
            JM = JM1(J) 
            IPJMK = IP_OF(IJMK) 
            IJKS = SOUTH_OF(IJK) 
            IJKSE = EAST_OF(IJKS) 
	    IF(V(IJMK) >= ZERO)THEN
	       MOM_LO = U_S(IJMK,M)
               IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJMK,M), & 
                            U_S(JM_OF(IJMK),M), TMP4(JMMM4))
	    ELSE
	       MOM_LO = U_S(IJK,M)
               IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_S(IJMK,M), U_S(IJK,M), & 
                            U_S(IJPK,M), U_S(JP_OF(IJPK),M))
	    ENDIF
            IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_N(IJMK)*U_S(IJK,M)+ &
                               (1.0-XSI_N(IJMK))*U_S(IJMK,M)
	    Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IPJMK,M))
	    SOUTH_DC = Flux * (MOM_LO - MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i+1/2, j, k-1/2)
!
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               KM = KM1(K) 
               IPJKM = IP_OF(IJKM) 
               IJKB = BOTTOM_OF(IJK) 
               IJKBE = EAST_OF(IJKB) 
	       IF(WW(IJK) >= ZERO)THEN
		    MOM_LO = U_S(IJKM,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJKM,M), & 
                            U_S(KM_OF(IJKM),M), TMP4(KMMM4))
	       ELSE
		    MOM_LO = U_S(IJK,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(U_S(IJKM,M), U_S(IJK,M), & 
                            U_S(IJKP,M), U_S(KP_OF(IJKP),M))
	       ENDIF
               IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_T(IJKM)*U_S(IJK,M)+ &
                               (1.0-XSI_T(IJKM))*U_S(IJKM,M)
	       Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IPJKM,M))
	       BOTTOM_DC = Flux * (MOM_LO - MOM_HO)
            ELSE
	       BOTTOM_DC = ZERO
            ENDIF
!
!		CONTRIBUTION DUE TO DEFERRED CORRECTION
!
		B_M(IJK,M) = B_M(IJK,M)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC&
				+BOTTOM_DC-TOP_DC
!
         ENDIF 
      END DO 
     
      call unlock_tmp4_array
      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_U_SDC 


!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_U_s1(A_U_s, M, IER)
!  Purpose: Determine convection diffusion terms for U_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; Higher order C
!  See source_u_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_U_S1(A_U_S, M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc 
      USE fldvar
      USE output
      USE vshear
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar  
      USE mflux 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
! 
!                      Indices 
      INTEGER          I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,& 
                       IJKNE, IJKP, IJKT, IJKTE 
      INTEGER          IMJK, IM, IJKW 
      INTEGER          IJMK, JM, IPJMK, IJKS, IJKSE 
      INTEGER          IJKM, KM, IPJKM, IJKB, IJKBE 
! 
!                      Solids phase 
      INTEGER          M 
! 
! loezos                     
      INTEGER incr    
!loezos

!                      Face mass flux 
      DOUBLE PRECISION Flux 
      
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_U_s 
      DOUBLE PRECISION A_U_s(DIMENSION_3, -3:3, M:M) 
! 
!                      Convection weighting factors 
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),& 
!                       XSI_t(DIMENSION_3) 
!      DOUBLE PRECISION U(DIMENSION_3),& 
!                       V(DIMENSION_3), WW(DIMENSION_3) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'


      call lock_tmp_array
      call lock_xsi_array

!
!  Calculate convection factors
!
!$omp parallel do private(IJK,I,IP,IPJK,IJKE)
      DO IJK = ijkstart3, ijkend3
!
         I = I_OF(IJK) 
         IP = IP1(I) 
         IPJK = IP_OF(IJK) 
         IJKE = EAST_OF(IJK) 


!
!
!           East face (i+1, j, k)
         U(IJK) = AVG_X_E(U_S(IJK,M),U_S(IPJK,M),IP) 
!
!
!           North face (i+1/2, j+1/2, k)
         V(IJK) = AVG_X(V_S(IJK,M),V_S(IPJK,M),I) 
!
!
!           Top face (i+1/2, j, k+1/2)
         IF (DO_K) WW(IJK) = AVG_X(W_S(IJK,M),W_S(IPJK,M),I) 
      END DO 

! loezos
	incr=1		
! loezos

      CALL CALC_XSI (DISCRETIZE(3), U_S(1,M), U, V, WW, XSI_E, XSI_N, XSI_T,&
	incr) 

! loezos      
! update to true velocity
      IF (SHEAR) THEN
!$omp  parallel do private(IJK)
	 DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
	   V(IJK)=V(IJK)+VSHE(IJK)		
          END IF
        END DO
      END IF
! loezos

!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!$omp      parallel do 	&
!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,	&
!$omp&                    IJKNE, IJKP, IJKT, IJKTE,   D_f,	&
!$omp&                    IMJK, IM, IJKW,	&
!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,	&
!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE)
      DO IJK = ijkstart3, ijkend3 
!
         IF (FLOW_AT_E(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKE 
            ELSE 
               IJKC = IJK 
            ENDIF 
            IP = IP1(I) 
            IJKN = NORTH_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           East face (i+1, j, k)
	    Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IPJK,M))
            D_F = MU_S(IJKE,M)*ODX(IP)*AYZ_U(IJK) 
!
            A_U_S(IJK,E,M) = D_F - XSI_E(IJK) * Flux 
!
            A_U_S(IPJK,W,M) = D_F + (ONE - XSI_E(IJK)) * Flux 
!
!
!           North face (i+1/2, j+1/2, k)
	    Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IPJK,M))
            D_F = AVG_X_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(&
               IJKE,M),MU_S(IJKNE,M),J),I)*ODY_N(J)*AXZ_U(IJK) 
!
            A_U_S(IJK,N,M) = D_F - XSI_N(IJK) * Flux
!
            A_U_S(IJPK,S,M) = D_F + (ONE - XSI_N(IJK)) * Flux
!
!
!           Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 
!
	       Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IPJK,M))
               D_F = AVG_X_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),AVG_Z_H(MU_S(&
                  IJKE,M),MU_S(IJKTE,M),K),I)*OX_E(I)*ODZ_T(K)*AXY_U(IJK) 
!
               A_U_S(IJK,T,M) = D_F - XSI_T(IJK) * Flux 
!
               A_U_S(IJKP,B,M) = D_F + (ONE - XSI_T(IJK)) * Flux 
            ENDIF 
!
!           West face (i, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
!
	       Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IJK,M))
               D_F = MU_S(IJKC,M)*ODX(I)*AYZ_U(IMJK) 
!
               A_U_S(IJK,W,M) = D_F + (ONE - XSI_E(IMJK)) * Flux 
            ENDIF 
!
!           South face (i+1/2, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IJMK)) THEN 
               JM = JM1(J) 
               IPJMK = IP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKSE = EAST_OF(IJKS) 
!
	       Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IPJMK,M))
               D_F = AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),AVG_Y_H(MU_S&
                  (IJKSE,M),MU_S(IJKE,M),JM),I)*ODY_N(JM)*AXZ_U(IJMK) 
!
               A_U_S(IJK,S,M) = D_F + (ONE - XSI_N(IJMK)) * Flux 
            ENDIF 
!
!           Bottom face (i+1/2, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLOW_AT_E(IJKM)) THEN 
                  KM = KM1(K) 
                  IPJKM = IP_OF(IJKM) 
                  IJKB = BOTTOM_OF(IJK) 
                  IJKBE = EAST_OF(IJKB) 
!
	          Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IPJKM,M))
                  D_F = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),AVG_Z_H(&
                     MU_S(IJKBE,M),MU_S(IJKE,M),KM),I)*OX_E(I)*ODZ_T(KM)*AXY_U(&
                     IJKM) 
!
                  A_U_S(IJK,B,M) = D_F + (ONE - XSI_T(IJKM)) * Flux 
               ENDIF 
            ENDIF 
!
         ENDIF 
      END DO 
      
      call unlock_tmp_array
      call unlock_xsi_array

      
      RETURN  
      END SUBROUTINE STORE_A_U_S1 
!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3

                                                                         conv_dif_v_g.f                                                                                      0100644 0002444 0000146 00000075266 10250052576 012034  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_V_g(A_m, B_m, IER)                            C
!  Purpose: Determine convection diffusion terms for V_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_v_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_DIF_V_G(A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE compar 
      
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!
      IF (.NOT.MOMENTUM_Y_EQ(0)) RETURN
!	IF DEFERRED CORRECTION IS TO BE USED TO SOLVE V_G
!
      IF (DEF_COR) THEN
	CALL STORE_A_V_G0 (A_M(1,-3,0), IER) 
        IF (DISCRETIZE(4) > 1)CALL STORE_A_V_GDC (A_M(1,-3,0), B_M(1,0), IER) 
      ELSE  
!
        IF (DISCRETIZE(4) == 0) THEN               ! 0 & 1 => FOUP 
          CALL STORE_A_V_G0 (A_M(1,-3,0), IER) 
        ELSE 
          CALL STORE_A_V_G1 (A_M(1,-3,0), IER) 
        ENDIF 
      ENDIF

      CALL DIF_V_IS (MU_GT, A_M, B_M, 0, IER) 

      RETURN  
      END SUBROUTINE CONV_DIF_V_G 
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_V_g0(A_V_g, IER)                               C
!  Purpose: Determine convection diffusion terms for V_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; FOUP         C
!  See source_v_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-JUN-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_V_G0(A_V_G, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE toleranc 
      USE physprop
      USE fldvar
      USE output
      USE compar
      USE mflux  
      
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, JP, IJKE,& 
                       IJKNE, IJKP, IJKT, IJKTN, IJK 
! 
!                      Indices 
      INTEGER          IMJK, IM, IJKW, IJKWN, IMJPK 
      INTEGER          IJMK, JM, IJKS 
      INTEGER          IJKM, KM, IJKB, IJKBN, IJPKM 
! 
!                      Solids phase 
      INTEGER          M 
! 
!                      Face mass flux 
      DOUBLE PRECISION Flux 
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_V_g 
      DOUBLE PRECISION A_V_g(DIMENSION_3, -3:3) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'


!
!  Calculate convection-diffusion fluxes through each of the faces
!
!     Fluid phase
      M = 0 
      
!$omp      parallel do                                            &
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,  &
!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK, V_f, D_f, &
!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,                  &
!$omp&             IJMK, JM, IJKS,                                &
!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM )
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_N(IJK)) THEN 
!
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK)	     
            IJKN = NORTH_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKN 
            ELSE 
               IJKC = IJK 
            ENDIF 
            JP = JP1(J) 
            IJKE = EAST_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!
!           East face (i+1/2, j+1/2, k)
	    Flux = HALF * (Flux_gE(IJK) + Flux_gE(IJPK))
            D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKN&
               ),MU_GT(IJKNE),I),J)*ODX_E(I)*AYZ_V(IJK) 
            IF (Flux >= ZERO) THEN 
               A_V_G(IJK,E) = D_F 
               A_V_G(IPJK,W) = D_F + Flux
            ELSE 
               A_V_G(IJK,E) = D_F - Flux
               A_V_G(IPJK,W) = D_F 
            ENDIF 
!
!           North face (i, j+1, k)
	    Flux = HALF * (Flux_gN(IJK) + Flux_gN(IJPK))
            D_F = MU_GT(IJKN)*ODY(JP)*AXZ_V(IJK) 
            IF (Flux >= ZERO) THEN 
               A_V_G(IJK,N) = D_F 
               A_V_G(IJPK,S) = D_F + Flux
            ELSE 
               A_V_G(IJK,N) = D_F - Flux
               A_V_G(IJPK,S) = D_F 
            ENDIF 
!
!           Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN 
               IJKT = TOP_OF(IJK) 
               IJKTN = NORTH_OF(IJKT) 
	       Flux = HALF * (Flux_gT(IJK) + Flux_gT(IJPK))
               D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                  IJKN),MU_GT(IJKTN),K),J)*OX(I)*ODZ_T(K)*AXY_V(IJK) 
               IF (Flux >= ZERO) THEN 
                  A_V_G(IJK,T) = D_F 
                  A_V_G(IJKP,B) = D_F + Flux
               ELSE 
                  A_V_G(IJK,T) = D_F - Flux
                  A_V_G(IJKP,B) = D_F 
               ENDIF 
            ENDIF 
!
!           West face (i-1/2, j+1/2, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_N(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
               IJKWN = NORTH_OF(IJKW) 
               IMJPK = JP_OF(IMJK) 
	       Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IMJPK))
               D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                  IJKWN),MU_GT(IJKN),IM),J)*ODX_E(IM)*AYZ_V(IMJK) 
               IF (Flux >= ZERO) THEN 
                  A_V_G(IJK,W) = D_F + Flux
               ELSE 
                  A_V_G(IJK,W) = D_F 
               ENDIF 
            ENDIF 
!
!           South face (i, j, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_N(IJMK)) THEN 
               JM = JM1(J) 
               IJKS = SOUTH_OF(IJK) 
	       Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IJK))
               D_F = MU_GT(IJKC)*ODY(J)*AXZ_V(IJMK) 
               IF (Flux >= ZERO) THEN 
                  A_V_G(IJK,S) = D_F + Flux
               ELSE 
                  A_V_G(IJK,S) = D_F 
               ENDIF 
            ENDIF 
!
!           Bottom face (i, j+1/2, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLOW_AT_N(IJKM)) THEN 
                  KM = KM1(K) 
                  IJKB = BOTTOM_OF(IJK) 
                  IJKBN = NORTH_OF(IJKB) 
                  IJPKM = JP_OF(IJKM) 
	          Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IJPKM))
                  D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                     MU_GT(IJKBN),MU_GT(IJKN),KM),J)*OX(I)*ODZ_T(KM)*AXY_V(IJKM&
                     ) 
                  IF (Flux >= ZERO) THEN 
                     A_V_G(IJK,B) = D_F + Flux
                  ELSE 
                     A_V_G(IJK,B) = D_F 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE STORE_A_V_G0

!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_V_gdc(A_V_g, B_M, IER)                         C
!  Purpose: TO USE DEFERRED CORRECTION METHOD TO SOLVE THE V-MOMENTUM  C
!  EQUATION. THIS METHOD COMBINES FIRST ORDER UPWIND AND A USER        C
!  SPECIFIED HIGH ORDER METHOD                                         C
!                                                                      C
!  Author: C. GUENTHER                                 Date:8-APR-99   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_V_GDC(A_V_G, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE toleranc 
      USE physprop
      USE fldvar
      USE output
      Use xsi_array
      USE vshear
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar    
      USE sendrecv
      USE sendrecv3
      USE mflux

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, JP, IJKE,& 
                       IJKNE, IJKP, IJKT, IJKTN, IJK 
      INTEGER          IMJK, IM, IJKW, IJKWN, IMJPK 
      INTEGER          IJMK, JM, IJKS 
      INTEGER          IJKM, KM, IJKB, IJKBN, IJPKM 
      INTEGER          IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER          IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
! 
! loezos
	INTEGER incr
! loezos

!                      Diffusion parameter 
      DOUBLE PRECISION D_f
! 
!                      Septadiagonal matrix A_V_g 
      DOUBLE PRECISION A_V_g(DIMENSION_3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3) 
!
!	DEFERRED CORRCTION CONTRIBUTION FORM HIGH ORDER METHOD
	DOUBLE PRECISION MOM_HO
!
!	LOW ORDER APPROXIMATION 
	DOUBLE PRECISION MOM_LO
!
!	CONVECTION FACTOR AT THE FACE
	DOUBLE PRECISION Flux
!
!	DEFERRED CORRECTION CONTRIBUTIONS FROM EACH FACE
	DOUBLE PRECISION 	EAST_DC
	DOUBLE PRECISION 	WEST_DC
	DOUBLE PRECISION 	NORTH_DC
	DOUBLE PRECISION 	SOUTH_DC
        DOUBLE PRECISION  TOP_DC
        DOUBLE PRECISION  BOTTOM_DC
!
! 
!-----------------------------------------------
!
!---------------------------------------------------------------
!	EXTERNAL FUNCTIONS
!---------------------------------------------------------------
	DOUBLE PRECISION , EXTERNAL :: FPFOI_OF
!---------------------------------------------------------------
!
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'function3.inc'


      call lock_tmp4_array


      call lock_tmp_array
      call lock_xsi_array

!
!  Calculate convection factors
!
!
! Send recv the third ghost layer
      IF ( FPFOI ) THEN
         Do IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = V_G(IJK)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

!$omp parallel do private(IJK,J,IJPK,IJKN)
      DO IJK = ijkstart3, ijkend3 
         J = J_OF(IJK) 
         IJPK = JP_OF(IJK) 
         IJKN = NORTH_OF(IJK) 
!
!
!           East face (i+1/2, j+1/2, k)
         U(IJK) = AVG_Y(U_G(IJK),U_G(IJPK),J) 
!
!
!           North face (i, j+1, k)
         V(IJK) = AVG_Y_N(V_G(IJK),V_G(IJPK)) 
!
!
!           Top face (i, j+1/2, k+1/2)
         IF (DO_K) WW(IJK) = AVG_Y(W_G(IJK),W_G(IJPK),J) 
      END DO 

! loezos
	incr=2		
! loezos

      CALL CALC_XSI (DISCRETIZE(4), V_G, U, V, WW, XSI_E, XSI_N, XSI_T,incr) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!

!$omp      parallel do                                             &
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,   &
!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK,  D_f,      &
!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,                   &
!$omp&             IJMK, JM, IJKS,                                 &
!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM , &
!$omp&              MOM_HO, MOM_LO, CONV_FAC,EAST_DC,WEST_DC,NORTH_DC,&
!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_N(IJK)) THEN 
!
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)
!
!           Third Ghost layer information
            IPPP  = IP_OF(IP_OF(IPJK))
            IPPP4 = funijk3(I_OF(IPPP), J_OF(IPPP), K_OF(IPPP))
            IMMM  = IM_OF(IM_OF(IMJK))
            IMMM4 = funijk3(I_OF(IMMM), J_OF(IMMM), K_OF(IMMM))
            JPPP  = JP_OF(JP_OF(IJPK))
            JPPP4 = funijk3(I_OF(JPPP), J_OF(JPPP), K_OF(JPPP))
            JMMM  = JM_OF(JM_OF(IJMK))
            JMMM4 = funijk3(I_OF(JMMM), J_OF(JMMM), K_OF(JMMM))
            KPPP  = KP_OF(KP_OF(IJKP))
            KPPP4 = funijk3(K_OF(IPPP), J_OF(KPPP), K_OF(KPPP))
            KMMM  = KM_OF(KM_OF(IJKM))
            KMMM4 = funijk3(I_OF(KMMM), J_OF(KMMM), K_OF(KMMM))
!
!

            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK)
            IJKN = NORTH_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKN 
            ELSE 
               IJKC = IJK 
            ENDIF 
            JP = JP1(J) 
            IJKE = EAST_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1/2, j+1/2, k)
!            
		IF(U(IJK) >= ZERO)THEN
		    MOM_LO = V_G(IJK)
                     IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IPJK), V_G(IJK), & 
                            V_G(IMJK), V_G(IM_OF(IMJK)))
		ELSE
		    MOM_LO = V_G(IPJK)
                     IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IPJK), & 
                            V_G(IP_OF(IPJK)), TMP4(IPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_E(IJK)*V_G(IPJK)+(1.0-XSI_E(IJK))*V_G(IJK)
	        Flux = HALF * (Flux_gE(IJK) + Flux_gE(IJPK))
		EAST_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1, k)
!
		IF(V(IJK) >= ZERO)THEN
		    MOM_LO = V_G(IJK)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJPK), V_G(IJK), & 
                            V_G(IJMK), V_G(JM_OF(IJMK)))
		ELSE
		    MOM_LO = V_G(IJPK)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IJPK), & 
                            V_G(JP_OF(IJPK)), TMP4(JPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_N(IJK)*V_G(IJPK)+(1.0-XSI_N(IJK))*V_G(IJK)
	        Flux = HALF * (Flux_gN(IJK) + Flux_gN(IJPK))
		NORTH_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j+1/2, k+1/2)
!
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTN = NORTH_OF(IJKT) 
	       IF(WW(IJK) >= ZERO)THEN
		    MOM_LO = V_G(IJK)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJKP), V_G(IJK), & 
                            V_G(IJKM), V_G(KM_OF(IJKM)))
		ELSE
		    MOM_LO = V_G(IJKP)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IJKP), & 
                            V_G(KP_OF(IJKP)), TMP4(KPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_T(IJK)*V_G(IJKP)+(1.0-XSI_T(IJK))*V_G(IJK)
	        Flux = HALF * (Flux_gT(IJK) + Flux_gT(IJPK))
		TOP_DC = Flux*(MOM_LO-MOM_HO)
	    ELSE
		TOP_DC = ZERO
	    
            ENDIF
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE West face (i-1/2, j+1/2, k)
!
            IMJK = IM_OF(IJK) 
            IM = IM1(I) 
            IJKW = WEST_OF(IJK) 
            IJKWN = NORTH_OF(IJKW) 
            IMJPK = JP_OF(IMJK) 
	    IF(U(IMJK) >= ZERO)THEN
	      MOM_LO = V_G(IMJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IMJK), & 
                            V_G(IM_OF(IMJK)), TMP4(IMMM4))
	    ELSE
	      MOM_LO = V_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IMJK), V_G(IJK), & 
                            V_G(IPJK), V_G(IP_OF(IPJK)))
	    ENDIF
            IF (.NOT. FPFOI ) &
              MOM_HO = XSI_E(IMJK)*V_G(IJK)+(1.0-XSI_E(IMJK))*V_G(IMJK)
	    Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IMJPK))
	    WEST_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i, j, k)
!
            IJMK = JM_OF(IJK) 
            JM = JM1(J) 
            IJKS = SOUTH_OF(IJK) 
	    IF(V(IJMK) >= ZERO)THEN
	      MOM_LO = V_G(IJMK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IJMK), & 
                            V_G(JM_OF(IJMK)), TMP4(JMMM4))
	    ELSE
	      MOM_LO = V_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJMK), V_G(IJK), & 
                            V_G(IJPK), V_G(JP_OF(IJPK)))
	    ENDIF
            IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_N(IJMK)*V_G(IJK)+(1.0-XSI_N(IJMK))*V_G(IJMK)
	    Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IJK))
	    SOUTH_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i, j+1/2, k-1/2)
!
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               KM = KM1(K) 
               IJKB = BOTTOM_OF(IJK) 
               IJKBN = NORTH_OF(IJKB) 
               IJPKM = JP_OF(IJKM) 
               IF(WW(IJK) >= ZERO)THEN
		 MOM_LO = V_G(IJKM)
                 IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJK), V_G(IJKM), & 
                            V_G(KM_OF(IJKM)), TMP4(KMMM4))
	       ELSE
		 MOM_LO = V_G(IJK)
                 IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_G(IJKM), V_G(IJK), & 
                            V_G(IJKP), V_G(KP_OF(IJKP)))
	       ENDIF
               IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_T(IJKM)*V_G(IJK)+(1.0-XSI_T(IJKM))*V_G(IJKM)
	       Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IJPKM))
	       BOTTOM_DC = Flux*(MOM_LO-MOM_HO)
            ELSE
	       BOTTOM_DC = ZERO
            ENDIF
!
!		CONTRIBUTION DUE TO DEFERRED CORRECTION
!
	    B_M(IJK) = B_M(IJK)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC&
				+BOTTOM_DC-TOP_DC
! 
         ENDIF 
      END DO 

      call unlock_tmp4_array
      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_V_GDC 

 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_V_g1(A_V_g, IER)                               C
!  Purpose: Determine convection diffusion terms for V_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; Higher order C
!  See source_v_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date:20-MAR-97   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_V_G1(A_V_G, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE toleranc 
      USE physprop
      USE fldvar
      USE output
      Use xsi_array
      USE vshear
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar 
      USE mflux   
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, JP, IJKE,& 
                       IJKNE, IJKP, IJKT, IJKTN, IJK 
      INTEGER          IMJK, IM, IJKW, IJKWN, IMJPK 
      INTEGER          IJMK, JM, IJKS 
      INTEGER          IJKM, KM, IJKB, IJKBN, IJPKM 
!
! start loezos
      INTEGER incr   
! end loezos
 
!                      Face mass flux 
      DOUBLE PRECISION Flux 
 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_V_g 
      DOUBLE PRECISION A_V_g(DIMENSION_3, -3:3) 
! 
!                      Convection weighting factors 
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),& 
!                       XSI_t(DIMENSION_3) 
!      DOUBLE PRECISION U(DIMENSION_3),& 
!                       V(DIMENSION_3), WW(DIMENSION_3) 
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'

      call lock_tmp_array
      call lock_xsi_array

!
!  Calculate convection factors
!
!$omp parallel do private(IJK,J,IJPK,IJKN)
      DO IJK = ijkstart3, ijkend3 
         J = J_OF(IJK) 
         IJPK = JP_OF(IJK) 
         IJKN = NORTH_OF(IJK) 
!
!
!           East face (i+1/2, j+1/2, k)
         U(IJK) = AVG_Y(U_G(IJK),U_G(IJPK),J) 
!
!
!           North face (i, j+1, k)
         V(IJK) = AVG_Y_N(V_G(IJK),V_G(IJPK)) 
!
!
!           Top face (i, j+1/2, k+1/2)
         IF (DO_K) WW(IJK) = AVG_Y(W_G(IJK),W_G(IJPK),J) 
      END DO 

! loezos
	incr=2		
! loezos

      CALL CALC_XSI (DISCRETIZE(4), V_G, U, V, WW, XSI_E, XSI_N, XSI_T,incr) 

! loezos    
! update to true velocity
      IF (SHEAR) THEN
!$omp parallel do private(IJK)  
	 DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
	   V(IJK)=V(IJK)+VSH(IJK)	
          END IF
        END DO
      END IF
! loezos
	

!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do                                             &
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,   &
!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK,  D_f,      &
!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,                   &
!$omp&             IJMK, JM, IJKS,                                 &
!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM )
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_N(IJK)) THEN 
!
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKN 
            ELSE 
               IJKC = IJK 
            ENDIF 
            JP = JP1(J) 
            IJKE = EAST_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           East face (i+1/2, j+1/2, k)
	    Flux = HALF * (Flux_gE(IJK) + Flux_gE(IJPK))
            D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKN&
               ),MU_GT(IJKNE),I),J)*ODX_E(I)*AYZ_V(IJK) 
!
            A_V_G(IJK,E) = D_F - XSI_E(IJK)*Flux
!
            A_V_G(IPJK,W) = D_F + (ONE - XSI_E(IJK))*Flux
!
!
!           North face (i, j+1, k)
	    Flux = HALF * (Flux_gN(IJK) + Flux_gN(IJPK))
            D_F = MU_GT(IJKN)*ODY(JP)*AXZ_V(IJK) 
            A_V_G(IJK,N) = D_F - XSI_N(IJK)*Flux
!
            A_V_G(IJPK,S) = D_F + (ONE - XSI_N(IJK))*Flux
!
!
!           Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTN = NORTH_OF(IJKT) 
	       Flux = HALF * (Flux_gT(IJK) + Flux_gT(IJPK))
               D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKC),MU_GT(IJKT),K),AVG_Z_H(MU_GT(&
                  IJKN),MU_GT(IJKTN),K),J)*OX(I)*ODZ_T(K)*AXY_V(IJK) 
!
               A_V_G(IJK,T) = D_F - XSI_T(IJK)*Flux
!
               A_V_G(IJKP,B) = D_F + (ONE - XSI_T(IJK))*Flux
            ENDIF 
!
!           West face (i-1/2, j+1/2, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_N(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
               IJKWN = NORTH_OF(IJKW) 
               IMJPK = JP_OF(IMJK) 
!
	       Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IMJPK))
               D_F = AVG_Y_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                  IJKWN),MU_GT(IJKN),IM),J)*ODX_E(IM)*AYZ_V(IMJK) 
!
               A_V_G(IJK,W) = D_F + (ONE - XSI_E(IMJK))*Flux
            ENDIF 
!
!           South face (i, j, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_N(IJMK)) THEN 
               JM = JM1(J) 
               IJKS = SOUTH_OF(IJK) 
!
	       Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IJK))
               D_F = MU_GT(IJKC)*ODY(J)*AXZ_V(IJMK) 
!
               A_V_G(IJK,S) = D_F + (ONE - XSI_N(IJMK))*Flux
            ENDIF 
!
!           Bottom face (i, j+1/2, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLOW_AT_N(IJKM)) THEN 
                  KM = KM1(K) 
                  IJKB = BOTTOM_OF(IJK) 
                  IJKBN = NORTH_OF(IJKB) 
                  IJPKM = JP_OF(IJKM) 
!
	          Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IJPKM))
                  D_F = AVG_Y_H(AVG_Z_H(MU_GT(IJKB),MU_GT(IJKC),KM),AVG_Z_H(&
                     MU_GT(IJKBN),MU_GT(IJKN),KM),J)*OX(I)*ODZ_T(KM)*AXY_V(IJKM&
                     ) 
!
                  A_V_G(IJK,B) = D_F + (ONE - XSI_T(IJKM))*Flux
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 



      call unlock_tmp_array
      call unlock_xsi_array

      
      RETURN  
      END SUBROUTINE STORE_A_V_G1 
      
!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
                                                                                                                                                                                                                                                                                                                                          conv_dif_v_s.f                                                                                      0100644 0002444 0000146 00000076616 10250057104 012040  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_V_s(A_m, B_m, IER)                            C
!  Purpose: Determine convection diffusion terms for V_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_v_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_DIF_V_S(A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE compar   
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Solids phase index 
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
  
 
!-----------------------------------------------
!
!



      DO M = 1, MMAX 
        IF  (MOMENTUM_Y_EQ(M)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	IF DEFERRED CORRECTION IS USED TO SOLVE V_S
           IF (DEF_COR) THEN
	     CALL STORE_A_V_S0 (A_M(1,-3,M), M, IER)
	     IF (DISCRETIZE(4) > 1)CALL STORE_A_V_SDC (A_M(1,-3,M), M, B_M, IER)
           ELSE    
!
             IF (DISCRETIZE(4) == 0) THEN         ! 0 & 1 => FOUP 
               CALL STORE_A_V_S0 (A_M(1,-3,M), M, IER) 
             ELSE 
               CALL STORE_A_V_S1 (A_M(1,-3,M), M, IER) 
             ENDIF 
           ENDIF
!
           CALL DIF_V_IS (MU_S(1,M), A_M, B_M, M, IER) 
	   
!
        ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE CONV_DIF_V_S 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_V_s0(A_V_s, M, IER)                            C
!  Purpose: Determine convection diffusion terms for V_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; FOUP         C
!  See source_v_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-JUN-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_V_S0(A_V_S, M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc 
      USE fldvar
      USE output
      USE compar 
      USE mflux  
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, JP, IJKE,& 
                       IJKNE, IJKP, IJKT, IJKTN, IJK 
      INTEGER          IMJK, IM, IJKW, IJKWN, IMJPK 
      INTEGER          IJMK, JM, IJKS 
      INTEGER          IJKM, KM, IJKB, IJKBN, IJPKM 
! 
!                      Solids phase 
      INTEGER          M 
! 
!                      Face mass flux
      DOUBLE PRECISION Flux 
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_V_s 
      DOUBLE PRECISION A_V_s(DIMENSION_3, -3:3, M:M) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do                                                  &
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,	&
!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK, V_f, D_f,	&
!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,	&
!$omp&             IJMK, JM, IJKS,	&
!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM ) 
      DO IJK = ijkstart3, ijkend3 
!
         IF (FLOW_AT_N(IJK)) THEN 
!
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 	    	    
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK)
            IJKN = NORTH_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKN 
            ELSE 
               IJKC = IJK 
            ENDIF 
            JP = JP1(J) 
            IJKE = EAST_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           East face (i+1/2, j+1/2, k)
	    Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IJPK,M))
            D_F = AVG_Y_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),AVG_X_H(MU_S(&
               IJKN,M),MU_S(IJKNE,M),I),J)*ODX_E(I)*AYZ_V(IJK) 
            IF (Flux >= ZERO) THEN 
               A_V_S(IJK,E,M) = D_F 
               A_V_S(IPJK,W,M) = D_F + Flux
            ELSE 
               A_V_S(IJK,E,M) = D_F - Flux
               A_V_S(IPJK,W,M) = D_F 
            ENDIF 
!
!           North face (i, j+1, k)
	    Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJPK,M))
            D_F = MU_S(IJKN,M)*ODY(JP)*AXZ_V(IJK) 
            IF (Flux >= ZERO) THEN 
               A_V_S(IJK,N,M) = D_F 
               A_V_S(IJPK,S,M) = D_F + Flux
            ELSE 
               A_V_S(IJK,N,M) = D_F - Flux 
               A_V_S(IJPK,S,M) = D_F 
            ENDIF 
!
!           Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTN = NORTH_OF(IJKT) 
	       Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJPK,M))
               D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),AVG_Z_H(MU_S(&
                  IJKN,M),MU_S(IJKTN,M),K),J)*OX(I)*ODZ_T(K)*AXY_V(IJK) 
               IF (Flux >= ZERO) THEN 
                  A_V_S(IJK,T,M) = D_F 
                  A_V_S(IJKP,B,M) = D_F + Flux
               ELSE 
                  A_V_S(IJK,T,M) = D_F - Flux
                  A_V_S(IJKP,B,M) = D_F 
               ENDIF 
            ENDIF 
!
!           West face (i-1/2, j+1/2, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_N(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
               IJKWN = NORTH_OF(IJKW) 
               IMJPK = JP_OF(IMJK) 
	       Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IMJPK,M))
               D_F = AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),AVG_X_H(MU_S&
                  (IJKWN,M),MU_S(IJKN,M),IM),J)*ODX_E(IM)*AYZ_V(IMJK) 
               IF (Flux >= ZERO) THEN 
                  A_V_S(IJK,W,M) = D_F + Flux
               ELSE 
                  A_V_S(IJK,W,M) = D_F 
               ENDIF 
            ENDIF 
!
!           South face (i, j, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_N(IJMK)) THEN 
               JM = JM1(J) 
               IJKS = SOUTH_OF(IJK) 
	       Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IJK,M))
               D_F = MU_S(IJKC,M)*ODY(J)*AXZ_V(IJMK) 
               IF (Flux >= ZERO) THEN 
                  A_V_S(IJK,S,M) = D_F + Flux
               ELSE 
                  A_V_S(IJK,S,M) = D_F 
               ENDIF 
            ENDIF 
!
!           Bottom face (i, j+1/2, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLOW_AT_N(IJKM)) THEN 
                  KM = KM1(K) 
                  IJKB = BOTTOM_OF(IJK) 
                  IJKBN = NORTH_OF(IJKB) 
                  IJPKM = JP_OF(IJKM) 
	          Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IJPKM,M))
                  D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),AVG_Z_H(&
                     MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)*OX(I)*ODZ_T(KM)*AXY_V(&
                     IJKM) 
                  IF (Flux >= ZERO) THEN 
                     A_V_S(IJK,B,M) = D_F + Flux
                  ELSE 
                     A_V_S(IJK,B,M) = D_F 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE STORE_A_V_S0

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_V_sdc(A_V_s, M, B_M, IER)                      C
!  Purpose: TO USE DEFERRED CORRECTION METHOD TO SOLVE THE U-MOMENTUM  C
!  EQUATION. THIS METHOD COMBINES FIRST ORDER UPWIND AND A USER        C
!  SPECIFIED HIGH ORDER METHOD                                         C
!                                                                      C
!  Author: C. GUENTHER                                 Date:8-APR-99   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_V_SDC(A_V_S, M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc 
      USE fldvar
      USE output
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar     
      USE sendrecv
      USE sendrecv3
      USE mflux
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, JP, IJKE,& 
                       IJKNE, IJKP, IJKT, IJKTN, IJK 
      INTEGER          IMJK, IM, IJKW, IJKWN, IMJPK 
      INTEGER          IJMK, JM, IJKS 
      INTEGER          IJKM, KM, IJKB, IJKBN, IJPKM 
      INTEGER          IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER          IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
! 
!                      Solids phase 
      INTEGER          M 
!
! loezos
	INTEGER  incr
! loezos
 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_V_s 
      DOUBLE PRECISION A_V_s(DIMENSION_3, -3:3, M:M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Convection weighting factors 
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),& 
!                       XSI_t(DIMENSION_3) 
!      DOUBLE PRECISION U(DIMENSION_3),& 
!                       V(DIMENSION_3), WW(DIMENSION_3)
!
!	DEFERRED CORRCTION CONTRIBUTION FORM HIGH ORDER METHOD
	DOUBLE PRECISION MOM_HO
!
!	LOW ORDER APPROXIMATION 
	DOUBLE PRECISION MOM_LO
!
!	CONVECTION FACTOR AT THE FACE
	DOUBLE PRECISION Flux
!
!	DEFERRED CORRECTION CONTRIBUTIONS FROM EACH FACE
	DOUBLE PRECISION 	EAST_DC
	DOUBLE PRECISION 	WEST_DC
	DOUBLE PRECISION 	NORTH_DC
	DOUBLE PRECISION 	SOUTH_DC
        DOUBLE PRECISION  TOP_DC
        DOUBLE PRECISION  BOTTOM_DC
!
! 
!-----------------------------------------------

!
!---------------------------------------------------------------
!	EXTERNAL FUNCTIONS
!---------------------------------------------------------------
	DOUBLE PRECISION , EXTERNAL :: FPFOI_OF
!---------------------------------------------------------------
!
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'function3.inc'


      call lock_tmp4_array

      call lock_tmp_array
      call lock_xsi_array
!
!  Calculate convection factors
!
!
! Send recv the third ghost layer
      IF ( FPFOI ) THEN
         Do IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = V_S(IJK,M)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

!$omp parallel do private(IJK,J,IJPK,IJKN)
      DO IJK = ijkstart3, ijkend3
         J = J_OF(IJK) 
         IJPK = JP_OF(IJK) 
         IJKN = NORTH_OF(IJK) 
!
!
!           East face (i+1/2, j+1/2, k)
         U(IJK) = AVG_Y(U_S(IJK,M),U_S(IJPK,M),J) 
!
!
!           North face (i, j+1, k)
         V(IJK) = AVG_Y_N(V_S(IJK,M),V_S(IJPK,M)) 
!
!
!           Top face (i, j+1/2, k+1/2)
         IF (DO_K) WW(IJK) = AVG_Y(W_S(IJK,M),W_S(IJPK,M),J) 
      END DO 

! loezos
	incr=2		
! loezos


      CALL CALC_XSI (DISCRETIZE(4), V_S(1,M), U, V, WW, XSI_E, XSI_N, XSI_T,&
			incr) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do 	&
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,	&
!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK,  D_f,	&
!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,	&
!$omp&             IJMK, JM, IJKS,	&
!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM, &
!$omp&              MOM_HO, MOM_LO, CONV_FAC,EAST_DC,WEST_DC,NORTH_DC,&
!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC )
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_N(IJK)) THEN 
!
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK)
            IJKN = NORTH_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKN 
            ELSE 
               IJKC = IJK 
            ENDIF 
            JP = JP1(J) 
            IJKE = EAST_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           Third Ghost layer information
            IPPP  = IP_OF(IP_OF(IPJK))
            IPPP4 = funijk3(I_OF(IPPP), J_OF(IPPP), K_OF(IPPP))
            IMMM  = IM_OF(IM_OF(IMJK))
            IMMM4 = funijk3(I_OF(IMMM), J_OF(IMMM), K_OF(IMMM))
            JPPP  = JP_OF(JP_OF(IJPK))
            JPPP4 = funijk3(I_OF(JPPP), J_OF(JPPP), K_OF(JPPP))
            JMMM  = JM_OF(JM_OF(IJMK))
            JMMM4 = funijk3(I_OF(JMMM), J_OF(JMMM), K_OF(JMMM))
            KPPP  = KP_OF(KP_OF(IJKP))
            KPPP4 = funijk3(K_OF(IPPP), J_OF(KPPP), K_OF(KPPP))
            KMMM  = KM_OF(KM_OF(IJKM))
            KMMM4 = funijk3(I_OF(KMMM), J_OF(KMMM), K_OF(KMMM))
!
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1/2, j+1/2, k)
!            
		IF(U(IJK) >= ZERO)THEN
		    MOM_LO = V_S(IJK,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_S(IPJK,M), V_S(IJK,M), & 
                            V_S(IMJK,M), V_S(IM_OF(IMJK),M))
		ELSE
		    MOM_LO = V_S(IPJK,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_S(IJK,M), V_S(IPJK,M), & 
                            V_S(IP_OF(IPJK),M), TMP4(IPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_E(IJK)*V_S(IPJK,M)+ &
                               (1.0-XSI_E(IJK))*V_S(IJK,M)
	        Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IJPK,M))
		EAST_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1, k)
!
		IF(V(IJK) >= ZERO)THEN
		    MOM_LO = V_S(IJK,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_S(IJPK,M), V_S(IJK,M), & 
                            V_S(IJMK,M), V_S(JM_OF(IJMK),M))
		ELSE
		    MOM_LO = V_S(IJPK,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_S(IJK,M), V_S(IJPK,M), & 
                            V_S(JP_OF(IJPK),M), TMP4(JPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_N(IJK)*V_S(IJPK,M)+ &
                               (1.0-XSI_N(IJK))*V_S(IJK,M)
	        Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJPK,M))
		NORTH_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j+1/2, k+1/2)
!
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTN = NORTH_OF(IJKT) 
	       IF(WW(IJK) >= ZERO)THEN
		    MOM_LO = V_S(IJK,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_S(IJKP,M), V_S(IJK,M), & 
                            V_S(IJKM,M), V_S(KM_OF(IJKM),M))
		ELSE
		    MOM_LO = V_S(IJKP,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_S(IJK,M), V_S(IJKP,M), & 
                            V_S(KP_OF(IJKP),M), TMP4(KPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		      MOM_HO = XSI_T(IJK)*V_S(IJKP,M)+ &
                               (1.0-XSI_T(IJK))*V_S(IJK,M)
	        Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJPK,M))
		TOP_DC = Flux*(MOM_LO-MOM_HO)
	    ELSE
		TOP_DC = ZERO
            ENDIF
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE West face (i-1/2, j+1/2, k)
!
            IMJK = IM_OF(IJK) 
            IM = IM1(I) 
            IJKW = WEST_OF(IJK) 
            IJKWN = NORTH_OF(IJKW) 
            IMJPK = JP_OF(IMJK) 
	    IF(U(IMJK) >= ZERO)THEN
	      MOM_LO = V_S(IMJK,M)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_S(IJK,M), V_S(IMJK,M), & 
                            V_S(IM_OF(IMJK),M), TMP4(IMMM4))
	    ELSE
	      MOM_LO = V_S(IJK,M)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_S(IMJK,M), V_S(IJK,M), & 
                            V_S(IPJK,M), V_S(IP_OF(IPJK),M))
	    ENDIF
            IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_E(IMJK)*V_S(IJK,M)+ &
                               (1.0-XSI_E(IMJK))*V_S(IMJK,M)
	    Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IMJPK,M))
	    WEST_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i, j, k)
!
            IJMK = JM_OF(IJK) 
            JM = JM1(J) 
            IJKS = SOUTH_OF(IJK) 
	    IF(V(IJMK) >= ZERO)THEN
	       MOM_LO = V_S(IJMK,M)
               IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_S(IJK,M), V_S(IJMK,M), & 
                            V_S(JM_OF(IJMK),M), TMP4(JMMM4))
	    ELSE
	       MOM_LO = V_S(IJK,M)
               IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_S(IJMK,M), V_S(IJK,M), & 
                            V_S(IJPK,M), V_S(JP_OF(IJPK),M))
	    ENDIF
            IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_N(IJMK)*V_S(IJK,M)+ &
                               (1.0-XSI_N(IJMK))*V_S(IJMK,M)
	    Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IJK,M))
	    SOUTH_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i, j+1/2, k-1/2)
!
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               KM = KM1(K) 
               IJKB = BOTTOM_OF(IJK) 
               IJKBN = NORTH_OF(IJKB) 
               IJPKM = JP_OF(IJKM) 
               IF(WW(IJK) >= ZERO)THEN
		 MOM_LO = V_S(IJKM,M)
                 IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_S(IJK,M), V_S(IJKM,M), & 
                            V_S(KM_OF(IJKM),M), TMP4(KMMM4))
	       ELSE
		 MOM_LO = V_S(IJK,M)
                 IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(V_S(IJKM,M), V_S(IJK,M), & 
                            V_S(IJKP,M), V_S(KP_OF(IJKP),M))
	       ENDIF
               IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_T(IJKM)*V_S(IJK,M)+ &
                               (1.0-XSI_T(IJKM))*V_S(IJKM,M)
	       Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IJPKM,M))
	       BOTTOM_DC = Flux*(MOM_LO-MOM_HO)
            ELSE
	       BOTTOM_DC = ZERO
            ENDIF
!
!	    CONTRIBUTION DUE TO DEFERRED CORRECTION
!
	    B_M(IJK,M) = B_M(IJK,M)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC&
				+BOTTOM_DC-TOP_DC
! 
         ENDIF 
      END DO 

      call unlock_tmp4_array
      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_V_SDC 

 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_V_s1(A_V_s, M, IER)                            C
!  Purpose: Determine convection diffusion terms for V_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; Higher order C
!  See source_v_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date:20-MAR-97   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_V_S1(A_V_S, M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc 
      USE fldvar
      USE output
      USE vshear
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar 
      USE mflux   
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, JP, IJKE,& 
                       IJKNE, IJKP, IJKT, IJKTN, IJK 
      INTEGER          IMJK, IM, IJKW, IJKWN, IMJPK 
      INTEGER          IJMK, JM, IJKS 
      INTEGER          IJKM, KM, IJKB, IJKBN, IJPKM 
! 
!                      Solids phase 
      INTEGER          M 
!
! start loezos
      INTEGER incr   
! end loezos
 
!                      Face mass flux 
      DOUBLE PRECISION Flux 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_V_s 
      DOUBLE PRECISION A_V_s(DIMENSION_3, -3:3, M:M) 
! 
!                      Convection weighting factors 
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),& 
!                       XSI_t(DIMENSION_3) 
!      DOUBLE PRECISION U(DIMENSION_3),& 
!                       V(DIMENSION_3), WW(DIMENSION_3) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'

      call lock_tmp_array
      call lock_xsi_array

!
!  Calculate convection factors
!

!$omp parallel do private(IJK,J,IJPK,IJKN)
      DO IJK = ijkstart3, ijkend3 
         J = J_OF(IJK) 
	
         IJPK = JP_OF(IJK) 
         IJKN = NORTH_OF(IJK) 
!
!
!           East face (i+1/2, j+1/2, k)
         U(IJK) = AVG_Y(U_S(IJK,M),U_S(IJPK,M),J) 
!
!
!           North face (i, j+1, k)
         V(IJK) = AVG_Y_N(V_S(IJK,M),V_S(IJPK,M)) 
!
!
!           Top face (i, j+1/2, k+1/2)
         IF (DO_K) WW(IJK) = AVG_Y(W_S(IJK,M),W_S(IJPK,M),J) 
      END DO 

! loezos
	incr=2		
! loezos

      CALL CALC_XSI (DISCRETIZE(4), V_S(1,M), U, V, WW, XSI_E, XSI_N, XSI_T,incr) 


! loezos      
! update to true velocity
      IF (SHEAR) THEN
!$omp      parallel do private(IJK)
	 DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
	   V(IJK)=V(IJK)+VSH(IJK)	
          END IF
        END DO

      END IF
! loezos


!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do 	&
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, JP,	&
!$omp&             IJKE, IJKNE, IJKP, IJKT, IJKTN, IJK,  D_f,	&
!$omp&             IMJK, IM, IJKW, IJKWN, IMJPK,	&
!$omp&             IJMK, JM, IJKS,	&
!$omp&             IJKM, KM, IJKB, IJKBN, IJPKM )
      DO IJK = ijkstart3, ijkend3 
!
         IF (FLOW_AT_N(IJK)) THEN 
!
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKN 
            ELSE 
               IJKC = IJK 
            ENDIF 
            JP = JP1(J) 
            IJKE = EAST_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           East face (i+1/2, j+1/2, k)
	    Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IJPK,M))
            D_F = AVG_Y_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),AVG_X_H(MU_S(&
               IJKN,M),MU_S(IJKNE,M),I),J)*ODX_E(I)*AYZ_V(IJK) 
!
            A_V_S(IJK,E,M) = D_F - XSI_E(IJK)*Flux
!
            A_V_S(IPJK,W,M) = D_F + (ONE - XSI_E(IJK))*Flux
!
!
!           North face (i, j+1, k)
	    Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJPK,M))
            D_F = MU_S(IJKN,M)*ODY(JP)*AXZ_V(IJK) 
            A_V_S(IJK,N,M) = D_F - XSI_N(IJK)*Flux 
!
            A_V_S(IJPK,S,M) = D_F + (ONE - XSI_N(IJK))*Flux
!
!
!           Top face (i, j+1/2, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTN = NORTH_OF(IJKT) 
	       Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJPK,M))
               D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),AVG_Z_H(MU_S(&
                  IJKN,M),MU_S(IJKTN,M),K),J)*OX(I)*ODZ_T(K)*AXY_V(IJK) 
!
               A_V_S(IJK,T,M) = D_F - XSI_T(IJK)*Flux
!
               A_V_S(IJKP,B,M) = D_F + (ONE - XSI_T(IJK))*Flux
            ENDIF 
!
!           West face (i-1/2, j+1/2, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_N(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
               IJKWN = NORTH_OF(IJKW) 
               IMJPK = JP_OF(IMJK) 
!
	       Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IMJPK,M))
               D_F = AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),AVG_X_H(MU_S&
                  (IJKWN,M),MU_S(IJKN,M),IM),J)*ODX_E(IM)*AYZ_V(IMJK) 
!
               A_V_S(IJK,W,M) = D_F + (ONE - XSI_E(IMJK))*Flux
            ENDIF 
!
!           South face (i, j, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_N(IJMK)) THEN 
               JM = JM1(J) 
               IJKS = SOUTH_OF(IJK) 
!
	       Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IJK,M))
               D_F = MU_S(IJKC,M)*ODY(J)*AXZ_V(IJMK) 
!
               A_V_S(IJK,S,M) = D_F + (ONE - XSI_N(IJMK))*Flux
            ENDIF 
!
!           Bottom face (i, j+1/2, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLOW_AT_N(IJKM)) THEN 
                  KM = KM1(K) 
                  IJKB = BOTTOM_OF(IJK) 
                  IJKBN = NORTH_OF(IJKB) 
                  IJPKM = JP_OF(IJKM) 
!
	          Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IJPKM,M))
                  D_F = AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),AVG_Z_H(&
                     MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)*OX(I)*ODZ_T(KM)*AXY_V(&
                     IJKM) 
!
                  A_V_S(IJK,B,M) = D_F + (ONE - XSI_T(IJKM))*Flux
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 


      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_V_S1 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
                                                                                                                  conv_dif_w_g.f                                                                                      0100644 0002444 0000146 00000074140 10250064554 012022  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_W_g(A_m, B_m, IER)                            C
!  Purpose: Determine convection diffusion terms for W_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_w_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_DIF_W_G(A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE compar   
      
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
 
!
      IF (.NOT.MOMENTUM_Z_EQ(0)) RETURN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	IF DEFERRED CORRECTION IS USED TO SOLVE W_G
      IF (DEF_COR) THEN
	CALL STORE_A_W_G0 (A_M(1,-3,0), IER) 
	IF (DISCRETIZE(5) > 1)CALL STORE_A_W_GDC (A_M(1,-3,0), B_M(1,0), IER)
      ELSE  
!
        IF (DISCRETIZE(5) == 0) THEN               ! 0 & 1 => FOUP 
          CALL STORE_A_W_G0 (A_M(1,-3,0), IER) 
        ELSE 
          CALL STORE_A_W_G1 (A_M(1,-3,0), IER) 
        ENDIF 
      ENDIF
!      
      CALL DIF_W_IS (MU_GT, A_M, B_M, 0, IER) 
!

      RETURN  
      END SUBROUTINE CONV_DIF_W_G 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_W_g0(A_W_g, IER)                               C
!  Purpose: Determine convection diffusion terms for W_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; FOUP         C
!  See source_w_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-JUN-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_W_G0(A_W_G, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE toleranc 
      USE physprop
      USE fldvar
      USE output
      USE compar 
      USE mflux    
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,& 
                       IJKTE, IJKP, IJKT, IJKTN, IJK 
      INTEGER          IMJK, IM, IJKW, IJKWT, IMJKP 
      INTEGER          IJMK, JM, IJMKP, IJKS, IJKST 
      INTEGER          IJKM, KM, IJKB 
! 
!                      Solids phase 
      INTEGER          M 
! 
!                      Face mass flux 
      DOUBLE PRECISION Flux
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_W_g 
      DOUBLE PRECISION A_W_g(DIMENSION_3, -3:3) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'

!
!  Calculate convection-diffusion fluxes through each of the faces
!
!     Fluid phase
      M = 0 

!$omp      parallel do                                               &
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,     &
!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK, V_f, D_f,    &
!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,                     &
!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,                     &
!$omp&             IJKM, KM, IJKB)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_T(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 	    
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKN = NORTH_OF(IJK ) 
            IJKT = TOP_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKT 
            ELSE 
               IJKC = IJK 
            ENDIF 
            KP = KP1(K) 
            IJKE = EAST_OF(IJK) 
            IJKP = KP_OF(IJK) 
            IJKTN = NORTH_OF(IJKT) 
            IJKTE = EAST_OF(IJKT) 
!
!           East face (i+1/2, j, k+1/2)
	    Flux = HALF * (Flux_gE(IJK) + Flux_gE(IJKP))
            D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKT&
               ),MU_GT(IJKTE),I),K)*ODX_E(I)*AYZ_W(IJK) 
            IF (Flux >= ZERO) THEN 
               A_W_G(IJK,E) = D_F 
               A_W_G(IPJK,W) = D_F + Flux
            ELSE 
               A_W_G(IJK,E) = D_F - Flux
               A_W_G(IPJK,W) = D_F 
            ENDIF 
!
!           North face (i, j+1/2, k+1/2)
	    Flux = HALF * (Flux_gN(IJK) + Flux_gN(IJKP))
            D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKT&
               ),MU_GT(IJKTN),J),K)*ODY_N(J)*AXZ_W(IJK) 
            IF (Flux >= ZERO) THEN 
               A_W_G(IJK,N) = D_F 
               A_W_G(IJPK,S) = D_F + Flux
            ELSE 
               A_W_G(IJK,N) = D_F - Flux
               A_W_G(IJPK,S) = D_F 
            ENDIF 
!
!           Top face (i, j, k+1)
	    Flux = HALF * (Flux_gT(IJK) + Flux_gT(IJKP))
            D_F = MU_GT(IJKT)*OX(I)*ODZ(KP)*AXY_W(IJK) 
            IF (Flux >= ZERO) THEN 
               A_W_G(IJK,T) = D_F 
               A_W_G(IJKP,B) = D_F + Flux
            ELSE 
               A_W_G(IJK,T) = D_F - Flux
               A_W_G(IJKP,B) = D_F 
            ENDIF 
!
!           West face (i-1/2, j, k+1/2)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
               IJKWT = TOP_OF(IJKW) 
               IMJKP = KP_OF(IMJK) 
	       Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IMJKP))
               D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                  IJKWT),MU_GT(IJKT),IM),K)*ODX_E(IM)*AYZ_W(IMJK) 
               IF (Flux >= ZERO) THEN 
                  A_W_G(IJK,W) = D_F + Flux
               ELSE 
                  A_W_G(IJK,W) = D_F 
               ENDIF 
            ENDIF 
!
!           South face (i, j-1/2, k+1/2)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IJMK)) THEN 
               JM = JM1(J) 
               IJMKP = KP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKST = TOP_OF(IJKS) 
	       Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IJMKP))
               D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                  IJKST),MU_GT(IJKT),JM),K)*ODY_N(JM)*AXZ_W(IJMK) 
               IF (Flux >= ZERO) THEN 
                  A_W_G(IJK,S) = D_F + Flux
               ELSE 
                  A_W_G(IJK,S) = D_F 
               ENDIF 
            ENDIF 
!
!           Bottom face (i, j, k)
            IJKM = KM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IJKM)) THEN 
               KM = KM1(K) 
               IJKB = BOTTOM_OF(IJK) 
	       Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IJK))
               D_F = MU_GT(IJK)*OX(I)*ODZ(K)*AXY_W(IJKM) 
               IF (Flux >= ZERO) THEN 
                  A_W_G(IJK,B) = D_F + Flux
               ELSE 
                  A_W_G(IJK,B) = D_F 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE STORE_A_W_G0 

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_W_gdc(A_W_g, B_M, IER)                         C
!  Purpose: TO USE DEFERRED CORRECTION METHOD TO SOLVE THE W-MOMENTUM  C
!  EQUATION. THIS METHOD COMBINES FIRST ORDER UPWIND AND A USER        C
!  SPECIFIED HIGH ORDER METHOD                                         C
!                                                                      C
!  Author: C. GUENTHER                                 Date:8-APR-99   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_W_GDC(A_W_G, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE toleranc 
      USE physprop
      USE fldvar
      USE output
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar    
      USE sendrecv
      USE sendrecv3
      USE mflux
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,& 
                       IJKTE, IJKP, IJKT, IJKTN, IJK 
      INTEGER          IMJK, IM, IJKW, IJKWT, IMJKP 
      INTEGER          IJMK, JM, IJMKP, IJKS, IJKST 
      INTEGER          IJKM, KM, IJKB 
      INTEGER          IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER          IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
!
! loezos
	INTEGER incr
! loezos
 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
!
!                      Septadiagonal matrix A_W_g 
      DOUBLE PRECISION A_W_g(DIMENSION_3, -3:3)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3)
!
!	FACE VELOCITY
	DOUBLE PRECISION V_F
!
!	DEFERRED CORRCTION CONTRIBUTION FORM HIGH ORDER METHOD
	DOUBLE PRECISION MOM_HO
!
!	LOW ORDER APPROXIMATION 
	DOUBLE PRECISION MOM_LO
!
!	CONVECTION FACTOR AT THE FACE
	DOUBLE PRECISION Flux
!
!	DEFERRED CORRECTION CONTRIBUTIONS FROM EACH FACE
	DOUBLE PRECISION 	EAST_DC
	DOUBLE PRECISION 	WEST_DC
	DOUBLE PRECISION 	NORTH_DC
	DOUBLE PRECISION 	SOUTH_DC
        DOUBLE PRECISION  TOP_DC
        DOUBLE PRECISION  BOTTOM_DC
!
! 
!-----------------------------------------------
!
!---------------------------------------------------------------
!	EXTERNAL FUNCTIONS
!---------------------------------------------------------------
	DOUBLE PRECISION , EXTERNAL :: FPFOI_OF
!---------------------------------------------------------------
!
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'function3.inc'


      call lock_tmp4_array
 
     
      call lock_tmp_array
      call lock_xsi_array

!
!  Calculate convection factors
!
!
! Send recv the third ghost layer
      IF ( FPFOI ) THEN
         Do IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = W_G(IJK)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF


!$omp parallel do private(IJK,K,IJKT,IJKP)
      DO IJK = ijkstart3, ijkend3
         K = K_OF(IJK) 
         IJKT = TOP_OF(IJK) 
         IJKP = KP_OF(IJK) 
!
!           East face (i+1/2, j, k+1/2)
         U(IJK) = AVG_Z(U_G(IJK),U_G(IJKP),K) 
!
!
!           North face (i, j+1/2, k+1/2)
         V(IJK) = AVG_Z(V_G(IJK),V_G(IJKP),K) 
!
!
!           Top face (i, j, k+1)
         WW(IJK) = AVG_Z_T(W_G(IJK),W_G(IJKP)) 
      END DO 

! loezos
	incr=0		
! loezos

      CALL CALC_XSI (DISCRETIZE(5), W_G, U, V, WW, XSI_E, XSI_N,&
	 XSI_T,incr) 
!
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!

!$omp      parallel do                                               &
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,     &
!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK,  D_f,        &
!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,                     &
!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,                     &
!$omp&             IJKM, KM, IJKB, &
!$omp&              MOM_HO, MOM_LO, CONV_FAC,EAST_DC,WEST_DC,NORTH_DC,&
!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3 
!
         IF (FLOW_AT_T(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKN = NORTH_OF(IJK) 
            IJKT = TOP_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKT 
            ELSE 
               IJKC = IJK 
            ENDIF 
            KP = KP1(K) 
            IJKE = EAST_OF(IJK)
            IJKP = KP_OF(IJK) 
            IJKTN = NORTH_OF(IJKT) 
            IJKTE = EAST_OF(IJKT) 
!
!           Third Ghost layer information
            IPPP  = IP_OF(IP_OF(IPJK))
            IPPP4 = funijk3(I_OF(IPPP), J_OF(IPPP), K_OF(IPPP))
            IMMM  = IM_OF(IM_OF(IMJK))
            IMMM4 = funijk3(I_OF(IMMM), J_OF(IMMM), K_OF(IMMM))
            JPPP  = JP_OF(JP_OF(IJPK))
            JPPP4 = funijk3(I_OF(JPPP), J_OF(JPPP), K_OF(JPPP))
            JMMM  = JM_OF(JM_OF(IJMK))
            JMMM4 = funijk3(I_OF(JMMM), J_OF(JMMM), K_OF(JMMM))
            KPPP  = KP_OF(KP_OF(IJKP))
            KPPP4 = funijk3(K_OF(IPPP), J_OF(KPPP), K_OF(KPPP))
            KMMM  = KM_OF(KM_OF(IJKM))
            KMMM4 = funijk3(I_OF(KMMM), J_OF(KMMM), K_OF(KMMM))
!
!

!
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1/2, j, k+1/2)
!            
	    IF(U(IJK) >= ZERO)THEN
	      MOM_LO = W_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IPJK), W_G(IJK), & 
                            W_G(IMJK), W_G(IM_OF(IMJK)))
	    ELSE
	      MOM_LO = W_G(IPJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJK), W_G(IPJK), & 
                            W_G(IP_OF(IPJK)), TMP4(IPPP4))
	    ENDIF
            IF (.NOT. FPFOI ) &
	              MOM_HO = XSI_E(IJK)*W_G(IPJK)+ &
                               (1.0-XSI_E(IJK))*W_G(IJK)
	    Flux = HALF * (Flux_gE(IJK) + Flux_gE(IJKP))
	    EAST_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1/2, k+1/2)
!            
	    IF(V(IJK) >= ZERO)THEN
	      MOM_LO = W_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJPK), W_G(IJK), & 
                            W_G(IJMK), W_G(JM_OF(IJMK)))
	    ELSE
	      MOM_LO = W_G(IJPK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJK), W_G(IJPK), & 
                            W_G(JP_OF(IJPK)), TMP4(JPPP4))
	    ENDIF
            IF (.NOT. FPFOI ) &
	               MOM_HO = XSI_N(IJK)*W_G(IJPK)+ &
                                (1.0-XSI_N(IJK))*W_G(IJK)
	    Flux = HALF * (Flux_gN(IJK) + Flux_gN(IJKP))
	    NORTH_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j, k+1)
!
	    IF(WW(IJK) >= ZERO)THEN
	      MOM_LO = W_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJKP), W_G(IJK), & 
                            W_G(IJKM), W_G(KM_OF(IJKM)))
	    ELSE
	      MOM_LO = W_G(IJKP)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJK), W_G(IJKP), & 
                            W_G(KP_OF(IJKP)), TMP4(KPPP4))
	    ENDIF
            IF (.NOT. FPFOI ) &
	               MOM_HO = XSI_T(IJK)*W_G(IJKP)+ &
                                (1.0-XSI_T(IJK))*W_G(IJK)
	    Flux = HALF * (Flux_gT(IJK) + Flux_gT(IJKP))
	    TOP_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE West face (i-1/2, j, k+1/2)
!
            IMJK = IM_OF(IJK) 
            IM = IM1(I) 
            IJKW = WEST_OF(IJK) 
            IJKWT = TOP_OF(IJKW) 
            IMJKP = KP_OF(IMJK) 
	    IF(U(IMJK) >= ZERO)THEN
	      MOM_LO = W_G(IMJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJK), W_G(IMJK), & 
                            W_G(IM_OF(IMJK)), TMP4(IMMM4))
	    ELSE
	      MOM_LO = W_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IMJK), W_G(IJK), & 
                            W_G(IPJK), W_G(IP_OF(IPJK)))
	    ENDIF
            IF (.NOT. FPFOI ) &
	               MOM_HO = XSI_E(IMJK)*W_G(IJK)+ &
                                (1.0-XSI_E(IMJK))*W_G(IMJK)
	    Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IMJKP))
	    WEST_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i, j-1/2, k+1/2)
!
            IJMK = JM_OF(IJK) 
            JM = JM1(J) 
            IJMKP = KP_OF(IJMK) 
            IJKS = SOUTH_OF(IJK) 
            IJKST = TOP_OF(IJKS) 
            IF(V(IJMK) >= ZERO)THEN
	      MOM_LO = W_G(IJMK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJK), W_G(IJMK), & 
                            W_G(JM_OF(IJMK)), TMP4(JMMM4))
	    ELSE
	      MOM_LO = W_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJMK), W_G(IJK), & 
                            W_G(IJPK), W_G(JP_OF(IJPK)))
	    ENDIF
            IF (.NOT. FPFOI ) &
	                MOM_HO = XSI_N(IJMK)*W_G(IJK)+ &
                                 (1.0-XSI_N(IJMK))*W_G(IJMK)
	    Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IJMKP))
	    SOUTH_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i, j, k)
!
            IJKM = KM_OF(IJK) 
            KM = KM1(K) 
            IJKB = BOTTOM_OF(IJK) 
	    IF(WW(IJK) >= ZERO)THEN
	      MOM_LO = W_G(IJKM)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJK), W_G(IJKM), & 
                            W_G(KM_OF(IJKM)), TMP4(KMMM4))
	    ELSE
	      MOM_LO = W_G(IJK)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_G(IJKM), W_G(IJK), & 
                            W_G(IJKP), W_G(KP_OF(IJKP)))
	    ENDIF
            IF (.NOT. FPFOI ) &
	               MOM_HO = XSI_T(IJKM)*W_G(IJK)+ &
                                (1.0-XSI_T(IJKM))*W_G(IJKM)
	    Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IJK))
	    BOTTOM_DC = Flux*(MOM_LO-MOM_HO)
!
!		CONTRIBUTION DUE TO DEFERRED CORRECTION
!
            B_M(IJK) = B_M(IJK)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC&
				+BOTTOM_DC-TOP_DC
! 
         ENDIF 
      END DO 
      
      call unlock_tmp4_array
      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_W_GDC 


!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_W_g1(A_W_g, IER)                               C
!  Purpose: Determine convection diffusion terms for W_g momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative. Higher order C
!  discretization.                                                     C
!  See source_w_g                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date:19-MAR-97   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_W_G1(A_W_G, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE visc_g
      USE toleranc 
      USE physprop
      USE fldvar
      USE output
      USE vshear
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar
      USE mflux   
      
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,& 
                       IJKTE, IJKP, IJKT, IJKTN, IJK 
      INTEGER          IMJK, IM, IJKW, IJKWT, IMJKP 
      INTEGER          IJMK, JM, IJMKP, IJKS, IJKST 
      INTEGER          IJKM, KM, IJKB 
!
! loezos
	INTEGER incr
! loezos     
 
!                      Face mass flux 
      DOUBLE PRECISION Flux 
 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_W_g 
      DOUBLE PRECISION A_W_g(DIMENSION_3, -3:3) 
! 
!                      Convection weighting factors 
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),& 
!                       XSI_t(DIMENSION_3) 
!      DOUBLE PRECISION U(DIMENSION_3),& 
!                       V(DIMENSION_3), WW(DIMENSION_3) 
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'


      call lock_tmp_array
      call lock_xsi_array


!
!  Calculate convection factors
!

!$omp parallel do private(IJK,K,IJKT,IJKP)
      DO IJK = ijkstart3, ijkend3
         K = K_OF(IJK) 
         IJKT = TOP_OF(IJK) 
         IJKP = KP_OF(IJK) 
!
!
!           East face (i+1/2, j, k+1/2)
         U(IJK) = AVG_Z(U_G(IJK),U_G(IJKP),K) 
!
!
!           North face (i, j+1/2, k+1/2)
         V(IJK) = AVG_Z(V_G(IJK),V_G(IJKP),K) 
!
!
!           Top face (i, j, k+1)
         WW(IJK) = AVG_Z_T(W_G(IJK),W_G(IJKP)) 
      END DO 

! loezos
	incr=0
! loezos

      CALL CALC_XSI (DISCRETIZE(5), W_G, U, V, WW, XSI_E, XSI_N,&
               XSI_T,incr) 

! loezos    
! update to true velocity
      IF (SHEAR) THEN
!$omp parallel do private(IJK)  
	 DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
	   V(IJK)=V(IJK)+VSH(IJK)	
          END IF
        END DO
      END IF
! loezos

!
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!

!!$omp      parallel do                                               &
!!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,     &
!!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK,  D_f,        &
!!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,                     &
!!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,                     &
!!$omp&             IJKM, KM, IJKB)
      DO IJK = ijkstart3, ijkend3 
!
         IF (FLOW_AT_T(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            IJKT = TOP_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKT 
            ELSE 
               IJKC = IJK 
            ENDIF 
            KP = KP1(K) 
            IJKE = EAST_OF(IJK) 
            IJKP = KP_OF(IJK) 
            IJKTN = NORTH_OF(IJKT) 
            IJKTE = EAST_OF(IJKT) 
!
!           East face (i+1/2, j, k+1/2)
	    Flux = HALF * (Flux_gE(IJK) + Flux_gE(IJKP))
            D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),AVG_X_H(MU_GT(IJKT&
               ),MU_GT(IJKTE),I),K)*ODX_E(I)*AYZ_W(IJK) 
!
            A_W_G(IJK,E) = D_F - XSI_E(IJK)*Flux
!
            A_W_G(IPJK,W) = D_F + (ONE - XSI_E(IJK))*Flux
!
!           North face (i, j+1/2, k+1/2)
	    Flux = HALF * (Flux_gN(IJK) + Flux_gN(IJKP))
            D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),AVG_Y_H(MU_GT(IJKT&
               ),MU_GT(IJKTN),J),K)*ODY_N(J)*AXZ_W(IJK) 
!
            A_W_G(IJK,N) = D_F - XSI_N(IJK)*Flux
!
            A_W_G(IJPK,S) = D_F + (ONE - XSI_N(IJK))*Flux
!
!           Top face (i, j, k+1)
	    Flux = HALF * (Flux_gT(IJK) + Flux_gT(IJKP))
            D_F = MU_GT(IJKT)*OX(I)*ODZ(KP)*AXY_W(IJK) 
            A_W_G(IJK,T) = D_F - XSI_T(IJK)*Flux
            A_W_G(IJKP,B) = D_F + (ONE - XSI_T(IJK))*Flux
!
!           West face (i-1/2, j, k+1/2)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
               IJKWT = TOP_OF(IJKW) 
               IMJKP = KP_OF(IMJK) 
	       Flux = HALF * (Flux_gE(IMJK) + Flux_gE(IMJKP))
               D_F = AVG_Z_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),AVG_X_H(MU_GT(&
                  IJKWT),MU_GT(IJKT),IM),K)*ODX_E(IM)*AYZ_W(IMJK) 
               A_W_G(IJK,W) = D_F + (ONE - XSI_E(IMJK))*Flux
            ENDIF 
!
!           South face (i, j-1/2, k+1/2)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IJMK)) THEN 
               JM = JM1(J) 
               IJMKP = KP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKST = TOP_OF(IJKS) 
	       Flux = HALF * (Flux_gN(IJMK) + Flux_gN(IJMKP))
               D_F = AVG_Z_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),AVG_Y_H(MU_GT(&
                  IJKST),MU_GT(IJKT),JM),K)*ODY_N(JM)*AXZ_W(IJMK) 
               A_W_G(IJK,S) = D_F + (ONE - XSI_N(IJMK))*Flux
            ENDIF 
!
!           Bottom face (i, j, k)
            IJKM = KM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IJKM)) THEN 
               KM = KM1(K) 
               IJKB = BOTTOM_OF(IJK) 
	       Flux = HALF * (Flux_gT(IJKM) + Flux_gT(IJK))
               D_F = MU_GT(IJK)*OX(I)*ODZ(K)*AXY_W(IJKM) 
               A_W_G(IJK,B) = D_F + (ONE - XSI_T(IJKM))*Flux
            ENDIF 
         ENDIF 
      END DO 
      
      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_W_G1 
      
!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3

                                                                                                                                                                                                                                                                                                                                                                                                                                conv_dif_w_s.f                                                                                      0100644 0002444 0000146 00000075112 10250070310 012021  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_W_s(A_m, B_m, IER)                            C
!  Purpose: Determine convection diffusion terms for W_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_W_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_DIF_W_S(A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE compar    
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Solids phase index 
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
!
      DO M = 1, MMAX 
        IF (MOMENTUM_Z_EQ(M)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	IF DEFERRED CORRECTION IS TO BE USED TO SOLVE W_S
!
          IF (DEF_COR) THEN
	    CALL STORE_A_W_S0 (A_M(1,-3,M), M, IER) 
	    IF (DISCRETIZE(5) > 1)CALL STORE_A_W_SDC (A_M(1,-3,M), M, B_M, IER)
          ELSE
!	   NO DEFERRED CORRECTION IS TO BE USED TO SOLVE FOR W_S
!
            IF (DISCRETIZE(5) == 0) THEN         ! 0 & 1 => FOUP 
              CALL STORE_A_W_S0 (A_M(1,-3,M), M, IER) 
            ELSE 
              CALL STORE_A_W_S1 (A_M(1,-3,M), M, IER) 
            ENDIF 
          ENDIF
!

          CALL DIF_W_IS (MU_S(1,M), A_M, B_M, M, IER) 
!
        ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE CONV_DIF_W_S 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_W_s0(A_W_s, M, IER)                            C
!  Purpose: Determine convection diffusion terms for W_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative. FOUP         C
!  See source_w_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-JUN-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_W_S0(A_W_S, M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc 
      USE fldvar
      USE output
      USE compar 
      USE mflux  
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,& 
                       IJKTE, IJKP, IJKT, IJKTN, IJK 
      INTEGER          IMJK, IM, IJKW, IJKWT, IMJKP 
      INTEGER          IJMK, JM, IJMKP, IJKS, IJKST 
      INTEGER          IJKM, KM, IJKB 
! 
!                      Solids phase 
      INTEGER          M 
! 
!                      Face mass flux 
      DOUBLE PRECISION Flux
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_W_s 
      DOUBLE PRECISION A_W_s(DIMENSION_3, -3:3, M:M) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'

!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do 	&
!$omp&     private(IJK,  I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,	&
!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, V_f, D_f,	&
!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,	&
!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,	&
!$omp&             IJKM, KM, IJKB)
      DO IJK = ijkstart3, ijkend3 
!
         IF (FLOW_AT_T(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            IJKT = TOP_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKT 
            ELSE 
               IJKC = IJK 
            ENDIF 
            KP = KP1(K) 
            IJKE = EAST_OF(IJK) 
            IJKP = KP_OF(IJK) 
            IJKTN = NORTH_OF(IJKT) 
            IJKTE = EAST_OF(IJKT) 
!
!           East face (i+1/2, j, k+1/2)
	    Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IJKP,M))
            D_F = AVG_Z_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),AVG_X_H(MU_S(&
               IJKT,M),MU_S(IJKTE,M),I),K)*ODX_E(I)*AYZ_W(IJK) 
            IF (Flux >= ZERO) THEN 
               A_W_S(IJK,E,M) = D_F 
               A_W_S(IPJK,W,M) = D_F + Flux
            ELSE 
               A_W_S(IJK,E,M) = D_F - Flux
               A_W_S(IPJK,W,M) = D_F 
            ENDIF 
!
!           North face (i, j+1/2, k+1/2)
	    Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IJKP,M))
            D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(&
               IJKT,M),MU_S(IJKTN,M),J),K)*ODY_N(J)*AXZ_W(IJK) 
            IF (Flux >= ZERO) THEN 
               A_W_S(IJK,N,M) = D_F 
               A_W_S(IJPK,S,M) = D_F + Flux
            ELSE 
               A_W_S(IJK,N,M) = D_F - Flux
               A_W_S(IJPK,S,M) = D_F 
            ENDIF 
!
!           Top face (i, j, k+1)
	    Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJKP,M))
            D_F = MU_S(IJKT,M)*OX(I)*ODZ(KP)*AXY_W(IJK) 
            IF (Flux >= ZERO) THEN 
               A_W_S(IJK,T,M) = D_F 
               A_W_S(IJKP,B,M) = D_F + Flux
            ELSE 
               A_W_S(IJK,T,M) = D_F - Flux
               A_W_S(IJKP,B,M) = D_F 
            ENDIF 
!
!           West face (i-1/2, j, k+1/2)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
               IJKWT = TOP_OF(IJKW) 
               IMJKP = KP_OF(IMJK) 
	       Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IMJKP,M))
               D_F = AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),AVG_X_H(MU_S&
                  (IJKWT,M),MU_S(IJKT,M),IM),K)*ODX_E(IM)*AYZ_W(IMJK) 
               IF (Flux >= ZERO) THEN 
                  A_W_S(IJK,W,M) = D_F + Flux
               ELSE 
                  A_W_S(IJK,W,M) = D_F 
               ENDIF 
            ENDIF 
!
!           South face (i, j-1/2, k+1/2)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IJMK)) THEN 
               JM = JM1(J) 
               IJMKP = KP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKST = TOP_OF(IJKS) 
	       Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IJMKP,M))
               D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),AVG_Y_H(MU_S&
                  (IJKST,M),MU_S(IJKT,M),JM),K)*ODY_N(JM)*AXZ_W(IJMK) 
               IF (Flux >= ZERO) THEN 
                  A_W_S(IJK,S,M) = D_F + Flux
               ELSE 
                  A_W_S(IJK,S,M) = D_F 
               ENDIF 
            ENDIF 
!
!           Bottom face (i, j, k)
            IJKM = KM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IJKM)) THEN 
               KM = KM1(K) 
               IJKB = BOTTOM_OF(IJK) 
	       Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IJK,M))
               D_F = MU_S(IJK,M)*OX(I)*ODZ(K)*AXY_W(IJKM) 
               IF (Flux >= ZERO) THEN 
                  A_W_S(IJK,B,M) = D_F + Flux
               ELSE 
                  A_W_S(IJK,B,M) = D_F 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 
      
      RETURN  
      END SUBROUTINE STORE_A_W_S0

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_W_sdc(A_W_s, M, B_M, IER)                      C
!  Purpose: TO USE DEFERRED CORRECTION METHOD TO SOLVE THE U-MOMENTUM  C
!  EQUATION. THIS METHOD COMBINES FIRST ORDER UPWIND AND A USER        C
!  SPECIFIED HIGH ORDER METHOD                                         C
!                                                                      C
!  Author: C. GUENTHER                                 Date:8-APR-99   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_W_SDC(A_W_S, M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc 
      USE fldvar
      USE output
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar      
      USE sendrecv
      USE sendrecv3
      USE mflux
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,& 
                       IJKTE, IJKP, IJKT, IJKTN, IJK 
      INTEGER          IMJK, IM, IJKW, IJKWT, IMJKP 
      INTEGER          IJMK, JM, IJMKP, IJKS, IJKST 
      INTEGER          IJKM, KM, IJKB 
      INTEGER          IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER          IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
! 
!                      Solids phase 
      INTEGER          M 
!
! loezos
	INTEGER incr
! loezos
 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_W_s 
      DOUBLE PRECISION A_W_s(DIMENSION_3, -3:3, M:M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
! 
!
!	FACE VELOCITY
	DOUBLE PRECISION V_F
!
!	DEFERRED CORRCTION CONTRIBUTION FORM HIGH ORDER METHOD
	DOUBLE PRECISION MOM_HO
!
!	LOW ORDER APPROXIMATION 
	DOUBLE PRECISION MOM_LO
!
!	CONVECTION FACTOR AT THE FACE
	DOUBLE PRECISION Flux
!
!	DEFERRED CORRECTION CONTRIBUTIONS FROM EACH FACE
	DOUBLE PRECISION 	EAST_DC
	DOUBLE PRECISION 	WEST_DC
	DOUBLE PRECISION 	NORTH_DC
	DOUBLE PRECISION 	SOUTH_DC
        DOUBLE PRECISION  TOP_DC
        DOUBLE PRECISION  BOTTOM_DC
!
! 
!-----------------------------------------------
!
!---------------------------------------------------------------
!	EXTERNAL FUNCTIONS
!---------------------------------------------------------------
	DOUBLE PRECISION , EXTERNAL :: FPFOI_OF
!---------------------------------------------------------------
! 
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'function3.inc'


      call lock_tmp4_array



      call lock_tmp_array
      call lock_xsi_array

!
!  Calculate convection factors
!
!
! Send recv the third ghost layer
      IF ( FPFOI ) THEN
         Do IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = W_S(IJK,M)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF


!$omp parallel do private(IJK,K,IJKT,IJKP )
      DO IJK = ijkstart3, ijkend3
         K = K_OF(IJK) 
         IJKT = TOP_OF(IJK) 
         IJKP = KP_OF(IJK) 
!
!
!           East face (i+1/2, j, k+1/2)
         U(IJK) = AVG_Z(U_S(IJK,M),U_S(IJKP,M),K) 
!
!
!           North face (i, j+1/2, k+1/2)
         V(IJK) = AVG_Z(V_S(IJK,M),V_S(IJKP,M),K) 
!
!
!           Top face (i, j, k+1)
         WW(IJK) = AVG_Z_T(W_S(IJK,M),W_S(IJKP,M)) 
      END DO 

! loezos
	 incr=0
! loezos

      CALL CALC_XSI (DISCRETIZE(5), W_S(1,M), U, V, WW, XSI_E, XSI_N,& 
	XSI_T,incr) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!

!$omp      parallel do 	&
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,	&
!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK,  D_f,	&
!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,	&
!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,	&
!$omp&             IJKM, KM, IJKB, &
!$omp&              MOM_HO, MOM_LO, CONV_FAC,EAST_DC,WEST_DC,NORTH_DC,&
!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3 
!
         IF (FLOW_AT_T(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKN = NORTH_OF(IJK) 
            IJKT = TOP_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKT 
            ELSE 
               IJKC = IJK 
            ENDIF 
            KP = KP1(K) 
            IJKE = EAST_OF(IJK)
            IJKTN = NORTH_OF(IJKT) 
            IJKTE = EAST_OF(IJKT) 
!
!           Third Ghost layer information
            IPPP  = IP_OF(IP_OF(IPJK))
            IPPP4 = funijk3(I_OF(IPPP), J_OF(IPPP), K_OF(IPPP))
            IMMM  = IM_OF(IM_OF(IMJK))
            IMMM4 = funijk3(I_OF(IMMM), J_OF(IMMM), K_OF(IMMM))
            JPPP  = JP_OF(JP_OF(IJPK))
            JPPP4 = funijk3(I_OF(JPPP), J_OF(JPPP), K_OF(JPPP))
            JMMM  = JM_OF(JM_OF(IJMK))
            JMMM4 = funijk3(I_OF(JMMM), J_OF(JMMM), K_OF(JMMM))
            KPPP  = KP_OF(KP_OF(IJKP))
            KPPP4 = funijk3(K_OF(IPPP), J_OF(KPPP), K_OF(KPPP))
            KMMM  = KM_OF(KM_OF(IJKM))
            KMMM4 = funijk3(I_OF(KMMM), J_OF(KMMM), K_OF(KMMM))
!
!
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1/2, j, k+1/2)
!           
		IF(U(IJK) >= ZERO)THEN
		    MOM_LO = W_S(IJK,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_S(IPJK,M), W_S(IJK,M), & 
                            W_S(IMJK,M), W_S(IM_OF(IMJK),M))
		ELSE
		    MOM_LO = W_S(IPJK,M)
                     IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_S(IJK,M), W_S(IPJK,M), & 
                            W_S(IP_OF(IPJK),M), TMP4(IPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		       MOM_HO = XSI_E(IJK)*W_S(IPJK,M)+ &
                               (1.0-XSI_E(IJK))*W_S(IJK,M)
	        Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IJKP,M))
		EAST_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i, j+1/2, k+1/2)
!            
		IF(V(IJK) >= ZERO)THEN
		    MOM_LO = W_S(IJK,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_S(IJPK,M), W_S(IJK,M), & 
                            W_S(IJMK,M), W_S(JM_OF(IJMK),M))
		ELSE
		    MOM_LO = W_S(IJPK,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_S(IJK,M), W_S(IJPK,M), & 
                            W_S(JP_OF(IJPK),M), TMP4(JPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		       MOM_HO = XSI_N(IJK)*W_S(IJPK,M)+ &
                                (1.0-XSI_N(IJK))*W_S(IJK,M)
	        Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IJKP,M))
		NORTH_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i, j, k+1)
!            
		IF(WW(IJK) >= ZERO)THEN
		    MOM_LO = W_S(IJK,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_S(IJKP,M), W_S(IJK,M), & 
                            W_S(IJKM,M), W_S(KM_OF(IJKM),M))
		ELSE
		    MOM_LO = W_S(IJKP,M)
                    IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_S(IJK,M), W_S(IJKP,M), & 
                            W_S(KP_OF(IJKP),M), TMP4(KPPP4))
		ENDIF
                IF (.NOT. FPFOI ) &
		       MOM_HO = XSI_T(IJK)*W_S(IJKP,M)+ &
                               (1.0-XSI_T(IJK))*W_S(IJK,M)
	        Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJKP,M))
		TOP_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE West face (i-1/2, j, k+1/2)
!
            IMJK = IM_OF(IJK) 
            IM = IM1(I) 
            IJKW = WEST_OF(IJK) 
            IJKWT = TOP_OF(IJKW) 
            IMJKP = KP_OF(IMJK) 
            IF(U(IMJK) >= ZERO)THEN
	      MOM_LO = W_S(IMJK,M)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_S(IJK,M), W_S(IMJK,M), & 
                            W_S(IM_OF(IMJK),M), TMP4(IMMM4))
	    ELSE
	      MOM_LO = W_S(IJK,M)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_S(IMJK,M), W_S(IJK,M), & 
                            W_S(IPJK,M), W_S(IP_OF(IPJK),M))
	    ENDIF
            IF (.NOT. FPFOI ) &
	               MOM_HO = XSI_E(IMJK)*W_S(IJK,M)+ &
                                (1.0-XSI_E(IMJK))*W_S(IMJK,M)
	    Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IMJKP,M))
	    WEST_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i, j-1/2, k+1/2)
!
            IJMK = JM_OF(IJK) 
            JM = JM1(J) 
            IJMKP = KP_OF(IJMK) 
            IJKS = SOUTH_OF(IJK) 
            IJKST = TOP_OF(IJKS) 
	    IF(V(IJMK) >= ZERO)THEN
	      MOM_LO = W_S(IJMK,M)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_S(IJK,M), W_S(IJMK,M), & 
                            W_S(JM_OF(IJMK),M), TMP4(JMMM4))
	    ELSE
	      MOM_LO = W_S(IJK,M)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_S(IJMK,M), W_S(IJK,M), & 
                            W_S(IJPK,M), W_S(JP_OF(IJPK),M))
	    ENDIF
            IF (.NOT. FPFOI ) &
	               MOM_HO = XSI_N(IJMK)*W_S(IJK,M)+ &
                                (1.0-XSI_N(IJMK))*W_S(IJMK,M)
	    Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IJMKP,M))
	    SOUTH_DC = Flux*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i, j, k)
!
            IJKM = KM_OF(IJK) 
            KM = KM1(K) 
            IJKB = BOTTOM_OF(IJK) 
            IF(WW(IJK) >= ZERO)THEN
	      MOM_LO = W_S(IJKM,M)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_S(IJK,M), W_S(IJKM,M), & 
                            W_S(KM_OF(IJKM),M), TMP4(KMMM4))
	    ELSE
	      MOM_LO = W_S(IJK,M)
              IF ( FPFOI ) &
                      MOM_HO = FPFOI_OF(W_S(IJKM,M), W_S(IJK,M), & 
                            W_S(IJKP,M), W_S(KP_OF(IJKP),M))
	    ENDIF
            IF (.NOT. FPFOI ) &
	               MOM_HO = XSI_T(IJKM)*W_S(IJK,M)+ &
                                (1.0-XSI_T(IJKM))*W_S(IJKM,M)
	    Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IJK,M))
	    BOTTOM_DC = Flux*(MOM_LO-MOM_HO)
!
!
!		CONTRIBUTION DUE TO DEFERRED CORRECTION
!
            B_M(IJK,M) = B_M(IJK,M)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC&
				+BOTTOM_DC-TOP_DC
! 
         ENDIF 
      END DO 

      call unlock_tmp4_array
      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_W_SDC 

 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_W_s1(A_W_s, M, IER)                            C
!  Purpose: Determine convection diffusion terms for W_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative. Higher order C
!  discretization.                                                     C
!  See source_w_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date:19-MAR-97   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE STORE_A_W_S1(A_W_S, M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc 
      USE fldvar
      USE output
      USE vshear
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar 
      USE mflux  
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I,  J, K, IPJK, IJPK, IJKN, IJKC, KP, IJKE,& 
                       IJKTE, IJKP, IJKT, IJKTN, IJK 
      INTEGER          IMJK, IM, IJKW, IJKWT, IMJKP 
      INTEGER          IJMK, JM, IJMKP, IJKS, IJKST 
      INTEGER          IJKM, KM, IJKB 
! 
!                      Solids phase 
      INTEGER          M 
!
! loezos
	INTEGER incr
! loezos
 
!                      Face mass flux 
      DOUBLE PRECISION Flux 
 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_W_s 
      DOUBLE PRECISION A_W_s(DIMENSION_3, -3:3, M:M) 
! 
!                      Convection weighting factors 
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),& 
!                       XSI_t(DIMENSION_3) 
!      DOUBLE PRECISION U(DIMENSION_3),& 
!                       V(DIMENSION_3), WW(DIMENSION_3) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'



      call lock_tmp_array
      call lock_xsi_array
      
!
!  Calculate convection factors
!

!$omp parallel do private(IJK,K,IJKT,IJKP )
      DO IJK = ijkstart3, ijkend3 
         K = K_OF(IJK) 
         IJKT = TOP_OF(IJK) 
         IJKP = KP_OF(IJK) 
!
!
!           East face (i+1/2, j, k+1/2)
         U(IJK) = AVG_Z(U_S(IJK,M),U_S(IJKP,M),K) 
!
!
!           North face (i, j+1/2, k+1/2)
         V(IJK) = AVG_Z(V_S(IJK,M),V_S(IJKP,M),K) 
!
!
!           Top face (i, j, k+1)
         WW(IJK) = AVG_Z_T(W_S(IJK,M),W_S(IJKP,M)) 
      END DO 

! loezos
	 incr=0
! loezos

      CALL CALC_XSI (DISCRETIZE(5), W_S(1,M), U, V, WW, XSI_E, XSI_N,&
	 XSI_T,incr) 


! loezos      
! ! update to true velocity
      IF (SHEAR) THEN
!$omp      parallel do private(IJK)
	 DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
	   V(IJK)=V(IJK)+VSH(IJK)	
          END IF
        END DO

      END IF
! loezos

!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do 	&
!$omp&     private( I,  J, K, IPJK, IJPK, IJKN, IJKC, KP,	&
!$omp&             IJKE, IJKTE, IJKP, IJKT, IJKTN, IJK,  D_f,	&
!$omp&             IMJK, IM, IJKW, IJKWT, IMJKP,	&
!$omp&             IJMK, JM, IJMKP, IJKS, IJKST,	&
!$omp&             IJKM, KM, IJKB)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_T(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            IJKT = TOP_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKT 
            ELSE 
               IJKC = IJK 
            ENDIF 
            KP = KP1(K) 
            IJKE = EAST_OF(IJK) 
            IJKP = KP_OF(IJK) 
            IJKTN = NORTH_OF(IJKT) 
            IJKTE = EAST_OF(IJKT) 
!
!           East face (i+1/2, j, k+1/2)
	    Flux = HALF * (Flux_sE(IJK,M) + Flux_sE(IJKP,M))
            D_F = AVG_Z_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),AVG_X_H(MU_S(&
               IJKT,M),MU_S(IJKTE,M),I),K)*ODX_E(I)*AYZ_W(IJK) 
!
            A_W_S(IJK,E,M) = D_F - XSI_E(IJK)*Flux
!
            A_W_S(IPJK,W,M) = D_F + (ONE - XSI_E(IJK))*Flux
!
!           North face (i, j+1/2, k+1/2)
	    Flux = HALF * (Flux_sN(IJK,M) + Flux_sN(IJKP,M))
            D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(&
               IJKT,M),MU_S(IJKTN,M),J),K)*ODY_N(J)*AXZ_W(IJK) 
!
            A_W_S(IJK,N,M) = D_F - XSI_N(IJK)*Flux
!
            A_W_S(IJPK,S,M) = D_F + (ONE - XSI_N(IJK))*Flux
!
!           Top face (i, j, k+1)
	    Flux = HALF * (Flux_sT(IJK,M) + Flux_sT(IJKP,M))
            D_F = MU_S(IJKT,M)*OX(I)*ODZ(KP)*AXY_W(IJK) 
!
            A_W_S(IJK,T,M) = D_F - XSI_T(IJK)*Flux 
!
            A_W_S(IJKP,B,M) = D_F + (ONE - XSI_T(IJK))*Flux
!
!           West face (i-1/2, j, k+1/2)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK)
               IJKWT = TOP_OF(IJKW) 
               IMJKP = KP_OF(IMJK) 
	       Flux = HALF * (Flux_sE(IMJK,M) + Flux_sE(IMJKP,M))
               D_F = AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),AVG_X_H(MU_S&
                  (IJKWT,M),MU_S(IJKT,M),IM),K)*ODX_E(IM)*AYZ_W(IMJK) 
               A_W_S(IJK,W,M) = D_F + (ONE - XSI_E(IMJK))*Flux
            ENDIF 
!
!           South face (i, j-1/2, k+1/2)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IJMK)) THEN 
               JM = JM1(J) 
               IJMKP = KP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKST = TOP_OF(IJKS) 
	       Flux = HALF * (Flux_sN(IJMK,M) + Flux_sN(IJMKP,M))
               D_F = AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),AVG_Y_H(MU_S&
                  (IJKST,M),MU_S(IJKT,M),JM),K)*ODY_N(JM)*AXZ_W(IJMK) 
               A_W_S(IJK,S,M) = D_F + (ONE - XSI_N(IJMK))*Flux
            ENDIF 
!
!           Bottom face (i, j, k)
            IJKM = KM_OF(IJK) 
            IF (.NOT.FLOW_AT_T(IJKM)) THEN 
               KM = KM1(K) 
               IJKB = BOTTOM_OF(IJK) 
	       Flux = HALF * (Flux_sT(IJKM,M) + Flux_sT(IJK,M))
               D_F = MU_S(IJK,M)*OX(I)*ODZ(K)*AXY_W(IJKM) 
               A_W_S(IJK,B,M) = D_F + (ONE - XSI_T(IJKM))*Flux
            ENDIF 
         ENDIF 
      END DO 


      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_W_S1 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
                                                                                                                                                                                                                                                                                                                                                                                                                                                      conv_pp_g.f                                                                                         0100644 0002444 0000146 00000014466 10247636227 011366  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_Pp_g(A_m, B_m, IER)
!  Purpose: Determine convection terms for Pressure                    C
!  correction equation.  The off-diagonal coefficients calculated here C
!  must be positive. The center coefficient and the source vector are  C
!  negative. Multiplication with factors d_e, d_n, and d_t are carried C
!  out in source_pp_g.  Constant pressure boundaries are handled by
!  holding the Pp_g at the boundaries zero.  For specified mass flow   C
!  boundaries (part of) a's are calculated here since b is calculated  C
!  from a's in source_pp_g.  After calculating b, a's are multiplied by
!  d and at the flow boundaries get are set to zero.                   C                                !  Author: M. Syamlal                                 Date: 20-JUN-96  C
!  Reviewer:                                          Date:            C
!   Revision to use face densities calculated in CONV_ROP              C
!  Author: M. Syamlal                                 Date: 1-JUN-05  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_PP_G(A_M, B_M, IER) 
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE run 
      USE parallel 
      USE matrix 
      USE physprop
      USE geometry
      USE indices
      USE pgcor
      USE compar
      USE mflux   
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Indices 
      INTEGER          IJK, IPJK, IJPK, IJKP 
      INTEGER          IMJK, IJMK, IJKM 
      INTEGER          M 
! 
!                      local value of A_m 
      DOUBLE PRECISION am 
!-----------------------------------------------
      INCLUDE 'function.inc'
      

!
!  Calculate convection fluxes through each of the faces

!$omp  parallel do private( IJK, IPJK, IJPK, M, AM, IMJK, IJMK, IJKM) &
!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLUID_AT(IJK)) THEN 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKP = KP_OF(IJK) 
!
!         East face (i+1/2, j, k)
            AM = ROP_GE(IJK)*AYZ(IJK) 
            A_M(IJK,E,0) = AM 
            A_M(IPJK,W,0) = AM 
!
!         North face (i, j+1/2, k)
            AM = ROP_GN(IJK)*AXZ(IJK) 
            A_M(IJK,N,0) = AM 
            A_M(IJPK,S,0) = AM 
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               AM = ROP_GT(IJK)*AXY(IJK) 
               A_M(IJK,T,0) = AM 
               A_M(IJKP,B,0) = AM 
            ENDIF 
!
!         West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               AM = ROP_GE(IMJK)*AYZ(IMJK) 
               A_M(IJK,W,0) = AM 
            ENDIF 
!
!         South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               AM = ROP_GN(IJMK)*AXZ(IJMK) 
               A_M(IJK,S,0) = AM 
            ENDIF 
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  AM = ROP_GT(IJKM)*AXY(IJKM) 
                  A_M(IJK,B,0) = AM 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 
      DO M = 1, MMAX 
         IF (.NOT.CLOSE_PACKED(M)) THEN
            DO IJK = ijkstart3, ijkend3
               IF (FLUID_AT(IJK)) THEN 
                  IPJK = IP_OF(IJK) 
                  IJPK = JP_OF(IJK) 
                  IJKP = KP_OF(IJK) 
!
!             East face (i+1/2, j, k)
                  AM = ROP_SE(IJK,M)*AYZ(IJK) 
                  A_M(IJK,E,M) = AM 
                  A_M(IPJK,W,M) = AM 
!
!             North face (i, j+1/2, k)
                  AM = ROP_SN(IJK,M)*AXZ(IJK) 
                  A_M(IJK,N,M) = AM 
                  A_M(IJPK,S,M) = AM 
!
!             Top face (i, j, k+1/2)
                  IF (DO_K) THEN 
                     AM = ROP_ST(IJK,M)*AXY(IJK) 
                     A_M(IJK,T,M) = AM 
                     A_M(IJKP,B,M) = AM 
                  ENDIF 
!
!             West face (i-1/2, j, k)
                  IMJK = IM_OF(IJK) 
                  IF (.NOT.FLUID_AT(IMJK)) THEN 
                     AM = ROP_SE(IMJK,M)*AYZ(IMJK) 
                     A_M(IJK,W,M) = AM 
                  ENDIF 
!
!             South face (i, j-1/2, k)
                  IJMK = JM_OF(IJK) 
                  IF (.NOT.FLUID_AT(IJMK)) THEN 
                     AM = ROP_SN(IJMK,M)*AXZ(IJMK) 
                     A_M(IJK,S,M) = AM 
                  ENDIF 
!
!             Bottom face (i, j, k-1/2)
                  IF (DO_K) THEN 
                     IJKM = KM_OF(IJK) 
                     IF (.NOT.FLUID_AT(IJKM)) THEN 
                        AM = ROP_ST(IJKM,M)*AXY(IJKM) 
                        A_M(IJK,B,M) = AM 
                     ENDIF 
                  ENDIF 
               ENDIF 
            END DO 
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE CONV_PP_G 
                                                                                                                                                                                                          conv_rop.f                                                                                          0100644 0002444 0000146 00000030146 10247142224 011217  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_ROP(IER)                                          C
!  Purpose: Calculate the face value of density used for calculating   C
!           convection fluxes. Master routine.                         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_ROP(IER) 
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE mflux
      USE physprop
      USE run
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      solids phase index 
      INTEGER          M 
! 
! 
!
!
      IF (DISCRETIZE(1) == 0) THEN               ! 0 & 1 => first order upwinding 
         CALL CONV_ROP0 (ROP_g, U_g, V_g, W_g, ROP_gE, ROP_gN, ROP_gT, IER) 
      ELSE 
         CALL CONV_ROP1 (DISCRETIZE(1), ROP_g, U_g, V_g, W_g, ROP_gE, ROP_gN, ROP_gT, IER) 
      ENDIF 
      
      IF (DISCRETIZE(2) == 0) THEN               ! 0 & 1 => first order upwinding 
        DO M = 1, MMAX
          CALL CONV_ROP0 (ROP_s(1, M), U_s(1, M), V_s(1, M), W_s(1, M), &
	                  ROP_sE(1, M), ROP_sN(1, M), ROP_sT(1, M), IER) 
        ENDDO
      ELSE 
        DO M = 1, MMAX
          CALL CONV_ROP1 (DISCRETIZE(2), ROP_s(1, M), U_s(1, M), V_s(1, M), W_s(1, M), &
	                  ROP_sE(1, M), ROP_sN(1, M), ROP_sT(1, M), IER) 
        ENDDO
      ENDIF 
      
      RETURN  
      END SUBROUTINE CONV_ROP 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!  Module name: CONV_ROP0(ROP, U, V, W, ROP_E, ROP_N, ROP_T, IER)      C
!  Purpose: Calculate the face value of density used for calculating   C
!           convection fluxes. FOU routine.                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_ROP0(ROP, U, V, W, ROP_E, ROP_N, ROP_T, IER) 
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run 
      USE parallel 
      USE physprop
      USE geometry
      USE indices
      USE compar   
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
!
!                      macroscopic density (rho_prime)
      DOUBLE PRECISION ROP(DIMENSION_3) 
!
!                      Velocity components
      DOUBLE PRECISION U(DIMENSION_3), V(DIMENSION_3), W(DIMENSION_3) 
!
!                      Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION ROP_E(DIMENSION_3), ROP_N(DIMENSION_3), ROP_T(DIMENSION_3) 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, IJKE, IJKN, IJKT
      INTEGER          IJKW, IJKS, IJKB, IMJK, IJMK, IJKM 
!-----------------------------------------------
      INCLUDE 'function.inc'
      

!
!  Interpolate the face value of density for calculating the convection fluxes 
!$omp  parallel do private( IJK, IJKE, IJKN, IJKT, IJKW, IJKS, IJKB, IMJK, IJMK, IJKM) &
!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLUID_AT(IJK)) THEN 
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            IJKT = TOP_OF(IJK) 
!
!         East face (i+1/2, j, k)
            IF (U(IJK) >= ZERO) THEN 
               ROP_E(IJK) = ROP(IJK) 
            ELSE 
               ROP_E(IJK) = ROP(IJKE) 
            ENDIF 
!
!         North face (i, j+1/2, k)
            IF (V(IJK) >= ZERO) THEN 
               ROP_N(IJK) = ROP(IJK) 
            ELSE 
               ROP_N(IJK) = ROP(IJKN) 
            ENDIF 
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               IF (W(IJK) >= ZERO) THEN 
                  ROP_T(IJK) = ROP(IJK) 
               ELSE 
                  ROP_T(IJK) = ROP(IJKT) 
               ENDIF 
            ENDIF 
!
!         West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               IJKW = WEST_OF(IJK) 
               IF (U(IMJK) >= ZERO) THEN 
                  ROP_E(IMJK) = ROP(IJKW) 
               ELSE 
                  ROP_E(IMJK) = ROP(IJK) 
               ENDIF 
            ENDIF 
!
!         South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               IJKS = SOUTH_OF(IJK) 
               IF (V(IJMK) >= ZERO) THEN 
                 ROP_N(IJMK) = ROP(IJKS) 
               ELSE 
                 ROP_N(IJMK) = ROP(IJK) 
               ENDIF 
            ENDIF 
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  IJKB = BOTTOM_OF(IJK) 
                  IF (W(IJKM) >= ZERO) THEN 
                     ROP_T(IJKM) = ROP(IJKB) 
                  ELSE 
                     ROP_T(IJKM) = ROP(IJK) 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE CONV_ROP0 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!  Module name: CONV_ROP1(DISC, ROP, U, V, W, ROP_E, ROP_N, ROP_T, IER)C
!  Purpose: Calculate the face value of density used for calculating   C
!           convection fluxes. HR routine.                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!  Reviewer:                                          Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CONV_ROP1(DISC, ROP, U, V, W, ROP_E, ROP_N, ROP_T, IER) 
!
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run
      USE parallel 
      USE physprop
      USE geometry
      USE indices
      Use xsi_array
      USE compar    
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!
!                      Discretization scheme 
      INTEGER          DISC
!
!                      macroscopic density (rho_prime)
      DOUBLE PRECISION ROP(DIMENSION_3) 
!
!                      Velocity components
      DOUBLE PRECISION U(DIMENSION_3), V(DIMENSION_3), W(DIMENSION_3) 
!
!                      Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION ROP_E(DIMENSION_3), ROP_N(DIMENSION_3), ROP_T(DIMENSION_3) 
!
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, IJKE, IJKN, IJKT 
      INTEGER          IJKW, IJKS, IJKB, IMJK, IJMK, IJKM 
      
      Integer          incr
!-----------------------------------------------
      INCLUDE 'function.inc'
      

      call lock_xsi_array
!
!  Calculate convection factors
!
   
       incr=0	
       CALL CALC_XSI (DISC, ROP, U, V, W, XSI_E, XSI_N,	XSI_T,incr) 
!
!
!  Calculate convection fluxes through each of the faces
!
!$omp  parallel do private(IJK, IJKE, IJKN, IJKT, IJKW, IJKS, IJKB, IMJK, IJMK, IJKM) &
!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3 
!
         IF (FLUID_AT(IJK)) THEN 
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            IJKT = TOP_OF(IJK) 
!
!         East face (i+1/2, j, k)
            ROP_E(IJK) = ((ONE-XSI_E(IJK))*ROP(IJK)+XSI_E(IJK)*ROP(IJKE)) 
!
!         North face (i, j+1/2, k)
            ROP_N(IJK) = ((ONE-XSI_N(IJK))*ROP(IJK)+XSI_N(IJK)*ROP(IJKN)) 
!
!         Top face (i, j, k+1/2)
            IF (DO_K) THEN 
               ROP_T(IJK) = ((ONE - XSI_T(IJK))*ROP(IJK)+XSI_T(IJK)*ROP(IJKT))
            ENDIF 
!
!         West face (i-1/2, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLUID_AT(IMJK)) THEN 
               IJKW = WEST_OF(IJK) 
               ROP_E(IMJK) = ((ONE - XSI_E(IMJK))*ROP(IJKW)+XSI_E(IMJK)*ROP(IJK))
            ENDIF 
!
!         South face (i, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLUID_AT(IJMK)) THEN 
               IJKS = SOUTH_OF(IJK) 
               ROP_N(IJMK) = ((ONE - XSI_N(IJMK))*ROP(IJKS)+XSI_N(IJMK)*ROP(IJK))
            ENDIF 
!
!         Bottom face (i, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLUID_AT(IJKM)) THEN 
                  IJKB = BOTTOM_OF(IJK) 
                  ROP_T(IJKM) = ((ONE - XSI_T(IJKM))*ROP(IJKB)+XSI_T(IJKM)*ROP(IJK))
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 
      
      call unlock_xsi_array
 
           
      RETURN  
      END SUBROUTINE CONV_ROP1 
                                                                                                                                                                                                                                                                                                                                                                                                                          correct_0.f                                                                                         0100644 0002444 0000146 00000023422 10247644574 011270  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !
!
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CORRECT_0(IER)                                         C
!  Purpose: Correct the fluid pressure and gas and solids velocities   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CORRECT_0(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE fldvar
      USE pgcor
      USE ur_facs 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER IER 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!
!                      error index
!
      CALL CORRECT_0G (PP_G, UR_FAC(1), D_E, D_N, D_T, P_G, U_G, V_G, W_G, IER) 
!      CALL CORRECT_0S (PP_G, D_E, D_N, D_T, U_S, V_S, W_S, IER) 
      RETURN  
      END SUBROUTINE CORRECT_0 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CORRECT_0g(Pp_g, UR_fac, d_e, d_n, d_t,                C
!     &                P_g, U_g, V_g, W_g, IER)                        C
!  Purpose: Correct the fluid pressure and velocities.                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CORRECT_0G(PP_G,UR_FAC,D_E,D_N,D_T,P_G,U_G,V_G,W_G,IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE indices
      USE physprop
      USE compar 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Pressure correction 
      DOUBLE PRECISION Pp_g(DIMENSION_3) 
! 
!                      Under relaxation factor for Pressure correction 
      DOUBLE PRECISION UR_fac 
! 
!                      Pressure correction coefficient -- East 
      DOUBLE PRECISION d_e(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Pressure correction coefficient -- North 
      DOUBLE PRECISION d_n(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Pressure correction coefficient -- Top 
      DOUBLE PRECISION d_t(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Pressure 
      DOUBLE PRECISION P_g(DIMENSION_3) 
! 
!                      Velocity components 
      DOUBLE PRECISION U_g(DIMENSION_3), V_g(DIMENSION_3),& 
                       W_g(DIMENSION_3) 
! 
!                      error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, IJKE, IJKN, IJKT, M 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!  Underrelax pressure correction.  Velocity corrections should not be
!  underrelaxed, so that the continuity eq. is satisfied.
!
!$omp    parallel do private(IJK,IJKE,IJKN,IJKT)
      DO IJK = ijkstart3, ijkend3 
         IF (FLUIDORP_FLOW_AT(IJK)) THEN 
            P_G(IJK) = P_G(IJK) + UR_FAC*PP_G(IJK) 
!
            IJKE = EAST_OF(IJK) 
            IJKN = NORTH_OF(IJK) 
            U_G(IJK) = U_G(IJK) - D_E(IJK,0)*(PP_G(IJKE)-PP_G(IJK)) 
            V_G(IJK) = V_G(IJK) - D_N(IJK,0)*(PP_G(IJKN)-PP_G(IJK)) 
            IF (DO_K) THEN 
               IJKT = TOP_OF(IJK) 
               W_G(IJK) = W_G(IJK) - D_T(IJK,0)*(PP_G(IJKT)-PP_G(IJK)) 
            ENDIF 
         ENDIF 
      END DO 

      RETURN  
      END SUBROUTINE CORRECT_0G 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CORRECT_0s(Pp_g, d_e, d_n, d_t, U_s, V_s, W_s, IER)    C
!  Purpose: Correct the solids velocities.                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CORRECT_0S(PP_G, D_E, D_N, D_T, U_S, V_S, W_S, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE indices
      USE physprop
      USE compar 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Pressure correction 
      DOUBLE PRECISION Pp_g(DIMENSION_3) 
! 
!                      Pressure correction coefficient -- East 
      DOUBLE PRECISION d_e(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Pressure correction coefficient -- North 
      DOUBLE PRECISION d_n(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Pressure correction coefficient -- Top 
      DOUBLE PRECISION d_t(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Velocity components 
      DOUBLE PRECISION U_s(DIMENSION_3, DIMENSION_M),& 
                       V_s(DIMENSION_3, DIMENSION_M),& 
                       W_s(DIMENSION_3, DIMENSION_M) 
! 
!                      error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          IJK, IJKE, IJKN, IJKT, M 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!  Velocity corrections should not be
!  underrelaxed, so that the continuity eq. is satisfied.
!
      DO M = 1, MMAX 
!$omp    parallel do private(IJK,IJKE,IJKN,IJKT)
         DO IJK = ijkstart3, ijkend3
            IF (FLUIDORP_FLOW_AT(IJK)) THEN 
!
               IJKE = EAST_OF(IJK) 
               IJKN = NORTH_OF(IJK) 
               U_S(IJK,M) = U_S(IJK,M) - D_E(IJK,M)*(PP_G(IJKE)-PP_G(IJK)) 
               V_S(IJK,M) = V_S(IJK,M) - D_N(IJK,M)*(PP_G(IJKN)-PP_G(IJK)) 
               IF (DO_K) THEN 
                  IJKT = TOP_OF(IJK) 
                  W_S(IJK,M) = W_S(IJK,M) - D_T(IJK,M)*(PP_G(IJKT)-PP_G(IJK)) 
               ENDIF 
            ENDIF 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE CORRECT_0S 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
                                                                                                                                                                                                                                              iterate.f                                                                                           0100644 0002444 0000146 00000046746 10250122776 011050  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ITERATE(IER, NIT)                                      C
!  Purpose: This module controls the iterations for solving equations  C
!           Version 2.0                                                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 12-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate the "Kinetic" flag so taht the solids       C
!  calculations are not done using kinetic when doing DES              C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE ITERATE(IER, NIT) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
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
      USE funits 
      USE time_cpu 
      USE pscor
      USE coeff
      USE leqsol 
      USE visc_g
      USE pgcor
      USE cont
      USE scalars
      USE compar   
      USE mpi_utility 
      USE discretelement

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Number of iterations 
      INTEGER          NIT 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
!                      current cpu time used 
      DOUBLE PRECISION CPU_NOW 
! 
!                      MUSTIT = 0 implies complete convergence. 
      INTEGER          MUSTIT 
! 
!                      Sum of solids densities 
      DOUBLE PRECISION SUM 
! 
!                      Weight of solids in the reactor 
      DOUBLE PRECISION SMASS 
! 
!                      Heat loss from the reactor 
      DOUBLE PRECISION HLOSS 
! 
!                      phase index 
      INTEGER          M 
! 
!                      Normalization factor for gas pressure residual 
      DOUBLE PRECISION NORMg 
! 
!                      Normalization factor for solids pressure residual 
      DOUBLE PRECISION NORMs 
! 
!                      Set normalization factor for gas pressure residual 
      LOGICAL          SETg 
! 
!                      Set normalization factor for solids pressure residual 
      LOGICAL          SETs 
! 
!                      gas pressure residual 
      DOUBLE PRECISION RESg 
! 
!                      solids pressure residual 
      DOUBLE PRECISION RESs 
  
      DOUBLE PRECISION TLEFT 

      DOUBLE PRECISION errorpercent(0:MMAX)
       
      LOGICAL          ABORT_IER
      CHARACTER*4 TUNIT 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: VAVG_U_G, VAVG_V_G, VAVG_W_G, VAVG_U_S, &
         VAVG_V_S, VAVG_W_S 
	 external cpu_time   !use the subroutine from machine.f
!-----------------------------------------------
!
!
      DT_prev = DT
      NIT = 0 
      RESG = ZERO 
      RESS = ZERO 
!
      IF (NORM_G == UNDEFINED) THEN 
         NORMG = ONE 
         SETG = .FALSE. 
      ELSE 
         NORMG = NORM_G 
         SETG = .TRUE. 
      ENDIF 
!
      IF (NORM_S == UNDEFINED) THEN 
         NORMS = ONE 
         SETS = .FALSE. 
      ELSE 
         NORMS = NORM_S 
         SETS = .TRUE. 
      ENDIF 
!
      LEQ_ADJUST = .FALSE. 

!
!     Initialize residuals
!
      CALL INIT_RESID (IER) 

!	Initialize the routine for holding gas mass flux constant with cyclic bc
      IF(CYCLIC) CALL GoalSeekMassFlux(0, 0, .false.)
!
!
!     CPU time left
!
      IF (FULL_LOG) THEN      !//
         TLEFT = (TSTOP - TIME)*CPUOS 
         CALL GET_TUNIT (TLEFT, TUNIT) 
!
         IF (DT == UNDEFINED) THEN 
            CALL GET_SMASS (SMASS) 

	    IF(myPE.eq.PE_IO) THEN
              WRITE (*, '(/A,G10.5, A,F9.3,1X,A)') ' Starting solids mass = ', &
                 SMASS, '    CPU time left = ', TLEFT, TUNIT 
	    ENDIF
         ELSE 
            IF(myPE.eq.PE_IO) THEN
!AE TIME 050801 add CN_ON check to print original timestep size
!             IF (CN_ON) THEN
             IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. & 
               (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
	        WRITE (*, '(/A,G12.5, A,G12.5, A,F9.3,1X,A)') ' Time = ', TIME, &
                 '  Dt = ', 2.*DT, '    CPU time left = ', TLEFT, TUNIT 
             ELSE
              WRITE (*, '(/A,G12.5, A,G12.5, A,F9.3,1X,A)') ' Time = ', TIME, &
                 '  Dt = ', DT, '    CPU time left = ', TLEFT, TUNIT 
             ENDIF
	    ENDIF
!
         ENDIF 
      ENDIF 
!
!

      call CALC_RESID_MB(0, errorpercent)

!
!     Calculate the face values of densities and mass fluxes for the first solve_vel_star call.
      CALL CONV_ROP(IER)
      CALL CALC_MFLUX (IER)
      
!     Begin iterations
!
!------------------------------------------------------------------------------
!
!
   50 CONTINUE 
      MUSTIT = 0 
      NIT = NIT + 1 
!
      IF (.NOT.SETG) THEN 
         IF (RESG > SMALL_NUMBER) THEN 
            NORMG = RESG 
            SETG = .TRUE. 
         ENDIF 
      ENDIF 
!
      IF (.NOT.SETS) THEN 
         IF (RESS > SMALL_NUMBER) THEN 
            NORMS = RESS 
            SETS = .TRUE. 
         ENDIF 
      ENDIF 
!
!
!     Call user-defined subroutine to set quantities that need to be updated
!     every iteration
!
      IF (CALL_USR) CALL USR2       
!
!
!     Calculate coefficients.  Explicitly set flags for all the quantities
!     that need to be calculated before calling CALC_COEFF.
!
      VISC(0) = RECALC_VISC_G 
! 	The IF (GRANULAR_ENERGY) statement was commented out to allow calling calc_mu_s even
! 	when the algebraic granular equation is solved (not only the PDE form).
! 	This may enhance convergence. sof, March-10-2005.
       VISC(1:MMAX) = .TRUE. 
!
      CALL CALC_COEFF (DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, RRATE, DRAGCOEF, &
         HEAT_TR, WALL_TR, IER) 

!
!     DIffusion coefficient and source terms for user-defined scalars
      IF(NScalar /= 0)CALL SCALAR_PROP(IER) 

!
!     DIffusion coefficient and source terms for K & Epsilon Eq.
      IF(K_Epsilon) CALL K_Epsilon_PROP(IER)

!      
!
!     Solve strarred velocitiy components
!

      CALL SOLVE_VEL_STAR (IER) 
!
!     Calculate density and reaction rates. Do not change density or reaction rate before the call to
!     solve_vel_star.
      IF (RO_G0 == UNDEFINED) DENSITY(0) = .TRUE. 
      IF (ANY_SPECIES_EQ) RRATE = .TRUE. 
      CALL CALC_COEFF (DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, RRATE, DRAGCOEF, &
         HEAT_TR, WALL_TR, IER) 

!
!     Solve solids volume fraction correction equation for close-packed
!     solids phases
!
      IF(.NOT.DISCRETE_ELEMENT) THEN
        IF (MMAX > 0) THEN
          IF(MMAX == 1)THEN
            CALL CALC_K_CP (K_CP, IER)
	    CALL SOLVE_EPP (NORMS, RESS, IER)
            CALL CORRECT_1 (IER) 
          ELSE
            DO M=1,MMAX 
!   	      IF (M .EQ. MCp) THEN !Volume fraction correction technique for multiparticle types is 
   	      IF (.FALSE.) THEN    !not implemented.  This will only slow down convergence.
                CALL CALC_K_CP (K_CP, IER)
	        CALL SOLVE_EPP (NORMS, RESS, IER)
                CALL CORRECT_1 (IER) 

      	      ELSE
	        CALL SOLVE_CONTINUITY(M,IER)
	             
	      ENDIF
	    END DO
          ENDIF

          CALL CALC_VOL_FR (P_STAR, RO_G, ROP_G, EP_G, ROP_S, IER) 

	  abort_ier = ier.eq.1
	  call global_all_or(abort_ier)
          IF (abort_ier) THEN 
	      ier = 1
              MUSTIT = 2                           !indicates divergence 
              IF(DT/=UNDEFINED)GO TO 1000 
          ENDIF 
 
        ENDIF
      END IF
!
!
!     Calculate P_star in cells where solids continuity equation is
!     solved
!
      IF(.NOT.DISCRETE_ELEMENT) THEN
        IF (MMAX > 0) CALL CALC_P_STAR (EP_G, P_STAR, IER) 
      END IF
!
!     Calculate the face values of densities.
      CALL CONV_ROP(IER)
!
!     Solve fluid pressure correction equation
!
      IF (RO_G0 /= ZERO) CALL SOLVE_PP_G (NORMG, RESG, IER) 
!
!
!     Correct pressure, velocities, and density
!
      IF (RO_G0 /= ZERO) CALL CORRECT_0 (IER) 
     
      IF (RO_G0 == UNDEFINED) THEN
        DENSITY(0) = .TRUE. 
        CALL CALC_COEFF (DENSITY, PSIZE, SP_HEAT, VISC, COND, DIFF, RRATE, DRAGCOEF, &
         HEAT_TR, WALL_TR, IER) 
      ENDIF 

!
!  Update wall velocities
!  modified by sof to force wall functions even when NSW or FSW are declared
!  default wall BC will still be treated as NSW and no wall functions will be used
 
      IF(.NOT. K_EPSILON) CALL SET_WALL_BC (IER) 
!
!     Calculate the face values of densities and mass fluxes.
      CALL CONV_ROP(IER)
      CALL CALC_MFLUX (IER)
!
!     Solve energy equations
!
      IF (ENERGY_EQ) CALL SOLVE_ENERGY_EQ (IER) 
!
!     Solve granular energy equation
!
      IER = 0
      IF (GRANULAR_ENERGY) CALL SOLVE_GRANULAR_ENERGY (IER) 

      abort_ier = ier.eq.1
      call global_all_or(abort_ier)
      IF (abort_ier) THEN
         ier = 1
         MUSTIT = 2                              !indicates divergence 
         IF(DT/=UNDEFINED)GO TO 1000 
!
      ENDIF 

      
!
!     Solve species mass balance equations
!
      CALL SOLVE_SPECIES_EQ (IER) 
!
!     Solve other scalar transport equations
!
      IF(NScalar /= 0) CALL SOLVE_Scalar_EQ (IER) 
!
!     Solve K & Epsilon transport equations
!
      IF(K_Epsilon) CALL SOLVE_K_Epsilon_EQ (IER) 
      
!
!    User-defined linear equation solver parameters may be adjusted after
!    the first iteration
!
      IF (.NOT.CYCLIC) LEQ_ADJUST = .TRUE. 
!
!
!     Check for convergence
!
      call CALC_RESID_MB(1, errorpercent)
      CALL CHECK_CONVERGENCE (NIT, errorpercent(0), MUSTIT, IER) 
      
      IF(CYCLIC)THEN
        IF(MUSTIT==0 .OR. NIT >= MAX_NIT) CALL GoalSeekMassFlux(NIT, MUSTIT, .true.)
        IF(AUTOMATIC_RESTART) RETURN
      ENDIF

!
!      If not converged continue iterations; else exit subroutine.
!
 1000 CONTINUE 
!
!     Display residuals
!
      IF (FULL_LOG) CALL DISPLAY_RESID (NIT, IER) 
      
      IF (MUSTIT == 0) THEN 
         IF (DT==UNDEFINED .AND. NIT==1) GO TO 50!Iterations converged 
         IF (MOD(NSTEP,NLOG) == 0) THEN 
            CALL CPU_TIME (CPU_NOW)
!
            CPUOS = (CPU_NOW - CPU_NLOG)/(TIME - TIME_NLOG) 
            CPU_NLOG = CPU_NOW 
            TIME_NLOG = TIME 
!
            CPU_NOW = CPU_NOW - CPU0 
            call CALC_RESID_MB(1, errorpercent)
            CALL GET_SMASS (SMASS) 
            IF (ENERGY_EQ) CALL GET_HLOSS (HLOSS) 
!
!
            IF (ENERGY_EQ) THEN 
               WRITE (UNIT_LOG, 5000) TIME, DT, NIT, SMASS, errorpercent(0), HLOSS, CPU_NOW 
               IF(FULL_LOG.and.myPE.eq.PE_IO) &          
                       WRITE(*,5000)TIME,DT,NIT,SMASS,errorpercent(0), HLOSS,CPU_NOW        !//
            ELSE 
               WRITE (UNIT_LOG, 5001) TIME, DT, NIT, SMASS, errorpercent(0), CPU_NOW 
               IF (FULL_LOG .and. myPE.eq.PE_IO) &
                       WRITE (*, 5001) TIME, DT, NIT, SMASS, errorpercent(0), CPU_NOW       !//
            ENDIF 
            CALL START_LOG 
            IF (.NOT.FULL_LOG) THEN 
               TLEFT = (TSTOP - TIME)*CPUOS 
               CALL GET_TUNIT (TLEFT, TUNIT) 
               WRITE (UNIT_LOG, '(46X, A, F9.3, 1X, A)') '    CPU time left = '&
                  , TLEFT, TUNIT 
            ENDIF 
!
            IF (CYCLIC_X .OR. CYCLIC_Y .OR. CYCLIC_Z) THEN 
               IF (DO_I) WRITE (UNIT_LOG, 5050) 'U_g = ', VAVG_U_G() 
               IF (DO_J) WRITE (UNIT_LOG, 5050) 'V_g = ', VAVG_V_G() 
               IF (DO_K) WRITE (UNIT_LOG, 5050) 'W_g = ', VAVG_W_G() 
               DO M = 1, MMAX 
                  IF (DO_I) WRITE (UNIT_LOG, 5060) 'U_s(', M, ') = ', VAVG_U_S(&
                     M) 
                  IF (DO_J) WRITE (UNIT_LOG, 5060) 'V_s(', M, ') = ', VAVG_V_S(&
                     M) 
                  IF (DO_K) WRITE (UNIT_LOG, 5060) 'W_s(', M, ') = ', VAVG_W_S(&
                     M) 
               END DO 
            ENDIF 
!
            CALL END_LOG 
         ENDIF 
         IER = 0 
         RETURN  
!                                                ! diverged or
      ELSE IF (MUSTIT==2 .AND. DT/=UNDEFINED) THEN 
         IF (FULL_LOG) THEN 
            CALL START_LOG 
            call CALC_RESID_MB(1, errorpercent)
            WRITE (UNIT_LOG, 5200) TIME, DT, NIT, errorpercent(0) 
            CALL END_LOG 

            if (myPE.eq.PE_IO) WRITE (*, 5200) TIME, DT, NIT, errorpercent(0)   !//
         ENDIF 
         IER = 1 
         RETURN  
      ENDIF 
!
!
      IF (NIT < MAX_NIT) THEN 
         MUSTIT = 0 
         GO TO 50 
      ENDIF 
!
!------------------------------------------------------------------------------
!     End iterations
!
      CALL GET_SMASS (SMASS) 
      if (myPE.eq.PE_IO) WRITE (UNIT_OUT, 5100) TIME, DT, NIT, SMASS    !//
      CALL START_LOG 
      WRITE (UNIT_LOG, 5100) TIME, DT, NIT, SMASS 
      CALL END_LOG 
!
      IER = 0 
      RETURN  
 5000 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT=',I3,' Sm=',G10.5,'MbErr%=', G10.4,' Hl=',G12.5,&
         T84,'CPU=',F8.0,' s') 
 5001 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT=',I3,' Sm=',G10.5,'MbErr%=', G10.4,T84,'CPU=',F8.0&
         ,' s') 
 5050 FORMAT(5X,'Average ',A,G12.5) 
 5060 FORMAT(5X,'Average ',A,I2,A,G12.5) 
 5100 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT>',I3,' Sm= ',G10.5, 'MbErr%=', G10.4) 
 5200 FORMAT(1X,'t=',F10.4,' Dt=',G10.4,' NIT=',&
      I3,'MbErr%=', G10.4, ': Run diverged/stalled :-(') 

      END SUBROUTINE ITERATE 
!
!
      SUBROUTINE GET_TUNIT(TLEFT, TUNIT) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION TLEFT 
      CHARACTER TUNIT*4 
!-----------------------------------------------
!
!
      IF (TLEFT < 3600.0d0) THEN 
         TUNIT = 's' 
      ELSE 
         TLEFT = TLEFT/3600.0d0 
         TUNIT = 'h' 
         IF (TLEFT >= 24.) THEN 
            TLEFT = TLEFT/24.0d0 
            TUNIT = 'days' 
         ENDIF 
      ENDIF 
!
      RETURN  
      END SUBROUTINE GET_TUNIT 



!  Purpose:  In the following module the mass flux across a periodic domain
!            with pressure drop is held constant at a user-specified value.  
!            This module is activated only if the user specifies a value for
!            the keyword Flux_g in the datafile.
      subroutine GoalSeekMassFlux(NIT, MUSTIT, doit)

      USE bc
      USE geometry
      USE constant
      USE compar 
      USE run
      USE time_cpu 
      IMPLICIT NONE

      INTEGER, PARAMETER :: MAXOUTIT = 500
      DOUBLE PRECISION, PARAMETER          :: omega = 0.9
      DOUBLE PRECISION, PARAMETER          :: TOL = 1E-03
      
      INTEGER :: NIT, MUSTIT
      INTEGER, SAVE :: OUTIT
      
      LOGICAl :: doit
      LOGICAl, SAVE :: firstPass = .true.
      
      DOUBLE PRECISION, Save  :: mdot_n, mdot_nm1, delp_n, delp_nm1, err
      DOUBLE PRECISION          :: mdot_0, delp_xyz 
      
      CHARACTER, Save  :: Direction
      
      DOUBLE PRECISION , EXTERNAL :: VAVG_Flux_U_G, VAVG_Flux_V_G, VAVG_Flux_W_G
!
!                      To check when mdot_n or delp_n becomes NaN to auto_restart MFIX
      CHARACTER *80 notnumber
!
      
      IF(CYCLIC_X_MF)THEN
        delp_n = delp_x
      ELSEIF(CYCLIC_Y_MF)THEN
        delp_n = delp_y
      ELSEIF(CYCLIC_Z_MF)THEN
        delp_n = delp_z
      ELSE
        return
      ENDIF
     
      if(.not.doit) then
        OUTIT = 0
	return
      endif
      
      OUTIT = OUTIT + 1
      if(OUTIT > MAXOUTIT) then
        Write(*,5400) MAXOUTIT
        call mfix_exit(0)
      endif
      
      mdot_0 = Flux_g
      
      
      ! calculate the average gas mass flux and error
      IF(CYCLIC_X_MF)THEN
        mdot_n = VAVG_Flux_U_G()
      ELSEIF(CYCLIC_Y_MF)THEN
        mdot_n = VAVG_Flux_V_G()
      ELSEIF(CYCLIC_Z_MF)THEN
        mdot_n = VAVG_Flux_W_G()
      ENDIF
      
      WRITE(notnumber,*) mdot_n
! Check for NaN's in mdot_n
! See if velocity (a real number) contains a letter "n" or symbol "?"
! in which case it's a NaN (Not a Number)
!
      IF(INDEX(notnumber,'?') > 0 .OR.     &
         INDEX(notnumber,'n') > 0 .OR.     &
         INDEX(notnumber,'N') > 0 ) THEN
        write(*,*) mdot_n, ' NaN being caught in GoalSeekMassFlux '
        AUTOMATIC_RESTART = .TRUE.
	RETURN
      ENDIF

      WRITE(notnumber,*) delp_n
! Check for NaN's in delp_n
!
      IF(INDEX(notnumber,'?') > 0 .OR.     &
         INDEX(notnumber,'n') > 0 .OR.     &
         INDEX(notnumber,'N') > 0 ) THEN
        write(*,*) delp_n, ' NaN being caught in GoalSeekMassFlux '
        AUTOMATIC_RESTART = .TRUE.
	RETURN
      ENDIF
! end of NaN's checking...
! 
      err = abs((mdot_n - mdot_0)/mdot_0)
      if( err < TOL)then
        MUSTIT = 0
      else
        MUSTIT = 1
	NIT = 1
      endif
      
      ! correct delp
      if(.not.firstPass)then
!        delp_xyz = delp_n - omega * (delp_n - delp_nm1) * (mdot_n - mdot_0) &
!	                   / (mdot_n - mdot_nm1)
! Fail-Safe Newton's method (below) works better than the regular Newton method (above)
!
        delp_xyz = delp_n - omega * (delp_n - delp_nm1) * ((mdot_n - mdot_0)/(mdot_nm1 - mdot_0)) &
	                   / ((mdot_n - mdot_0)/(mdot_nm1 - mdot_0) - ONE)
      else
        firstPass=.false.
        delp_xyz = delp_n*0.99
      endif
      IF(myPE.eq.PE_IO) Write(*,5500) OUTIT, err, delp_xyz, mdot_n
     
      mdot_nm1 = mdot_n
      delp_nm1 = delp_n
      
      IF(CYCLIC_X_MF)THEN
        delp_x = delp_xyz
      ELSEIF(CYCLIC_Y_MF)THEN
        delp_y = delp_xyz
      ELSEIF(CYCLIC_Z_MF)THEN
        delp_z = delp_xyz
      ENDIF
!
      
      return
5400 FORMAT(/1X,70('*')//' From: GoalSeekMassFlux',/&
      ' Message: Number of outer iterations exceeded ', I4,/1X,70('*')/) 
5500  Format('  MassFluxIteration:', I4, ' Err=', G12.5, ' DelP=', G12.5, ' Flux=', G12.5)
     
    
      end subroutine GoalSeekMassFlux
                          mflux_mod.f                                                                                         0100644 0002444 0000146 00000006044 10247137011 011362  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: mflux_mod                                              C
!  Purpose: Module for mass fluxes and densities at faces               C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: None                                C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
 
 
      MODULE mflux
 
 
      Use param
      Use param1
 
 
!
!                      x-component of gas mass flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_gE 
!
!                      x-component of solids mass flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Flux_sE 
!
!                      y-component of gas mass flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_gN 
!
!                      y-component of solids mass flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Flux_sN 
!
!                      z-component of gas mass flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_gT 
!
!                      z-component of solids mass flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Flux_sT 
!
!
!                      macroscopic gas density at east face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gE 
!
!                      macroscopic solids density at east face
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_sE 
!
!                      macroscopic gas density at north face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gN 
!
!                      macroscopic solids density at north face
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_sN 
!
!                      macroscopic gas density at top face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gT 
!
!                      macroscopic solids density at top face
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_sT 
 

      END MODULE mflux
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            rrates.f                                                                                            0100644 0002444 0000146 00000021735 10240127240 010670  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RRATES(IER)                                            C
!  Purpose: Calculate reaction rates for various reactions in cell ijk C
!                                                                      C
!  Author:                                            Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: MMAX, IJK, T_g, T_s1, D_p, X_g, X_s, EP_g,    C
!            P_g, HOR_g, HOR_s                                         C
!                                                                      C
!                                                                      C
!  Variables modified: M, N, R_gp, R_sp, RoX_gc, RoX_sc, SUM_R_g,      C
!                      SUM_R_s                                         C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
!
      SUBROUTINE RRATES(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE rxns
      USE energy
      USE geometry
      USE run
      USE indices
      USE physprop
      USE constant
      USE funits 
      USE compar        !//d
      USE sendrecv      !// 400
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Local phase and species indices
      INTEGER          L, LM, M, N

!                      cell index
      INTEGER          IJK
      
      DOUBLE PRECISION R_tmp(0:MMAX, 0:MMAX)
!
!-----------------------------------------------
      INCLUDE 'function.inc'



!************************************************************************
      R_tmp = UNDEFINED
!
!  ---  Remember to include all the local variables here for parallel
!  ---- processing
!$omp  parallel do firstprivate(R_tmp), &
!$omp  parallel do private(ijk, L, LM, M, N)

      DO IJK = IJKSTART3, IJKEND3 
      
         IF (FLUID_AT(IJK)) THEN 
!
!
!  User input is required in sections 1 through 4.
!
!1111111111111111111111111111111111111111111111111111111111111111111111111111111
!
! 1. Write the rates of various reactions:
!    Write the reaction rates for each of the reactions as RXNxF and RXNxB (both
!    quantities >= 0), where x identifies the reaction, F stands for forward
!    rate, and B stands for the backward rate.  The rates can be in
!    g-mole/(cm^3.s) or g/(cm^3.s).  For the sake of clarity, give the reaction
!    scheme and the units in a comment statement above the rate expression.
!    The volume (cm^3) is that of the computational cell.  Therefore, for
!    example, the rate term of a gas phase reaction will have a multiplicative
!    factor of epsilon. Note that X_g and X_s are mass fractions
!
!
!
!2222222222222222222222222222222222222222222222222222222222222222222222222222222
!
! 2. Write the formation and consumption rates of various species:
!    Obtain the rates of formation and consumption of various species
!    in g/(cm^3.s) from the rate expressions RXNxF and RXNxB obtained in the
!    previous section.  Pay attention to the units of RXNxF and RXNxB.
!    the formation rates for gas species n are added to get R_gp (IJK, n).
!    All the consumption rates are added and then divided by X_g(IJK, n) to
!    get RoX_gc(IJK, n).  If X_g(IJK, n) is zero and species n is likely
!    to be consumed in a reaction then it is recommended that RoX_gc (IJK, n)
!    be initialized to the derivative of the consumption rate with respect to
!    X_g at X_g=0.
!    If the slope is not known analytically a small value such as 1.0e-9 may
!    instead be used.  A similar procedure is used for all the species in the
!    solids phases also.
!
!  GAS SPECIES
!
!
!  SOLIDS SPECIES
!
!
!3333333333333333333333333333333333333333333333333333333333333333333333333333333
!
! 3.  Determine the g/(cm^3.s) transferred from one phase to the other.
!          R_tmp(To phase #, From phase #)
!     e.g. R_tmp(0,1) -  mass generation of gas phase from solids-1,
!          R_tmp(0,2) -  mass generation of gas phase from solids-2,
!          R_tmp(1,0) -  mass generation of solid-1 from gas = -R_tmp(0,1)
!          R_tmp(1,2) -  mass generation of solid-1 from solids-2.
!     Note, for example, that if gas is generated from solids-1 then
!     R_tmp(0,1) > 0.
!     The R-phase matrix is skew-symmetric and diagonal elements are not needed.
!     Only one of the two skew-symmetric elements -- e.g., R_tmp(0,1) or
!     R_tmp(1,0) -- needs to be specified.
!
!
      if(MMAX > 0) R_tmp(0,1) =  ZERO
!
!4444444444444444444444444444444444444444444444444444444444444444444444444444444
!
! 4.  Determine the heat of reactions in cal/(cm^3.s) at the
!     temperature T_g or T_s1.  Note that for exothermic reactions
!     HOR_g (or HOR_s) will be negative. The assignment of heat of reaction
!     is user defined as it depends upon the microphysics near the interface,
!     which is averaged out in the multiphase flow equations.  For example,
!     heat of Reaction for the C + O2 reaction is split into parts;
!     CO formation is assigned to the solid phase and CO2 formation from CO to
!     the gas phase.
!
!
!==============================================================================
!
!     No user input is required below this line
!-----------------------------------------------------------------------------
!   Determine g/(cm^3.s) of mass generation for each of the phases by adding
!   the reaction rates of all the individual species.

            SUM_R_G(IJK) = ZERO 
            IF (SPECIES_EQ(0)) THEN 
               IF (NMAX(0) > 0) THEN 
                  SUM_R_G(IJK) = SUM_R_G(IJK) + SUM(R_GP(IJK,:NMAX(0))-ROX_GC(&
                     IJK,:NMAX(0))*X_G(IJK,:NMAX(0))) 
               ENDIF
	    ELSE
	      DO M = 1, MMAX
	        IF(R_tmp(0,M) .NE. UNDEFINED)THEN
		  SUM_R_G(IJK) = SUM_R_G(IJK) + R_tmp(0,M)
		ELSEIF(R_tmp(M,0) .NE. UNDEFINED)THEN
		  SUM_R_G(IJK) = SUM_R_G(IJK) - R_tmp(M,0)
		ENDIF
	      ENDDO 
            ENDIF 
!
            DO M = 1, MMAX 
               SUM_R_S(IJK,M) = ZERO 
               IF (SPECIES_EQ(M)) THEN 
                  IF (NMAX(M) > 0) THEN 
                     SUM_R_S(IJK,M) = SUM_R_S(IJK,M) + SUM(R_SP(IJK,M,:NMAX(M))&
                        -ROX_SC(IJK,M,:NMAX(M))*X_S(IJK,M,:NMAX(M))) 
                  ENDIF 
	       ELSE
 	         DO L = 0, MMAX
	           IF(R_tmp(M,L) .NE. UNDEFINED)THEN
		     SUM_R_s(IJK,M) = SUM_R_s(IJK,M) + R_tmp(M,L)
		   ELSEIF(R_tmp(L,M) .NE. UNDEFINED)THEN
		     SUM_R_s(IJK,M) = SUM_R_s(IJK,M) - R_tmp(L,M)
		   ENDIF
	         ENDDO 
               ENDIF 
            END DO 
	    
!
!
!     Store R_tmp values in an array.  Only store the upper triangle without
!     the diagonal of R_tmp array.
!
            DO L = 0, MMAX 
               DO M = L + 1, MMAX 
                  LM = L + 1 + (M - 1)*M/2 
                  IF (R_TMP(L,M) /= UNDEFINED) THEN 
                     R_PHASE(IJK,LM) = R_TMP(L,M) 
                  ELSE IF (R_TMP(M,L) /= UNDEFINED) THEN 
                     R_PHASE(IJK,LM) = -R_TMP(M,L) 
                  ELSE 
                     CALL START_LOG 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) L, M 
                     CALL END_LOG 
                     call mfix_exit(myPE)  
                  ENDIF 
               END DO 
            END DO 
	   
         ENDIF 
      END DO 
      
 1000 FORMAT(/1X,70('*')//' From: RRATES',/&
         ' Message: Mass transfer between phases ',I2,' and ',I2,&
         ' (R_tmp) not specified',/1X,70('*')/) 
      RETURN  
      END SUBROUTINE RRATES 

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
                                   solve_energy_eq.f                                                                                   0100644 0002444 0000146 00000020416 10247405442 012563  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_ENERGY_EQ(IER)                                   C
!  Purpose: Solve energy equations                                     C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-APR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To eliminate kinetic solids calculations when doing DES    C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SOLVE_ENERGY_EQ(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
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
      Use tmp_array, S_p => ARRAY1, S_C => ARRAY2, EPs => ARRAY3, &
                     TxCp => ARRAY4
      Use tmp_array1, VxGama => ARRAYm1
      USE compar   
      USE discretelement 
      USE mflux     
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
!                      phase index 
      INTEGER          m 
      INTEGER          TEMP_MMAX
! 
!                      Septadiagonal matrix A_m 
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Source term on LHS.  Must be positive. 
!      DOUBLE PRECISION S_p(DIMENSION_3) 
! 
!                      Source term on RHS 
!      DOUBLE PRECISION S_C(DIMENSION_3) 
! 
!                      Solids volume fraction 
!      DOUBLE PRECISION EPs(DIMENSION_3) 
! 
!                      ROP * Cp 
!      DOUBLE PRECISION ROPxCp(DIMENSION_3) 
! 
!                      Volume x average gama at cell centers 
!      DOUBLE PRECISION VxGama(DIMENSION_3, DIMENSION_M) 
! 
! 
      DOUBLE PRECISION apo 
! 
!                      Indices 
      INTEGER          IJK 
! 
!                      linear equation solver method and iterations 
      INTEGER          LEQM, LEQI 
!-----------------------------------------------
      INCLUDE 'radtn1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'radtn2.inc'

      call lock_ambm
      call lock_tmp_array
      call lock_tmp_array1

      TEMP_MMAX = MMAX
      IF(DISCRETE_ELEMENT) THEN
         MMAX = 0   ! Only the gas calculations are needed
      END IF            
!
      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      END DO 

      DO IJK = IJKSTART3, IJKEND3
!
         IF(.NOT.WALL_AT(IJK))THEN
            TxCP(IJK) = T_G(IJK)*C_PG(IJK) 
	 ELSE
            TxCP(IJK) = ZERO 
	 ENDIF
	 
         IF (FLUID_AT(IJK)) THEN 
            APO = ROP_GO(IJK)*C_PG(IJK)*VOL(IJK)*ODT 
            S_P(IJK) = APO + S_RPG(IJK)*VOL(IJK) 
            S_C(IJK)=APO*T_GO(IJK)-HOR_G(IJK)*VOL(IJK)+S_RCG(IJK)*VOL(IJK) 
         ELSE 
!
            S_P(IJK) = ZERO 
            S_C(IJK) = ZERO 
!
         ENDIF 
      END DO 
      CALL CONV_DIF_PHI (TxCP, K_G, DISCRETIZE(6), U_G, V_G, W_G, Flux_gE, Flux_gN, Flux_gT, 0, A_M&
         , B_M, IER) 
!
      CALL BC_PHI (T_g, BC_T_G, BC_TW_G, BC_HW_T_G, BC_C_T_G, 0, A_M, B_M, IER) 
!
      CALL SOURCE_PHI (S_P, S_C, EP_G, T_G, 0, A_M, B_M, IER) 
!
      DO M = 1, MMAX 
!
         DO IJK = IJKSTART3, IJKEND3
!
            IF(.NOT.WALL_AT(IJK))THEN
               TxCP(IJK) = T_S(IJK,M)*C_PS(IJK,M) 
	    ELSE
               TxCP(IJK) = ZERO 
	    ENDIF
	 
            IF (FLUID_AT(IJK)) THEN 
               APO = ROP_SO(IJK,M)*C_PS(IJK,M)*VOL(IJK)*ODT 
               S_P(IJK) = APO + S_RPS(IJK,M)*VOL(IJK) 
               S_C(IJK) = APO*T_SO(IJK,M) - HOR_S(IJK,M)*VOL(IJK) + S_RCS(IJK,M&
                  )*VOL(IJK) 
!
               VXGAMA(IJK,M) = GAMA_GS(IJK,M)*VOL(IJK) 
               EPS(IJK) = EP_S(IJK,M) 
!
            ELSE 
!
               S_P(IJK) = ZERO 
               S_C(IJK) = ZERO 
               VXGAMA(IJK,M) = ZERO 
               EPS(IJK) = ZERO 
!
            ENDIF 
         END DO 
         CALL CONV_DIF_PHI (TxCP, K_S(1,M), DISCRETIZE(6), U_S(1,M), V_S(1,&
            M), W_S(1,M), Flux_sE(1,M), Flux_sN(1,M), Flux_sT(1,M), M, A_M, B_M, IER) 
!
         CALL BC_PHI (T_s(1,M), BC_T_S(1,M), BC_TW_S(1,M), BC_HW_T_S(1,M), BC_C_T_S(1,M)&
            , M, A_M, B_M, IER) 
!
         CALL SOURCE_PHI (S_P, S_C, EPS, T_S(1,M), M, A_M, B_M, IER) 
      END DO 
      IF (MMAX > 0) CALL PARTIAL_ELIM_S (T_G, T_S, VXGAMA, A_M, B_M, IER) 
!
      CALL CALC_RESID_S (T_G, A_M, B_M, 0, RESID(RESID_T,0), MAX_RESID(RESID_T,&
         0), IJK_RESID(RESID_T,0), ZERO, IER) 
!
      CALL UNDER_RELAX_S (T_G, A_M, B_M, 0, UR_FAC(6), IER) 
!
!        call check_ab_m(a_m, b_m, 0, .false., ier)
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!        write(*,*)
!     &    resid(resid_t, 0), max_resid(resid_t, 0),
!     &    ijk_resid(resid_t, 0)
!
!
      DO M = 1, MMAX 
!
         CALL CALC_RESID_S (T_S(1,M), A_M, B_M, M, RESID(RESID_T,M), MAX_RESID(&
            RESID_T,M), IJK_RESID(RESID_T,M), ZERO, IER) 
!
         CALL UNDER_RELAX_S (T_S(1,M), A_M, B_M, M, UR_FAC(6), IER) 
      END DO 
      CALL ADJUST_LEQ(RESID(RESID_T,0),LEQ_IT(6),LEQ_METHOD(6),LEQI,LEQM,IER) 
!         call test_lin_eq(a_m(1, -3, 0), LEQI, LEQM, LEQ_SWEEP(6), LEQ_TOL(6), 0, ier)
!
      CALL SOLVE_LIN_EQ ('T_g', T_G, A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(6), LEQ_TOL(6),IER)  

!       bound the temperature
         DO IJK = IJKSTART3, IJKEND3
            IF(.NOT.WALL_AT(IJK))&
              T_g(IJK) = MIN(TMAX, MAX(TMIN, T_g(IJK)))
	 ENDDO

!        call out_array(T_g, 'T_g')
!
      DO M = 1, MMAX 
!
         CALL ADJUST_LEQ (RESID(RESID_T,M), LEQ_IT(6), LEQ_METHOD(6), LEQI, &
            LEQM, IER) 
!         call test_lin_eq(a_m(1, -3, M), LEQI, LEQM, LEQ_SWEEP(6), LEQ_TOL(6), 0, ier)
!
         CALL SOLVE_LIN_EQ ('T_s', T_S(1,M), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(6), LEQ_TOL(6),IER) 

!       bound the temperature
        DO IJK = IJKSTART3, IJKEND3
          IF(.NOT.WALL_AT(IJK))&
            T_s(IJK, M) = MIN(TMAX, MAX(TMIN, T_s(IJK, M))) 
        ENDDO

      END DO 
      
      call unlock_ambm
      call unlock_tmp_array
      call unlock_tmp_array1
      
      MMAX = TEMP_MMAX
      
      RETURN  
      END SUBROUTINE SOLVE_ENERGY_EQ 


!// Comments on the modifications for DMP version implementation      
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
                                                                                                                                                                                                                                                  solve_granular_energy.f                                                                             0100644 0002444 0000146 00000015510 10247407216 013771  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_GRANULAR_ENERGY(IER)                             C
!  Purpose: Solve granular energy equations                            C
!                                                                      C
!                                                                      C
!  Author: K. Agrawal                                 Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SOLVE_GRANULAR_ENERGY(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      USE run
      USE physprop
      USE geometry
      USE fldvar
      USE constant
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
      Use tmp_array, S_p => Array1, S_c => Array2, EPs => Array3, TxCp => Array4
      USE compar      
      USE mflux     
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      phase index 
      INTEGER          m 
! 
!                      Septadiagonal matrix A_m 
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Source term on LHS.  Must be positive. 
!      DOUBLE PRECISION S_p(DIMENSION_3) 
! 
!                      Source term on RHS 
!      DOUBLE PRECISION S_C(DIMENSION_3) 
! 
!                      Solids volume fraction 
!      DOUBLE PRECISION EPs(DIMENSION_3) 
! 
!                      ROP * Cp 
!      DOUBLE PRECISION ROPxCp(DIMENSION_3) 
! 
! 
      DOUBLE PRECISION apo, sourcelhs, sourcerhs 
! 
!                      Indices 
      INTEGER          IJK 
! 
!                      linear equation solver method and iterations 
      INTEGER          LEQM, LEQI 
!-----------------------------------------------
      INCLUDE 'radtn1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'radtn2.inc'
      
      call lock_ambm
      call lock_tmp_array

!
      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      END DO 
      DO M = 1, MMAX 
!
!
         DO IJK = ijkstart3, ijkend3
!
            IF (FLUID_AT(IJK)) THEN 
!
               CALL SOURCE_GRANULAR_ENERGY (SOURCELHS, SOURCERHS, IJK, M, IER) 
               TXCP(IJK) = 1.5D0*THETA_M(IJK,M) 
               APO = 1.5D0*ROP_SO(IJK,M)*VOL(IJK)*ODT 
               S_P(IJK) = APO + SOURCELHS + ZMAX(SUM_R_S(IJK,M)) * VOL(IJK) 
               S_C(IJK) = APO*THETA_MO(IJK,M) + SOURCERHS + &
	                  THETA_M(IJK,M)*ZMAX((-SUM_R_S(IJK,M))) * VOL(IJK)
               EPS(IJK) = EP_S(IJK,M) 
!
            ELSE 
!
               EPS(IJK) = ZERO 
               TXCP(IJK) = ZERO 
               S_P(IJK) = ZERO 
               S_C(IJK) = ZERO 
!
            ENDIF 
         END DO 
         CALL CONV_DIF_PHI (TXCP, KTH_S(1,M), DISCRETIZE(8), U_S(1,M), &
            V_S(1,M), W_S(1,M), Flux_sE(1,M), Flux_sN(1,M), Flux_sT(1,M), M, A_M, B_M, IER) 
!
         CALL BC_PHI (THETA_M(1,M), BC_THETA_M(1,M), BC_THETAW_M(1,M), BC_HW_THETA_M(1,M), &
            BC_C_THETA_M(1,M), M, A_M, B_M, IER) 
!
         CALL BC_THETA (M, A_M, B_M, IER)        !override bc settings if 
!                                        Johnson-Jackson bcs are specified
!
         CALL SOURCE_PHI (S_P, S_C, EPS, THETA_M(1,M), M, A_M, B_M, IER) 
!
!
! Adjusting the values of theta_m to zero when Ep_g < EP_star (Shaeffer, 1987)
! This is done here instead of calc_mu_s.f to avoid convergence problems. (sof)
!
         IF (SCHAEFFER) THEN
           DO IJK = ijkstart3, ijkend3
!
              IF (FLUID_AT(IJK) .AND. EP_g(IJK) .LT. EP_star) THEN 
!

                 A_M(IJK,1,M) = ZERO 
                 A_M(IJK,-1,M) = ZERO 
                 A_M(IJK,2,M) = ZERO 
                 A_M(IJK,-2,M) = ZERO 
                 A_M(IJK,3,M) = ZERO 
                 A_M(IJK,-3,M) = ZERO 
                 A_M(IJK,0,M) = -ONE 		  
                 B_M(IJK,M) = ZERO
	      ENDIF
	   END DO
	 ENDIF	 
! End of Shaeffer adjustments, sof.
!
         CALL CALC_RESID_S (THETA_M(1,M), A_M, B_M, M, RESID(RESID_TH,M), &
            MAX_RESID(RESID_TH,M), IJK_RESID(RESID_TH,M), ZERO, IER) 
!
         CALL UNDER_RELAX_S (THETA_M(1,M), A_M, B_M, M, UR_FAC(8), IER) 
!
!        call check_ab_m(a_m, b_m, m, .true., ier)
!          write(*,*)
!     &      resid(resid_th, m), max_resid(resid_th, m),
!     &      ijk_resid(resid_th, m)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!
!
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &    ier)
!
         CALL ADJUST_LEQ (RESID(RESID_TH,M), LEQ_IT(8), LEQ_METHOD(8), LEQI, &
            LEQM, IER) 
!
         CALL SOLVE_LIN_EQ ('Theta_m', THETA_M(1,M), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(8), LEQ_TOL(8),IER) 
!          call out_array(Theta_m(1,m), 'Theta_m')
!
!
!        Remove very small negative values of theta caused by leq solvers
         CALL ADJUST_THETA (M, IER) 
         IF (IER /= 0) RETURN                    !large negative granular temp -> divergence 
      END DO 
      
      call unlock_ambm
      call unlock_tmp_array

      RETURN  
      END SUBROUTINE SOLVE_GRANULAR_ENERGY 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
                                                                                                                                                                                        solve_k_epsilon_eq.f                                                                                0100644 0002444 0000146 00000030406 10247407021 013250  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_K_Epsilon_EQ(IER)                                C
!  Purpose: Solve K & Epsilon equations for a turbulent flow           C
!                                                                      C
!                                                                      C
!  Author: S. Benyahia                                Date: MAY-13-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SOLVE_K_Epsilon_EQ(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      USE run
      USE physprop
      USE geometry
      !USE matrix
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
      USE turb
      USE usr
      Use ambm
      Use tmp_array, S_p => Array1, S_c => Array2, EPs => Array3, VxGama => Array4
      USE compar      
      USE mflux     
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Error index 
      INTEGER          IER 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
! 
!                      phase index 
      INTEGER          m , I,IM, I1, J, K, LC
! 
!                      species index 
      INTEGER          n 
! 
      DOUBLE PRECISION apo
! 

! 
!                      temporary variables in residual computation 
      DOUBLE PRECISION res1, mres1
      INTEGER          ires1 
! 
!                      Indices 
      INTEGER          IJK, IPJK, IMJK,IJMK,IJPK
! 
!                      linear equation solver method and iterations 
      INTEGER          LEQM, LEQI 
      
!                      A default zero flux will be defined for both K & Epsilon at walls
      DOUBLE PRECISION BC_hw_K_Turb_G (DIMENSION_BC),  BC_hw_E_Turb_G (DIMENSION_BC)    
      DOUBLE PRECISION BC_K_Turb_GW (DIMENSION_BC),   BC_E_Turb_GW (DIMENSION_BC) 
      DOUBLE PRECISION BC_C_K_Turb_G (DIMENSION_BC),  BC_C_E_Turb_G (DIMENSION_BC)    
!
      character*8      Vname
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: IS_SMALL 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'

      IF( .NOT. K_Epsilon) RETURN
      
      call lock_ambm
      call lock_tmp_array
      
      RESID(RESID_ke,0) = ZERO
      MAX_RESID(RESID_ke,0) = ZERO
      IJK_RESID(RESID_ke,0) = 0

! Setting default zero flux for K & Epsilon since we use wall functions.
! If an expert user want to use Low Re K-Epilon model and needs to set
! the turbulence quatities to zero at walls, then set the hw's to UNDEFINE will
! do it. All the variables below can be changed in the same way as in the 
! MFIX data file in the boundary conditions section.
!
      DO LC = 1, DIMENSION_BC 
        BC_hw_K_Turb_G (LC) = ZERO
	BC_hw_E_Turb_G (LC) = ZERO
	BC_K_Turb_GW (LC) = ZERO
	BC_E_Turb_GW (LC) = ZERO
	BC_C_K_Turb_G (LC) = ZERO
	BC_C_E_Turb_G (LC) = ZERO
      ENDDO
! End of setting default zero flux for K & Epsilon wall boundary conditions
! 
! Equations solved for gas phase, thus M = 0
      M = 0  
          CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 

! Solve first fot the K_Turb_G Equation
  
          DO IJK = IJKSTART3, IJKEND3
!
             I = I_OF(IJK)
             J = J_OF(IJK)
             K = K_OF(IJK)
               IF (FLUID_AT(IJK)) THEN
   
                  APO = ROP_G(IJK)*VOL(IJK)*ODT
                  S_P(IJK) = APO + (  ZMAX(SUM_R_G(IJK)) &
                                    + K_Turb_G_p(IJK) )*VOL(IJK) 
                  S_C(IJK) =   APO*K_Turb_GO(IJK) &
                            + K_Turb_G(IJK)*ZMAX((-SUM_R_G(IJK)))*VOL(IJK) &
                            + K_Turb_G_c(IJK) *VOL(IJK)
               ELSE 
!
                  S_P(IJK) = ZERO 
                  S_C(IJK) = ZERO 
!
               ENDIF 

            END DO   
    
            CALL CONV_DIF_PHI (K_Turb_G, DIF_K_Turb_G, DISCRETIZE(9), &
                               U_G, V_G, W_G, Flux_gE, Flux_gN, Flux_gT, M, A_M, B_M, IER) 
!
!
            CALL BC_PHI (K_Turb_G, BC_K_Turb_G, BC_K_Turb_GW, BC_HW_K_Turb_G, &
                         BC_C_K_Turb_G, M, A_M, B_M, IER)
!
!
            CALL SOURCE_PHI (S_P, S_C, EP_G, K_Turb_G, M, A_M, B_M, IER) 
!           
	    CALL CALC_RESID_S (K_Turb_G, A_M, B_M, M, res1, &
               mres1, ires1, ZERO, IER) 
               RESID(RESID_ke,0) = RESID(RESID_ke,0)+res1
               if(mres1 .gt. MAX_RESID(RESID_ke,0))then
                 MAX_RESID(RESID_ke,0) = mres1
                 IJK_RESID(RESID_ke,0) = ires1
               endif
!
            CALL UNDER_RELAX_S (K_Turb_G, A_M, B_M, M, UR_FAC(9), IER) 
!
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!          write(*,*) res1, mres1, &
!           ires1
!
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K, &
!          ier)
!
            CALL ADJUST_LEQ (res1, LEQ_IT(9), LEQ_METHOD(9), &
               LEQI, LEQM, IER) 
!
            write(Vname, '(A,I2)')'K_Turb_G'
            CALL SOLVE_LIN_EQ (Vname, K_Turb_G, A_M, B_M, M, LEQI, LEQM, &
                             LEQ_SWEEP(9), LEQ_TOL(9),IER) 
!          call out_array(K_Turb_G, Vname)

!!!!!!!!!!!!!
! Now, we'll solve for the E_Turb_G (dissipation) Equation. 
!!!!!!!!!!!!!
!
! Initiate (again) the Am Bm matrix. This has to be done for every scalar equation.
          CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
  
          DO IJK = IJKSTART3, IJKEND3
!
             I = I_OF(IJK)
             J = J_OF(IJK)
             K = K_OF(IJK)
               IF (FLUID_AT(IJK)) THEN
   
                  APO = ROP_G(IJK)*VOL(IJK)*ODT
                  S_P(IJK) = APO + (  ZMAX(SUM_R_G(IJK)) &
                                    + E_Turb_G_p(IJK)  )*VOL(IJK)
                  S_C(IJK) =   APO*E_Turb_GO(IJK) &
                            + E_Turb_G(IJK)*ZMAX((-SUM_R_G(IJK)))*VOL(IJK) &
                            + E_Turb_G_c(IJK) *VOL(IJK)
               ELSE 
!
                  S_P(IJK) = ZERO 
                  S_C(IJK) = ZERO 
!
               ENDIF 

            END DO   
    
            CALL CONV_DIF_PHI (E_Turb_G, DIF_E_Turb_G, DISCRETIZE(9), &
                               U_G, V_G, W_G, Flux_gE, Flux_gN, Flux_gT, M, A_M, B_M, IER) 
!
!
            CALL BC_PHI (E_Turb_G, BC_E_Turb_G, BC_E_Turb_GW, BC_HW_E_Turb_G, &
                         BC_C_E_Turb_G, M, A_M, B_M, IER)
!
!
            CALL SOURCE_PHI (S_P, S_C, EP_G, E_Turb_G, M, A_M, B_M, IER) 
!
! When implementing the wall functions, The Epsilon (dissipation) value at the fluid cell
! near the walls needs to be set.

            DO IJK = IJKSTART3, IJKEND3
             I = I_OF(IJK)
             J = J_OF(IJK)
             K = K_OF(IJK)
!
               IF (FLUID_AT(IJK)) THEN  
!
                   IF(WALL_AT(JP_OF(IJK)).OR.WALL_AT(JM_OF(IJK))) THEN 
                     A_M(IJK,1,M) = ZERO 
                     A_M(IJK,-1,M) = ZERO 
                     A_M(IJK,2,M) = ZERO 
                     A_M(IJK,-2,M) = ZERO 
                     A_M(IJK,3,M) = ZERO 
                     A_M(IJK,-3,M) = ZERO 
                     A_M(IJK,0,M) = -ONE 
                     B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)/DY(J) &
                                 *2.0D+0/0.42D+0

                   ELSE IF(WALL_AT(KP_OF(IJK)).OR.WALL_AT(KM_OF(IJK))) THEN 
                     A_M(IJK,1,M) = ZERO 
                     A_M(IJK,-1,M) = ZERO 
                     A_M(IJK,2,M) = ZERO 
                     A_M(IJK,-2,M) = ZERO 
                     A_M(IJK,3,M) = ZERO 
                     A_M(IJK,-3,M) = ZERO 
                     A_M(IJK,0,M) = -ONE 
                     B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)* &
                                  (ODZ(K)*OX(I)*2.0D+0)/0.42D+0
                   ENDIF  !for identifying wall cells in J or K direction
   
                 IF(CYLINDRICAL) THEN
                  IF (WALL_AT(IP_OF(IJK)))  THEN
                   A_M(IJK,1,M) = ZERO 
                   A_M(IJK,-1,M) = ZERO 
                   A_M(IJK,2,M) = ZERO 
                   A_M(IJK,-2,M) = ZERO 
                   A_M(IJK,3,M) = ZERO 
                   A_M(IJK,-3,M) = ZERO 
                   A_M(IJK,0,M) = -ONE 
                   B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)/DX(I) &
                               *2.0D+0/0.42D+0

                     ENDIF! for wall cells in I direction

                 ELSE IF (WALL_AT(IP_OF(IJK)).OR.WALL_AT(IM_OF(IJK))) THEN
                     A_M(IJK,1,M) = ZERO 
                     A_M(IJK,-1,M) = ZERO 
                     A_M(IJK,2,M) = ZERO 
                     A_M(IJK,-2,M) = ZERO 
                     A_M(IJK,3,M) = ZERO 
                     A_M(IJK,-3,M) = ZERO 
                     A_M(IJK,0,M) = -ONE 
                     B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)/DX(I) &
                                   *2.0D+0/0.42D+0

                  ENDIF ! for cylindrical 

!			
               ENDIF  !for fluid at ijk
            ENDDO   	    
           
	    CALL CALC_RESID_S (E_Turb_G, A_M, B_M, M, res1, &
               mres1, ires1, ZERO, IER) 
               RESID(RESID_ke,0) = RESID(RESID_ke,0)+res1
               if(mres1 .gt. MAX_RESID(RESID_ke,0))then
                 MAX_RESID(RESID_ke,0) = mres1
                 IJK_RESID(RESID_ke,0) = ires1
               endif
!
            CALL UNDER_RELAX_S (E_Turb_G, A_M, B_M, M, UR_FAC(9), IER) 
!
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!          write(*,*) res1, mres1, &
!           ires1
!
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K, &
!          ier)
!
            CALL ADJUST_LEQ (res1, LEQ_IT(9), LEQ_METHOD(9), &
               LEQI, LEQM, IER) 
!
            write(Vname, '(A,I2)')'E_Turb_G'
            CALL SOLVE_LIN_EQ (Vname, E_Turb_G, A_M, B_M, M, LEQI, LEQM, &
                             LEQ_SWEEP(9), LEQ_TOL(9),IER) 
!          call out_array(E_Turb_G, Vname)
!
! remove small negative K and Epsilon values generated by linear solver 
! same as adjust_theta.f
!
           DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN 
	     IF(K_Turb_G(IJK) < ZERO_EP_S) K_Turb_G(IJK) = ZERO_EP_S 
	     IF(E_Turb_G(IJK) < ZERO_EP_S) E_Turb_G(IJK) = ZERO_EP_S 
!
            ENDIF 
           END DO 

      call unlock_ambm
      call unlock_tmp_array

      RETURN  
      END SUBROUTINE SOLVE_K_Epsilon_EQ 


!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
                                                                                                                                                                                                                                                          solve_scalar_eq.f                                                                                   0100644 0002444 0000146 00000017772 10247406516 012555  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_SCALAR_EQ(IER)                                   C
!  Purpose: Solve scalar transport equations                           C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 4-12-99    C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SOLVE_Scalar_EQ(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
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
      USE scalars
      Use ambm
      Use tmp_array, S_p => Array1, S_c => Array2, EPs => Array3, VxGama => Array4
      USE compar      
      USE mflux     
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Error index 
      INTEGER          IER 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
! 
!                      phase index 
      INTEGER          m 
! 
!                      species index 
      INTEGER          n 
! 
      DOUBLE PRECISION apo 
! 

! 
!                      temporary variables in residual computation 
      DOUBLE PRECISION res1, mres1
      INTEGER          ires1 
! 
!                      Indices 
      INTEGER          IJK 
! 
!                      linear equation solver method and iterations 
      INTEGER          LEQM, LEQI 
      
      character*8      Vname
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: IS_SMALL 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
      
      call lock_ambm
      call lock_tmp_array
      
      RESID(RESID_sc,0) = ZERO
      MAX_RESID(RESID_sc,0) = ZERO
      IJK_RESID(RESID_sc,0) = 0
!
!     Fluid phase species mass balance equations
!
       DO N = 1, NScalar
	    
	  M = Phase4Scalar(N)
	    
          CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
	    
	  IF(M == 0) THEN
	  
          DO IJK = IJKSTART3, IJKEND3
!
               IF (FLUID_AT(IJK)) THEN 
                  APO = ROP_GO(IJK)*VOL(IJK)*ODT 
                  S_P(IJK) = APO + (  ZMAX(SUM_R_G(IJK)) &
		                    + Scalar_p(IJK, N)      )*VOL(IJK) 
                  S_C(IJK) =   APO*ScalarO(IJK,N) &
		             + Scalar(IJK,N)*ZMAX((-SUM_R_G(IJK)))*VOL(IJK) &
			     + Scalar_c(IJK, N)*VOL(IJK)
               ELSE 
!
                  S_P(IJK) = ZERO 
                  S_C(IJK) = ZERO 
!
               ENDIF 
            END DO 
	    
    
            CALL CONV_DIF_PHI (Scalar(1,N), DIF_Scalar(1,N), DISCRETIZE(9), &
	                       U_G, V_G, W_G, Flux_gE, Flux_gN, Flux_gT, M, A_M, B_M, IER) 
!
!
            CALL BC_PHI (Scalar(1,N), BC_Scalar(1,N), BC_ScalarW(1,N), BC_HW_Scalar(1,N), &
	                 BC_C_Scalar(1,N), M, A_M, B_M, IER) 
!
!
            CALL SOURCE_PHI (S_P, S_C, EP_G, Scalar(1,N), M, A_M, B_M, IER) 
!
            CALL CALC_RESID_S (Scalar(1,N), A_M, B_M, M, res1, &
               mres1, ires1, ZERO, IER) 
               RESID(RESID_sc,0) = RESID(RESID_sc,0)+res1
	       if(mres1 .gt. MAX_RESID(RESID_sc,0))then
	         MAX_RESID(RESID_sc,0) = mres1
		 IJK_RESID(RESID_sc,0) = ires1
	       endif
!
            CALL UNDER_RELAX_S (Scalar(1,N), A_M, B_M, M, UR_FAC(9), IER) 
!
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!          write(*,*) res1, mres1, &
!           ires1
!
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K, &
!          ier)
!
            CALL ADJUST_LEQ (res1, LEQ_IT(9), LEQ_METHOD(9), &
               LEQI, LEQM, IER) 
!
            write(Vname, '(A,I2)')'Scalar',N
            CALL SOLVE_LIN_EQ (Vname, Scalar(1,N), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(9), LEQ_TOL(9),IER) 
!          call out_array(Scalar(1, N), Vname)
!
         ELSE
	 
            DO IJK = IJKSTART3, IJKEND3
!
               IF (FLUID_AT(IJK)) THEN 
                  APO = ROP_sO(IJK, M)*VOL(IJK)*ODT 
                  S_P(IJK) = APO + (  ZMAX(SUM_R_s(IJK, M)) &
		                    + Scalar_p(IJK, N)       )*VOL(IJK) 
                  S_C(IJK) =   APO*ScalarO(IJK,N) &
		             + Scalar(IJK,N)*ZMAX((-SUM_R_s(IJK, M)))*VOL(IJK)&
			     + Scalar_c(IJK, N)*VOL(IJK)
                  EPs(IJK) = EP_s(IJK, M)
               ELSE 
!
                  S_P(IJK) = ZERO 
                  S_C(IJK) = ZERO 
                  EPS(IJK) = ZERO 
!
               ENDIF 
            END DO 
	    
	    
            CALL CONV_DIF_PHI (Scalar(1,N), DIF_Scalar(1,N), DISCRETIZE(9), &
	                       U_s(1,m), V_s(1,m), W_s(1,m), Flux_sE(1,M), Flux_sN(1,M), Flux_sT(1,M), M, &
			       A_M, B_M, IER) 
!
!
            CALL BC_PHI (Scalar(1,N), BC_Scalar(1,N), BC_ScalarW(1,N), BC_HW_Scalar(1,N), &
	                 BC_C_Scalar(1,N), M, A_M, B_M, IER) 
!
!
            CALL SOURCE_PHI (S_P, S_C, EPs, Scalar(1,N), M, A_M, B_M, IER) 
!
            CALL CALC_RESID_S (Scalar(1,N), A_M, B_M, M, res1, &
               mres1, ires1, ZERO, IER) 
               RESID(RESID_sc,0) = RESID(RESID_sc,0)+res1
	       if(mres1 .gt. MAX_RESID(RESID_sc,0))then
	         MAX_RESID(RESID_sc,0) = mres1
		 IJK_RESID(RESID_sc,0) = ires1
	       endif
!
            CALL UNDER_RELAX_S (Scalar(1,N), A_M, B_M, M, UR_FAC(9), IER) 
!
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!          write(*,*) res1, mres1, ires1
!
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K, &
!          ier)
!
            CALL ADJUST_LEQ (res1, LEQ_IT(9), LEQ_METHOD(9), &
               LEQI, LEQM, IER) 
!
            write(Vname, '(A,I2)')'Scalar',N
            CALL SOLVE_LIN_EQ (Vname, Scalar(1,N), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(9), LEQ_TOL(9),IER) 
!          call out_array(Scalar(1, N), Vname)
!
         END IF 
      END DO 
      
      call unlock_ambm
      call unlock_tmp_array

      RETURN  
      END SUBROUTINE SOLVE_Scalar_EQ 


!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
      solve_species_eq.f                                                                                  0100644 0002444 0000146 00000017317 10247142360 012730  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_SPECIES_EQ(IER)                                  C
!  Purpose: Solve species mass balance equations                       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 11-FEB-98  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SOLVE_SPECIES_EQ(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
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
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Error index 
      INTEGER          IER 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
! 
!                      phase index 
      INTEGER          m 
! 
!                      species index 
      INTEGER          ln 
! 
! 
      DOUBLE PRECISION apo
      DOUBLE PRECISION errorpercent(0:MMAX)
! 
!                      Indices 
      INTEGER          IJK, IMJK, IJMK, IJKM 
 
! 
!                      linear equation solver method and iterations 
      INTEGER          LEQM, LEQI 
      
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: IS_SMALL 
      DOUBLE PRECISION , EXTERNAL :: Check_conservation 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
     
     
      call lock_ambm
      call lock_tmp_array
 
!
!     Fluid phase species mass balance equations
!
      IF (SPECIES_EQ(0)) THEN 
          if(chi_scheme) call set_chi(DISCRETIZE(7), X_g, NMAX(0), U_g, V_g, W_g, IER)
         DO LN = 1, NMAX(0) 
            CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0, IER) 
!$omp    parallel do private(IJK, APO)
            DO IJK = ijkstart3, ijkend3 
!
               IF (FLUID_AT(IJK)) THEN
                   APO = ROP_GO(IJK)*VOL(IJK)*ODT 
                   S_P(IJK) = APO + (ZMAX(SUM_R_G(IJK))+ROX_GC(IJK,LN))*VOL(IJK) 
                   S_C(IJK) = APO*X_GO(IJK,LN) + X_G(IJK,LN)*ZMAX((-SUM_R_G(IJK)))&
                    *VOL(IJK) + R_GP(IJK,LN)*VOL(IJK)
               ELSE 
!
                  S_P(IJK) = ZERO 
                  S_C(IJK) = ZERO 
!
               ENDIF 
            END DO
	    
	     
            CALL CONV_DIF_PHI (X_G(1,LN), DIF_G(1,LN), DISCRETIZE(7), U_G, V_G, &
               W_G, Flux_gE, Flux_gN, Flux_gT, 0, A_M, B_M, IER)

!
            CALL BC_PHI (X_G(1,LN), BC_X_G(1,LN), BC_XW_G(1,LN), BC_HW_X_G(1,LN), BC_C_X_G(1,&
               LN), 0, A_M, B_M, IER) 
!
            CALL SOURCE_PHI (S_P, S_C, EP_G, X_G(1,LN), 0, A_M, B_M, IER)
!
            CALL CALC_RESID_S (X_G(1,LN), A_M, B_M, 0, RESID(RESID_X+(LN-1),0), &
               MAX_RESID(RESID_X+(LN-1),0), IJK_RESID(RESID_X+(LN-1),0), &
	       ZERO_X_GS, IER) 

            CALL UNDER_RELAX_S (X_G(1,LN), A_M, B_M, 0, UR_FAC(7), IER) 

            CALL ADJUST_LEQ (RESID(RESID_X+(LN-1),0), LEQ_IT(7), LEQ_METHOD(7), &
                LEQI, LEQM, IER)
!
            CALL SOLVE_LIN_EQ ('X_g', X_G(1,LN), A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(7), LEQ_TOL(7),IER) 

            CALL BOUND_X (X_G(1,LN), IJKMAX2, IER) 

         END DO 
         if(chi_scheme) call unset_chi(IER)
      ENDIF 
!
!
!     Granular phase species balance equations
!
      DO M = 1, MMAX 
!
         IF (SPECIES_EQ(M)) THEN 
!
!
            DO LN = 1, NMAX(M) 
!
               CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
!
!$omp    parallel do private(IJK, APO)
               DO IJK = ijkstart3, ijkend3 
!
                  IF (FLUID_AT(IJK)) THEN 
!
                    APO = ROP_SO(IJK,M)*VOL(IJK)*ODT 
                    S_P(IJK) = APO + (ZMAX(SUM_R_S(IJK,M))+ROX_SC(IJK,M,LN))*&
                        VOL(IJK) 
                    S_C(IJK) = APO*X_SO(IJK,M,LN) + X_S(IJK,M,LN)*ZMAX((-SUM_R_S&
                        (IJK,M)))*VOL(IJK) + R_SP(IJK,M,LN)*VOL(IJK)
!
                    EPS(IJK) = EP_S(IJK,M) 
!
                  ELSE 
!
                     S_P(IJK) = ZERO 
                     S_C(IJK) = ZERO 
                     EPS(IJK) = ZERO 
!
                  ENDIF 
               END DO 
               CALL CONV_DIF_PHI (X_S(1,M,LN), DIF_S(1,M,LN), DISCRETIZE(7), U_S(&
                  1,M), V_S(1,M), W_S(1,M), Flux_sE(1,M), Flux_sN(1,M), Flux_sT(1,M), M, A_M, B_M, IER)
!
               CALL BC_PHI (X_S(1,M,LN), BC_X_S(1,M,LN), BC_XW_S(1,M,LN), BC_HW_X_S(1,M,LN), &
                  BC_C_X_S(1,M,LN), M, A_M, B_M, IER) 
!
!
               CALL SOURCE_PHI (S_P, S_C, EPS, X_S(1,M,LN), M, A_M, B_M, IER)

               CALL CALC_RESID_S (X_S(1,M,LN), A_M, B_M, M, RESID(RESID_X+(LN-1),&
                  M), MAX_RESID(RESID_X+(LN-1),M), IJK_RESID(RESID_X+(LN-1),M), &
                  ZERO_X_GS, IER) 
!
!
               CALL UNDER_RELAX_S (X_S(1,M,LN), A_M, B_M, M, UR_FAC(7), IER) 
!           call check_ab_m(a_m, b_m, m, .false., ier)
!           write(*,*)
!     &      resid(resid_x+(LN-1), m), max_resid(resid_x+(LN-1), m),
!     &      ijk_resid(resid_x+(LN-1), m)
!           call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!
!           call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &     ier)
!
               CALL ADJUST_LEQ (RESID(RESID_X+(LN-1),M), LEQ_IT(7), LEQ_METHOD(7&
                  ), LEQI, LEQM, IER) 
!
               CALL SOLVE_LIN_EQ ('X_s', X_S(1,M,LN), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(7), LEQ_TOL(7),IER) 
!            call out_array(X_s(1,m,LN), 'X_s')
!               CALL BOUND_X (X_S(1,M,LN), IJKMAX2, IER) 
!
            END DO 
         ENDIF 
      END DO 
      
      call unlock_ambm
      call unlock_tmp_array

      RETURN  
      END SUBROUTINE SOLVE_SPECIES_EQ 
      

                                                                                                                                                                                                                                                                                                                 solve_vel_star.f                                                                                    0100644 0002444 0000146 00000036731 10250130232 012415  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_VEL_STAR(IER)                                    C
!  Purpose: Solve starred velocity components                          C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow multiparticle in D_E, D_N and D_T claculations       C
!           and account for the averaged Solid-Solid drag              C
!                                                                      C
!  Author: S. Dartevelle, LANL                        Date: 28-FEb-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: To flag the solids calculations whether Kinetic or DES.    C
!  And to change the gas arrays incorporating the drag, when doing DES C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SOLVE_VEL_STAR(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
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
      Use ambm
      Use tmp_array1,  VxF_gs => Arraym1
      Use tmp_array,  VxF_ss => ArrayLM     !S. Dartevelle, LANL, MARCH 2004
      USE compar
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
!                      Error index 
      INTEGER          IER 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
! 
! 
!                      phase index 
      INTEGER          m, UVEL, VVEL, WVEL 
      INTEGER          IJK 
! 
!                      temporary velocity arrays 
      DOUBLE PRECISION U_gtmp(DIMENSION_3),  V_gtmp(DIMENSION_3), W_gtmp(DIMENSION_3)
      DOUBLE PRECISION U_stmp(DIMENSION_3, DIMENSION_M),  V_stmp(DIMENSION_3, DIMENSION_M), W_stmp(DIMENSION_3, DIMENSION_M)
! 
!                      Vector b_m 
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Volume x average at momentum cell centers 
!      DOUBLE PRECISION VxF_gs(DIMENSION_3, DIMENSION_M) 
! 
!                      linear equation solver method and iterations 
      INTEGER          LEQM, LEQI 
!-----------------------------------------------
!
      call lock_ambm
      call lock_tmp_array1
      call lock_tmp_array2
      
!     Store the velocities so that the order of solving the momentum equations does not matter
      DO IJK = ijkstart3, ijkend3
        U_gtmp(IJK) = U_g(IJK)
        V_gtmp(IJK) = V_g(IJK)
        W_gtmp(IJK) = W_g(IJK)
      ENDDO
      DO M = 1, MMAX
        DO IJK = ijkstart3, ijkend3
          U_stmp(IJK, M) = U_s(IJK, M)
          V_stmp(IJK, M) = V_s(IJK, M)
          W_stmp(IJK, M) = W_s(IJK, M)
        ENDDO
      ENDDO
      
      IF(DISCRETE_ELEMENT) THEN
         UVEL = 0
         VVEL = 0
         WVEL = 0
      END IF
!
      IF (MMAX == 0) CALL ZERO_ARRAY (VXF_GS(1,1), IER)
	  IF (MMAX == 1) CALL ZERO_ARRAY (VXF_SS(1,1), IER)
!
!  2.1 Calculate U_m_star and residuals
!
      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      END DO 
      
      CALL CONV_DIF_U_G (A_M, B_M, IER) 
      
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL CONV_DIF_U_S (A_M, B_M, IER) 
      END IF
!
      CALL SOURCE_U_G (A_M, B_M, IER) 
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL SOURCE_U_S (A_M, B_M, IER) 
      END IF
!
      IF (MMAX > 0) CALL VF_GS_X (VXF_GS, IER)
      IF(.NOT.DISCRETE_ELEMENT) THEN
         IF (MMAX > 0) CALL VF_SS_X (VXF_SS, IER)   !S. Dartevelle, LANL, Feb.2004
      END IF
!
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL CALC_D_E (A_M, VXF_GS, VXF_SS, D_E, IER)  !S. Dartevelle, LANL, Feb.2004
      ELSE IF(DISCRETE_ELEMENT) THEN
         CALL DES_CALC_D_E (A_M, VXF_GS, D_E, IER)
      END IF
!  
      IF(.NOT.DISCRETE_ELEMENT) THEN
         IF (MMAX > 0) CALL CALC_E_E (A_M, MCP, E_E, IER) 
         IF (MMAX > 0) CALL PARTIAL_ELIM_U (U_G, U_S, VXF_GS, A_M, B_M, IER) 
      END IF
!
      CALL ADJUST_A_U_G (A_M, B_M, IER) 
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL ADJUST_A_U_S (A_M, B_M, IER)
      END IF 
      
      
      IF((CALLED.GT.3).AND.(DES_CONTINUUM_COUPLED)) THEN
         UVEL = 1
         CALL GAS_DRAG(A_M, B_M, VXF_GS, IER, UVEL, VVEL, WVEL)
         UVEL = 0
      END IF
!
      IF (MOMENTUM_X_EQ(0)) THEN 
         CALL CALC_RESID_U (U_G, V_G, W_G, A_M, B_M, 0, RESID(RESID_U,0), &
            MAX_RESID(RESID_U,0), IJK_RESID(RESID_U,0), IER) 
         CALL UNDER_RELAX_U (U_G, A_M, B_M, 0, UR_FAC(3), IER) 
!
!        call check_ab_m(a_m, b_m, 0, .false., ier)
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!        write(*,*)
!     &    resid(resid_u, 0), max_resid(resid_u, 0),
!     &    ijk_resid(resid_u, 0)
!
      ENDIF 
      
!
      DO M = 1, MMAX 
         IF (MOMENTUM_X_EQ(M)) THEN 
            CALL CALC_RESID_U (U_S(1,M), V_S(1,M), W_S(1,M), A_M, B_M, M, RESID&
               (RESID_U,M), MAX_RESID(RESID_U,M), IJK_RESID(RESID_U,M), IER) 
            CALL UNDER_RELAX_U (U_S(1,M), A_M, B_M, M, UR_FAC(3), IER) 
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          write(*,*)
!     &      resid(resid_u, m), max_resid(resid_u, m),
!     &      ijk_resid(resid_u, m)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
         ENDIF 
      END DO 

      
      IF (MOMENTUM_X_EQ(0)) THEN 
!        call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K,
!     &    ier)
!
         CALL ADJUST_LEQ (RESID(RESID_U,0), LEQ_IT(3), LEQ_METHOD(3), LEQI, &
            LEQM, IER) 
!
         CALL SOLVE_LIN_EQ ('U_g', U_Gtmp, A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(3), LEQ_TOL(3),IER) 
!        call out_array(u_g, 'u_g')
      ENDIF 
!
      DO M = 1, MMAX 
         IF (MOMENTUM_X_EQ(M)) THEN 
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &    ier)
            CALL ADJUST_LEQ (RESID(RESID_U,M), LEQ_IT(3), LEQ_METHOD(3), LEQI, &
               LEQM, IER) 
!
            CALL SOLVE_LIN_EQ ('U_s', U_Stmp(1,M), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(3), LEQ_TOL(3),IER) 
!          call out_array(u_s(1,m), 'u_s')
         ENDIF 
      END DO 
      
      DO M = 0, MMAX 
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
      END DO 
     
      CALL CONV_DIF_V_G (A_M, B_M, IER) 
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL CONV_DIF_V_S (A_M, B_M, IER) 
      END IF
!
      CALL SOURCE_V_G (A_M, B_M, IER) 
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL SOURCE_V_S (A_M, B_M, IER)
      END IF 

!
      IF (MMAX > 0) CALL VF_GS_Y (VXF_GS, IER)
      IF(.NOT.DISCRETE_ELEMENT) THEN
         IF (MMAX > 0) CALL VF_SS_Y (VXF_SS, IER)    !S. Dartevelle, LANL, Feb.2004
      END IF
!
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL CALC_D_N (A_M, VXF_GS, VXF_SS, D_N, IER)   !S. Dartevelle, LANL, Feb.2004
      ELSE IF(DISCRETE_ELEMENT) THEN
         CALL DES_CALC_D_N (A_M, VXF_GS, D_N, IER)
      END IF
!
      IF(.NOT.DISCRETE_ELEMENT) THEN
         IF (MMAX > 0) CALL CALC_E_N (A_M, MCP, E_N, IER) 
         IF (MMAX > 0) CALL PARTIAL_ELIM_V (V_G, V_S, VXF_GS, A_M, B_M, IER) 
      END IF
!
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
      CALL ADJUST_A_V_G (A_M, B_M, IER)
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
      IF(.NOT.DISCRETE_ELEMENT) THEN
         CALL ADJUST_A_V_S (A_M, B_M, IER) 
      END IF
!
      IF((CALLED.GT.3).AND.(DES_CONTINUUM_COUPLED)) THEN
         VVEL = 1     
         CALL GAS_DRAG(A_M, B_M, VXF_GS, IER, UVEL, VVEL, WVEL)
         VVEL = 0
      END IF
      
      IF (MOMENTUM_Y_EQ(0)) THEN 
         CALL CALC_RESID_V (U_G, V_G, W_G, A_M, B_M, 0, RESID(RESID_V,0), &
            MAX_RESID(RESID_V,0), IJK_RESID(RESID_V,0), IER) 
         CALL UNDER_RELAX_V (V_G, A_M, B_M, 0, UR_FAC(4), IER) 
!
!        call check_ab_m(a_m, b_m, 0, .false., ier)
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!        write(*,*)
!     &    resid(resid_v, 0), max_resid(resid_v, 0),
!     &    ijk_resid(resid_v, 0)
      ENDIF 
!

      DO M = 1, MMAX 
         IF (MOMENTUM_Y_EQ(M)) THEN 
            CALL CALC_RESID_V (U_S(1,M), V_S(1,M), W_S(1,M), A_M, B_M, M, RESID&
               (RESID_V,M), MAX_RESID(RESID_V,M), IJK_RESID(RESID_V,M), IER) 
            CALL UNDER_RELAX_V (V_S(1,M), A_M, B_M, M, UR_FAC(4), IER) 
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          write(*,*)
!     &      resid(resid_v, m), max_resid(resid_v, m),
!     &      ijk_resid(resid_v, m)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
         ENDIF 
      END DO 
      IF (MOMENTUM_Y_EQ(0)) THEN 
!        call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K,
!     &    ier)
!
         CALL ADJUST_LEQ (RESID(RESID_V,0), LEQ_IT(4), LEQ_METHOD(4), LEQI, &
            LEQM, IER) 
!
         CALL SOLVE_LIN_EQ ('V_g', V_Gtmp, A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(4), LEQ_TOL(4),IER) 
!        call out_array(v_g, 'v_g')
      ENDIF 
!
      DO M = 1, MMAX 
         IF (MOMENTUM_Y_EQ(M)) THEN 
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &    ier)
!
            CALL ADJUST_LEQ (RESID(RESID_V,M), LEQ_IT(4), LEQ_METHOD(4), LEQI, &
               LEQM, IER) 
!
            CALL SOLVE_LIN_EQ ('V_s', V_Stmp(1,M), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(4), LEQ_TOL(4),IER) 
!          call out_array(v_s(1,m), 'v_s')
         ENDIF 
      END DO 

!
      IF (DO_K)THEN

        DO M = 0, MMAX 
          CALL INIT_AB_M (A_M, B_M, IJKMAX2, M, IER) 
        END DO 
        CALL CONV_DIF_W_G (A_M, B_M, IER) 
!
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
        IF(.NOT.DISCRETE_ELEMENT) THEN
          CALL CONV_DIF_W_S (A_M, B_M, IER) 
        END IF
!
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!         
        CALL SOURCE_W_G (A_M, B_M, IER) 

!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!     
        IF(.NOT.DISCRETE_ELEMENT) THEN    
          CALL SOURCE_W_S (A_M, B_M, IER) 
        END IF
!
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
! call mfix_exit(myPE)
!
        IF (MMAX > 0) CALL VF_GS_Z (VXF_GS, IER)
        IF(.NOT.DISCRETE_ELEMENT) THEN
          IF (MMAX > 0) CALL VF_SS_Z (VXF_SS, IER)   !S. Dartevelle, LANL, Feb.2004
        END IF
!
        IF(.NOT.DISCRETE_ELEMENT) THEN
          CALL CALC_D_T (A_M, VXF_GS, VXF_SS, D_T, IER)  !S. Dartevelle, LANL, Feb.2004
        ELSE IF(DISCRETE_ELEMENT) THEN
          IF(DIMN.EQ.3) THEN
            CALL DES_CALC_D_T (A_M, VXF_GS, D_T, IER)
          END IF
        END IF
!
        IF(.NOT.DISCRETE_ELEMENT) THEN
          IF (MMAX > 0) CALL CALC_E_T (A_M, MCP, E_T, IER) 
          IF (MMAX > 0) CALL PARTIAL_ELIM_W (W_G, W_S, VXF_GS, A_M, B_M, IER)
        END IF 
!
        CALL ADJUST_A_W_G (A_M, B_M, IER)
        IF(.NOT.DISCRETE_ELEMENT) THEN 
          CALL ADJUST_A_W_S (A_M, B_M, IER)
        END IF 
!
        IF(DIMN.EQ.3) THEN
          IF((CALLED.GT.3).AND.(DES_CONTINUUM_COUPLED)) THEN
            WVEL = 1
            CALL GAS_DRAG(A_M, B_M, VXF_GS, IER, UVEL, VVEL, WVEL)
            WVEL = 0
          END IF
        END IF
                                              
        IF (MOMENTUM_Z_EQ(0)) THEN 
          CALL CALC_RESID_W (U_G, V_G, W_G, A_M, B_M, 0, RESID(RESID_W,0), &
            MAX_RESID(RESID_W,0), IJK_RESID(RESID_W,0), IER) 
          CALL UNDER_RELAX_W (W_G, A_M, B_M, 0, UR_FAC(5), IER) 
!
!        call check_ab_m(a_m, b_m, 0, .false., ier)
!     &      resid(resid_w, 0), max_resid(resid_w, 0),
!     &      ijk_resid(resid_w, 0)
        ENDIF 
!
        DO M = 1, MMAX 
          IF (MOMENTUM_Z_EQ(M)) THEN 
            CALL CALC_RESID_W (U_S(1,M), V_S(1,M), W_S(1,M), A_M, B_M, M, RESID&
               (RESID_W,M), MAX_RESID(RESID_W,M), IJK_RESID(RESID_W,M), IER) 
            CALL UNDER_RELAX_W (W_S(1,M), A_M, B_M, M, UR_FAC(5), IER) 
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          write(*,*)
!     &      resid(resid_w, m), max_resid(resid_w, m),
!     &      ijk_resid(resid_w, m)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
          ENDIF 
        END DO 
        IF (MOMENTUM_Z_EQ(0)) THEN 
!        call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K,
!     &    ier)
!
          CALL ADJUST_LEQ (RESID(RESID_W,0), LEQ_IT(5), LEQ_METHOD(5), LEQI, &
            LEQM, IER) 
!
          CALL SOLVE_LIN_EQ ('W_g', W_Gtmp, A_M, B_M, 0, LEQI, LEQM, &
	                     LEQ_SWEEP(5), LEQ_TOL(5),IER) 
!        call out_array(w_g, 'w_g')
        ENDIF 
!
        DO M = 1, MMAX 
          IF (MOMENTUM_Z_EQ(M)) THEN 
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1, DO_K,
!     &    ier)
!
            CALL ADJUST_LEQ (RESID(RESID_W,M), LEQ_IT(5), LEQ_METHOD(5), LEQI, &
               LEQM, IER) 
!
            CALL SOLVE_LIN_EQ ('W_s', W_Stmp(1,M), A_M, B_M, M, LEQI, LEQM, &
	                     LEQ_SWEEP(5), LEQ_TOL(5),IER) 
!          call out_array(w_s(1,m), 'w_s')
          ENDIF 
        END DO 
      ENDIF
      
!     Now update all velocity components
      DO IJK = ijkstart3, ijkend3
        U_g(IJK) = U_gtmp(IJK)
        V_g(IJK) = V_gtmp(IJK)
        W_g(IJK) = W_gtmp(IJK)
      ENDDO
      DO M = 1, MMAX
        DO IJK = ijkstart3, ijkend3
        U_s(IJK, M) = U_stmp(IJK, M)
        V_s(IJK, M) = V_stmp(IJK, M)
        W_s(IJK, M) = W_stmp(IJK, M)
        ENDDO
      ENDDO

      call unlock_ambm
      call unlock_tmp_array1
      call unlock_tmp_array2

      RETURN  
      END SUBROUTINE SOLVE_VEL_STAR 
      
!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
                                       source_pp_g.f                                                                                       0100644 0002444 0000146 00000025045 10250123040 011666  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_Pp_g(A_m, B_m, IER)
!  Purpose: Determine source terms for Pressure                        C
!  correction equation.  The off-diagonal coefficients are             C
!   positive. The center coefficient and the source vector are         C
!  negative.                                                           C
!  See conv_Pp_g
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SOURCE_PP_G(A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE physprop
      USE fldvar
      USE rxns
      USE run
      USE geometry
      USE indices
      USE pgcor
      USE bc
      USE vshear
      Use xsi_array
      USE compar    
      USE ur_facs 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Convection weighting factors 
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),& 
!                       XSI_t(DIMENSION_3) 
! 
! 
!                      under relaxation factor for pressure
      DOUBLE PRECISION fac 

!                      Indices 
      INTEGER          IJK, IMJK, IJMK, IJKM, I, J, K 
      INTEGER          M, IJKE, IJKW, IJKN, IJKS, IJKT, IJKB 
!
! loezos
      integer incr
 
!                      error message 
      CHARACTER*80     LINE(1) 
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: DROODP_G 
!-----------------------------------------------
 
!// Need to extract i, j, k from ijk_p_g to determine the processor which
!     acts on ijk_p_g to fix the value of pressure
!       ----------------------------------------------------------
!       inline functions for determining i, j, k for global ijk_p_g
!       -----------------------------------------------------------
        integer i_of_g,j_of_g,k_of_g
        integer ijk_p_g_local

      INCLUDE 'function.inc'

        k_of_g(ijk) = int( (ijk-1)/( (imax3-imin3+1)*(jmax3-jmin3+1) ) ) + kmin3
        i_of_g(ijk) = int( ( (ijk-  (k_of_g(ijk)-kmin3)*((imax3-imin3+1)*(jmax3-jmin3+1))) &
                      - 1)/(jmax3-jmin3+1)) + imin3
        j_of_g(ijk) = ijk - (i_of_g(ijk)-imin3)*(jmax3-jmin3+1) - &
                      (k_of_g(ijk)-kmin3)*((imax3-imin3+1)*(jmax3-jmin3+1)) - 1 + jmin3

! loezos
! update to true velocity
      IF (SHEAR) THEN
!$omp parallel do private(IJK) 
	 DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN  
	   V_G(IJK)=V_G(IJK)+VSH(IJK)	
         END IF
        END DO 
 
      END IF

      call lock_xsi_array

!
!  Calculate convection-diffusion fluxes through each of the faces
!
!$omp    parallel do private(IJK, IMJK, IJMK, IJKM, M)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN 
            IMJK = IM_OF(IJK) 
            IJMK = JM_OF(IJK) 
            IJKM = KM_OF(IJK) 
            B_M(IJK,0) = -((-((ROP_G(IJK)-ROP_GO(IJK))*VOL(IJK)*ODT+A_M(IJK,E,0&
               )*U_G(IJK)-A_M(IJK,W,0)*U_G(IMJK)+A_M(IJK,N,0)*V_G(IJK)-A_M(IJK,&
               S,0)*V_G(IJMK)+A_M(IJK,T,0)*W_G(IJK)-A_M(IJK,B,0)*W_G(IJKM)))+&
               SUM_R_G(IJK)*VOL(IJK)) 
            A_M(IJK,E,0) = A_M(IJK,E,0)*D_E(IJK,0) 
            A_M(IJK,W,0) = A_M(IJK,W,0)*D_E(IMJK,0) 
            A_M(IJK,N,0) = A_M(IJK,N,0)*D_N(IJK,0) 
            A_M(IJK,S,0) = A_M(IJK,S,0)*D_N(IJMK,0) 
            A_M(IJK,T,0) = A_M(IJK,T,0)*D_T(IJK,0) 
            A_M(IJK,B,0) = A_M(IJK,B,0)*D_T(IJKM,0) 
            DO M = 1, MMAX 
               IF (.NOT.CLOSE_PACKED(M)) THEN 
                  B_M(IJK,0) = B_M(IJK,0) - ((-((ROP_S(IJK,M)-ROP_SO(IJK,M))*&
                     VOL(IJK)*ODT+A_M(IJK,E,M)*U_S(IJK,M)-A_M(IJK,W,M)*U_S(IMJK&
                     ,M)+A_M(IJK,N,M)*V_S(IJK,M)-A_M(IJK,S,M)*V_S(IJMK,M)+A_M(&
                     IJK,T,M)*W_S(IJK,M)-A_M(IJK,B,M)*W_S(IJKM,M)))+SUM_R_S(IJK&
                     ,M)*VOL(IJK)) 
                  A_M(IJK,E,0) = A_M(IJK,E,0) + A_M(IJK,E,M)*D_E(IJK,M) 
                  A_M(IJK,W,0) = A_M(IJK,W,0) + A_M(IJK,W,M)*D_E(IMJK,M) 
                  A_M(IJK,N,0) = A_M(IJK,N,0) + A_M(IJK,N,M)*D_N(IJK,M) 
                  A_M(IJK,S,0) = A_M(IJK,S,0) + A_M(IJK,S,M)*D_N(IJMK,M) 
                  A_M(IJK,T,0) = A_M(IJK,T,0) + A_M(IJK,T,M)*D_T(IJK,M) 
                  A_M(IJK,B,0) = A_M(IJK,B,0) + A_M(IJK,B,M)*D_T(IJKM,M) 
               ENDIF 
            END DO 
            A_M(IJK,0,0) = -(A_M(IJK,E,0)+A_M(IJK,W,0)+A_M(IJK,N,0)+A_M(IJK,S,0&
               )+A_M(IJK,T,0)+A_M(IJK,B,0)) 
!
            IF (ABS(A_M(IJK,0,0)) < SMALL_NUMBER) THEN 
               IF (ABS(B_M(IJK,0)) < SMALL_NUMBER) THEN 
                  A_M(IJK,0,0) = -ONE 
                  B_M(IJK,0) = ZERO 
               ELSE 
!$omp             critical
                  WRITE (LINE, '(A,I6,A,I1,A,G12.5)') 'Error: At IJK = ', IJK, &
                     ' M = ', 0, ' A = 0 and b = ', B_M(IJK,0) 
                  CALL WRITE_ERROR ('SOURCE_Pp_g', LINE, 1) 
!$omp             end critical
               ENDIF 
            ENDIF 
!
         ELSE 
            A_M(IJK,E,0) = ZERO 
            A_M(IJK,W,0) = ZERO 
            A_M(IJK,N,0) = ZERO 
            A_M(IJK,S,0) = ZERO 
            A_M(IJK,T,0) = ZERO 
            A_M(IJK,B,0) = ZERO 
            A_M(IJK,0,0) = -ONE 
            B_M(IJK,0) = ZERO 
         ENDIF 
      END DO 

! loezos
      IF (SHEAR) THEN
!$omp parallel do private(IJK) 
	 DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN  
	   V_G(IJK)=V_G(IJK)-VSH(IJK)	
         END IF
       END DO 
      END IF
! loezos

      IF (RO_G0 == UNDEFINED) THEN 
        fac = UR_FAC(1)  !since p_g = p_g* + ur_fac * pp_g

! loezos
	incr=0		
! loezos
         CALL CALC_XSI(DISCRETIZE(1),ROP_G,U_G,V_G,W_G,XSI_E,XSI_N,XSI_T,incr) 
	
!$omp    parallel do                                                     &
!$omp&   private(IJK,I,J,K,                                       &
!$omp&            IMJK,IJMK,IJKM,IJKE,IJKW,IJKN,IJKS,IJKT,IJKB)
         DO IJK = ijkstart3, ijkend3 
            IF (FLUID_AT(IJK)) THEN 
               I = I_OF(IJK) 
               J = J_OF(IJK) 
               K = K_OF(IJK) 
               IMJK = IM_OF(IJK) 
               IJMK = JM_OF(IJK) 
               IJKM = KM_OF(IJK) 
               IJKE = EAST_OF(IJK) 
               IJKW = WEST_OF(IJK) 
               IJKN = NORTH_OF(IJK) 
               IJKS = SOUTH_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKB = BOTTOM_OF(IJK) 
!
               A_M(IJK,0,0) = A_M(IJK,0,0) - fac*DROODP_G(RO_G(IJK),P_G(IJK))*EP_G(&
                  IJK)*((ONE - XSI_E(IJK))*U_G(IJK)*AYZ(IJK)-XSI_E(IMJK)*U_G(&
                  IMJK)*AYZ(IMJK)+(ONE-XSI_N(IJK))*V_G(IJK)*AXZ(IJK)-XSI_N(IJMK&
                  )*V_G(IJMK)*AXZ(IJMK)+VOL(IJK)*ODT) 
!
               A_M(IJK,E,0) = A_M(IJK,E,0) - EP_G(IJKE)*fac*DROODP_G(RO_G(IJKE),P_G&
                  (IJKE))*XSI_E(IJK)*U_G(IJK)*AYZ(IJK) 
               A_M(IJK,W,0) = A_M(IJK,W,0) + EP_G(IJKW)*fac*DROODP_G(RO_G(IJKW),P_G&
                  (IJKW))*(ONE - XSI_E(IMJK))*U_G(IMJK)*AYZ(IMJK) 
               A_M(IJK,N,0) = A_M(IJK,N,0) - EP_G(IJKN)*fac*DROODP_G(RO_G(IJKN),P_G&
                  (IJKN))*XSI_N(IJK)*V_G(IJK)*AXZ(IJK) 
               A_M(IJK,S,0) = A_M(IJK,S,0) + EP_G(IJKS)*fac*DROODP_G(RO_G(IJKS),P_G&
                  (IJKS))*(ONE - XSI_N(IJMK))*V_G(IJMK)*AXZ(IJMK) 
               IF (DO_K) THEN 
                  A_M(IJK,0,0) = A_M(IJK,0,0) - fac*DROODP_G(RO_G(IJK),P_G(IJK))*&
                     EP_G(IJK)*((ONE - XSI_T(IJK))*W_G(IJK)*AXY(IJK)-XSI_T(IJKM&
                     )*W_G(IJKM)*AXY(IJKM)) 
                  A_M(IJK,T,0) = A_M(IJK,T,0) - EP_G(IJKT)*fac*DROODP_G(RO_G(IJKT),&
                     P_G(IJKT))*XSI_T(IJK)*W_G(IJK)*AXY(IJK) 
                  A_M(IJK,B,0) = A_M(IJK,B,0) + EP_G(IJKB)*fac*DROODP_G(RO_G(IJKB),&
                     P_G(IJKB))*(ONE - XSI_T(IJKM))*W_G(IJKM)*AXY(IJKM) 
               ENDIF 
!
            ENDIF 
         END DO 
      ENDIF 
!
!
!  Specify P' to zero at a certain location for incompressible flows and
!  cyclic boundary conditions.
!
!// Parallel implementation of fixing a pressure at a point
   I = I_OF_G(IJK_P_G)
   J = J_OF_G(IJK_P_G)
   K = K_OF_G(IJK_P_G)

   IF(IS_ON_myPE_OWNS(I,J,K)) THEN
      IF (IJK_P_G /= UNDEFINED_I) THEN 
         IJK_P_G_LOCAL = FUNIJK(I,J,K)
         A_M(IJK_P_G_LOCAL,E,0) = ZERO 
         A_M(IJK_P_G_LOCAL,W,0) = ZERO 
         A_M(IJK_P_G_LOCAL,N,0) = ZERO 
         A_M(IJK_P_G_LOCAL,S,0) = ZERO 
         A_M(IJK_P_G_LOCAL,T,0) = ZERO 
         A_M(IJK_P_G_LOCAL,B,0) = ZERO 
         A_M(IJK_P_G_LOCAL,0,0) = -ONE 
         B_M(IJK_P_G_LOCAL,0) = ZERO 
      ENDIF 
   ENDIF
!
      call unlock_xsi_array

      RETURN  
      END SUBROUTINE SOURCE_PP_G 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           source_u_g.f                                                                                        0100644 0002444 0000146 00000076717 10250045301 011532  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_U_g(A_m, B_m, IER)                              C
!  Purpose: Determine source terms for U_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.    The center      C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SOURCE_U_G(A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_g
      USE rxns
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is
      USE tau_g
      USE bc
      USE compar    
      USE sendrecv  
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I, IJK, IJKE, IPJK, IJKM, IPJKM 
! 
!                      Phase index 
      INTEGER          M 
! 
!                      Internal surface 
      INTEGER          ISV 
! 
!                      Pressure at east cell 
      DOUBLE PRECISION PgE 
! 
!                      Average volume fraction 
      DOUBLE PRECISION EPGA 
! 
!                      Average density 
      DOUBLE PRECISION ROPGA, ROGA 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Average viscosity 
      DOUBLE PRECISION MUGA 
! 
!                      Average viscosity 
      DOUBLE PRECISION EPMU_gte, EPMU_gbe, EPMUGA 
! 
!                      Average W_g 
      DOUBLE PRECISION Wge 
! 
!                      Average dW/Xdz 
      DOUBLE PRECISION dWoXdz 
! 
!                      Source terms (Surface) 
      DOUBLE PRECISION Sdp 
! 
!                      Source terms (Volumetric) 
      DOUBLE PRECISION V0, Vpm, Vmt, Vbf, Vcf, Vtza 
! 
!                      error message 
      CHARACTER*80     LINE 
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!
      M = 0 
      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN  
!
!$omp    parallel do private(I, IJK, IJKE, IJKM, IPJK, IPJKM,     &
!$omp&                  ISV, Sdp, V0, Vpm, Vmt, Vbf,              &
!$omp&                  Vcf, EPMUGA, VTZA, WGE, PGE, ROGA,        &
!$omp&                  MUGA, ROPGA, EPGA )
      DO IJK = ijkstart3, ijkend3 
         I = I_OF(IJK) 
         IJKE = EAST_OF(IJK) 
         IJKM = KM_OF(IJK) 
         IPJK = IP_OF(IJK) 
         IPJKM = IP_OF(IJKM) 
         EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I) 
         IF (IP_AT_E(IJK)) THEN 
            A_M(IJK,E,M) = ZERO 
            A_M(IJK,W,M) = ZERO 
            A_M(IJK,N,M) = ZERO 
            A_M(IJK,S,M) = ZERO 
            A_M(IJK,T,M) = ZERO 
            A_M(IJK,B,M) = ZERO 
            A_M(IJK,0,M) = -ONE 
            B_M(IJK,M) = ZERO 
!
!       dilute flow
         ELSE IF (EPGA <= DIL_EP_S) THEN 
            A_M(IJK,E,M) = ZERO 
            A_M(IJK,W,M) = ZERO 
            A_M(IJK,N,M) = ZERO 
            A_M(IJK,S,M) = ZERO 
            A_M(IJK,T,M) = ZERO 
            A_M(IJK,B,M) = ZERO 
            A_M(IJK,0,M) = -ONE 
            B_M(IJK,M) = ZERO 
!
            IF (EP_G(WEST_OF(IJK)) > DIL_EP_S) THEN 
               A_M(IJK,W,M) = ONE 
            ELSE IF (EP_G(EAST_OF(IJK)) > DIL_EP_S) THEN 
               A_M(IJK,E,M) = ONE 
            ELSE 
               B_M(IJK,M) = -U_G(IJK) 
            ENDIF 
         ELSE 
!
!       Surface forces
!
!         Pressure term
            PGE = P_G(IJKE) 
            IF (CYCLIC_X_PD) THEN 
               IF (IMAP(I_OF(IJK)).EQ.IMAX1) PGE = P_G(IJKE) - DELP_X 
            ENDIF 
            IF (MODEL_B) THEN 
               SDP = -P_SCALE*(PGE - P_G(IJK))*AYZ(IJK) 
!
            ELSE 
               SDP = -P_SCALE*EPGA*(PGE - P_G(IJK))*AYZ(IJK) 
!
            ENDIF 
!
!       Volumetric forces
            ROPGA = HALF * (VOL(IJK)*ROP_G(IJK) + VOL(IPJK)*ROP_G(IJKE))/VOL_U(IJK)
            ROGA  = HALF * (VOL(IJK)*RO_G(IJK) + VOL(IPJK)*RO_G(IJKE))/VOL_U(IJK) 
!
!         Previous time step
            V0 = HALF * (VOL(IJK)*ROP_GO(IJK) + VOL(IPJK)*ROP_GO(IJKE))*ODT/VOL_U(IJK) 
!
!         pressure drop through porous media
            IF (SIP_AT_E(IJK)) THEN 
               ISV = IS_ID_AT_E(IJK) 
               MUGA = AVG_X(MU_G(IJK),MU_G(IJKE),I) 
               VPM = MUGA/IS_PC(ISV,1) 
               IF (IS_PC(ISV,2) /= ZERO) VPM = VPM + HALF*IS_PC(ISV,2)*ROPGA*ABS(&
                  U_G(IJK)) 
            ELSE 
               VPM = ZERO 
            ENDIF 
!
!         Interphase mass transfer
            VMT = HALF * (VOL(IJK)*SUM_R_G(IJK) + VOL(IPJK)*SUM_R_G(IJKE))/VOL_U(IJK) 
!
!         Body force
            IF (MODEL_B) THEN 
               VBF = ROGA*BFX_G(IJK) 
!
            ELSE                                 !Model A 
               VBF = ROPGA*BFX_G(IJK) 
!
            ENDIF 
!
!         Special terms for cylindrical coordinates
            IF (CYLINDRICAL) THEN 
!
!           centrifugal force
               WGE = AVG_X(HALF*(W_G(IJK)+W_G(IJKM)),HALF*(W_G(IPJK)+W_G(IPJKM)&
                  ),I) 
               VCF = ROPGA*WGE**2*OX_E(I) 
!
!           -(2mu/x)*(u/x) part of Tau_zz/X
               EPMUGA = AVG_X(MU_GT(IJK),MU_GT(IJKE),I) 
               VTZA = 2.*EPMUGA*OX_E(I)*OX_E(I) 
            ELSE 
               VCF = ZERO 
               VTZA = ZERO 
            ENDIF 
!
!         Collect the terms
            A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(IJK,S,M&
               )+A_M(IJK,T,M)+A_M(IJK,B,M)+(V0+VPM+ZMAX(VMT)+VTZA)*VOL_U(IJK)) 
            B_M(IJK,M) = -(SDP + TAU_U_G(IJK)+((V0+ZMAX((-VMT)))*U_GO(IJK)+VBF+&
               VCF)*VOL_U(IJK))+B_M(IJK,M)
	ENDIF 
      END DO 
      CALL SOURCE_U_G_BC (A_M, B_M, IER) 
!

      RETURN  
      END SUBROUTINE SOURCE_U_G 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_U_g_BC(A_m, B_m, IER)                           C
!  Purpose: Determine source terms for U_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.    The center      C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SOURCE_U_G_BC(A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_g 
      USE rxns 
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is 
      USE tau_g 
      USE bc
      USE output
      USE compar    
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Boundary condition 
      INTEGER          L 
! 
!                      Indices 
      INTEGER          I,  J, K, IM, I1, I2, J1, J2, K1, K2, IJK,& 
                       JM, KM, IJKW, IMJK, IP, IPJK 
! 
!                      Solids phase 
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Turbulent shear stress
      DOUBLE PRECISION  W_F_Slip
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!
      M = 0 
!
!
!  Set the default boundary conditions
!
      IF (DO_K) THEN 
         K1 = 1 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    
               IJK = FUNIJK(I1,J1,K1) 	       
               IF (NS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = -ONE 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSE IF (FS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ONE 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            END DO 
         END DO 
	 
         K1 = KMAX2 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    	    
               IJK = FUNIJK(I1,J1,K1) 
               IF (NS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = -ONE 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSE IF (FS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ONE 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            END DO 
         END DO 
      ENDIF 
!
      J1 = 1 
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3
  	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    	    	 
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = -ONE 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSE IF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ONE 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      
      J1 = JMAX2 
      DO K1 = kmin3, kmax3      
         DO I1 = imin3, imax3 
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	 
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = -ONE 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSE IF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ONE 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
            IF (BC_TYPE(L) == 'NO_SLIP_WALL' .AND. .NOT. K_Epsilon) THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     
                        IJK = FUNIJK(I,J,K) 
                        IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                           A_M(IJK,N,M) = -ONE 
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                           A_M(IJK,S,M) = -ONE 
                        ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                           A_M(IJK,T,M) = -ONE 
                        ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                           A_M(IJK,B,M) = -ONE 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ELSE IF (BC_TYPE(L) == 'FREE_SLIP_WALL' .AND. .NOT. K_Epsilon) THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IJK = FUNIJK(I,J,K) 
                        IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                           A_M(IJK,N,M) = ONE 
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                           A_M(IJK,S,M) = ONE 
                        ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                           A_M(IJK,T,M) = ONE 
                        ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                           A_M(IJK,B,M) = ONE 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ELSE IF (BC_TYPE(L) == 'PAR_SLIP_WALL' .AND. .NOT. K_Epsilon) THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IJK = FUNIJK(I,J,K) 
                        IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                        JM = JM1(J) 
                        KM = KM1(K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN
			   IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,N,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_UW_G(L) 
                           ELSE 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODY_N(J)) 
                              A_M(IJK,N,M) = -(HALF*BC_HW_G(L)-ODY_N(J)) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L) 
                           ENDIF 
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,S,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_UW_G(L) 
                           ELSE 
                              A_M(IJK,S,M) = -(HALF*BC_HW_G(L)-ODY_N(JM)) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODY_N(JM)) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L) 
                           ENDIF 
                        ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN  
                           IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,T,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_UW_G(L) 
                           ELSE 
                              A_M(IJK,0,M)=-(HALF*BC_HW_G(L)+ODZ_T(K)*OX_E(I)) 
                              A_M(IJK,T,M)=-(HALF*BC_HW_G(L)-ODZ_T(K)*OX_E(I)) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L) 
                           ENDIF 
                        ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN   
                           IF (BC_HW_G(L) == UNDEFINED) THEN 
                              A_M(IJK,B,M) = -HALF 
                              A_M(IJK,0,M) = -HALF 
                              B_M(IJK,M) = -BC_UW_G(L) 
                           ELSE 
                              A_M(IJK,B,M) = -(HALF*BC_HW_G(L)-ODZ_T(KM)*OX_E(I&
                                 )) 
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODZ_T(KM)*OX_E(I&
                                 )) 
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L) 
                           ENDIF 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
! wall functions for U-momentum are specify in this section of the code
            ELSE IF (BC_TYPE(L) == 'PAR_SLIP_WALL'   .OR.  &
	             BC_TYPE(L) == 'NO_SLIP_WALL'    .OR.  &
		     BC_TYPE(L) == 'FREE_SLIP_WALL'  .AND. &
		     K_Epsilon                            )THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IJK = FUNIJK(I,J,K) 
                        IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                        JM = JM1(J) 
                        KM = KM1(K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = ZERO 
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN  
			     CALL Wall_Function(IJK,NORTH_OF(IJK),ODY_N(J),W_F_Slip)
                             A_M(IJK,N,M) = W_F_Slip
                             A_M(IJK,0,M) = -ONE 
                             IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_UW_G(L)
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN
			     CALL Wall_Function(IJK,SOUTH_OF(IJK),ODY_N(JM),W_F_Slip)
                             A_M(IJK,S,M) = W_F_Slip
                             A_M(IJK,0,M) = -ONE 
                             IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_UW_G(L)
                        ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
			     CALL Wall_Function(IJK,TOP_OF(IJK),ODZ_T(K)*OX_E(I),W_F_Slip)
                             A_M(IJK,T,M) = W_F_Slip
                             A_M(IJK,0,M) = -ONE
                             IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_UW_G(L)
                        ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN
			     CALL Wall_Function(IJK,BOTTOM_OF(IJK),ODZ_T(KM)*OX_E(I),W_F_Slip)
                             A_M(IJK,B,M) = W_F_Slip
                             A_M(IJK,0,M) = -ONE
                             IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_UW_G(L) 
                        ENDIF 
                     END DO 
                  END DO 
               END DO
! end of wall functions
            ELSE IF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'W') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ONE 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
               ENDIF 
            ELSE IF (BC_TYPE(L) == 'OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'W') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ONE 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
!
                           IM = IM1(I) 
                           IMJK = IM_OF(IJK) 
                           A_M(IMJK,E,M) = ZERO 
                           A_M(IMJK,W,M) = X_E(IM)/X_E(IM1(IM)) 
                           A_M(IMJK,N,M) = ZERO 
                           A_M(IMJK,S,M) = ZERO 
                           A_M(IMJK,T,M) = ZERO 
                           A_M(IMJK,B,M) = ZERO 
                           A_M(IMJK,0,M) = -ONE 
                           B_M(IMJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
               ELSE IF (BC_PLANE(L) == 'E') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IJK = FUNIJK(I,J,K) 
!
                           IP = IP1(I) 
                           IPJK = IP_OF(IJK) 
                           A_M(IPJK,E,M) = X_E(IP)/X_E(I) 
                           A_M(IPJK,W,M) = ZERO 
                           A_M(IPJK,N,M) = ZERO 
                           A_M(IPJK,S,M) = ZERO 
                           A_M(IPJK,T,M) = ZERO 
                           A_M(IPJK,B,M) = ZERO 
                           A_M(IPJK,0,M) = -ONE 
                           B_M(IPJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
               ENDIF 
            ELSE 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IJK = FUNIJK(I,J,K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = -U_G(IJK) 
                        IF (BC_PLANE(L) == 'W') THEN 
                           IJKW = WEST_OF(IJK) 
                           A_M(IJKW,E,M) = ZERO 
                           A_M(IJKW,W,M) = ZERO 
                           A_M(IJKW,N,M) = ZERO 
                           A_M(IJKW,S,M) = ZERO 
                           A_M(IJKW,T,M) = ZERO 
                           A_M(IJKW,B,M) = ZERO 
                           A_M(IJKW,0,M) = -ONE 
                           B_M(IJKW,M) = -U_G(IJKW) 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE SOURCE_U_G_BC 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Wall_Function(IJK1,IJK2,ODX_WF,W_F_Slip)               C
!  Purpose: Calculate Slip velocity using wall functions               C
!                                                                      C
!  Author: S. Benyahia                                Date: MAY-13-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE Wall_Function(IJK1,IJK2,ODX_WF,W_F_Slip)
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE physprop 
      USE fldvar
      USE visc_g  
      USE geometry 
      USE indices 
      USE bc
      USE compar 
      USE turb        
      USE mpi_utility 
      IMPLICIT NONE
!
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

!                      IJK indices for wall cell and fluid cell
      INTEGER          IJK1, IJK2

!                      ODX_WF: 1/dx, and W_F_Slip: value of turb. shear stress at walls
      DOUBLE PRECISION ODX_WF, W_F_Slip

!                      C_mu and Kappa are constants in turb. viscosity and Von Karmen const.
      DOUBLE PRECISION C_mu, Kappa
!-----------------------------------------------
!
!
	C_mu = 0.09D+0
	Kappa = 0.42D+0
	
		W_F_Slip = (ONE - ONE/ODX_WF* RO_g(IJK2)*C_mu**0.25   &
			   *SQRT(K_Turb_G(IJK2))/MU_gT(IJK2)	      &
			   *Kappa/LOG(9.81D+0/(ODX_WF*2.D+0)*         &
			    RO_g(IJK2)*C_mu**0.25*                    &
			   SQRT(K_Turb_G(IJK2))/MU_g(IJK2)))
!
      RETURN  
      END SUBROUTINE Wall_Function

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3      
!// 360 Check if i,j,k resides on current processor
                                                 source_u_s.f                                                                                        0100644 0002444 0000146 00000113042 10250046017 011532  0                                                                                                    ustar   msyaml                                                                                                                                                                                                                                                 !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_U_s(A_m, B_m, IER)                              C
!  Purpose: Determine source terms for U_s momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.  The center        C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose: Allow for partial-slip boundary conditions proposed by     C
!           by Johnson & Jackson (1987) if the Granular Temperature    C
!           equation is used.                                          C
!  Author: K. Agrawal, Princeton University           Date: 24-JAN-98  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SOURCE_U_S(A_M, B_M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is
      USE tau_s
      USE bc
      USE compar    
      USE sendrecv  
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I, IJK, IJKE, IJKM, IPJK, IPJKM 
! 
!                      Phase index 
      INTEGER          M,MM 
      DOUBLE PRECISION   SUM_EPS_CP
! 
!                      Internal surface number 
      INTEGER          ISV 
! 
!                      Pressure at east cell 
      DOUBLE PRECISION PgE 
! 
!                      average volume fraction 
      DOUBLE PRECISION EPSA 
! 
!                      Average density 
      DOUBLE PRECISION ROPSA 
! 
!                      Average density difference 
      DOUBLE PRECISION dro1, dro2, droa 
! 
!                      Average quantities 
      DOUBLE PRECISION wse, EPMUGA 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Source terms (Surface) 
      DOUBLE PRECISION Sdp, Sdps 
! 
!                      Source terms (Volumetric) 
      DOUBLE PRECISION V0, Vmt, Vbf, Vcf, Vtza 
! 
!                      error message 
      CHARACTER*80     LINE(2) 
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!
      DO M = 1, MMAX 
         IF (MOMENTUM_X_EQ(M)) THEN 
!
!$omp  parallel do private( IJK, IJKE, ISV, Sdp, Sdps, V0, Vmt, Vbf, &
!$omp&  I,PGE,DRO1,DRO2,DROA, IJKM,IPJK,IPJKM,  WSE,VCF,EPMUGA,VTZA, &
!$omp&  EPSA, ROPSA, LINE,SUM_EPS_CP,MM) &
!$omp&  schedule(static)
            DO IJK = ijkstart3, ijkend3 
!
!           Wall or impermeable internal surface
               I = I_OF(IJK) 
               IJKE = EAST_OF(IJK) 
               IJKM = KM_OF(IJK) 
               IPJK = IP_OF(IJK) 
               IPJKM = IP_OF(IJKM) 
               EPSA = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I) 
               IF (IP_AT_E(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSE IF (SIP_AT_E(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  ISV = IS_ID_AT_E(IJK) 
                  B_M(IJK,M) = -IS_VEL_S(ISV,M) 
!
!           dilute flow
               ELSE IF (EPSA <= DIL_EP_S) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
!
                  IF (EP_S(WEST_OF(IJK),M) > DIL_EP_S) THEN 
                     A_M(IJK,W,M) = ONE 
                  ELSE IF (EP_S(EAST_OF(IJK),M) > DIL_EP_S) THEN 
                     A_M(IJK,E,M) = ONE 
                  ELSE 
                     B_M(IJK,M) = -U_S(IJK,M) 
                  ENDIF 
!
!           Normal case
               ELSE 
!
!           Surface forces
!
!             Pressure terms
                  PGE = P_G(IJKE) 
                  IF (CYCLIC_X_PD) THEN 
                     IF (CYCLIC_AT_E(IJK)) PGE = P_G(IJKE) - DELP_X 
                  ENDIF 
                  IF (MODEL_B) THEN 
                     SDP = ZERO 
!
                  ELSE 
                     SDP = -P_SCALE*EPSA*(PGE - P_G(IJK))*AYZ(IJK) 
!
                  ENDIF 
!
                  IF (CLOSE_PACKED(M)) THEN 
		     SUM_EPS_CP=0.0 
		     DO MM=1,MMAX
		       IF (CLOSE_PACKED(MM))&
			     SUM_EPS_CP=SUM_EPS_CP+EP_S(IJK,MM)
		     END DO
		     SUM_EPS_CP = Max(SUM_EPS_CP, small_number)
                     SDPS = -((P_S(IJKE,M)-P_S(IJK,M))+(EP_S(IJK,M)/SUM_EPS_CP)*&
			(P_STAR(IJKE)-P_STAR(IJK&
                        )))*AYZ(IJK) 
                  ELSE 
                     SDPS = -(P_S(IJKE,M)-P_S(IJK,M))*AYZ(IJK) 
                  ENDIF 
!
!             Shear stress terms
!
!           Volumetric forces
                  ROPSA = HALF * (VOL(IJK)*ROP_S(IJK,M) + VOL(IPJK)*ROP_S(IJKE,M))/VOL_U(IJK) 
!
!             Previous time step
                  V0 = HALF * (VOL(IJK)*ROP_SO(IJK,M) + VOL(IPJK)*ROP_SO(IJKE,M))*ODT/VOL_U(IJK) 
!
!             Interphase mass transfer
                  VMT = HALF * (VOL(IJK)*SUM_R_S(IJK,M) + VOL(IPJK)*SUM_R_S(IJKE,M))/VOL_U(IJK) 
!
!             Body force
                  IF (MODEL_B) THEN 
                     DRO1 = (RO_S(M)-RO_G(IJK))*EP_S(IJK,M) 
                     DRO2 = (RO_S(M)-RO_G(IJKE))*EP_S(IJKE,M) 
                     DROA = AVG_X(DRO1,DRO2,I) 
!
                     VBF = DROA*BFX_S(IJK,M) 
!
                  ELSE 
                     VBF = ROPSA*BFX_S(IJK,M) 
                  ENDIF 
!
!
!           Special terms for cylindrical coordinates
                  IF (CYLINDRICAL) THEN 
!
!             centrifugal force
                     WSE = AVG_X(HALF*(W_S(IJK,M)+W_S(IJKM,M)),HALF*(W_S(IPJK,M&
                        )+W_S(IPJKM,M)),I) 
                     VCF = ROPSA*WSE**2*OX_E(I) 
!
!             -(2mu/x)*(u/x) part of Tau_zz/X
                     EPMUGA = AVG_X(MU_S(IJK,M),MU_S(IJKE,M),I) 
                     VTZA = 2.*EPMUGA*OX_E(I)*OX_E(I) 
                  ELSE 
                     VCF = ZERO 
                     VTZA = ZERO 
                  ENDIF 
!
!             Collect the terms
                  A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(&
                     IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M)+(V0+ZMAX(VMT)+VTZA)*&
                     VOL_U(IJK)) 
                  B_M(IJK,M) = -(SDP + SDPS + TAU_U_S(IJK,M)+((V0+ZMAX((-VMT)))&
                     *U_SO(IJK,M)+VBF+VCF)*VOL_U(IJK))+B_M(IJK,M) 
               ENDIF 
            END DO 
            CALL SOURCE_U_S_BC (A_M, B_M, M, IER) 
         ENDIF 

      END DO 

      
      RETURN  
      END SUBROUTINE SOURCE_U_S 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_U_s_BC(A_m, B_m, M, IER)                        C
!  Purpose: Determine source terms for U_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.    The center      C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SOURCE_U_S_BC(A_M, B_M, M, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns 
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is 
      USE tau_s 
      USE bc
      USE output
      USE compar   
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Boundary condition 
      INTEGER          L 
! 
!                      Indices 
      INTEGER          I,  J, K, IM, I1, I2, J1, J2, K1, K2, IJK,& 
                       JM, KM, IJKW, IMJK, IPJK, IP 
! 
!                      Solids phase 
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!
!
!  Set the default boundary conditions
!
      IF (DO_K) THEN 
         K1 = 1 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 	 
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    
               IJK = FUNIJK(I1,J1,K1) 
               IF (NS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = -ONE 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSE IF (FS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ONE 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            END DO 
         END DO 
	 
         K1 = KMAX2 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 	 
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	 
               IJK = FUNIJK(I1,J1,K1) 
               IF (NS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = -ONE 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSE IF (FS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ONE 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            END DO 
         END DO 
      ENDIF 
      J1 = 1 
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    	    	 
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = -ONE 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSE IF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ONE 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      J1 = JMAX2 
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	    	    	 
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = -ONE 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSE IF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ONE 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
            IF (BC_TYPE(L) == 'NO_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
!
               IF (BC_JJ_PS(L) == 0) THEN 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IJK = FUNIJK(I,J,K) 
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ZERO 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                           IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                              A_M(IJK,N,M) = -ONE 
                           ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                              A_M(IJK,S,M) = -ONE 
                           ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                              A_M(IJK,T,M) = -ONE 
                           ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                              A_M(IJK,B,M) = -ONE 
                           ENDIF 
                        END DO 
                     END DO 
                  END DO 
		  
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_U_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!
               ENDIF 
!
            ELSE IF (BC_TYPE(L) == 'FREE_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
!
               IF (BC_JJ_PS(L) == 0) THEN 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IJK = FUNIJK(I,J,K) 
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ZERO 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                           IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                              A_M(IJK,N,M) = ONE 
                           ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                              A_M(IJK,S,M) = ONE 
                           ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                              A_M(IJK,T,M) = ONE 
                           ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                              A_M(IJK,B,M) = ONE 
                           ENDIF 
                        END DO 
                     END DO 
                  END DO 
		  
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_U_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!
               ENDIF 
!
            ELSE IF (BC_TYPE(L) == 'PAR_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
!
               IF (BC_JJ_PS(L) == 0) THEN 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     			
                           IJK = FUNIJK(I,J,K) 
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                           JM = JM1(J) 
                           KM = KM1(K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ZERO 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                           IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,N,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODY_N(J)) 
                                 A_M(IJK,N,M) = -(HALF*BC_HW_S(L,M)-ODY_N(J)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,S,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,S,M) = -(HALF*BC_HW_S(L,M)-ODY_N(JM)) 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODY_N(JM)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,T,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODZ_T(K)*&
                                    OX_E(I)) 
                                 A_M(IJK,T,M) = -(HALF*BC_HW_S(L,M)-ODZ_T(K)*&
                                    OX_E(I)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,B,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,B,M) = -(HALF*BC_HW_S(L,M)-ODZ_T(KM)*&
                                    OX_E(I)) 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODZ_T(KM)*&
                                    OX_E(I)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ENDIF 
                        END DO 
                     END DO 
                  END DO 
		  
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_U_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!
               ENDIF 
            ELSE IF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'W') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE			
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ONE 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
		  
               ENDIF 
            ELSE IF (BC_TYPE(L) == 'OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'W') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ONE 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
!
                           IM = IM1(I) 
                           IMJK = IM_OF(IJK) 
                           A_M(IMJK,E,M) = ZERO 
                           A_M(IMJK,W,M) = X_E(IM)/X_E(IM1(IM)) 
                           A_M(IMJK,N,M) = ZERO 
                           A_M(IMJK,S,M) = ZERO 
                           A_M(IMJK,T,M) = ZERO 
                           A_M(IMJK,B,M) = ZERO 
                           A_M(IMJK,0,M) = -ONE 
                           B_M(IMJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
		  
               ELSE IF (BC_PLANE(L) == 'E') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     

                           IJK = FUNIJK(I,J,K) 
!
                           IP = IP1(I) 
                           IPJK = IP_OF(IJK) 
                           A_M(IPJK,E,M) = X_E(IP)/X_E(I) 
                           A_M(IPJK,W,M) = ZERO 
                           A_M(IPJK,N,M) = ZERO 
                           A_M(IPJK,S,M) = ZERO 
                           A_M(IPJK,T,M) = ZERO 
                           A_M(IPJK,B,M) = ZERO 
                           A_M(IPJK,0,M) = -ONE 
                           B_M(IPJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
		  
               ENDIF 
            ELSE 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
               	        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     
		     
                        IJK = FUNIJK(I,J,K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = -U_S(IJK,M) 
                        IF (BC_PLANE(L) == 'W') THEN 
                           IJKW = WEST_OF(IJK) 
                           A_M(IJKW,E,M) = ZERO 
                           A_M(IJKW,W,M) = ZERO 
                           A_M(IJKW,N,M) = ZERO 
                           A_M(IJKW,S,M) = ZERO 
                           A_M(IJKW,T,M) = ZERO 
                           A_M(IJKW,B,M) = ZERO 
                           A_M(IJKW,0,M) = -ONE 
                           B_M(IJKW,M) = -U_S(IJKW,M) 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
	       
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE SOURCE_U_S_BC 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: JJ_BC_U_s(I1, I2, J1, J2, K1, K2, L, M, A_m, b_m)      C
!  Purpose: Implement Johnson and Jackson boundary condition           C
!                                                                      C
!  Author: K. Agrawal & A. Srivastava,                Date: 14-APR-98  C
!          Princeton University                                        C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!  Modified: S. Benyahia, Fluent Inc.                 Date: 02-FEB-05  C
!      Added the argument L to calc_grbdry                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE JJ_BC_U_S(I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s 
      USE rxns 
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is 
      USE tau_s 
      USE bc
      USE output 
      USE compar   
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Boundary condition
      INTEGER          L
!
!                      Indices
      INTEGER          I,  J, K, I1, I2, J1, J2, K1, K2, IJK, &
                       JM, KM
!
!                      Solids phase
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      coefficients for granular bc
      DOUBLE PRECISION hw, gw, cw
!      
 !-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
!
      DO K = K1, K2 
         DO J = J1, J2 
            DO I = I1, I2 
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     

               IJK = FUNIJK(I,J,K) 
               IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
               JM = JM1(J) 
               KM = KM1(K) 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
!
               IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                  IF (EP_S(NORTH_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,N,M) = -ONE 
                  ELSE 
! start anuj 4/20
                     IF (FRICTION .AND. EP_S(IJK,M)>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, NORTH_OF(IJK), 'N', 'U', L, &
                           M, GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, NORTH_OF(IJK), 'N', 'U', M, &
                              L, HW) 
                           GW = 1D0 
                           CW = HW*BC_UW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_UW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,N,M) = -(HALF*HW - ODY_N(J)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODY_N(J)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
!
               ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                  IF (EP_S(SOUTH_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,S,M) = -ONE 
                  ELSE 
                     IF (FRICTION .AND. EP_S(IJK,M)>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, SOUTH_OF(IJK), 'S', 'U', L, &
                           M, GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, SOUTH_OF(IJK), 'S', 'U', M, &
                              L, HW) 
                           GW = 1D0 
                           CW = HW*BC_UW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_UW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,S,M) = -(HALF*HW - ODY_N(JM)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODY_N(JM)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
!
               ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                  IF (EP_S(TOP_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,T,M) = -ONE 
                  ELSE 
                     IF (FRICTION .AND. EP_S(IJK,M)>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, TOP_OF(IJK), 'T', 'U', L, M&
                           , GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, TOP_OF(IJK), 'T', 'U', M, L, HW) 
                           GW = 1D0 
                           CW = HW*BC_UW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_UW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,T,M) = -(HALF*HW - ODZ_T(K)*OX_E(I)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODZ_T(K)*OX_E(I)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
!
               ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                  IF (EP_S(BOTTOM_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,B,M) = -ONE 
                  ELSE 
                     IF (FRICTION .AND. EP_S(IJK,M)>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, BOTTOM_OF(IJK), 'B', 'U', L&
                           , M, GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, BOTTOM_OF(IJK), 'B', 'U', M, &
                              L, HW) 
                           GW = 1D0 
                           CW = HW*BC_UW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_UW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,B,M) = -(HALF*HW - ODZ_T(KM)*OX_E(I)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODZ_T(KM)*OX_E(I)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
               ENDIF 
            END DO 
         END DO 
      END DO 
     
      RETURN  
      END SUBROUTINE JJ_BC_U_S 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3      
!// 360 Check if i,j,k resides on current processor
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              