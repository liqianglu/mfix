!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
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
