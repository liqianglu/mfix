!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
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
	                       U_G, V_G, W_G, ROP_G, M, A_M, B_M, IER) 
!
!
            CALL BC_PHI (BC_Scalar(1,N), BC_ScalarW(1,N), BC_HW_Scalar(1,N), &
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
	                       U_s(1,m), V_s(1,m), W_s(1,m), ROP_s(1,m), M, &
			       A_M, B_M, IER) 
!
!
            CALL BC_PHI (BC_Scalar(1,N), BC_ScalarW(1,N), BC_HW_Scalar(1,N), &
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
