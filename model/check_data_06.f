!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_06                                          C
!  Purpose: check the initial conditions input section                 C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 27-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Check the specification of physical quantities             C
!  Author: M. Syamlal                                 Date: 24-JUL-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b        C
!                        IC_Z_t, DX, IMAX, DY, JMAX, DZ, KMAX          C
!                         IMIN1, IMAX1, JMIN1, JMAX1, KMIN1, KMAX1     C
!  Variables modified: IC_DEFINED, IC_I_w, IC_I_e, IC_J_s, IC_J_n      C
!                      IC_K_b, IC_k_t, ICBC_FLAG                       C
!                                                                      C
!  Local variables:  ICV, I_w, I_e, J_s, J_n, K_b, K_t, L4              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_06 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE ic
      USE fldvar
      USE physprop
      USE run
      USE indices
      USE funits 
      USE compar      !//AIKEPARDBG      
      USE mpi_utility !//AIKEPARDBG      
      
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!    ICV,I,J,k,L4 - loop/variable indices
!    I_w,I_e    -  cell indices in X direction from IC_X_e,IC_X_w
!    J_s,J_n    -  cell indices in Y direction from IC_Y_s,IC_Y_n
!    K_b,K_t    -  cell indices in Z direction from IC_Z_b,IC_Z_t
!    IC2 - Last two digits of IC
!
      INTEGER I_w , I_e , J_s , J_n , K_b , K_t , ICV
      INTEGER I, J, k, IJK, IC2, M, N
      DOUBLE PRECISION SUM, SUM_EP

      character*3, allocatable :: array1(:)   !//AIKEPARDBG      
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
! Initialize the icbc_flag array.  If not a NEW run then do not
! check the initial conditions.
!

!//AIKEPARDBGSTOP 0922
      write(UNIT_LOG,"('(PE ',I2,'): beginning of check_data_06')") myPE       !//AIKEPARDBG
!     write(UNIT_LOG,"('(PE ',I2,'): from chk_data_06.f ',&                    !//AIKEPARDBG
!                /,9X,'Kmin3 = ',I6,'  Kmax3 = ',I6,'  Kmax = ',I6, &   !//AIKEPARDBG
!                /,9X,'Jmin3 = ',I6,'  Jmax3 = ',I6,'  Jmax = ',I6,&    !//AIKEPARDBG
!                /,9X,'Imin3 = ',I6,'  Imax3 = ',I6,'  Imax = ',I6)") &  !//AIKEPARDBG
!                myPE,Kmin3,Kmax3,Kmax,Jmin3,Jmax3,Jmax,Imin3,Imax3,Imax !//AIKEPARDBG
!     write(UNIT_LOG,"('(PE ',I2,'): from chk_data_06.f ' &                    !//AIKEPARDBG
!                /,9X,'Kstart3 = ',I6,'  Kstart2 = ',I6,'  Kend3 = ',I6 &   !//AIKEPARDBG
!                /,9X,'Jstart3 = ',I6,'  Jstart2 = ',I6,'  Jend3 = ',I6 &    !//AIKEPARDBG
!	 /,9X,'Istart3 = ',I6,'  Istart2 = ',I6,'  Iend3 = ',I6 &
!         )") &  !//AIKEPARDBG
!                myPE,Kstart3,Kstart2,Kend3,Jstart3,Jstart2,Jend3, & !//AIKEPARDBG
!             Istart3,Istart2,Iend3 !//AIKEPARDBG
!     write(UNIT_LOG,"('(PE ',I2,'): ',/, &
!                /,9X,'IJKstart3 = ',I6,'  IJKend3 = ',I6 &   !//AIKEPARDBG      
!                 ,/,9X,'IJKstart2 = ',I6,'  IJKend2 = ',I6, &   !//AIKEPARDBG      		 
!                )") myPE, ijkstart3,ijkend3
		 
      
!      call exitMPI(myPE)   !//AIKEPARDBGSTOP

!// 200 1008 Changed the limits for the DO loop KMIN3-->kstart3, KMAX3-->Kend3
!     DO K = KMIN3, KMAX3 
       DO K = Kstart3, Kend3       
         DO J = Jstart3, Jend3 
           DO I = Istart3, Iend3 
!// 220 1004 Replaced with global FUNIJK		     
               IJK = FUNIJK_GL(I,J,K)

               IF (RUN_TYPE == 'NEW') THEN 
                  ICBC_FLAG(IJK) = '   ' 
               ELSE 
                  ICBC_FLAG(IJK) = '.--' 
               ENDIF 
	       
               IF (DO_K) THEN 
!// 200 1008 Changed the limits in K BC type label assignment in 2nd layer of ghost cells
!                  IF (K==1 .OR. K==KMAX2) THEN 
                  IF (K==KMIN3 .OR. K==KMIN2 .OR. &
		      K==KMAX2 .OR. K==KMAX3) THEN 		  
                     IF (CYCLIC_Z_PD) THEN 
                        ICBC_FLAG(IJK) = 'C--' 
                     ELSE IF (CYCLIC_Z) THEN 
                        ICBC_FLAG(IJK) = 'c--' 
                     ELSE 
                        ICBC_FLAG(IJK) = 'W--' 
!//AIKEPARDBG	       			
!	                write(UNIT_LOG,"('(PE ',I2,'): ICBC_FLAG(',I6,',',I6,',',I6, &
!	                      ') = ',A3)") myPE,I,J,K,ICBC_FLAG(IJK)   !//AIKEPARDBG
                     ENDIF 
                  ENDIF 
               ENDIF 

!      call mfix_exit(myPE) !//AIKEPARDBG

!//S 1008 do similar limit modifications for 2D/3D decomposition in following
!//S      by replacing with JMIN3,JMIN2,JMAX3	       
               IF (DO_J) THEN 
!// 200 1008 Changed the limits in J, BC type label assignment in 2nd layer of ghost cells
                  IF (J==JMIN3 .OR. J==JMIN2 .OR. &
		      J==JMAX2 .OR. J==JMAX3) THEN 		  		  
                     IF (CYCLIC_Y_PD) THEN 
                        ICBC_FLAG(IJK) = 'C--' 
                     ELSE IF (CYCLIC_Y) THEN 
                        ICBC_FLAG(IJK) = 'c--' 
                     ELSE 
                        ICBC_FLAG(IJK) = 'W--' 
                     ENDIF 
                  ENDIF 
               ENDIF 
!//S 1008 do similar limit modifications for 2D/3D decomposition in following	       
!//S      by replacing with IMIN3,IMIN2,IMAX3	       
               IF (DO_I) THEN 
!// 200 1008 Changed the limits in I, BC type label assignment in 2nd layer of ghost cells	       
                  IF (I==IMIN3 .OR. I==IMIN2 .OR. &
		      I==IMAX2 .OR. I==IMAX3) THEN 		  
		  
                     IF (CYCLIC_X_PD) THEN 
                        ICBC_FLAG(IJK) = 'C--' 
                     ELSE IF (CYCLIC_X) THEN 
                        ICBC_FLAG(IJK) = 'c--' 
                     ELSE 
                        ICBC_FLAG(IJK) = 'W--' 
                     ENDIF 
                  ENDIF 
                  IF (I==1 .AND. CYLINDRICAL .AND. XMIN==ZERO) ICBC_FLAG(IJK)&
                      = 'S--' 
               ENDIF 
!
!              corner cells are wall cells
!// 200 1008 modified the limits for additional ghost cells when detecting corners
!               IF ((I==1 .OR. I==IMAX2) .AND. (J==1 .OR. J==JMAX2) .AND. (K==1&
!                   .OR. K==KMAX2)) THEN 
               IF ((I==IMIN3 .OR. I==IMIN2 .OR. I==IMAX2 .OR. I==IMAX3) .AND. &
	           (J==JMIN3 .OR. J==JMIN2 .OR. J==JMAX2 .OR. J==JMIN3) .AND. &
		   (K==KMIN3 .OR. K==KMIN2 .OR. K==KMAX2 .OR. K==KMAX3)) THEN 		   
                  IF (ICBC_FLAG(IJK) /= 'S--') ICBC_FLAG(IJK) = 'W--' 
		  
               ENDIF 
            END DO 
         END DO 	 
      END DO 



!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): aft 1st DO loop in chk_data_06')") myPE !//AIKEPARDBG
!      if (myPE == 0) &
!      CALL OUT_ARRAY_C (ICBC_FLAG, 'BC/IC condition flags')  
!//AIKEPARDBG dump the ICBC_FLAG in matrix form to verify with serial version
!       DO K = kstart3, kend3                               !//AIKEPARDBG
!          write(UNIT_LOG,"('K = ',I5)") K                !//AIKEPARDBG 
! 	 write(UNIT_LOG,"(7X,14(I3,2X))") (I,i=IMIN3,IMAX3)  !//AIKEPARDBG
!          DO J = jstart3, Jend3                            !//AIKEPARDBG
!            write(UNIT_LOG,"(I5,')',$)") J               !//AIKEPARDBG	 
!            DO I = istart3, Iend3                          !//AIKEPARDBG
!              IJK = FUNIJK_GL(I,J,K)                     !//AIKEPARDBG
!              write(UNIT_LOG,"(2X,A3,$)") ICBC_FLAG(IJK) !//AIKEPARDBG
!            END DO                                       !//AIKEPARDBG
!            write(UNIT_LOG,"(/)")                        !//AIKEPARDBG
!          END DO                                         !//AIKEPARDBG
!       END DO                                            !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG
      
      DO ICV = 1, DIMENSION_IC 
         IC_DEFINED(ICV) = .FALSE. 
         IF (IC_X_W(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_X_E(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_Y_S(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_Y_N(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_Z_B(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_Z_T(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_I_W(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_I_E(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_J_S(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_J_N(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_K_B(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_K_T(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_DEFINED(ICV)) THEN 
            IF (IC_X_W(ICV)==UNDEFINED .AND. IC_I_W(ICV)==UNDEFINED_I) THEN 
               IF (NO_I) THEN 
                  IC_X_W(ICV) = ZERO 
               ELSE 
                     WRITE (UNIT_LOG, 1000) 'IC_X_w and IC_I_w ', ICV 
                     call mfix_exit(myPE) !// 990 0912 replaced STOP so that all PEs are aborted
               ENDIF 
            ENDIF 
            IF (IC_X_E(ICV)==UNDEFINED .AND. IC_I_E(ICV)==UNDEFINED_I) THEN 
               IF (NO_I) THEN 
                  IC_X_E(ICV) = XLENGTH 
               ELSE 
                 WRITE (UNIT_LOG, 1000) 'IC_X_e and IC_I_e ', ICV 
                 call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
            ENDIF 
            IF (IC_Y_S(ICV)==UNDEFINED .AND. IC_J_S(ICV)==UNDEFINED_I) THEN 
               IF (NO_J) THEN 
                  IC_Y_S(ICV) = ZERO 
               ELSE 
                  WRITE (UNIT_LOG, 1000) 'IC_Y_s and IC_J_s ', ICV 
                  call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
            ENDIF 
            IF (IC_Y_N(ICV)==UNDEFINED .AND. IC_J_N(ICV)==UNDEFINED_I) THEN 
               IF (NO_J) THEN 
                  IC_Y_N(ICV) = YLENGTH 
               ELSE 
                  WRITE (UNIT_LOG, 1000) 'IC_Y_n and IC_J_n ', ICV 
                  call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
            ENDIF 
            IF (IC_Z_B(ICV)==UNDEFINED .AND. IC_K_B(ICV)==UNDEFINED_I) THEN 
               IF (NO_K) THEN 
                  IC_Z_B(ICV) = ZERO 
               ELSE 
                  WRITE (UNIT_LOG, 1000) 'IC_Z_b and IC_K_b ', ICV 
                  call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
            ENDIF 
            IF (IC_Z_T(ICV)==UNDEFINED .AND. IC_K_T(ICV)==UNDEFINED_I) THEN 
               IF (NO_K) THEN 
                  IC_Z_T(ICV) = ZLENGTH 
               ELSE 
                  WRITE (UNIT_LOG, 1000) 'IC_Z_t and IC_K_t ', ICV 
                  call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
            ENDIF 
         ENDIF 
      END DO 

!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): aft 1st DO ICV loop in chk_data_06')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG

      DO ICV = 1, DIMENSION_IC 
!
         IF (IC_X_W(ICV)/=UNDEFINED .AND. IC_X_E(ICV)/=UNDEFINED) THEN 
            IF (NO_I) THEN 
               I_W = 1 
               I_E = 1 
            ELSE 
               CALL CALC_CELL (IC_X_W(ICV), DX, IMAX, I_W) 
               I_W = I_W + 1 
               CALL CALC_CELL (IC_X_E(ICV), DX, IMAX, I_E) 
            ENDIF 
            IF (IC_I_W(ICV)/=UNDEFINED_I .OR. IC_I_E(ICV)/=UNDEFINED_I) THEN 
               CALL LOCATION_CHECK (IC_I_W(ICV), I_W, ICV, 'IC - west') 
               CALL LOCATION_CHECK (IC_I_E(ICV), I_E, ICV, 'IC - east') 
            ELSE 
               IC_I_W(ICV) = I_W 
               IC_I_E(ICV) = I_E 
            ENDIF 
         ENDIF 
!
         IF (IC_Y_S(ICV)/=UNDEFINED .AND. IC_Y_N(ICV)/=UNDEFINED) THEN 
            IF (NO_J) THEN 
               J_S = 1 
               J_N = 1 
            ELSE 
               CALL CALC_CELL (IC_Y_S(ICV), DY, JMAX, J_S) 
               J_S = J_S + 1 
               CALL CALC_CELL (IC_Y_N(ICV), DY, JMAX, J_N) 
            ENDIF 
            IF (IC_J_S(ICV)/=UNDEFINED_I .OR. IC_J_N(ICV)/=UNDEFINED_I) THEN 
               CALL LOCATION_CHECK (IC_J_S(ICV), J_S, ICV, 'IC - south') 
               CALL LOCATION_CHECK (IC_J_N(ICV), J_N, ICV, 'IC - north') 
            ELSE 
               IC_J_S(ICV) = J_S 
               IC_J_N(ICV) = J_N 
            ENDIF 
         ENDIF 
!
!//D 1008 Calc_cell runs over KMAX so no modification due DDECOMP
!//? we may need to flag out the situations where IC_Z_T(ICV) or IC_K_T() may
!//? not be on this processors, flaggging out such circumstances may be worth.
!//? see the following loop where this check is done for the global domain!!!
         IF (IC_Z_B(ICV)/=UNDEFINED .AND. IC_Z_T(ICV)/=UNDEFINED) THEN 
            IF (NO_K) THEN 
               K_B = 1 
               K_T = 1 
            ELSE 
               CALL CALC_CELL (IC_Z_B(ICV), DZ, KMAX, K_B) 
               K_B = K_B + 1 
               CALL CALC_CELL (IC_Z_T(ICV), DZ, KMAX, K_T) 
            ENDIF 
            IF (IC_K_B(ICV)/=UNDEFINED_I .OR. IC_K_T(ICV)/=UNDEFINED_I) THEN 
               CALL LOCATION_CHECK (IC_K_B(ICV), K_B, ICV, 'IC - bottom') 
               CALL LOCATION_CHECK (IC_K_T(ICV), K_T, ICV, 'IC - top') 
            ELSE 
               IC_K_B(ICV) = K_B 
               IC_K_T(ICV) = K_T 
            ENDIF 
         ENDIF 
      END DO 

!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): aft 2nd DO ICV loop in chk_data_06')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG

      DO ICV = 1, DIMENSION_IC 
         IF (IC_DEFINED(ICV)) THEN 
!
!   For Restart runs IC is defined only if IC_TYPE='PATCH'
!
            IF (RUN_TYPE/='NEW' .AND. IC_TYPE(ICV)/='PATCH') THEN 
               IC_DEFINED(ICV) = .FALSE. 
               CYCLE  
            ENDIF 

!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): INTERCHK in chk_data_06')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG

!
            IF (IC_I_W(ICV) > IC_I_E(ICV)) GO TO 900 
            IF (IC_J_S(ICV) > IC_J_N(ICV)) GO TO 900 
            IF (IC_K_B(ICV) > IC_K_T(ICV)) GO TO 900 
            IF (IC_I_W(ICV)<IMIN1 .OR. IC_I_W(ICV)>IMAX1) GO TO 900 
            IF (IC_I_E(ICV)<IMIN1 .OR. IC_I_E(ICV)>IMAX1) GO TO 900 
            IF (IC_J_S(ICV)<JMIN1 .OR. IC_J_S(ICV)>JMAX1) GO TO 900 
            IF (IC_J_N(ICV)<JMIN1 .OR. IC_J_N(ICV)>JMAX1) GO TO 900 
!//D 1008 no changes for DDECOMP as this performs reality check w/ bndaries	    
!//? sanity check should be also based on kstart1 and kend1 so that the
!//? IC_Z_T() etc. not residing on current PE may be flagged easily here,i.e.
!//?        IF (IC_K_B(ICV)<kstart1 .OR. IC_K_B(ICV)>Kend1) GO TO 900 
            IF (IC_K_B(ICV)<KMIN1 .OR. IC_K_B(ICV)>KMAX1) GO TO 900 
            IF (IC_K_T(ICV)<KMIN1 .OR. IC_K_T(ICV)>KMAX1) GO TO 900 

!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): INTERCHK2 in chk_data_06')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG

!
!  If a 'PATCH' need not check whether all the variables are specified
!
            IF (IC_TYPE(ICV) /= 'PATCH') THEN 
!
!  Check the specification of physical quantities
!
               IF (IC_U_G(ICV) == UNDEFINED) THEN 
                  IF (NO_I) THEN 
                     IC_U_G(ICV) = ZERO 
                  ELSE 
                     WRITE (UNIT_LOG, 1000) 'IC_U_g', ICV 
                     call mfix_exit(myPE) !// 990 0912 replaced STOP 
                  ENDIF 
               ENDIF 
	       	       
               IF (IC_V_G(ICV) == UNDEFINED) THEN 
                  IF (NO_J) THEN 
                     IC_V_G(ICV) = ZERO 
                  ELSE 
                     WRITE (UNIT_LOG, 1000) 'IC_V_g', ICV 
                     call mfix_exit(myPE) !// 990 0912 replaced STOP 
                  ENDIF 
               ENDIF 
               IF (IC_W_G(ICV) == UNDEFINED) THEN 
                  IF (NO_K) THEN 
                     IC_W_G(ICV) = ZERO 
                  ELSE 
                     WRITE (UNIT_LOG, 1000) 'IC_W_g', ICV 
                     call mfix_exit(myPE) !// 990 0912 replaced STOP 
                  ENDIF 
               ENDIF 

!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): INTERCHK3 in chk_data_06')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG
	       
               IF (IC_EP_G(ICV) == UNDEFINED) THEN 
                  WRITE (UNIT_LOG, 1000) 'IC_EP_g', ICV 
                  call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
               IF (IC_P_G(ICV) /= UNDEFINED) THEN 
                  IF (RO_G0==UNDEFINED .AND. IC_P_G(ICV)<=ZERO) THEN 
                     WRITE (UNIT_LOG, 1010) ICV, IC_P_G(ICV) 
                     call mfix_exit(myPE) !// 990 0912 replaced STOP 
                  ENDIF 
               ENDIF 

!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): INTERCHK4 in chk_data_06')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG

!
               IF ((ENERGY_EQ .OR. RO_G0==UNDEFINED .OR. MU_G0==UNDEFINED)&
                   .AND. IC_T_G(ICV)==UNDEFINED) THEN 
                     WRITE (UNIT_LOG, 1000) 'IC_T_g', ICV 
                     call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
!
               IF (ENERGY_EQ) THEN 
                  IF (IC_GAMA_RG(ICV) < ZERO) THEN 
                     WRITE (UNIT_LOG, 1001) 'IC_GAMA_Rg', ICV 
                     call mfix_exit(myPE) !// 990 0912 replaced STOP 
                  ELSE IF (IC_GAMA_RG(ICV) > ZERO) THEN 
                     IF (IC_T_RG(ICV) == UNDEFINED) THEN 
                       WRITE (UNIT_LOG, 1000) 'IC_T_Rg', ICV 
                       call mfix_exit(myPE) !// 990 0912 replaced STOP 
                     ENDIF 
                  ENDIF 
               ENDIF 
!
!//AIKEPARDBGSTOP 0922
      write(*,"('(PE ',I2,'): INTERCHK5 in chk_data_06')") myPE !//AIKEPARDBG
      call mfix_exit(myPE) !//AIKEPARDBG

               SUM = ZERO 
               DO N = 1, NMAX(0) 
                  IF (IC_X_G(ICV,N) /= UNDEFINED) SUM = SUM + IC_X_G(ICV,N) 
               END DO 
               DO N = 1, NMAX(0) 
                  IF (IC_X_G(ICV,N) == UNDEFINED) THEN 
                     IF (.NOT.COMPARE(ONE,SUM)) WRITE (UNIT_LOG, 1050) ICV, N 
                     IC_X_G(ICV,N) = ZERO 
                  ENDIF 
               END DO 
               IF (.NOT.COMPARE(ONE,SUM)) THEN 
                  WRITE (UNIT_LOG, 1055) ICV 
                  IF (SPECIES_EQ(0) .OR. RO_G0==UNDEFINED .AND. MW_AVG==&
                     UNDEFINED) call mfix_exit(myPE)  
               ENDIF 
!
!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): INTERCHK6 in chk_data_06')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG

               SUM_EP = IC_EP_G(ICV) 
               DO M = 1, MMAX 
                  IF (IC_ROP_S(ICV,M) == UNDEFINED) THEN 
                     IF (IC_EP_G(ICV) == ONE) THEN 
                        IC_ROP_S(ICV,M) = ZERO 
                     ELSE IF (MMAX == 1) THEN 
                        IC_ROP_S(ICV,M) = (ONE - IC_EP_G(ICV))*RO_S(M) 
                     ELSE 
                           WRITE (UNIT_LOG, 1100) 'IC_ROP_s', ICV, M 
                           STOP  
                     ENDIF 
                  ENDIF 
                  SUM_EP = SUM_EP + IC_ROP_S(ICV,M)/RO_S(M) 
                  SUM = ZERO 
                  DO N = 1, NMAX(M) 
                     IF(IC_X_S(ICV,M,N)/=UNDEFINED)SUM=SUM+IC_X_S(ICV,M,N) 
                  END DO 
                  IF (IC_ROP_S(ICV,M)==ZERO .AND. SUM==ZERO) THEN 
                     IC_X_S(ICV,M,1) = ONE 
                     SUM = ONE 
                  ENDIF 
                  DO N = 1, NMAX(M) 
                     IF (IC_X_S(ICV,M,N) == UNDEFINED) THEN 
                        IF(.NOT.COMPARE(ONE,SUM))WRITE(UNIT_LOG,1110)ICV,M,N 
                        IC_X_S(ICV,M,N) = ZERO 
                     ENDIF 
                  END DO 
                  IF (.NOT.COMPARE(ONE,SUM)) THEN 
                        WRITE (UNIT_LOG, 1120) ICV, M 
                        IF (SPECIES_EQ(M)) STOP  
                  ENDIF 
                  IF (IC_U_S(ICV,M) == UNDEFINED) THEN 
                     IF (IC_ROP_S(ICV,M)==ZERO .OR. NO_I) THEN 
                        IC_U_S(ICV,M) = ZERO 
                     ELSE 
                        WRITE (UNIT_LOG, 1100) 'IC_U_s', ICV, M 
                        call mfix_exit(myPE) !// 990 0912 replaced STOP 
                     ENDIF 
                  ENDIF 
                  IF (IC_V_S(ICV,M) == UNDEFINED) THEN 
                     IF (IC_ROP_S(ICV,M)==ZERO .OR. NO_J) THEN 
                        IC_V_S(ICV,M) = ZERO 
                     ELSE 
                        WRITE (UNIT_LOG, 1100) 'IC_V_s', ICV, M 
                        call mfix_exit(myPE) !// 990 0912 replaced STOP 
                     ENDIF 
                  ENDIF 
                  IF (IC_W_S(ICV,M) == UNDEFINED) THEN 
                     IF (IC_ROP_S(ICV,M)==ZERO .OR. NO_K) THEN 
                        IC_W_S(ICV,M) = ZERO 
                     ELSE 
                        WRITE (UNIT_LOG, 1100) 'IC_W_s', ICV, M 
                        call mfix_exit(myPE) !// 990 0912 replaced STOP 
                     ENDIF 
                  ENDIF 
!
                  IF (ENERGY_EQ .AND. IC_T_S(ICV,M)==UNDEFINED) THEN 
                     IF (IC_ROP_S(ICV,M) == ZERO) THEN 
                        IC_T_S(ICV,M) = IC_T_G(ICV) 
                     ELSE 
                        WRITE (UNIT_LOG, 1100) 'IC_T_s', ICV, M 
                        call mfix_exit(myPE) !// 990 0912 replaced STOP 
                     ENDIF 
                  ENDIF 
!
!
                  IF (GRANULAR_ENERGY .AND. IC_THETA_M(ICV,M)==UNDEFINED) THEN 
                     IF (IC_ROP_S(ICV,M) == ZERO) THEN 
                        IC_THETA_M(ICV,M) = ZERO 
                     ELSE 
                        WRITE (UNIT_LOG, 1100) 'IC_Theta_m', ICV, M 
                        call mfix_exit(myPE) !// 990 0912 replaced STOP 
                     ENDIF 
                  ENDIF 
!
                  IF (ENERGY_EQ) THEN 
                     IF (IC_GAMA_RS(ICV,M) < ZERO) THEN 
                        WRITE (UNIT_LOG, 1101) 'IC_GAMA_Rs', ICV, M 
                        call mfix_exit(myPE) !// 990 0912 replaced STOP 
                     ELSE IF (IC_GAMA_RS(ICV,M) > ZERO) THEN 
                        IF (IC_T_RS(ICV,M) == UNDEFINED) THEN 
                           WRITE (UNIT_LOG, 1100) 'IC_T_Rs', ICV, M 
                           call mfix_exit(myPE) !// 990 0912 replaced STOP 
                        ENDIF 
                     ENDIF 
                  ENDIF 
               END DO 

!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): INTERCHK8 in chk_data_06')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG
	       
               IF (.NOT.COMPARE(ONE,SUM_EP)) THEN 
                     WRITE (UNIT_LOG, 1125) ICV 
                     call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
!
!  Set ICBC flag
!
            ENDIF 

!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): INTERCHK9 in chk_data_06')") myPE !//AIKEPARDBG
!      write(*,"('(PE ',I2,'): from chk_data_06.f at ICV = ',I6,&!//AIKEPARDBG
!                 /,9X,'IC_K_B = ',I6,'  IC_K_T = ',I6, &        !//AIKEPARDBG
!                 /,9X,'IC_J_S = ',I6,'  IC_J_N = ',I6, &        !//AIKEPARDBG
!		 /,9X,'IC_I_W = ',I6,'  IC_I_E = ',I6)") &       !//AIKEPARDBG
!                 myPE,ICV,IC_K_B(ICV), IC_K_T(ICV),&            !//AIKEPARDBG
!		 IC_J_S(ICV), IC_J_N(ICV),&                      !//AIKEPARDBG
!		 IC_I_W(ICV), IC_I_E(ICV)                        !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG

!//? WE NEED to make sure that if do loop limits yield to a IJK not residing
!//? on the current PE's subdomain, then this loop will attempt to assign 
!//? ijk values out of the current subdomain on each PE.
            DO I = IC_K_B(ICV), IC_K_T(ICV) 
               DO J = IC_J_S(ICV), IC_J_N(ICV) 
                  DO K = IC_I_W(ICV), IC_I_E(ICV) 
!// 220 1004 Replaced with global FUNIJK		     
!                     IJK = FUNIJK_GL(K,J,I) 
                     IJK = FUNIJK(K,J,I) 
        write(UNIT_LOG,"('(PE ',I2,'): i = ',I5,'  j = ',I5,'  k = ',I5,' ')") myPE,i,j,k
		     		     
		     if(IJK /= -99999) then
               write(*,"('(PE ',I2,'): IJK = ',I6)") myPE,IJK	 !//AIKEPARDBG	     
!               call mfix_exit(myPE) !//AIKEPARDBG	       
                     ICBC_FLAG(IJK)(1:1) = '.' 
                     IC2 = MOD(ICV,100) 
                     WRITE (ICBC_FLAG(IJK)(2:3), 1150) IC2 
		     endif
                  END DO 
               END DO 
            END DO 

!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): INTERCHK10 in chk_data_06')") myPE !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG
	    
         ELSE 
!
!  Check whether physical quantities are specified for undefined
!  initial conditions
!
            IF (IC_U_G(ICV) /= UNDEFINED) THEN 
                WRITE (UNIT_LOG, 1200) 'IC_U_g', ICV 
                call mfix_exit(myPE) !// 990 0912 replaced STOP 
            ENDIF 
            IF (IC_V_G(ICV) /= UNDEFINED) THEN 
                WRITE (UNIT_LOG, 1200) 'IC_V_g', ICV 
                call mfix_exit(myPE) !// 990 0912 replaced STOP 
            ENDIF 
            IF (IC_W_G(ICV) /= UNDEFINED) THEN 
                WRITE (UNIT_LOG, 1200) 'IC_W_g', ICV 
                call mfix_exit(myPE) !// 990 0912 replaced STOP 
            ENDIF 
            IF (IC_EP_G(ICV) /= UNDEFINED) THEN 
                WRITE (UNIT_LOG, 1200) 'IC_EP_g', ICV 
                call mfix_exit(myPE) !// 990 0912 replaced STOP 
            ENDIF 
            IF (IC_T_G(ICV) /= UNDEFINED) THEN 
                WRITE (UNIT_LOG, 1200) 'IC_T_g', ICV 
                call mfix_exit(myPE) !// 990 0912 replaced STOP 
            ENDIF 
            IF (IC_T_RG(ICV) /= UNDEFINED) THEN 
                WRITE (UNIT_LOG, 1200) 'IC_T_Rg', ICV 
                call mfix_exit(myPE) !// 990 0912 replaced STOP 
            ENDIF 
            DO N = 1, DIMENSION_N_G 
               IF (IC_X_G(ICV,N) /= UNDEFINED) THEN 
                  WRITE (UNIT_LOG, 1200) 'X_g', ICV 
                  call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
            END DO 
            DO M = 1, DIMENSION_M 
               IF (IC_ROP_S(ICV,M) /= UNDEFINED) THEN 
                   WRITE (UNIT_LOG, 1300) 'IC_ROP_s', ICV, M 
                   call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
               DO N = 1, DIMENSION_N_S 
                  IF (IC_X_S(ICV,M,N) /= UNDEFINED) THEN 
                      WRITE (UNIT_LOG, 1300) 'IC_X_s', ICV, M 
                      call mfix_exit(myPE) !// 990 0912 replaced STOP 
                  ENDIF 
               END DO 
               IF (IC_U_S(ICV,M) /= UNDEFINED) THEN 
                     WRITE (UNIT_LOG, 1300) 'IC_U_s', ICV, M 
                     call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
               IF (IC_V_S(ICV,M) /= UNDEFINED) THEN 
                     WRITE (UNIT_LOG, 1300) 'IC_V_s', ICV, M 
                     call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
               IF (IC_W_S(ICV,M) /= UNDEFINED) THEN 
                     WRITE (UNIT_LOG, 1300) 'IC_W_s', ICV, M 
                     call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
               IF (IC_T_S(ICV,M) /= UNDEFINED) THEN 
                     WRITE (UNIT_LOG, 1300) 'IC_T_s', ICV, M 
                     call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
               IF (IC_T_RS(ICV,M) /= UNDEFINED) THEN 
                     WRITE (UNIT_LOG, 1300) 'IC_T_Rs', ICV, M 
                     call mfix_exit(myPE) !// 990 0912 replaced STOP 
               ENDIF 
            END DO 
         ENDIF 
      END DO 

!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): end of chk_data_06')") myPE !//AIKEPARDBG
!//AIKEPARDBG dump the ICBC_FLAG in matrix form to verify with serial version
      DO K = Kstart3, Kend3                               !//AIKEPARDBG
         write(UNIT_LOG,"('K = ',I5)") K                !//AIKEPARDBG 
	 write(UNIT_LOG,"(7X,14(I3,2X))") (I,i=IMIN3,IMAX3)  !//AIKEPARDBG
         DO J = Jstart3, Jend3                            !//AIKEPARDBG
           write(UNIT_LOG,"(I5,')',$)") J               !//AIKEPARDBG	
           DO I = Istart3, Iend3                          !//AIKEPARDBG
             IJK = FUNIJK_GL(I,J,K)                     !//AIKEPARDBG
             write(UNIT_LOG,"(2X,A3,$)") ICBC_FLAG(IJK) !//AIKEPARDBG
           END DO                                       !//AIKEPARDBG
           write(UNIT_LOG,"(/)")                        !//AIKEPARDBG
         END DO                                         !//AIKEPARDBG
      END DO                                            !//AIKEPARDBG
      do icv=1,dimension_ic
         write(UNIT_LOG,"('(PE ',I2,'): IC_K_B(',I4,') = ',I5,'  IC_K_T= ',I5)") myPE,ICV,IC_K_B(ICV),IC_K_T(ICV)
         write(UNIT_LOG,"('(PE ',I2,'): IC_J_S       = ',I5,'  IC_J_N= ',I5)") myPE,IC_J_S(ICV),IC_J_N(ICV)
         write(UNIT_LOG,"('(PE ',I2,'): IC_I_W       = ',I5,'  IC_I_E= ',I5)") myPE,IC_I_W(ICV),IC_I_E(ICV)      
      end do      
      call mfix_exit(myPE) !//AIKEPARDBG


      RETURN  
!
! here if error in indices
!
  900 CONTINUE 
      CALL ERROR_ROUTINE ('check_data_06', 'Invalid IC region specified', 0, 2) 
         WRITE (UNIT_LOG, *) ' IC number = ', ICV
         WRITE (UNIT_LOG, *) ' IC_I_w(ICV) = ', IC_I_W(ICV) 
         WRITE (UNIT_LOG, *) ' IC_I_e(ICV) = ', IC_I_E(ICV) 
         WRITE (UNIT_LOG, *) ' IC_J_s(ICV) = ', IC_J_S(ICV) 
         WRITE (UNIT_LOG, *) ' IC_J_n(ICV) = ', IC_J_N(ICV) 
         WRITE (UNIT_LOG, *) ' IC_K_b(ICV) = ', IC_K_B(ICV) 
         WRITE (UNIT_LOG, *) ' IC_K_t(ICV) = ', IC_K_T(ICV) 
      CALL ERROR_ROUTINE (' ', ' ', 1, 3) 
!
 1000 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: ',A,'(',I2,&
         ') not specified',/1X,70('*')/) 
 1001 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: ',A,'(',I2,&
         ') value is unphysical',/1X,70('*')/) 
 1010 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: IC_P_g( ',I2,&
         ') = ',G12.5,/,&
         ' Pressure should be greater than 0.0 for compressible flow',/1X,70(&
         '*')/) 
 1050 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: IC_X_g(',I2,',',I2&
         ,') not specified',/1X,70('*')/) 
 1055 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: IC number:',I2,&
         ' - Sum of gas mass fractions is NOT equal to one',/1X,70('*')/) 
 1100 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: ',A,'(',I2,',',I1,&
         ') not specified',/1X,70('*')/) 
 1101 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: ',A,'(',I2,',',I1,&
         ') unphysical',/1X,70('*')/) 
 1110 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: IC_X_s(',I2,',',I2&
         ,',',I2,') not specified',/1X,70('*')/) 
 1120 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: IC number:',I2,&
         ' - Sum of solids-',I1,' mass fractions is NOT equal to one',/1X,70(&
         '*')/) 
 1125 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: IC number:',I2,&
         ' - Sum of volume fractions is NOT equal to one',/1X,70('*')/) 
 1150 FORMAT(I2.2) 
 1200 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: ',A,'(',I2,&
         ') specified',' for an undefined IC region',/1X,70('*')/) 
 1300 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: ',A,'(',I2,',',I1,&
         ') specified',' for an undefined IC region',/1X,70('*')/) 
      END SUBROUTINE CHECK_DATA_06 
