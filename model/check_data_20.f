!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CHECK_DATA_20                                           C
!  Purpose:                                                            C
!     - check whether field variables are initialized in all cells     C
!     - check whether the sum of void and volume fractions is 1.0      C
!       in all fluid and mass inflow cells                             C
!     - check whether mu_gmax is specified if k_epsilon or l_scale     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_DATA_20 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE toleranc 
      USE fldvar
      USE run
      USE geometry
      USE constant
      USE physprop
      USE indices
      USE funits 
      USE visc_g
      USE rxns
      USE scalars
      USE compar 
      USE sendrecv 
      USE discretelement
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP, &
                 IJKW, IJKE, IJKS, IJKN, IJKB, IJKT, &
                 IM, JM, KM
! Solids phase
      INTEGER :: M
! Species index
      INTEGER :: N
! Local index
      INTEGER :: L
! Logical variable to set, if there is an error
      LOGICAL :: ABORT
! Whether L_scale is nonzero
      LOGICAL :: NONZERO
! 1.0 - sum of all volume fractions
      DOUBLE PRECISION DIF
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

      call send_recv(p_g,2)
      call send_recv(ep_g,2)
      call send_recv(w_s,2)
      call send_recv(w_g,2)
      call send_recv(u_s,2)
      call send_recv(u_g,2)
      call send_recv(v_s,2)
      call send_recv(v_g,2)
      call send_recv(rop_s,2)
      call send_recv( P_STAR, 2 ) 
      call send_recv( ROP_G, 2 ) 
      call send_recv( ROP_S, 2 ) 
      call send_recv( RO_G, 2 ) 
      call send_recv( T_G, 2 ) 
      call send_recv( T_S, 2 ) 
      call send_recv( X_G, 2 ) 
      call send_recv( X_S, 2 ) 

      CALL START_LOG 
      ABORT = .FALSE. 
      NONZERO = .FALSE. 

! Check whether all field variables are initialized in all fluid cells
! and flow boundary cells
! ---------------------------------------------------------------->>>
      DO K = kstart2, kend2 
         DO J = jstart2, jend2 
            DO I = istart2, iend2 
               IJK = FUNIJK(I,J,K) 
               IF (.NOT.WALL_AT(IJK)) THEN 
                  CALL SET_INDEX1 (IJK, I, J, K, IMJK, IPJK, IJMK, &
                     IJPK, IJKM, IJKP, IJKW, IJKE, IJKS, IJKN, IJKB, &
                     IJKT, IM, JM, KM) 

! check gas phase fields                     
                  IF (EP_G(IJK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1010) I, J, K, 'EP_g' 
                  ENDIF 
                  IF (P_G(IJK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1010) I, J, K, 'P_g' 
                  ENDIF 
                  IF (P_STAR(IJK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1010) I, J, K, 'P_star'
                  ENDIF 
                  IF (RO_G(IJK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1010) I, J, K, 'RO_g' 
                  ENDIF 
                  IF (ROP_G(IJK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1010) I, J, K, 'ROP_g' 
                  ENDIF 
                  IF (T_G(IJK) == UNDEFINED) THEN 
                     IF (ENERGY_EQ .OR. RO_G0==UNDEFINED .OR.&
                         MU_G0==UNDEFINED) THEN 
                        IF (.NOT.ABORT) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                           ABORT = .TRUE. 
                        ENDIF 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1010) I, J, K, 'T_g'
                     ENDIF 
                  ENDIF 
                  IF (U_G(IMJK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1010) I-1, J, K, 'U_g' 
                  ENDIF 
                  IF (U_G(IJK)==UNDEFINED .AND. I/=IMAX2) THEN 
                     IF (.NOT.ABORT) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1010) I, J, K, 'U_g' 
                  ENDIF 
                  IF (V_G(IJMK) == UNDEFINED) THEN 
                     IF (.NOT.ABORT) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1010) I, J-1, K, 'V_g'
                  ENDIF 
                  IF (V_G(IJK)==UNDEFINED .AND. J/=JMAX2) THEN 
                     IF (.NOT.ABORT) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1010) I, J, K, 'V_g' 
                  ENDIF 

                  IF (W_G(IJKM) == UNDEFINED) THEN   
                     IF (.NOT.ABORT) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1010) I, J, K-1, 'W_g' 
                  ENDIF 
                  IF (W_G(IJK)==UNDEFINED .AND. K/=KMAX2) THEN 
                     IF (.NOT.ABORT) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                     ENDIF 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1010) I, J, K, 'W_g' 
                  ENDIF 
                  IF (SPECIES_EQ(0) .OR. RO_G0==UNDEFINED .AND.&
                      MW_AVG==UNDEFINED) THEN 
                     DO N = 1, NMAX(0) 
                        IF (X_G(IJK,N) == UNDEFINED) THEN 
                           IF (.NOT.ABORT) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                              ABORT = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1012) &
                              I, J, K, N, 'X_g' 
                        ENDIF 
                     ENDDO 
                  ENDIF 

                  DO N = 1, NScalar 
                    IF (Scalar(IJK,N) == UNDEFINED) THEN 
                      IF (.NOT.ABORT) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                        ABORT = .TRUE. 
                      ENDIF 
                      IF(DMP_LOG)WRITE (UNIT_LOG, 1012) &
                         I, J, K, N, 'Scalar' 
                    ENDIF 
                  ENDDO 

! check solids phase fields. these quantities are specified via the
! subroutines set_ic and set_bc0/set_bc1 that employ the initial and
! boundary conditions set in the mfix.dat.  
                  IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
! dont need the values at the ghost walls for setting pressure outflow 
! BC's for MPPIC. so not checking this for MPPIC
                     DO M = 1, SMAX 
                        IF (ROP_S(IJK,M) == UNDEFINED) THEN 
                           IF (.NOT.ABORT) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                              ABORT = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1011) &
                              I, J, K, M, 'ROP_s' 
                        ENDIF 
                        IF (T_S(IJK,M) == UNDEFINED) THEN 
                           IF (ENERGY_EQ) THEN 
                              IF (.NOT.ABORT) THEN 
                                 IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                                 ABORT = .TRUE. 
                              ENDIF 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1011) &
                                 I, J, K, M, 'T_s' 
                           ENDIF 
                        ENDIF 
                        IF (U_S(IMJK,M) == UNDEFINED) THEN 
                           IF (.NOT.ABORT) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                              ABORT = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1011) &
                              I-1, J, K, M, 'U_s' 
                        ENDIF 
                        IF (U_S(IJK,M)==UNDEFINED .AND. I/=IMAX2) THEN 
                           IF (.NOT.ABORT) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                              ABORT = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1011) &
                              I, J, K, M, 'U_s' 
                        ENDIF 
                        IF (V_S(IJMK,M) == UNDEFINED) THEN 
                           IF (.NOT.ABORT) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                              ABORT = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1011) &
                              I, J-1, K, M, 'V_s' 
                        ENDIF 
                        IF (V_S(IJK,M)==UNDEFINED .AND. J/=JMAX2) THEN 
                           IF (.NOT.ABORT) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                              ABORT = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1011) &
                              I, J, K, M, 'V_s' 
                        ENDIF 
                        IF (W_S(IJKM,M) == UNDEFINED) THEN 
                           IF (.NOT.ABORT) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                              ABORT = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1011) &
                              I, J, K - 1, M, 'W_s' 
                        ENDIF 
                        IF (W_S(IJK,M)==UNDEFINED .AND. K/=KMAX2) THEN 
                           IF (.NOT.ABORT) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                              ABORT = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1011) &
                              I, J, K, M, 'W_s' 
                        ENDIF 
                        IF (SPECIES_EQ(M)) THEN 
                           DO N = 1, NMAX(M) 
                              IF (X_S(IJK,M,N) == UNDEFINED) THEN 
                                 IF (.NOT.ABORT) THEN 
                                    IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 
                                    ABORT = .TRUE. 
                                 ENDIF 
                                 IF(DMP_LOG)WRITE (UNIT_LOG, 1013) &
                                    I, J, K, M, N, 'X_s'
                              ENDIF 
                           ENDDO 
                        ENDIF
                     ENDDO   ! end do m=1,smax
                  ENDIF      ! if (.not.discrete_element .or.
                             !     des_continuum_hybrid)

               ENDIF   ! IF (.NOT.WALL_AT(IJK)) THEN  
            ENDDO   ! end do I = istart2, iend2 
         ENDDO    ! end do J =j start2, jend2 
      ENDDO   ! end do K = kstart2, kend2 

! end check whether all field variables are initialized in all fluid
! cells and flow boundary cells
! ----------------------------------------------------------------<<<

      IF (ABORT) THEN 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 
         CALL MFIX_EXIT(myPE)
      ENDIF 


! Additional check for fluid or mass inflow cells
! ---------------------------------------------------------------->>>
      DO K = kstart2, kend2 
         DO J = jstart2, jend2 
            DO I = istart2, iend2 
               IJK = FUNIJK(I,J,K) 

               IF (FLAG(IJK)==1 .OR. FLAG(IJK)==20) THEN 

! Check whether L_scale is non-zero anywhere
                  IF (L_SCALE(IJK) /= ZERO) NONZERO = .TRUE. 


! check the sum of volume fractions. 
! this seems redundant with the existing checks in check_data_06 for 
! fluid cells (flag=1) and with check_data_07 for MI cells (flag=20).  
                  IF (.NOT.DISCRETE_ELEMENT) THEN
                     DIF = ONE - EP_G(IJK) 
                     IF (SMAX > 0) THEN 
!QX RO_S changed to RO_S0
                        SELECT CASE (TRIM(RUN_TYPE))  
                        CASE ('NEW')  
                           DIF = DIF - SUM(ROP_S(IJK,:MMAX)/RO_S(:MMAX))
                        CASE ('RESTART_1', 'RESTART_2')  
                           DIF = DIF - SUM(ROP_S(IJK,:MMAX)/RO_SV(IJK,:MMAX))
                        CASE DEFAULT
                     END SELECT 
!end
                     M = SMAX + 1 

                        IF (ABS(DIF) > SMALL_NUMBER) THEN 
                           IF (.NOT.ABORT) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1050)  
                              ABORT = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              I, J, K, (1.- dif) 
                        ENDIF  

! ep_g must have a value > 0 for all models, (sof)
                        IF (EP_G(IJK) < SMALL_NUMBER) THEN 
                           IF (.NOT.ABORT) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1070) 
                              ABORT = .TRUE. 
                           ENDIF 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1150) I, J, K
                        ENDIF

! ep_g cannot exceed unity for all models
                        IF (EP_G(IJK) > ONE) THEN 
                           IF (.NOT.ABORT) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1080) 
                              ABORT = .TRUE. 
                           ENDIF 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1150) I, J, K
                        ENDIF  
                     ENDIF ! for SMAX > 0

                  ELSE   ! else if discrete_element

! for discrete element simulations any solids bulk density is not yet
! assigned. that is, the particle configuration is not known until the
! discrete element portion of the simulation is invoked. (so if 
! run_type==new it will not be specified). so checking that the sum of
! void and volume fractions equal 1 is not meaningful yet. this holds
! true for the hybrid model.

! ep_g must have a value > 0 and < 1 
                     IF (EP_G(IJK) < ZERO .OR. EP_G(IJK) > ONE) THEN 
                        IF (.NOT.ABORT) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1075) 
                           ABORT = .TRUE. 
                        ENDIF 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1150) I, J, K
                     ENDIF
                  ENDIF ! end if/else (.not.discrete_element)

               ENDIF   ! IF (FLAG(IJK)==1 .OR. FLAG(IJK)==20) THEN 
            ENDDO   ! I = istart2, iend2 
         ENDDO    ! J =j start2, jend2 
      ENDDO   ! K = kstart2, kend2 
! end additional checks for fluid or mass inflow cells
! ----------------------------------------------------------------<<<


      IF (ABORT) THEN 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 
         CALL MFIX_EXIT (myPE)
      ENDIF 

!  Check whether MU_gmax is specified
      IF (NONZERO .AND. MU_GMAX==UNDEFINED) THEN 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1350) 
         CALL MFIX_EXIT(myPE) 
      ENDIF 

!  Check whether MU_gmax is specified for turbulence (sof)
      IF (K_Epsilon .AND. MU_GMAX==UNDEFINED) THEN 
         IF(DMP_LOG)WRITE (UNIT_LOG, 1360) 
         CALL MFIX_EXIT(myPE) 
      ENDIF 
      
      CALL END_LOG 
      RETURN  

 1000 FORMAT(/1X,70('*')//' From: CHECK_DATA_20',/&
         ' Message: The following field variables are undefined') 
 1010 FORMAT(1X,'I = ',I4,' J = ',I4,' K = ',I4,5X,A) 
 1011 FORMAT(1X,'I = ',I4,' J = ',I4,' K = ',I4,' M = ',I4,5X,A) 
 1012 FORMAT(1X,'I = ',I4,' J = ',I4,' K = ',I4,' N = ',I4,5X,A) 
 1013 FORMAT(1X,'I = ',I4,' J = ',I4,' K = ',I4,' M = ',I4,' N = ',I4,5X,A) 
 1050 FORMAT(/1X,70('*')//' From: CHECK_DATA_20',/&
         ' Message: The sum of volume fractions is not equal to 1',/&
         '          in the following cells:',/4X,'I',T14,'J',T24,'K')  
 1070 FORMAT(/1X,70('*')//' From: CHECK_DATA_20',/&
         ' Message: EP_g is less than SMALL_NUMBER ',/&
         '          in the following cells:',/4X,'I',T14,'J',T24,'K') 
 1075 FORMAT(/1X,70('*')//' From: CHECK_DATA_20',/&
         ' Message: EP_g is unphysical (>0 or <1) ',/&
         '          in the following cells:',/4X,'I',T14,'J',T24,'K') 
 1080 FORMAT(/1X,70('*')//' From: CHECK_DATA_20',/&
         ' Message: EP_g is greater than one ',/&
         '          in the following cells:',/4X,'I',T14,'J',T24,'K')
 1100 FORMAT(1X,I4,T11,I4,T21,I4,'  Sum of EP = ', G12.5, '.NE. 1')  
 1150 FORMAT(1X,I4,T11,I4,T21,I4) 
 1300 FORMAT(/1X,70('*')/) 
 1350 FORMAT(/1X,70('*')//' From: CHECK_DATA_20',/&
         ' Message: Turbulent length scale is nonzero. Specify MU_gmax.',/1X,70&
         ('*')/)  
 1360 FORMAT(/1X,70('*')//' From: CHECK_DATA_20',/&
         ' Message: K_Epsilon model is used. Specify MU_gmax in mfix.dat.',/1X,70&
         ('*')/)
      END SUBROUTINE CHECK_DATA_20 


