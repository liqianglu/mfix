!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: INTERP_RES                                             C
!  Purpose: Interpolate from old RES file to create a new RES file     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 03-DEC-93  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE INTERP_RES
!
      Use param
      Use param1
      Use geometry
      Use indices
      Use energy
      Use physprop
      Use fldvar
      Use post3d
      Use run
      Use scalars
      Use funits
      Use compar
      
      IMPLICIT NONE
      INCLUDE 'xforms.inc'
!
!  Function subroutines
!
      LOGICAL OPEN_FILEP
      INTEGER GET_INDEX, LC
!
!  Local variables
!
      DOUBLE PRECISION , DIMENSION(:), ALLOCATABLE ::  &
        EP_g_OLD, P_g_OLD, &
        P_star_OLD, RO_g_OLD, &
        ROP_g_OLD, T_g_OLD, &
	U_g_OLD, &
        V_g_OLD, W_g_OLD, &
        GAMA_RG_OLD, T_RG_OLD

      DOUBLE PRECISION , DIMENSION(:,:), ALLOCATABLE ::  &
        X_g_OLD

      DOUBLE PRECISION , DIMENSION(:,:), ALLOCATABLE ::  &
        Scalar_OLD
	  
      DOUBLE PRECISION , DIMENSION(:,:), ALLOCATABLE ::  &
        THETA_M_OLD, &
        T_s_OLD, &
        ROP_s_OLD, &
        U_s_OLD, &
        V_s_OLD,&
        W_s_OLD, &
        GAMA_RS_OLD, T_RS_OLD

      DOUBLE PRECISION , DIMENSION(:,:,:), ALLOCATABLE ::  &
        X_s_OLD
	
      DOUBLE PRECISION TIME_OLD
      REAL &
        XDIST_SC_OLD(DIM_I), XDIST_VEC_OLD(DIM_I), &
        YDIST_SC_OLD(DIM_J), YDIST_VEC_OLD(DIM_J),&
        ZDIST_SC_OLD(DIM_K), ZDIST_VEC_OLD(DIM_K)
      INTEGER &
        NMAX_OLD(0:DIMENSION_M), IMAX2_OLD, JMAX2_OLD, KMAX2_OLD,&
        IJMAX2_OLD, MMAX_OLD, FLAG_OLD(DIMENSION_3) 
      INTEGER &
        I_OLD, J_OLD, K_OLD, IV_OLD, JV_OLD, KV_OLD, IJK_OLD, &
        IM_OLD, JM_OLD, KM_OLD, IP_OLD, JP_OLD, KP_OLD, &
        IVJK_OLD, IJVK_OLD, IJKV_OLD, I1SAVE , J1SAVE , K1SAVE, &
        NSTEP_OLD, L
      LOGICAL EXT_I, EXT_J, EXT_K, DONE, SHIFT, INTERNAL

      INTEGER I, J, K, IJK, M, N
!
      INCLUDE 'function.inc'
!
      WRITE(*,*)' Processing data. Please wait. '



      Allocate(  EP_g_OLD(DIMENSION_3) )
      Allocate(  P_g_OLD(DIMENSION_3) )
      Allocate(  P_star_OLD(DIMENSION_3) )
      Allocate(  RO_g_OLD(DIMENSION_3) )
      Allocate(  ROP_g_OLD(DIMENSION_3) )
      Allocate(  U_g_OLD(DIMENSION_3) )
      Allocate(  V_g_OLD(DIMENSION_3) )
      Allocate(  W_g_OLD(DIMENSION_3) )
      Allocate(  T_g_OLD(DIMENSION_3) )
      Allocate(  GAMA_RG_OLD(DIMENSION_3) )
      Allocate(  T_RG_OLD(DIMENSION_3) )
      
      Allocate(  X_g_OLD(DIMENSION_3, DIMENSION_N_g) )
      
      if(Nscalar > 0) Allocate( Scalar_OLD(DIMENSION_3, NScalar) )
      
      Allocate(  T_s_OLD(DIMENSION_3, DIMENSION_M) )
      Allocate(  THETA_M_OLD(DIMENSION_3, DIMENSION_M) )
      Allocate(  ROP_s_OLD(DIMENSION_3, DIMENSION_M) )
      Allocate(  U_s_OLD(DIMENSION_3, DIMENSION_M) )
      Allocate(  V_s_OLD(DIMENSION_3, DIMENSION_M) )
      Allocate(  W_s_OLD(DIMENSION_3, DIMENSION_M) )
      Allocate(  GAMA_RS_OLD(DIMENSION_3, DIMENSION_M) )
      Allocate(  T_RS_OLD(DIMENSION_3, DIMENSION_M) )
      
      Allocate(  X_s_OLD(DIMENSION_3, DIMENSION_M, DIMENSION_N_s) )





      EXT_I = .FALSE.
      EXT_J = .FALSE.
      EXT_K = .FALSE.
      I1SAVE = I1
      J1SAVE = J1
      K1SAVE = K1
!
!  Read old RES file
!
      CALL READ_RES1
!
!  Save old values
!
      DO 100 K = 1, KMAX2
      DO 100 J = 1, JMAX2
      DO 100 I = 1, IMAX2
        IJK = FUNIJK(I, J, K)
        FLAG_OLD(IJK)    = FLAG(IJK)
        EP_g_OLD(IJK)    = EP_g(IJK)
        P_g_OLD(IJK)     = P_g(IJK)
        P_star_OLD(IJK)  = P_star(IJK)
        RO_g_OLD(IJK)    = RO_g(IJK)
        ROP_g_OLD(IJK)   = ROP_g(IJK)
        T_g_OLD(IJK)     = T_g(IJK)
        DO 80 N = 1, NMAX(0)
          X_g_OLD(IJK, N) = X_g(IJK, N)
80      CONTINUE
        U_g_OLD(IJK)     = U_g(IJK)
        V_g_OLD(IJK)     = V_g(IJK)
        W_g_OLD(IJK)     = W_g(IJK)
        DO 95 M = 1, MMAX
          ROP_s_OLD(IJK, M) = ROP_s(IJK, M)
          T_s_OLD(IJK, M)    = T_s(IJK, M)
          U_s_OLD(IJK, M)   = U_s(IJK, M)
          V_s_OLD(IJK, M)   = V_s(IJK, M)
          W_s_OLD(IJK, M)   = W_s(IJK, M)
	  Theta_m_OLD(IJK, M) = Theta_m(IJK, M)
          DO 85 N = 1, NMAX(M)
            X_s_OLD(IJK, M, N) = X_s(IJK, M, N)
85        CONTINUE
95      CONTINUE
!
!       Version 1.3

        DO LC = 1, NScalar 
          Scalar_OLD (IJK, LC) = Scalar (IJK, LC) 
        END DO
!
!       Version 1.4 -- write radiation variables in write_res1 
        GAMA_RG_OLD(IJK) = GAMA_RG(IJK) 
        T_RG_OLD(IJK) = T_RG(IJK) 

        DO LC = 1, MMAX 
          GAMA_RS_OLD(IJK, LC) = GAMA_RS(IJK, LC)
          T_RS_OLD(IJK, LC) = T_RS(IJK, LC)
        ENDDO 
100   CONTINUE

      DO 105 I = 1, IMAX2
        XDIST_SC_OLD(I)  = XDIST_SC(I)
        XDIST_VEC_OLD(I) = XDIST_VEC(I)
105   CONTINUE
      DO 110 J = 1, JMAX2
        YDIST_SC_OLD(J)  = YDIST_SC(J)
        YDIST_VEC_OLD(J) = YDIST_VEC(J)
110   CONTINUE
      DO 115 K = 1, KMAX2
        ZDIST_SC_OLD(K)  = ZDIST_SC(K)
        ZDIST_VEC_OLD(K) = ZDIST_VEC(K)
115   CONTINUE
      IMAX2_OLD  = IMAX2
      JMAX2_OLD  = JMAX2
      KMAX2_OLD  = KMAX2
      IJMAX2_OLD = IJMAX2
      MMAX_OLD   = MMAX
      DO 120 M = 0, MMAX
        NMAX_OLD(M) = NMAX(M)
120   CONTINUE
      TIME_OLD = TIME
      NSTEP_OLD = NSTEP
!
!  Read the new data file
!
      CALL DEALLOCATE_ARRAYS
      
      CALL INIT_NAMELIST
      CALL READ_NAMELIST(1)
      
      CALL ALLOCATE_ARRAYS
!
!  Do initial calculations
!
      SHIFT = .TRUE.
      CALL CHECK_DATA_03(SHIFT)
      CALL SET_GEOMETRY
      CALL CALC_DISTANCE (XMIN,DX,IMAX2,XDIST_SC,XDIST_VEC)
      CALL CALC_DISTANCE (ZERO,DY,JMAX2,YDIST_SC,YDIST_VEC)
      CALL CALC_DISTANCE (ZERO,DZ,KMAX2,ZDIST_SC,ZDIST_VEC)
      
      CALL CHECK_DATA_04                         ! solid phase section 
      CALL CHECK_DATA_05                         ! gas phase section 
!
!  Open new RES files
!
      CLOSE(UNIT_RES)
      IF (.NOT.DO_XFORMS) THEN
         WRITE(*,'(/A,$)') ' Enter a new RUN_NAME > '
         READ(*,'(A)')RUN_NAME
      ELSE
         RUN_NAME = TEMP_FILE(1:60)
         DO I = 1,60
            IF (RUN_NAME(I:I).EQ.CHAR(0)) RUN_NAME(I:I) = ' '
         END DO
      END IF
      CALL MAKE_UPPER_CASE(RUN_NAME, 60)
200   DONE = OPEN_FILEP(RUN_NAME, 'NEW', 0)
      IF(.NOT.DONE) THEN
        WRITE(*,'(/A,$)') &
          ' Unable to open RES file. Enter a new RUN_NAME > '
        READ(*,'(A)')RUN_NAME
        GOTO 200
      ENDIF
!
!  Interpolate values for the new grid
!
      DO 500 K = 1, KMAX2
      DO 500 J = 1, JMAX2
      DO 500 I = 1, IMAX2
        IJK = FUNIJK(I, J, K)
	INTERNAL = .TRUE.
!
!  compute I, J, and K for the old coordinate system
!
        IF(I .EQ. 1)THEN
          I_OLD = 1
          IV_OLD = 1
	  IF(.NOT.NO_I) INTERNAL = .FALSE.
	ELSEIF(I .EQ. IMAX2) THEN
          I_OLD  = IMAX2_OLD
          IV_OLD = IMAX2_OLD
	  IF(.NOT.NO_I) INTERNAL = .FALSE.
	ELSE
          I_OLD = GET_INDEX &
               (XDIST_SC(I), XDIST_SC_OLD, IMAX2_OLD, EXT_I, I1,'X')
          IV_OLD= GET_INDEX &
               (XDIST_VEC(I), XDIST_VEC_OLD, IMAX2_OLD, EXT_I, I1,'X_E')
	ENDIF
	
        IF(J .EQ. 1)THEN
          J_OLD = 1
          JV_OLD = 1
	  IF(.NOT.NO_J) INTERNAL = .FALSE.
	ELSEIF(J .EQ. JMAX2) THEN
          J_OLD = JMAX2_OLD
          JV_OLD = JMAX2_OLD
	  IF(.NOT.NO_J) INTERNAL = .FALSE.
	ELSE
          J_OLD = GET_INDEX &
               (YDIST_SC(J), YDIST_SC_OLD, JMAX2_OLD, EXT_J, J1,'Y')
          JV_OLD= GET_INDEX &
               (YDIST_VEC(J), YDIST_VEC_OLD, JMAX2_OLD, EXT_J, J1,'Y_N')
	ENDIF
	
        IF(K .EQ. 1)THEN
          K_OLD = 1
          KV_OLD = 1
	  IF(.NOT.NO_K) INTERNAL = .FALSE.
	ELSEIF(K .EQ. KMAX2) THEN
          K_OLD = KMAX2_OLD
          KV_OLD = KMAX2_OLD
	  IF(.NOT.NO_K) INTERNAL = .FALSE.
	ELSE
          K_OLD = GET_INDEX &
               (ZDIST_SC(K), ZDIST_SC_OLD, KMAX2_OLD, EXT_K, K1,'Z')
          KV_OLD= GET_INDEX &
               (ZDIST_VEC(K), ZDIST_VEC_OLD, KMAX2_OLD, EXT_K, K1,'Z_T')
	ENDIF
	
        IJK_OLD  = I_OLD + (J_OLD - 1) * IMAX2_OLD &
                  + (K_OLD - 1) * IJMAX2_OLD
!
!  If the old IJK location is a wall cell, search the near by cells
!  for a non-wall cell.  Although this interpolation will not be
!  accurate, it is essential for restarting a run, since non-zero values
!  are required for quantities such as void fraction, temperature etc.
!
        IF(FLAG_OLD(IJK_OLD) .GE. 100 .AND. INTERNAL) THEN
          DO 380 L = 1, 1000
            IM_OLD = MAX((I_OLD - L), 1)
            IJK_OLD  = IM_OLD + (J_OLD - 1) * IMAX2_OLD &
                  + (K_OLD - 1) * IJMAX2_OLD
            IF(FLAG_OLD(IJK_OLD) .LT. 100)EXIT
            IP_OLD = MIN((I_OLD + L), IMAX2_OLD)
            IJK_OLD  = IP_OLD + (J_OLD - 1) * IMAX2_OLD &
                  + (K_OLD - 1) * IJMAX2_OLD
            IF(FLAG_OLD(IJK_OLD) .LT. 100)EXIT
            JM_OLD = MAX((J_OLD - L), 1)
            IJK_OLD  = I_OLD + (JM_OLD - 1) * IMAX2_OLD &
                  + (K_OLD - 1) * IJMAX2_OLD
            IF(FLAG_OLD(IJK_OLD) .LT. 100)EXIT
            JP_OLD = MIN((J_OLD + L), JMAX2_OLD)
            IJK_OLD  = I_OLD + (JP_OLD - 1) * IMAX2_OLD &
                  + (K_OLD - 1) * IJMAX2_OLD
            IF(FLAG_OLD(IJK_OLD) .LT. 100)EXIT
            KM_OLD = MAX((K_OLD - L), 1)
            IJK_OLD  = I_OLD + (J_OLD - 1) * IMAX2_OLD &
                  + (KM_OLD - 1) * IJMAX2_OLD
            IF(FLAG_OLD(IJK_OLD) .LT. 100)EXIT
            KP_OLD = MIN((K_OLD + L), KMAX2_OLD)
            IJK_OLD  = I_OLD + (J_OLD - 1) * IMAX2_OLD &
                  + (KP_OLD - 1) * IJMAX2_OLD
            IF(FLAG_OLD(IJK_OLD) .LT. 100)EXIT
380       CONTINUE
        ENDIF
	
        IVJK_OLD = IV_OLD + (J_OLD - 1) * IMAX2_OLD &
                  + (K_OLD - 1) * IJMAX2_OLD
        IJVK_OLD = I_OLD + (JV_OLD - 1) * IMAX2_OLD &
                  + (K_OLD - 1) * IJMAX2_OLD
        IJKV_OLD = I_OLD + (J_OLD - 1) * IMAX2_OLD &
                  + (KV_OLD - 1) * IJMAX2_OLD
!
!  Set the values for the new arrays
!
        EP_g(IJK)    = EP_g_OLD(IJK_OLD)
        P_g(IJK)     = P_g_OLD(IJK_OLD)
        P_star(IJK)  = P_star_OLD(IJK_OLD)
        RO_g(IJK)    = RO_g_OLD(IJK_OLD)
        ROP_g(IJK)   = ROP_g_OLD(IJK_OLD)
        T_g(IJK)     = T_g_OLD(IJK_OLD)
        DO 400 N = 1, NMAX(0)
          IF(N .LE. NMAX_OLD(0))THEN
            X_g(IJK, N) = X_g_OLD(IJK_OLD, N)
          ELSE
            X_g(IJK, N) = ZERO
          ENDIF
400     CONTINUE
        IF(U_g_OLD(IVJK_OLD) .NE. UNDEFINED) THEN
          U_g(IJK)     = U_g_OLD(IVJK_OLD)
        ELSE
          U_g(IJK)     = ZERO
        ENDIF
        IF(V_g_OLD(IJVK_OLD) .NE. UNDEFINED) THEN
          V_g(IJK)     = V_g_OLD(IJVK_OLD)
        ELSE
          V_g(IJK)     = ZERO
        ENDIF
        IF(W_g_OLD(IJKV_OLD) .NE. UNDEFINED) THEN
          W_g(IJK)     = W_g_OLD(IJKV_OLD)
        ELSE
          W_g(IJK)     = ZERO
        ENDIF
        DO 450 M = 1, MMAX
          IF(M .LE. MMAX_OLD) THEN
            ROP_s(IJK, M) = ROP_s_OLD(IJK_OLD, M)
            T_s(IJK, M)    = T_s_OLD(IJK_OLD, M)
            IF(U_s_OLD(IVJK_OLD, M) .NE. UNDEFINED) THEN
              U_s(IJK, M)   = U_s_OLD(IVJK_OLD, M)
            ELSE
              U_s(IJK, M)   = ZERO
            ENDIF
            IF(V_s_OLD(IJVK_OLD, M) .NE. UNDEFINED) THEN
              V_s(IJK, M)   = V_s_OLD(IJVK_OLD, M)
            ELSE
              V_s(IJK, M)   = ZERO
            ENDIF
            IF(W_s_OLD(IJKV_OLD, M) .NE. UNDEFINED) THEN
              W_s(IJK, M)   = W_s_OLD(IJKV_OLD, M)
            ELSE
              W_s(IJK, M)   = ZERO
            ENDIF
            IF(Theta_m_OLD(IJKV_OLD, M) .NE. UNDEFINED) THEN
              Theta_m(IJK, M)   = Theta_m_OLD(IJKV_OLD, M)
            ELSE
              Theta_m(IJK, M)   = ZERO
            ENDIF
            DO 420 N = 1, NMAX(M)
              IF(N .LE. NMAX_OLD(M)) THEN
                X_s(IJK, M, N) = X_s_OLD(IJK_OLD, M, N)
              ELSE
                X_s(IJK, M, N) = ZERO
              ENDIF
420         CONTINUE
          ELSE
            ROP_s(IJK, M)  = ZERO
            U_s(IJK, M)    = ZERO
            V_s(IJK, M)    = ZERO
            W_s(IJK, M)    = ZERO
            Theta_m(IJK, M)   = ZERO
            X_s(IJK, M, 1) = ONE
            DO 440 N = 2, NMAX(M)
              X_s(IJK, M, N) = ZERO
440         CONTINUE
          ENDIF
450     CONTINUE
!
!       Version 1.3

        DO LC = 1, NScalar 
          Scalar (IJK, LC) = Scalar_OLD (IJK_OLD, LC) 
        END DO
!
!       Version 1.4 -- write radiation variables in write_res1 
        GAMA_RG(IJK) = GAMA_RG_OLD(IJK_OLD)
        T_RG(IJK) = T_RG_OLD(IJK_OLD)

        DO LC = 1, MMAX 
          GAMA_RS(IJK, LC) = GAMA_RS_OLD(IJK_OLD, LC) 
          T_RS(IJK, LC) = T_RS_OLD(IJK_OLD, LC)
        ENDDO 

500   CONTINUE
!
!  Write the new RES file
!
      IF (.NOT.DO_XFORMS) THEN
         WRITE(*,'(/A,$)') &
               ' Do you need time to be reset to 0.0 (Y/N) ?'
         READ(*,'(A)')RUN_NAME
      ELSE
         RUN_NAME = 'N'
         IF (RESET_TIME) RUN_NAME = 'Y'
      END IF
      IF(RUN_NAME(1:1) .EQ. 'Y' .OR. RUN_NAME(1:1) .EQ. 'y') THEN 
        TIME = ZERO
        NSTEP = 0
      ELSE
        TIME  = TIME_OLD
        NSTEP = NSTEP_OLD
       ENDIF
      CALL WRITE_RES0
      CALL WRITE_RES1
!
      WRITE(*,*)' New RES file written.  Start the new run by setting'
      WRITE(*,*)' RUN_TYPE = RESTART_2'
!
! SET IN CASE CODE IS CHANGED TO RETURN INSTEAD OF STOP
!
      I1 = I1SAVE
      J1 = J1SAVE
      K1 = K1SAVE
!
      STOP
      END
