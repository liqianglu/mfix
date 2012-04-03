!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_CARTESIAN                                   C
!  Purpose: check the data related to cartesian grid implementation    C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_DATA_CARTESIAN
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE constant 
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar             
      USE mpi_utility        
      USE bc
      USE DISCRETELEMENT

      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard
      USE stl


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
      INTEGER :: G,I,J,IJK,Q
      Character*80  Line(1)
      DOUBLE PRECISION :: norm, tan_half_angle
      CHARACTER(LEN=9) :: GR
!-----------------------------------------------
!

      IF(.NOT.CARTESIAN_GRID) RETURN

!      IF(GRANULAR_ENERGY) THEN
!         WRITE(*,*)'INPUT ERROR: CARTESIAN GRID OPTION NOT CURRENTLY'
!         WRITE(*,*)'AVALAIBLE WHEN SOLVING GRANULAR ENERGY EQUATION.'
!         WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
!         CALL MFIX_EXIT(MYPE)
!      ENDIF

      IF(DISCRETE_ELEMENT) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: CARTESIAN GRID OPTION NOT CURRENTLY'
            WRITE(*,*)'AVALAIBLE WITH DISCRETE ELEMENT MODEL.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(COORDINATES=='CYLINDRICAL') THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: CARTESIAN GRID OPTION NOT AVAILABLE'
            WRITE(*,*)'WITH CYLINDRICAL COORDINATE SYSTEM.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(USE_STL.AND.(.NOT.USE_MSH)) THEN
         IF(DO_K) THEN 
            CALL GET_STL_DATA
         ELSE
            IF(MyPE == PE_IO) WRITE(*,*) 'ERROR: STL METHOD VALID ONLY IN 3D.'
            CALL MFIX_EXIT(MYPE) 
         ENDIF
         IF(N_QUADRIC > 0) THEN
            IF(MyPE == PE_IO) WRITE(*,*) 'ERROR: BOTH QUADRIC(S) AND STL INPUT ARE SPECIFIED.'
            IF(MyPE == PE_IO) WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE OF SURFACE INPUT.'
            CALL MFIX_EXIT(MYPE) 
         ENDIF
!         IF(STL_BC_ID == UNDEFINED_I) THEN
!            IF(MyPE == PE_IO) WRITE(*,*) 'ERROR: STL_BC_ID NOT DEFINED.'
!            CALL MFIX_EXIT(MYPE) 
!         ENDIF
      ENDIF

      IF(USE_MSH.AND.(.NOT.USE_STL)) THEN
         IF(DO_K) THEN 
            CALL GET_MSH_DATA
         ELSE
            IF(MyPE == PE_IO) WRITE(*,*) 'ERROR: MSH METHOD VALID ONLY IN 3D.'
            CALL MFIX_EXIT(MYPE) 
         ENDIF
         IF(N_QUADRIC > 0) THEN
            IF(MyPE == PE_IO) WRITE(*,*) 'ERROR: BOTH QUADRIC(S) AND MSH INPUT ARE SPECIFIED.'
            IF(MyPE == PE_IO) WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE OF SURFACE INPUT.'
            CALL MFIX_EXIT(MYPE) 
         ENDIF
      ENDIF

      IF(USE_POLYGON) THEN
         IF(DO_K) THEN 
            IF(MyPE == PE_IO) WRITE(*,*) 'ERROR: POLYGON METHOD VALID ONLY IN 2D.'
            CALL MFIX_EXIT(MYPE) 
         ELSE
            CALL GET_POLY_DATA
         ENDIF
      ENDIF

      IF(N_QUADRIC > 0) THEN
         IF(N_POLYGON > 0) THEN 
            IF(MyPE == PE_IO) THEN
               WRITE(*,*) 'ERROR: BOTH QUADRIC(S) AND POLYGON(S) DEFINED.'
               WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE OF SURFACE INPUT.'
            ENDIF
            CALL MFIX_EXIT(MYPE) 
         ENDIF
         IF(N_USR_DEF > 0) THEN 
            IF(MyPE == PE_IO) THEN
               WRITE(*,*) 'ERROR: BOTH QUADRIC(S) AND USER-DEFINED FUNTION DEFINED.'
               WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE OF SURFACE.'
            ENDIF
            CALL MFIX_EXIT(MYPE) 
         ENDIF
         IF(QUADRIC_SCALE <= ZERO) THEN 
            IF(MyPE == PE_IO) THEN
               WRITE(*,*) 'ERROR: QUADRIC_SCALE MUST BE POSITIVE.'
            ENDIF
            CALL MFIX_EXIT(MYPE) 
         ELSEIF(QUADRIC_SCALE /= ONE) THEN
            DO Q = 1, N_QUADRIC
               lambda_x(Q)  = lambda_x(Q)  * quadric_scale**2
               lambda_y(Q)  = lambda_y(Q)  * quadric_scale**2
               lambda_z(Q)  = lambda_z(Q)  * quadric_scale**2
               Radius(Q)    = Radius(Q)    * quadric_scale
               t_x(Q)       = t_x(Q)       * quadric_scale
               t_y(Q)       = t_y(Q)       * quadric_scale
               t_z(Q)       = t_z(Q)       * quadric_scale
               clip_xmin(Q) = clip_xmin(Q) * quadric_scale
               clip_xmax(Q) = clip_xmax(Q) * quadric_scale
               clip_ymin(Q) = clip_ymin(Q) * quadric_scale
               clip_ymax(Q) = clip_ymax(Q) * quadric_scale
               clip_zmin(Q) = clip_zmin(Q) * quadric_scale
               clip_zmax(Q) = clip_zmax(Q) * quadric_scale
               piece_xmin(Q) = piece_xmin(Q) * quadric_scale
               piece_xmax(Q) = piece_xmax(Q) * quadric_scale
               piece_ymin(Q) = piece_ymin(Q) * quadric_scale
               piece_ymax(Q) = piece_ymax(Q) * quadric_scale
               piece_zmin(Q) = piece_zmin(Q) * quadric_scale
               piece_zmax(Q) = piece_zmax(Q) * quadric_scale
            ENDDO
         ENDIF
      ELSE
         IF((N_POLYGON > 0).AND.(N_USR_DEF > 0)) THEN 
            IF(MyPE == PE_IO) THEN
               WRITE(*,*) 'ERROR: POLYGON(S) AND USER-DEFINED FUNTION DEFINED.'
               WRITE(*,*) 'MFIX HANDLES ONLY ONE TYPE OF SURFACE.'
            ENDIF
            CALL MFIX_EXIT(MYPE) 
         ENDIF
      ENDIF

      
      IF(N_QUADRIC > DIM_QUADRIC) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF N_QUADRIC =', N_QUADRIC
            WRITE(*,*)'MAXIMUM ACCEPTABLE VALUE IS DIM_QUADRIC =', DIM_QUADRIC
            WRITE(*,*)'CHANGE MAXIMUM VALUE IN QUADRIC_MOD.F IF NECESSARY.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF


      DO Q = 1, N_QUADRIC



         SELECT CASE (TRIM(QUADRIC_FORM(Q)))
    
            CASE ('NORMAL')       

               lambda_x(Q) = lambda_x(Q)
               lambda_y(Q) = lambda_y(Q)
               lambda_z(Q) = lambda_z(Q)
               
               norm = dsqrt(lambda_x(Q)**2 + lambda_y(Q)**2 + lambda_z(Q)**2)

               IF(norm < TOL_F) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: QUADRIC:', Q, ' HAS ZERO COEFFICIENTS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)             
               ENDIF

            CASE ('PLANE')   ! The quadric is predefined as a plane

               lambda_x(Q) = n_x(Q)
               lambda_y(Q) = n_y(Q)
               lambda_z(Q) = n_z(Q)

               norm = dsqrt(lambda_x(Q)**2 + lambda_y(Q)**2 + lambda_z(Q)**2)

               IF( norm > TOL_F) THEN
                  lambda_x(Q) = lambda_x(Q) / norm
                  lambda_y(Q) = lambda_y(Q) / norm
                  lambda_z(Q) = lambda_z(Q) / norm
               ELSE
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: PLANE:', Q, ' HAS ZERO NORMAL VECTOR.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)             
               ENDIF

              dquadric(Q) = - (lambda_x(Q)*t_x(Q) + lambda_y(Q)*t_y(Q) + lambda_z(Q)*t_z(Q))

            CASE ('X_CYL_INT')   ! The quadric is predefined as a cylinder, along x-axis
                                 ! Internal flow

               IF( Radius(Q) <= ZERO) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, ' HAS ZERO RADIUS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = ZERO
                  lambda_y(Q) = ONE
                  lambda_z(Q) = ONE
                  dquadric(Q) = -Radius(Q)**2
               ENDIF

            CASE ('Y_CYL_INT')   ! The quadric is predefined as a cylinder, along y-axis
                                 ! Internal flow

               IF( Radius(Q) <= ZERO) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, ' HAS ZERO RADIUS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = ONE
                  lambda_y(Q) = ZERO
                  lambda_z(Q) = ONE
                  dquadric(Q) = -Radius(Q)**2
               ENDIF

            CASE ('Z_CYL_INT')   ! The quadric is predefined as a cylinder, along z-axis
                                 ! Internal flow

               IF( Radius(Q) <= ZERO) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, ' HAS ZERO RADIUS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = ONE
                  lambda_y(Q) = ONE
                  lambda_z(Q) = ZERO
                  dquadric(Q) = -Radius(Q)**2
               ENDIF


            CASE ('X_CYL_EXT')   ! The quadric is predefined as a cylinder, along x-axis
                                 ! External flow

               IF( Radius(Q) <= ZERO) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, ' HAS ZERO RADIUS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = ZERO
                  lambda_y(Q) = -ONE
                  lambda_z(Q) = -ONE
                  dquadric(Q) = Radius(Q)**2
               ENDIF

            CASE ('Y_CYL_EXT')   ! The quadric is predefined as a cylinder, along y-axis
                                 ! External flow

               IF( Radius(Q) <= ZERO) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, ' HAS ZERO RADIUS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = -ONE
                  lambda_y(Q) = ZERO
                  lambda_z(Q) = -ONE
                  dquadric(Q) = Radius(Q)**2
               ENDIF

            CASE ('Z_CYL_EXT')   ! The quadric is predefined as a cylinder, along z-axis
                                 ! External flow

               IF( Radius(Q) <= ZERO) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CYLINDER:', Q, ' HAS ZERO RADIUS.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = -ONE
                  lambda_y(Q) = -ONE
                  lambda_z(Q) = ZERO
                  dquadric(Q) = Radius(Q)**2
               ENDIF

            CASE ('SPHERE_INT')   ! The quadric is predefined as a sphere
                                  ! Internal flow

               IF( Radius(Q) <= ZERO) THEN
                  WRITE(*,*)'INPUT ERROR: SPHERE:', Q, ' HAS INVALID RADIUS.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = ONE
                  lambda_y(Q) = ONE
                  lambda_z(Q) = ONE
                  dquadric(Q) = -Radius(Q)**2
               ENDIF
 
           CASE ('SPHERE_EXT')   ! The quadric is predefined as a sphere
                                  ! External flow

               IF( Radius(Q) <= ZERO) THEN
                  WRITE(*,*)'INPUT ERROR: SPHERE:', Q, ' HAS INVALID RADIUS.'
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  lambda_x(Q) = -ONE
                  lambda_y(Q) = -ONE
                  lambda_z(Q) = -ONE
                  dquadric(Q) = Radius(Q)**2
               ENDIF
         

            CASE ('X_CONE')    ! The quadric is predefined as a cone, along x-axis
                               ! Internal flow

            IF(HALF_ANGLE(Q) <= ZERO .OR. HALF_ANGLE(Q) >= 90.0) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CONE:', Q, ' HAS INCORRECT HALF-ANGLE.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  tan_half_angle = DTAN(HALF_ANGLE(Q)/180.0*PI)
                  lambda_x(Q) = -ONE
                  lambda_y(Q) = ONE/(tan_half_angle)**2
                  lambda_z(Q) = ONE/(tan_half_angle)**2
                  dquadric(Q) = ZERO
               ENDIF

            CASE ('Y_CONE')    ! The quadric is predefined as a cone, along y-axis
                               ! Internal flow

            IF(HALF_ANGLE(Q) <= ZERO .OR. HALF_ANGLE(Q) >= 90.0) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CONE:', Q, ' HAS INCORRECT HALF-ANGLE.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  tan_half_angle = DTAN(HALF_ANGLE(Q)/180.0*PI)
                  lambda_x(Q) = ONE/(tan_half_angle)**2
                  lambda_y(Q) = -ONE
                  lambda_z(Q) = ONE/(tan_half_angle)**2
                  dquadric(Q) = ZERO
               ENDIF

            CASE ('Z_CONE')    ! The quadric is predefined as a cone, along z-axis
                               ! Internal flow

            IF(HALF_ANGLE(Q) <= ZERO .OR. HALF_ANGLE(Q) >= 90.0) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: CONE:', Q, ' HAS INCORRECT HALF-ANGLE.'
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)             
               ELSE
                  tan_half_angle = DTAN(HALF_ANGLE(Q)/180.0*PI)
                  lambda_x(Q) = ONE/(tan_half_angle)**2
                  lambda_y(Q) = ONE/(tan_half_angle)**2
                  lambda_z(Q) = -ONE
                  dquadric(Q) = ZERO
               ENDIF

            CASE ('C2C')        ! Cylinder to cylinder junction using cone
                                ! Internal flow

               CALL BUILD_CONE_FOR_C2C(Q)


            CASE DEFAULT
               IF(MyPE == PE_IO) THEN
                  WRITE(*,*)'INPUT ERROR: QUADRIC:', Q, ' HAS INCORRECT FORM: ',quadric_form(Q)
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
               ENDIF
               CALL MFIX_EXIT(MYPE)             

         END SELECT

         IF(BC_ID_Q(Q) == UNDEFINED_I) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)'INPUT ERROR: QUADRIC:', Q, ' HAS NO ASSIGNED BC ID.'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            CALL MFIX_EXIT(MYPE)
         ENDIF

      ENDDO
 

      IF(N_QUADRIC>0) THEN


         IF(N_GROUP > DIM_GROUP) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)'INPUT ERROR: INVALID VALUE OF N_GROUP =', N_GROUP
               WRITE(*,*)'MAXIMUM ACCEPTABLE VALUE IS DIM_GROUP =', DIM_GROUP
               WRITE(*,*)'CHANGE MAXIMUM VALUE IN QUADRIC_MOD.F IF NECESSARY.'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            CALL MFIX_EXIT(MYPE)
         ENDIF


         DO I = 1,N_GROUP

            IF(GROUP_SIZE(I) < 1 .OR. GROUP_SIZE(I) > N_QUADRIC) THEN
              IF(MyPE == PE_IO) THEN
                  WRITE(*,*)'INPUT ERROR: GROUP:', I, ' HAS INCORRECT SIZE:', GROUP_SIZE(I)
                  WRITE(*,*)'VALID GROUP SIZE RANGE IS:', 1, ' TO ', N_QUADRIC
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
               ENDIF
               CALL MFIX_EXIT(MYPE)
            ENDIF

            DO J = 1,GROUP_SIZE(I)
               IF(GROUP_Q(I,J) < 1 .OR. GROUP_Q(I,J) > N_QUADRIC) THEN
                  IF(MyPE == PE_IO) THEN
                     WRITE(*,*)'INPUT ERROR: GROUP_Q(', I,',',J, ') HAS INCORRECT VALUE:', GROUP_Q(I,J)
                     WRITE(*,*)'VALID GROUP_Q RANGE IS:', 1, ' TO ', N_QUADRIC
                     WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
                  ENDIF
                  CALL MFIX_EXIT(MYPE)
               ENDIF
            ENDDO
   
            GR = TRIM(GROUP_RELATION(I)) 

            IF(GR/='OR'.AND.GR/='AND'.AND.GR/='PIECEWISE') THEN
               IF(MyPE == PE_IO) THEN
                  WRITE(*,*)'INPUT ERROR: GROUP:', I, ' HAS INCORRECT GROUP RELATION: ', GR
                  WRITE(*,*)'VALID GROUP RELATIONS ARE ''OR'',''AND'', AND ''PIECEWISE''. '
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
               ENDIF
               CALL MFIX_EXIT(MYPE)
            ENDIF

         ENDDO

         DO I = 2,N_GROUP

            GR = TRIM(RELATION_WITH_PREVIOUS(I)) 

            IF(GR/='OR'.AND.GR/='AND') THEN
               IF(MyPE == PE_IO) THEN
                  WRITE(*,*)'INPUT ERROR: GROUP:', I, ' HAS INCORRECT RELATION WITH PREVIOUS: ', GR
                  WRITE(*,*)'VALID GROUP RELATIONS ARE ''OR'', AND ''AND''. '
                  WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
               ENDIF
               CALL MFIX_EXIT(MYPE)
            ENDIF
         
         ENDDO

      ENDIF


      IF(TOL_SNAP(1)<ZERO.OR.TOL_SNAP(1)>HALF) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SNAP IN X-DIRECTION =', TOL_SNAP(1)
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 0.5.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF   
             
      IF(TOL_SNAP(2)==UNDEFINED) TOL_SNAP(2)=TOL_SNAP(1)
   
      IF(TOL_SNAP(2)<ZERO.OR.TOL_SNAP(2)>HALF) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SNAP IN Y-DIRECTION =', TOL_SNAP(2)
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 0.5.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF   
                
      IF(TOL_SNAP(3)==UNDEFINED) TOL_SNAP(3)=TOL_SNAP(1)

      IF(TOL_SNAP(3)<ZERO.OR.TOL_SNAP(3)>HALF) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SNAP IN Z-DIRECTION =', TOL_SNAP(3)
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 0.5.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF   


      IF(TOL_DELH<ZERO.OR.TOL_DELH>ONE) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_DELH =', TOL_DELH
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 1.0.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(TOL_SMALL_CELL<ZERO.OR.TOL_SMALL_CELL>ONE) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SMALL_CELL =', TOL_SMALL_CELL
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 1.0.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(TOL_SMALL_AREA<ZERO.OR.TOL_SMALL_AREA>ONE) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF TOL_SMALL_AREA =', TOL_SMALL_AREA
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 1.0.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(ALPHA_MAX<ZERO) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF ALPHA_MAX =', ALPHA_MAX
            WRITE(*,*)'ACCEPTABLE VALUES ARE POSITIVE NUMBERS (E.G. 1.0).'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF


      IF(TOL_F<ZERO) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF TOL_F =', TOL_F
            WRITE(*,*)'ACCEPTABLE VALUES ARE SMALL POSITIVE NUMBERS (E.G. 1.0E-9).'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(TOL_POLY<ZERO) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF TOL_POLY =', TOL_POLY
            WRITE(*,*)'ACCEPTABLE VALUES ARE SMALL POSITIVE NUMBERS (E.G. 1.0E-9).'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(ITERMAX_INT<0) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF ITERMAX_INT =', ITERMAX_INT
            WRITE(*,*)'ACCEPTABLE VALUES ARE LARGE POSITIVE INTEGERS (E.G. 10000).'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(FAC_DIM_MAX_CUT_CELL<0.05.OR.FAC_DIM_MAX_CUT_CELL>ONE) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF FAC_DIM_MAX_CUT_CELL =', FAC_DIM_MAX_CUT_CELL
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.05 AND 1.0.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(.NOT.CARTESIAN_GRID) THEN
         IF(WRITE_VTK_FILES) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)'INPUT ERROR: VTK FILES CAN BE WRITTEN ONLY WHEN CARTESIAN GRID IS ACTIVATED.'
               WRITE(*,*)'PLEASE SET WRITE_VTK_FILES = .FALSE. IN MFIX.DAT AND TRY AGAIN.'
            ENDIF
            CALL MFIX_EXIT(MYPE) 
         ENDIF
      ENDIF


      IF(VTK_DT<ZERO) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: NEGATIVE VALUE OF VTK_DT =', VTK_DT
            WRITE(*,*)'ACCEPTABLE VALUES ARE POSITIVE NUMBERS (E.G. 0.1).'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(FRAME<-1) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INALID VALUE OF FRAME =', FRAME
            WRITE(*,*)'ACCEPTABLE VALUES ARE INTEGERS >= -1.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF((CG_SAFE_MODE(1)==1).AND.(PG_OPTION/=0)) THEN
         PG_OPTION = 0
         IF(MyPE == PE_IO) WRITE(*,*)'WARNING: SAFE_MODE ACTIVATED FOR GAS PRESSURE, REVERTING TO PG_OPTION = 0'
      ENDIF

      IF(PG_OPTION <0 .OR. PG_OPTION>2) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INALID VALUE OF PG_OPTION =', PG_OPTION
            WRITE(*,*)'ACCEPTABLE VALUES ARE 0,1,AND 2.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(CG_UR_FAC(2)<ZERO.OR.CG_UR_FAC(2)>ONE) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF CG_UR_FAC(2) =', CG_UR_FAC(2)
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 1.0.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(BAR_WIDTH<10.OR.BAR_WIDTH>80) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF BAR_WIDTH =', BAR_WIDTH
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 10 AND 80.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(BAR_RESOLUTION<ONE.OR.BAR_RESOLUTION>100.0) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INVALID VALUE OF BAR_RESOLUTION =', BAR_RESOLUTION
            WRITE(*,*)'ACCEPTABLE VALUES ARE BETWEEN 0.0 AND 100.0.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF

      IF(F_DASHBOARD<1) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,*)'INPUT ERROR: INALID VALUE OF F_DASHBOARD =', F_DASHBOARD
            WRITE(*,*)'ACCEPTABLE VALUES ARE INTEGERS >= 1.'
            WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
         ENDIF
         CALL MFIX_EXIT(MYPE)
      ENDIF


!======================================================================
! Data initialization for Dashboard
!======================================================================
      INIT_TIME = TIME
      SMMIN =  LARGE_NUMBER
      SMMAX = -LARGE_NUMBER

      DTMIN =  LARGE_NUMBER
      DTMAX = -LARGE_NUMBER

      NIT_MIN = MAX_NIT
      NIT_MAX = 0

      N_DASHBOARD = 0

      RETURN  
      END SUBROUTINE CHECK_DATA_CARTESIAN



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_BC_FLAGS                                         C
!  Purpose: check the boundary conditions flags                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_BC_FLAGS
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE constant 
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar             
      USE mpi_utility        
      USE bc
      
      USE fldvar
      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard


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
      INTEGER :: I,J,IJK,IJKW,IJKS,IJKB,M,N
      INTEGER :: IJKWW,IJKSS,IJKBB
      INTEGER :: BCV,BCV_U,BCV_V,BCV_W
      Character*80  Line(1)
!-----------------------------------------------
!
      DOUBLE PRECISION SUM, SUM_EP
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      LOGICAL , EXTERNAL :: COMPARE 
!-----------------------------------------------
      INCLUDE 'function.inc'

!======================================================================
! Boundary conditions
!======================================================================

      DO IJK = ijkstart3, ijkend3
         BCV = BC_ID(IJK)
         IF(BCV>0) THEN

            IF(BC_TYPE(BCV)  == 'CG_MI') THEN
 
               FLAG(IJK) = 20
               FLAG_E(IJK) = UNDEFINED_I
               FLAG_N(IJK) = UNDEFINED_I
               FLAG_T(IJK) = UNDEFINED_I

            ELSEIF(BC_TYPE(BCV)  == 'CG_PO') THEN
 
               FLAG(IJK) = 11
               FLAG_E(IJK) = UNDEFINED_I
               FLAG_N(IJK) = UNDEFINED_I
               FLAG_T(IJK) = UNDEFINED_I

               IJKW = WEST_OF(IJK)
               BCV_U = BC_U_ID(IJKW)
               IF(BCV_U>0) THEN
                  IF(BC_TYPE(BCV_U)  == 'CG_PO') THEN
                     FLAG(IJKW) = 11
                     FLAG_E(IJKW) = UNDEFINED_I
                     FLAG_N(IJKW) = UNDEFINED_I
                     FLAG_T(IJKW) = UNDEFINED_I
                  ENDIF
               ENDIF

               IJKS = SOUTH_OF(IJK)
               BCV_V = BC_V_ID(IJKS)
               IF(BCV_V>0) THEN
                  IF(BC_TYPE(BCV_V)  == 'CG_PO') THEN
                     FLAG(IJKS) = 11
                     FLAG_E(IJKS) = UNDEFINED_I
                     FLAG_N(IJKS) = UNDEFINED_I
                     FLAG_T(IJKS) = UNDEFINED_I
                  ENDIF
               ENDIF

               IF(DO_K) THEN
                  IJKB = BOTTOM_OF(IJK)
                  BCV_W = BC_W_ID(IJKB)
                  IF(BCV_W>0) THEN
                     IF(BC_TYPE(BCV_W)  == 'CG_PO') THEN
                        FLAG(IJKB) = 11
                        FLAG_E(IJKB) = UNDEFINED_I
                        FLAG_N(IJKB) = UNDEFINED_I
                        FLAG_T(IJKB) = UNDEFINED_I
                     ENDIF
                  ENDIF
               ENDIF

            ENDIF
         ENDIF
      ENDDO


      DO IJK = ijkstart3, ijkend3
         BCV = BC_ID(IJK)
         IF(BCV>0) THEN
            IF(BC_TYPE(BCV)  == 'CG_MI') THEN

               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN
                  FLAG_E(IJKW) = 2020
               ENDIF           

               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN
                  FLAG_N(IJKS) = 2020
               ENDIF           

               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN
                  FLAG_T(IJKB) = 2020
               ENDIF           

               IF (BC_U_G(BCV) == UNDEFINED) THEN 
                   IF (NO_I) THEN 
                       BC_U_G(BCV) = ZERO 
                   ELSEIF(BC_VOLFLOW_g(BCV)==UNDEFINED.AND. &
                          BC_MASSFLOW_g(BCV)==UNDEFINED.AND.&
                          BC_VELMAG_g(BCV)==UNDEFINED) THEN
                       IF(DMP_LOG)WRITE (UNIT_LOG, 900) 'BC_U_g', BCV 
                       call mfix_exit(myPE)
                   ENDIF 
               ENDIF 
               IF (BC_V_G(BCV) == UNDEFINED) THEN 
                   IF (NO_J) THEN 
                       BC_V_G(BCV) = ZERO 
                   ELSEIF(BC_VOLFLOW_g(BCV)==UNDEFINED.AND. &
                          BC_MASSFLOW_g(BCV)==UNDEFINED.AND.&
                          BC_VELMAG_g(BCV)==UNDEFINED) THEN
                       IF(DMP_LOG)WRITE (UNIT_LOG, 900) 'BC_V_g', BCV 
                       call mfix_exit(myPE)
                   ENDIF 
               ENDIF 
               IF (BC_W_G(BCV) == UNDEFINED) THEN 
                   IF (NO_K) THEN 
                       BC_W_G(BCV) = ZERO 
                   ELSEIF(BC_VOLFLOW_g(BCV)==UNDEFINED.AND. &
                          BC_MASSFLOW_g(BCV)==UNDEFINED.AND.&
                          BC_VELMAG_g(BCV)==UNDEFINED) THEN
                       IF(DMP_LOG)WRITE (UNIT_LOG, 900) 'BC_W_g', BCV 
                       call mfix_exit(myPE)
                   ENDIF 
               ENDIF  
               IF (K_Epsilon .AND. BC_K_Turb_G(BCV) == UNDEFINED) THEN
                   IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_K_Turb_G', BCV 
                   call mfix_exit(myPE) 
               ENDIF   
               IF (K_Epsilon .AND. BC_E_Turb_G(BCV) == UNDEFINED) THEN
                   IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_E_Turb_G', BCV 
                   call mfix_exit(myPE) 
               ENDIF 

!               Check whether the bc velocity components have the correct sign
!               SELECT CASE (BC_PLANE(BCV))  
!               CASE ('W')  
!                   IF (BC_U_G(BCV) > ZERO) THEN 
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_U_g', '<' 
!                       CALL MFIX_EXIT(myPE) 
!                   ENDIF 
!               CASE ('E')  
!                   IF (BC_U_G(BCV) < ZERO) THEN 
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_U_g', '>' 
!                       CALL MFIX_EXIT(myPE) 
!                   ENDIF 
!               CASE ('S')  
!                   IF (BC_V_G(BCV) > ZERO) THEN 
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_V_g', '<' 
!                       CALL MFIX_EXIT(myPE) 
!                   ENDIF 
!               CASE ('N')  
!                   IF (BC_V_G(BCV) < ZERO) THEN 
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_V_g', '>' 
!                       CALL MFIX_EXIT(myPE) 
!                   ENDIF 
!               CASE ('B')  
!                   IF (BC_W_G(BCV) > ZERO) THEN 
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_W_g', '<' 
!                       CALL MFIX_EXIT(myPE) 
!                   ENDIF 
!               CASE ('T')  
!                   IF (BC_W_G(BCV) < ZERO) THEN 
!                       IF(DMP_LOG)WRITE (UNIT_LOG, 1050) BCV, 'BC_W_g', '>' 
!                       CALL MFIX_EXIT(myPE) 
!                   ENDIF 
!               END SELECT 

               SUM_EP = BC_EP_G(BCV) 
               DO M = 1, MMAX 
                  IF (BC_ROP_S(BCV,M) == UNDEFINED) THEN 
                     IF (BC_EP_G(BCV) == ONE) THEN 
                        BC_ROP_S(BCV,M) = ZERO 
                     ELSEIF (MMAX == 1) THEN 
                         BC_ROP_S(BCV,M) = (ONE - BC_EP_G(BCV))*RO_S(M) 
                     ELSE 
                         IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_ROP_s', BCV, M 
                         call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 

                  SUM_EP = SUM_EP + BC_ROP_S(BCV,M)/RO_S(M) 
                  IF (SPECIES_EQ(M)) THEN 
                     SUM = ZERO 
                        DO N = 1, NMAX(M) 
                           IF(BC_X_S(BCV,M,N)/=UNDEFINED)SUM=SUM+BC_X_S(BCV,M,N) 
                        ENDDO 

                     IF (BC_ROP_S(BCV,M)==ZERO .AND. SUM==ZERO) THEN 
                        BC_X_S(BCV,M,1) = ONE 
                        SUM = ONE 
                     ENDIF 

                     DO N = 1, NMAX(M) 
                        IF (BC_X_S(BCV,M,N) == UNDEFINED) THEN 
                           IF(.NOT.COMPARE(ONE,SUM) .AND. DMP_LOG)WRITE (UNIT_LOG,1110)BCV,M,N 
                              BC_X_S(BCV,M,N) = ZERO 
                        ENDIF 
                     ENDDO 

                     IF (.NOT.COMPARE(ONE,SUM)) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1120) BCV, M 
                           call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 

                  IF (BC_U_S(BCV,M) == UNDEFINED) THEN 
                     IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_I) THEN 
                        BC_U_S(BCV,M) = ZERO 
                   ELSEIF(BC_VOLFLOW_s(BCV,M)==UNDEFINED.AND. &
                          BC_MASSFLOW_s(BCV,M)==UNDEFINED.AND.&
                          BC_VELMAG_s(BCV,M)==UNDEFINED) THEN
                        IF(DMP_LOG)WRITE (UNIT_LOG, 910) 'BC_U_s', BCV, M 
                            call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 
                  
                  IF (BC_V_S(BCV,M) == UNDEFINED) THEN 
                     IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_J) THEN 
                        BC_V_S(BCV,M) = ZERO 
                   ELSEIF(BC_VOLFLOW_s(BCV,M)==UNDEFINED.AND. &
                          BC_MASSFLOW_s(BCV,M)==UNDEFINED.AND.&
                          BC_VELMAG_s(BCV,M)==UNDEFINED) THEN
                        IF(DMP_LOG)WRITE (UNIT_LOG, 910) 'BC_V_s', BCV, M 
                            call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 
                  
                  IF (BC_W_S(BCV,M) == UNDEFINED) THEN 
                     IF (BC_ROP_S(BCV,M)==ZERO .OR. NO_K) THEN 
                        BC_W_S(BCV,M) = ZERO 
                   ELSEIF(BC_VOLFLOW_s(BCV,M)==UNDEFINED.AND. &
                          BC_MASSFLOW_s(BCV,M)==UNDEFINED.AND.&
                          BC_VELMAG_s(BCV,M)==UNDEFINED) THEN
                        IF(DMP_LOG)WRITE (UNIT_LOG, 910) 'BC_W_s', BCV, M 
                           call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 

                  IF (ENERGY_EQ .AND. BC_T_S(BCV,M)==UNDEFINED) THEN 
                     IF (BC_ROP_S(BCV,M) == ZERO) THEN 
                        BC_T_S(BCV,M) = BC_T_G(BCV) 
                     ELSE 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_T_s', BCV, M 
                           call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 

                  IF (GRANULAR_ENERGY .AND. BC_THETA_M(BCV,M)==UNDEFINED) THEN 
                     IF (BC_ROP_S(BCV,M) == ZERO) THEN 
                        BC_THETA_M(BCV,M) = ZERO 
                     ELSE 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1100) 'BC_Theta_m', BCV, M 
                          call mfix_exit(myPE)
                     ENDIF 
                  ENDIF 

!                   Check whether the bc velocity components have the correct sign
!                    SELECT CASE (TRIM(BC_PLANE(BCV)))  
!                    CASE ('W')  
!                        IF (BC_U_S(BCV,M) > ZERO) THEN 
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_U_s', M, '<' 
!                            CALL MFIX_EXIT(myPE) 
!                        ENDIF 
!                    CASE ('E')  
!                        IF (BC_U_S(BCV,M) < ZERO) THEN 
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_U_s', M, '>' 
!                            CALL MFIX_EXIT(myPE) 
!                        ENDIF 
!                    CASE ('S')  
!                        IF (BC_V_S(BCV,M) > ZERO) THEN 
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_V_s', M, '<' 
!                            CALL MFIX_EXIT(myPE) 
!                        ENDIF 
!                    CASE ('N')  
!                        IF (BC_V_S(BCV,M) < ZERO) THEN 
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_V_s', M, '>' 
!                            CALL MFIX_EXIT(myPE) 
!                        ENDIF 
!                    CASE ('B')  
!                        IF (BC_W_S(BCV,M) > ZERO) THEN 
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_W_s', M, '<' 
!                            CALL MFIX_EXIT(myPE) 
!                        ENDIF 
!                    CASE ('T')  
!                        IF (BC_W_S(BCV,M) < ZERO) THEN 
!                            IF(DMP_LOG)WRITE (UNIT_LOG, 1150) BCV, 'BC_W_s', M, '>' 
!                            CALL MFIX_EXIT(myPE) 
!                        ENDIF 
!                    END SELECT 


               ENDDO 

               IF (.NOT.COMPARE(ONE,SUM_EP)) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1125) BCV 
                     call mfix_exit(myPE)  
               ENDIF 
       
               DO N = 1, NScalar
                  IF (BC_Scalar(BCV,N) == UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1004) 'BC_Scalar', BCV, N 
                        CALL MFIX_EXIT(myPE)
                  ENDIF 
               ENDDO


            ELSEIF(BC_TYPE(BCV)  == 'CG_PO') THEN

               IJKW = WEST_OF(IJK)
               IF(FLUID_AT(IJKW)) THEN
                  FLAG_E(IJKW) = 2011
               ENDIF           

               BCV_U = BC_U_ID(IJKW)
               IF(BCV_U>0) THEN
                  IF(BC_TYPE(BCV_U)  == 'CG_PO') THEN
                    IJKWW = WEST_OF(IJKW)
                    IF(FLUID_AT(IJKWW)) THEN
                       FLAG_E(IJKWW) = 2011
                    ENDIF           
                  ENDIF
               ENDIF

               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN
                  FLAG_N(IJKS) = 2011
               ENDIF           

               BCV_V = BC_V_ID(IJKS)
               IF(BCV_V>0) THEN
                  IF(BC_TYPE(BCV_V)  == 'CG_PO') THEN
                    IJKSS = SOUTH_OF(IJKS)
                    IF(FLUID_AT(IJKSS)) THEN
                       FLAG_N(IJKSS) = 2011
                    ENDIF           
                  ENDIF
               ENDIF


               IF(DO_K) THEN
                  IJKB = BOTTOM_OF(IJK)
                  IF(FLUID_AT(IJKB)) THEN
                     FLAG_T(IJKB) = 2011
                  ENDIF 

                  BCV_W = BC_W_ID(IJKB)
                  IF(BCV_W>0) THEN
                     IF(BC_TYPE(BCV_W)  == 'CG_PO') THEN
                       IJKBB = BOTTOM_OF(IJKB)
                       IF(FLUID_AT(IJKBB)) THEN
                          FLAG_T(IJKBB) = 2011
                       ENDIF           
                     ENDIF
                  ENDIF

               ENDIF

               IF (BC_P_G(BCV) == UNDEFINED) THEN 
                   IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'BC_P_g', BCV 
                   call mfix_exit(myPE)  
               ELSEIF (BC_P_G(BCV)<=ZERO .AND. RO_G0==UNDEFINED) THEN 
                   IF(DMP_LOG)WRITE (UNIT_LOG, 1010) BCV, BC_P_G(BCV) 
                   call mfix_exit(myPE)  
               ENDIF 

            ENDIF

         ENDIF

      ENDDO

      RETURN  


 900 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,&
         ') not specified',/1X,'One of the following must be specified:',/1X,&
         'BC_VOLFLOW_g, BC_MASSFLOW_g or BC_VELMAG_g',/1X,70('*')/) 

 910 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I1,&
         ') not specified',/1X,'One of the following must be specified:',/1X,&
         'BC_VOLFLOW_g, BC_MASSFLOW_g or BC_VELMAG_g',/1X,70('*')/)

 1000 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,&
         ') not specified',/1X,70('*')/) 
 1001 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/&
         ' Message: Illegal BC_TYPE for BC # = ',I2,/'   BC_TYPE = ',A,/&
         '  Valid BC_TYPE are: ') 
 1002 FORMAT(5X,A16) 
 1003 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,&
         ') value is unphysical',/1X,70('*')/) 
 1004 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I2,&
         ') not specified',/1X,70('*')/) 
 1005 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I2,&
         ') value is unphysical',/1X,70('*')/) 
 1010 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC_P_g( ',I2,&
         ') = ',G12.5,/&
         ' Pressure should be greater than zero for compressible flow',/1X,70(&
         '*')/) 
 1050 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - ',A,' should be ',A,' zero.',/1X,70('*')/) 
 1060 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC_X_g(',I2,',',I2&
         ,') not specified',/1X,70('*')/) 
 1065 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - Sum of gas mass fractions is NOT equal to one',/1X,70('*')/) 
 1100 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I1,&
         ') not specified',/1X,70('*')/) 
 1103 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I1,&
         ') value is unphysical',/1X,70('*')/) 
 1104 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I2,&
         ',',I2,') not specified',/1X,70('*')/) 
 1105 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I2,&
         ',',I2,') value is unphysical',/1X,70('*')/) 
 1110 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC_X_s(',I2,',',I2&
         ,',',I2,') not specified',/1X,70('*')/) 
 1120 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - Sum of solids-',I1,' mass fractions is NOT equal to one',/1X,70(&
         '*')/) 
 1125 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - Sum of volume fractions is NOT equal to one',/1X,70('*')/) 
 1150 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: BC number:',I2,&
         ' - ',A,I1,' should be ',A,' zero.',/1X,70('*')/) 
 1160 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/&
         ' Message: Boundary condition no', &
         I2,' is a second outflow condition.',/1X,&
         '  Only one outflow is allowed.  Consider using P_OUTFLOW.',/1X, 70('*')/) 
 1200 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,&
         ') specified',' for an undefined BC location',/1X,70('*')/) 
 1300 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/' Message: ',A,'(',I2,',',I1,&
         ') specified',' for an undefined BC location',/1X,70('*')/) 
 1400 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/&
         ' Message: No initial or boundary condition specified',/&
         '    I       J       K') 
 1410 FORMAT(I5,3X,I5,3X,I5) 
 1420 FORMAT(/1X,70('*')/) 

 1500 FORMAT(/1X,70('*')//' From: CHECK_BC_FLAGS',/&
         ' Message: No initial or boundary condition specified',/&
         '    I       J       K') 


      END SUBROUTINE CHECK_BC_FLAGS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALL BUILD_CONE_FOR_C2C                                C
!  Purpose: Define cone parameters for Cylider to Cylinder junction.   C
!           The C2C quadric ID must be between the two cylinders ID    C
!           (e.g., if Quadric 4 is a C2C, then Quadrics 3 and 5        C
!            must be cylinders). The two cylinders must be aligned     C
!           in the same direction and be clipped to define the extent  C
!           of the conical junction.                                   C
!           This method is currentl available for internal flow only.  C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 02-Dec-10  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE BUILD_CONE_FOR_C2C(Q)
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE constant 
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar             
      USE mpi_utility        
      USE bc
      USE DISCRETELEMENT

      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard
      USE stl


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
      INTEGER :: Q,QM1,QP1
      Character*80  Line(1)
      DOUBLE PRECISION :: x1,x2,y1,y2,z1,z2,R1,R2
      DOUBLE PRECISION :: norm, tan_half_angle
      LOGICAL :: aligned
!-----------------------------------------------
!

      QM1 = Q-1
      QP1 = Q+1

      IF(MyPE == PE_IO) THEN
         WRITE(*,*)' INFO FOR QUADRIC', Q
         WRITE(*,*)' Defining Cone for Cylinder to Cylinder junction'
         WRITE(*,*)' Between Quadrics ',QM1,' AND ', QP1
      ENDIF


      IF((TRIM(QUADRIC_FORM(QM1))=='X_CYL_INT').AND.(TRIM(QUADRIC_FORM(QP1))=='X_CYL_INT')) THEN

         QUADRIC_FORM(Q) = 'X_CONE'

         aligned = (t_y(QM1)==t_y(QP1)).AND.(t_z(QM1)==t_z(QP1)) 
         IF(.NOT.aligned) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT ALIGNED'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         R1 = RADIUS(QM1)
         R2 = RADIUS(QP1)
         IF(R1==R2) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' HAVE THE SAME RADIUS'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         x1 = piece_xmax(QM1)
         x2 = piece_xmin(QP1)
         IF(x2<=x1) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT PIECED PROPERLY'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         tan_half_angle = (R2-R1)/(x2-x1)

         HALF_ANGLE(Q) = DATAN(tan_half_angle)/PI*180.0D0
         lambda_x(Q) = -ONE
         lambda_y(Q) = ONE/(tan_half_angle)**2
         lambda_z(Q) = ONE/(tan_half_angle)**2 
         dquadric(Q) = ZERO

         piece_xmin(Q) = x1
         piece_xmax(Q) = x2

         t_x(Q) = x1 - R1/tan_half_angle 
         t_y(Q) = t_y(QM1)
         t_z(Q) = t_z(QM1)

         IF(MyPE == PE_IO) THEN
            WRITE(*,*) ' QUADRIC:',Q, ' WAS DEFINED AS ',  TRIM(QUADRIC_FORM(Q))
            WRITE(*,*) ' WITH AN HALF-ANGLE OF ', HALF_ANGLE(Q), 'DEG.'
         ENDIF


      ELSEIF((TRIM(QUADRIC_FORM(QM1))=='Y_CYL_INT').AND.(TRIM(QUADRIC_FORM(QP1))=='Y_CYL_INT')) THEN

         QUADRIC_FORM(Q) = 'Y_CONE'

         aligned = (t_x(QM1)==t_x(QP1)).AND.(t_z(QM1)==t_z(QP1)) 
         IF(.NOT.aligned) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT ALIGNED'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         R1 = RADIUS(QM1)
         R2 = RADIUS(QP1)
         IF(R1==R2) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' HAVE THE SAME RADIUS'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         y1 = piece_ymax(QM1)
         y2 = piece_ymin(QP1)
         IF(y2<=y1) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT PIECED PROPERLY'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         tan_half_angle = (R2-R1)/(y2-y1)

         HALF_ANGLE(Q) = DATAN(tan_half_angle)/PI*180.0D0
         lambda_x(Q) = ONE/(tan_half_angle)**2
         lambda_y(Q) = -ONE
         lambda_z(Q) = ONE/(tan_half_angle)**2
         dquadric(Q) = ZERO

         piece_ymin(Q) = y1
         piece_ymax(Q) = y2

         t_x(Q) = t_x(QM1)
         t_y(Q) = y1 - R1/tan_half_angle
         t_z(Q) = t_z(QM1)

         IF(MyPE == PE_IO) THEN
            WRITE(*,*) ' QUADRIC:',Q, ' WAS DEFINED AS ',  TRIM(QUADRIC_FORM(Q))
            WRITE(*,*) ' WITH AN HALF-ANGLE OF ', HALF_ANGLE(Q), 'DEG.'
         ENDIF
       

      ELSEIF((TRIM(QUADRIC_FORM(QM1))=='Z_CYL_INT').AND.(TRIM(QUADRIC_FORM(QP1))=='Z_CYL_INT')) THEN

         QUADRIC_FORM(Q) = 'Z_CONE'

         aligned = (t_x(QM1)==t_x(QP1)).AND.(t_y(QM1)==t_y(QP1)) 
         IF(.NOT.aligned) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT ALIGNED'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         R1 = RADIUS(QM1)
         R2 = RADIUS(QP1)
         IF(R1==R2) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' HAVE THE SAME RADIUS'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         z1 = piece_zmax(QM1)
         z2 = piece_zmin(QP1)
         IF(z2<=z1) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)' ERROR: CYLINDERS ',QM1, ' AND ', QP1, ' ARE NOT PIECED PROPERLY'
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
            ENDIF
            call mfix_exit(myPE)
         ENDIF

         tan_half_angle = (R2-R1)/(z2-z1)

         HALF_ANGLE(Q) = DATAN(tan_half_angle)/PI*180.0D0
         lambda_x(Q) = ONE/(tan_half_angle)**2
         lambda_y(Q) = ONE/(tan_half_angle)**2
         lambda_z(Q) = -ONE
         dquadric(Q) = ZERO

         piece_zmin(Q) = z1
         piece_zmax(Q) = z2

         t_x(Q) = t_x(QM1)
         t_y(Q) = t_y(QM1)
         t_z(Q) = z1 - R1/tan_half_angle

         IF(MyPE == PE_IO) THEN
            WRITE(*,*) ' QUADRIC:',Q, ' WAS DEFINED AS ',  TRIM(QUADRIC_FORM(Q))
            WRITE(*,*) ' WITH AN HALF-ANGLE OF ', HALF_ANGLE(Q), 'DEG.'
         ENDIF

      ELSE
         IF(MyPE == PE_IO) THEN
            WRITE(*,*) ' ERROR: C2C MUST BE DEFINED BETWEEN 2 CYLINDERS'
            WRITE(*,*) ' QUADRIC:',QM1, ' IS ',  TRIM(QUADRIC_FORM(QM1))
            WRITE(*,*) ' QUADRIC:',QP1, ' IS ',  TRIM(QUADRIC_FORM(QP1))
         ENDIF
         call mfix_exit(myPE)

      ENDIF

      RETURN
      END


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_FLOW_TO_VEL                                         C
!  Purpose: Convert flow to velocity bc's                              C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CG_FLOW_TO_VEL
    


      USE physprop
      USE scales
      USE funits 

      USE param
      USE param1
      USE parallel
      USE constant
      USE bc
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      USE fldvar
      USE vtk

     
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
!     loop/variable indices 
      INTEGER :: IJK, M, BCV
      CHARACTER(LEN=9) :: BCT
!     Volumetric flow rate computed from mass flow rate 
      DOUBLE PRECISION :: VOLFLOW 
!     Solids phase volume fraction 
      DOUBLE PRECISION :: EPS 
!     Average molecular weight 
      DOUBLE PRECISION :: MW 
!
      INTEGER :: iproc,IERR
! 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      DOUBLE PRECISION , EXTERNAL :: EOSG, CALC_MW 
      LOGICAL , EXTERNAL :: COMPARE 
!-----------------------------------------------
!

      include "function.inc"

     
! Compute the Area of each boundary condition for cut cells

       DO BCV = 1, DIMENSION_BC 
          IF (BC_TYPE(BCV)(1:2)=='CG') BC_AREA(BCV) = ZERO
       ENDDO

      DO IJK = IJKSTART3, IJKEND3
         IF(CUT_CELL_AT(IJK)) THEN              
            BCV = BC_ID(IJK)
            IF(BCV > 0 ) THEN
               BCT = BC_TYPE(BCV)
               BC_AREA(BCV) = BC_AREA(BCV) + Area_CUT(IJK)
            ENDIF
         ENDIF
      END DO

!      IF (myPE == PE_IO) THEN
!          DO BCV = 1, DIMENSION_BC 
!             IF (BC_DEFINED(BCV).OR.BC_TYPE(BCV)(1:2)=='CG') THEN 
!                WRITE(*,100) 'BOUNDARY CONDITION ID   :',BCV
!                WRITE(*,110) 'BOUNDARY CONDITION TYPE :',BC_TYPE(BCV)
!                WRITE(*,120) 'BOUNDARY CONDITION AREA :',BC_AREA(BCV)
!             ENDIF
!          ENDDO
!       ENDIF


      DO BCV = 1, DIMENSION_BC 

         IF (BC_TYPE(BCV)=='CG_MI') THEN

            IF(BC_VELMAG_g(BCV)==UNDEFINED) THEN
!
!           If gas mass flow is defined convert it to volumetric flow
!
               IF (BC_MASSFLOW_G(BCV) /= UNDEFINED) THEN 
                  IF (RO_G0 /= UNDEFINED) THEN 
                     VOLFLOW = BC_MASSFLOW_G(BCV)/RO_G0 
                  ELSE 
                     IF (BC_P_G(BCV)/=UNDEFINED .AND. BC_T_G(BCV)/=UNDEFINED) &
                        THEN 
                        IF (MW_AVG == UNDEFINED) THEN 
                           MW = CALC_MW(BC_X_G,DIMENSION_BC,BCV,NMAX(0),MW_G) 
                        ELSE 
                           MW = MW_AVG 
                        ENDIF 
                        VOLFLOW = BC_MASSFLOW_G(BCV)/EOSG(MW,(BC_P_G(BCV)-P_REF), &
			                         BC_T_G(BCV))
                     ELSE 
                        IF (BC_TYPE(BCV) == 'CG_MO') THEN 
                           IF (BC_MASSFLOW_G(BCV) == ZERO) THEN 
                              VOLFLOW = ZERO 
                           ENDIF 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1020) BCV 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ENDIF 
                  ENDIF 
!
!             If volumetric flow is also specified compare both
!
                  IF (BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN 
                     IF (.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_G(BCV))) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1000) BCV, VOLFLOW, BC_VOLFLOW_G(BCV) 
                        call mfix_exit(myPE)  
                     ENDIF 
                  ELSE 
                     BC_VOLFLOW_G(BCV) = VOLFLOW 
                  ENDIF 
               ENDIF 
!
!           If gas volumetric flow is defined convert it to velocity
!
               IF (BC_VOLFLOW_G(BCV) /= UNDEFINED) THEN 
                  IF (BC_EP_G(BCV) /= UNDEFINED) THEN 
                     BC_VELMAG_g(BCV) = BC_VOLFLOW_G(BCV)/(BC_AREA(BCV)*BC_EP_G(BCV)) 
                  ELSE 
                     RETURN                      !Error caught in Check_data_07 
                  ENDIF 
               ENDIF 

            ENDIF


!
!  Do flow conversions for solids phases
!
            DO M = 1, MMAX 

               IF(BC_VELMAG_s(BCV,M)==UNDEFINED) THEN
!
!             If solids mass flow is defined convert it to volumetric flow
!
                  IF (BC_MASSFLOW_S(BCV,M) /= UNDEFINED) THEN 
                     IF (RO_S(M) /= UNDEFINED) THEN 
                        VOLFLOW = BC_MASSFLOW_S(BCV,M)/RO_S(M) 
                     ELSE 
                        RETURN                   !  This error will be caught in a previous routine 
                     ENDIF 
!
!               If volumetric flow is also specified compare both
!
                     IF (BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN 
                        IF (.NOT.COMPARE(VOLFLOW,BC_VOLFLOW_S(BCV,M))) THEN 
                           IF(DMP_LOG)WRITE(UNIT_LOG,1200)BCV,VOLFLOW,M,BC_VOLFLOW_S(BCV,M) 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ELSE 
                        BC_VOLFLOW_S(BCV,M) = VOLFLOW 
                     ENDIF 
                  ENDIF 

                  IF (BC_ROP_S(BCV,M)==UNDEFINED .AND. MMAX==1) BC_ROP_S(BCV,M)&
                        = (ONE - BC_EP_G(BCV))*RO_S(M) 
                  IF (BC_VOLFLOW_S(BCV,M) /= UNDEFINED) THEN 
                     IF (BC_ROP_S(BCV,M) /= UNDEFINED) THEN 
                        EPS = BC_ROP_S(BCV,M)/RO_S(M) 
                        IF (EPS /= ZERO) THEN 
                           BC_VELMAG_s(BCV,M) = BC_VOLFLOW_S(BCV,M)/(BC_AREA(BCV)*EPS) 
                        ELSE 
                           IF (BC_VOLFLOW_S(BCV,M) == ZERO) THEN 
                              BC_VELMAG_s(BCV,M) = ZERO 
                           ELSE 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1250) BCV, M 
                              call mfix_exit(myPE)  
                           ENDIF 
                        ENDIF 
                     ELSE 
                        IF (BC_VOLFLOW_S(BCV,M) == ZERO) THEN 
                           BC_VELMAG_s(BCV,M) = ZERO 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1260) BCV, M 
                           call mfix_exit(myPE)  
                        ENDIF 
                     ENDIF 
                  ENDIF 

               ENDIF
            END DO 
         ENDIF 
      END DO 



100         FORMAT(1X,A,I8)
110         FORMAT(1X,A,A)
120         FORMAT(1X,A,F14.8,/)
130         FORMAT(1X,A,I8,F14.8,/)


 1000 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed volumetric flow is not equal to specified value',/,&
         ' Value computed from mass flow  = ',G14.7,/,&
         ' Specified value (BC_VOLFLOW_g) = ',G14.7,/1X,70('*')/) 


 1020 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,&
         '  BC_P_g, BC_T_g, and BC_X_g',/' should be specified',/1X,70('*')/) 


 1200 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Computed volumetric flow is not equal to specified value',/,&
         ' Value computed from mass flow  = ',G14.7,/,&
         ' Specified value (BC_VOLFLOW_s',I1,') = ',G14.7,/1X,70('*')/) 

 1250 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' Non-zero vol. or mass flow specified with BC_ROP_s',&
         I1,' = 0.',/1X,70('*')/) 
 1260 FORMAT(/1X,70('*')//' From: FLOW_TO_VEL',/' Message: BC No:',I2,/,&
         ' BC_ROP_s',I1,' not specified',/1X,70('*')/) 
      RETURN

      
      END SUBROUTINE CG_FLOW_TO_VEL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_DXYZ_FROM_CONTROL_POINTS                           C
!  Purpose: Define DX, DY, and DZ using control points                 C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 02-Dec-10  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_DXYZ_FROM_CONTROL_POINTS
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE constant 
      USE run
      USE physprop
      USE indices
      USE scalars
      USE funits
      USE leqsol
      USE compar             
      USE mpi_utility        
      USE bc
      USE DISCRETELEMENT

      USE cutcell
      USE quadric
      USE vtk
      USE polygon
      USE dashboard
      USE stl


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

      INTEGER :: N,NX,NY,NZ
      INTEGER :: I,I1,I2,J,J1,J2,K,K1,K2
      DOUBLE PRECISION :: L,CELL_RATIO

      LOGICAL :: SOLUTION_FOUND

      LOGICAL,DIMENSION(MAX_CP) :: INDEPENDENT_SEGMENT

      DOUBLE PRECISION :: AA

      DOUBLE PRECISION, EXTERNAL :: F


!-----------------------------------------------
!

!======================================================================
! X-DIRECTION
!======================================================================

! Step 1.  Input verification
!      1.1 Shift control points arrays such that the user only needs to enter 
!          CPX(1) and above, and CPX(0) is automatically set to zero.

      DO N = MAX_CP,1,-1
         CPX(N) = CPX(N-1)
      ENDDO

      CPX(0) = ZERO

!      1.2. Last control point must match domain length.

      NX = 0
      DO N = 1,MAX_CP
         IF(CPX(N)>ZERO) NX = NX + 1
      ENDDO

      IF(NX>0) THEN
         IF(MyPE==0)  WRITE(*,*)' INFO: DEFINING GRID SPACING IN X-DIRECTION... '
         IF(MyPE==0)  WRITE(*,*)' INFO: NUMBER OF CONTROL POINTS IN X-DIRECTION = ',NX
         IF(CPX(NX)/=XLENGTH) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: LAST CONTROL POINT MUST BE EQUAL TO XLENGTH.'
            IF(MyPE==0)  WRITE(*,*)' XLENGTH = ',XLENGTH
            IF(MyPE==0)  WRITE(*,*)' LAST CONTROL POINT = ',CPX(NX)
            call mfix_exit(myPE)
         ENDIF
      ENDIF
 
!      1.3. Check for acceptable values, and identify independent segments. If 
!           the first or last cell dimension is given, it is converted into an
!           expansion ratio.

      INDEPENDENT_SEGMENT = .TRUE.

      DO N = 1,NX   ! For each segment

         IF(CPX(N) <= CPX(N-1)) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: CONTROL POINTS ALONG X MUST BE SORTED IN ASCENDING ORDER.'
            IF(MyPE==0)  WRITE(*,*)' CPX = ',CPX(0:NX)
            call mfix_exit(myPE)
         ENDIF

         IF(NCX(N) <= 1) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: NUMBER OF CELLS MUST BE LARGER THAN 1 IN X-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' NCX = ',NCX(N)
            call mfix_exit(myPE)
         ENDIF

         IF(ERX(N) <= ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: EXPANSION RATIO MUST BE POSITIVE IN X-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' ERX = ',ERX(N)
            call mfix_exit(myPE)
         ENDIF

      ENDDO

      DO N = 1,NX   ! For each segment

         IF(FIRST_DX(N)/=ZERO.AND.LAST_DX(N)/=ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST AND LAST DX ARE DEFINED, WHICH IS NOT ALLOWED IN X-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' FIRST DX = ',FIRST_DX(N)
            IF(MyPE==0)  WRITE(*,*)' LAST  DX = ',LAST_DX(N)
            call mfix_exit(myPE)
         ELSEIF(FIRST_DX(N)>ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' INFO: FIRST DX DEFINED IN X-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' FIRST DX = ',FIRST_DX(N)
            L = CPX(N) - CPX(N-1)  ! Size of the current segment
            IF(L<=FIRST_DX(N)+TOL_F) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST DX IS NOT SMALLER THAN SEGMENT LENGTH IN X-SEGMENT :',N
               IF(MyPE==0)  WRITE(*,*)' FIRST DX = ',FIRST_DX(N)
               IF(MyPE==0)  WRITE(*,*)' SEGMENT LENGTH = ',L
               call mfix_exit(myPE)
            ENDIF
            CALL FIND_CELL_RATIO('FIRST',FIRST_DX(N),L,NCX(N),CELL_RATIO)
            ERX(N) = CELL_RATIO**(NCX(N)-1)
            IF(MyPE==0)  WRITE(*,*)' CORRESPONDING EXPANSION RATIO = ',ERX(N)
         ELSEIF(LAST_DX(N)>ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' INFO: LAST DX DEFINED IN X-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' LAST DX = ',LAST_DX(N)
            L = CPX(N) - CPX(N-1)  ! Size of the current segment
            IF(L<=LAST_DX(N)+TOL_F) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: LAST DX IS NOT SMALLER THAN SEGMENT LENGTH IN X-SEGMENT :',N
               IF(MyPE==0)  WRITE(*,*)' LAST DX = ',LAST_DX(N)
               IF(MyPE==0)  WRITE(*,*)' SEGMENT LENGTH = ',L
               call mfix_exit(myPE)
            ENDIF
            CALL FIND_CELL_RATIO('LAST ',LAST_DX(N),L,NCX(N),CELL_RATIO)
            ERX(N) = CELL_RATIO**(NCX(N)-1)
            IF(MyPE==0)  WRITE(*,*)' CORRESPONDING EXPANSION RATIO = ',ERX(N)
         ELSEIF(FIRST_DX(N)<ZERO) THEN
            IF(N==1) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST DX CANNOT MATCH PREVIOUS DX FOR FIRST SEGMENT.'
               call mfix_exit(myPE)
            ELSE
               IF(MyPE==0)  WRITE(*,*)' INFO: FIRST DX WILL ATTEMPT TO MATCH PREVIOUS DX FOR SEGMENT :',N
               INDEPENDENT_SEGMENT(N) = .FALSE.
            ENDIF
         ELSEIF(LAST_DX(N)<ZERO) THEN
            IF(N==NX) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: LAST DX CANNOT MATCH NEXT DX FOR LAST SEGMENT.'
               call mfix_exit(myPE)
            ELSE
               IF(MyPE==0)  WRITE(*,*)' INFO: LAST DX WILL ATTEMPT TO MATCH NEXT DX FOR SEGMENT :',N
               INDEPENDENT_SEGMENT(N) = .FALSE.
            ENDIF
         ENDIF

      ENDDO

! Step 3.  Computation of cell sizes.

!      3.1 First pass: Set-up all independent segments


      I1 = 0  ! First index of segment
      I2 = 0  ! Last index of segment

      DO N = 1,NX   ! For each segment

         I2 = I1 + NCX(N) - 1 

         IF(INDEPENDENT_SEGMENT(N)) THEN

            L = CPX(N) - CPX(N-1)  ! Size of the current segment

            IF(ERX(N)/=ONE) THEN
               CELL_RATIO = ERX(N)**(ONE/DFLOAT(NCX(N)-1))                     ! Ratio between two consecutive cells
               DX(I1) = L * (ONE - CELL_RATIO) / (ONE - CELL_RATIO**NCX(N))     ! First cell size

               DO I = I1+1,I2                                                   ! All other cell sizes, geometric series
                 DX(I) = DX(I-1) * CELL_RATIO
               ENDDO

            ELSE
               DX(I1:I2) = L / NCX(N)                                           ! Uniform size if expansion ratio is unity.
            ENDIF

         ENDIF

         I1 = I2 + 1                                                            ! Prepare First index for next segment

      ENDDO

!      3.2 Second pass: Set-up all dependent segments


      I1 = 0  ! First index of segment
      I2 = 0  ! Last index of segment

      DO N = 1,NX   ! For each segment

         I2 = I1 + NCX(N) - 1 

         IF(.NOT.INDEPENDENT_SEGMENT(N)) THEN

            L = CPX(N) - CPX(N-1)  ! Size of the current segment

            IF(FIRST_DX(N)<ZERO) THEN
               DX(I1) = DX(I1-1)                                                ! First cell size
               CALL FIND_CELL_RATIO('FIRST',DX(I1),L,NCX(N),CELL_RATIO)
               DO I = I1+1,I2                                                   ! All other cell sizes, geometric series
                 DX(I) = DX(I-1) * CELL_RATIO
               ENDDO
            ELSEIF(LAST_DX(N)<ZERO) THEN
               DX(I2) = DX(I2+1)                                                ! Last cell size
               CALL FIND_CELL_RATIO('LAST ',DX(I2),L,NCX(N),CELL_RATIO)
               DO I = I2-1,I1,-1                                                ! All other cell sizes, geometric series
                 DX(I) = DX(I+1) / CELL_RATIO
               ENDDO
            ENDIF

         ENDIF

         I1 = I2 + 1                                                  ! Prepare First index for next segment

      ENDDO


! Step 4. Verify that the sum of cells among all segment matches the total number of cells

      IF(I1>0.AND.I1/=IMAX) THEN
         IF(MyPE==0)  WRITE(*,*)' ERROR: IMAX MUST BE EQUAL TO THE SUM OF NCX.'
         IF(MyPE==0)  WRITE(*,*)' IMAX = ', IMAX
         IF(MyPE==0)  WRITE(*,*)' SUM OF NCX = ', I1
         call mfix_exit(myPE)
      ENDIF


!======================================================================
! Y-DIRECTION
!======================================================================

! Step 1.  Input verification
!      1.1 Shift control points arrays such that the user only needs to enter 
!          CPY(1) and above, and CPY(0) is automatically set to zero.

      DO N = MAX_CP,1,-1
         CPY(N) = CPY(N-1)
      ENDDO

      CPY(0) = ZERO

!      1.2. Last control point must match domain length.

      NY = 0
      DO N = 1,MAX_CP
         IF(CPY(N)>ZERO) NY = NY + 1
      ENDDO

      IF(NY>0) THEN
         IF(MyPE==0)  WRITE(*,*)' INFO: DEFINING GRID SPACING IN Y-DIRECTION... '
         IF(MyPE==0)  WRITE(*,*)' INFO: NUMBER OF CONTROL POINTS IN Y-DIRECTION = ',NY
         IF(CPY(NY)/=YLENGTH) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: LAST CONTROL POINT MUST BE EQUAL TO YLENGTH.'
            IF(MyPE==0)  WRITE(*,*)' YLENGTH = ',YLENGTH
            IF(MyPE==0)  WRITE(*,*)' LAST CONTROL POINT = ',CPY(NY)
            call mfix_exit(myPE)
         ENDIF
      ENDIF
 
!      1.3. Check for acceptable values, and identify independent segments. If 
!           the first or last cell dimension is given, it is converted into an
!           expansion ratio.

      INDEPENDENT_SEGMENT = .TRUE.

      DO N = 1,NY   ! For each segment

         IF(CPY(N) <= CPY(N-1)) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: CONTROL POINTS ALONG Y MUST BE SORTED IN ASCENDING ORDER.'
            IF(MyPE==0)  WRITE(*,*)' CPY = ',CPY(0:NY)
            call mfix_exit(myPE)
         ENDIF

         IF(NCY(N) <= 1) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: NUMBER OF CELLS MUST BE LARGER THAN 1 IN Y-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' NCY = ',NCY(N)
            call mfix_exit(myPE)
         ENDIF

         IF(ERY(N) <= ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: EXPANSION RATIO MUST BE POSITIVE IN Y-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' ERY = ',ERY(N)
            call mfix_exit(myPE)
         ENDIF

      ENDDO

      DO N = 1,NY   ! For each segment

         IF(FIRST_DY(N)/=ZERO.AND.LAST_DY(N)/=ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST AND LAST DY ARE DEFINED, WHICH IS NOT ALLOWED IN Y-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' FIRST DY = ',FIRST_DY(N)
            IF(MyPE==0)  WRITE(*,*)' LAST  DY = ',LAST_DY(N)
            call mfix_exit(myPE)
         ELSEIF(FIRST_DY(N)>ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' INFO: FIRST DY DEFINED IN Y-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' FIRST DY = ',FIRST_DY(N)
            L = CPY(N) - CPY(N-1)  ! Size of the current segment
            IF(L<=FIRST_DY(N)+TOL_F) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST DY IS NOT SMALLER THAN SEGMENT LENGTH IN Y-SEGMENT :',N
               IF(MyPE==0)  WRITE(*,*)' FIRST DY = ',FIRST_DY(N)
               IF(MyPE==0)  WRITE(*,*)' SEGMENT LENGTH = ',L
               call mfix_exit(myPE)
            ENDIF
            CALL FIND_CELL_RATIO('FIRST',FIRST_DY(N),L,NCY(N),CELL_RATIO)
            ERY(N) = CELL_RATIO**(NCY(N)-1)
            IF(MyPE==0)  WRITE(*,*)' CORRESPONDING EXPANSION RATIO = ',ERY(N)
         ELSEIF(LAST_DY(N)>ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' INFO: LAST DY DEFINED IN Y-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' LAST DY = ',LAST_DY(N)
            L = CPY(N) - CPY(N-1)  ! Size of the current segment
            IF(L<=LAST_DY(N)+TOL_F) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: LAST DY IS NOT SMALLER THAN SEGMENT LENGTH IN Y-SEGMENT :',N
               IF(MyPE==0)  WRITE(*,*)' LAST DY = ',LAST_DY(N)
               IF(MyPE==0)  WRITE(*,*)' SEGMENT LENGTH = ',L
               call mfix_exit(myPE)
            ENDIF
            CALL FIND_CELL_RATIO('LAST ',LAST_DY(N),L,NCY(N),CELL_RATIO)
            ERY(N) = CELL_RATIO**(NCY(N)-1)
            IF(MyPE==0)  WRITE(*,*)' CORRESPONDING EXPANSION RATIO = ',ERY(N)
         ELSEIF(FIRST_DY(N)<ZERO) THEN
            IF(N==1) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST DY CANNOT MATCH PREVIOUS DY FOR FIRST SEGMENT.'
               call mfix_exit(myPE)
            ELSE
               IF(MyPE==0)  WRITE(*,*)' INFO: FIRST DY WILL ATTEMPT TO MATCH PREVIOUS DY FOR SEGMENT :',N
               INDEPENDENT_SEGMENT(N) = .FALSE.
            ENDIF
         ELSEIF(LAST_DY(N)<ZERO) THEN
            IF(N==NY) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: LAST DY CANNOT MATCH NEXT DY FOR LAST SEGMENT.'
               call mfix_exit(myPE)
            ELSE
               IF(MyPE==0)  WRITE(*,*)' INFO: LAST DY WILL ATTEMPT TO MATCH NEXT DY FOR SEGMENT :',N
               INDEPENDENT_SEGMENT(N) = .FALSE.
            ENDIF
         ENDIF

      ENDDO

! Step 3.  Computation of cell sizes.

!      3.1 First pass: Set-up all independent segments


      J1 = 0  ! First index of segment
      J2 = 0  ! Last index of segment

      DO N = 1,NY   ! For each segment

         J2 = J1 + NCY(N) - 1 

         IF(INDEPENDENT_SEGMENT(N)) THEN

            L = CPY(N) - CPY(N-1)  ! Size of the current segment

            IF(ERY(N)/=ONE) THEN
               CELL_RATIO = ERY(N)**(ONE/DFLOAT(NCY(N)-1))                     ! Ratio between two consecutive cells
               DY(J1) = L * (ONE - CELL_RATIO) / (ONE - CELL_RATIO**NCY(N))     ! First cell size

               DO J = J1+1,J2                                                   ! All other cell sizes, geometric series
                 DY(J) = DY(J-1) * CELL_RATIO
               ENDDO

            ELSE
               DY(J1:J2) = L / NCY(N)                                           ! Uniform size if expansion ratio is unity.
            ENDIF

         ENDIF

         J1 = J2 + 1                                                            ! Prepare First index for next segment

      ENDDO

!      3.2 Second pass: Set-up all dependent segments


      J1 = 0  ! First index of segment
      J2 = 0  ! Last index of segment

      DO N = 1,NY   ! For each segment

         J2 = J1 + NCY(N) - 1 

         IF(.NOT.INDEPENDENT_SEGMENT(N)) THEN

            L = CPY(N) - CPY(N-1)  ! Size of the current segment

            IF(FIRST_DY(N)<ZERO) THEN
               DY(J1) = DY(J1-1)                                                ! First cell size
               CALL FIND_CELL_RATIO('FIRST',DY(J1),L,NCY(N),CELL_RATIO)
               DO J = J1+1,J2                                                   ! All other cell sizes, geometric series
                 DY(J) = DY(J-1) * CELL_RATIO
               ENDDO
            ELSEIF(LAST_DY(N)<ZERO) THEN
               DY(J2) = DY(J2+1)                                                ! Last cell size
               CALL FIND_CELL_RATIO('LAST ',DY(J2),L,NCY(N),CELL_RATIO)
               DO J = J2-1,J1,-1                                                ! All other cell sizes, geometric series
                 DY(J) = DY(J+1) / CELL_RATIO
               ENDDO
            ENDIF

         ENDIF

         J1 = J2 + 1                                                  ! Prepare First index for next segment

      ENDDO


! Step 4. Verify that the sum of cells among all segment matches the total number of cells

      IF(J1>0.AND.J1/=JMAX) THEN
         IF(MyPE==0)  WRITE(*,*)' ERROR: JMAX MUST BE EQUAL TO THE SUM OF NCY.'
         IF(MyPE==0)  WRITE(*,*)' JMAX = ', JMAX
         IF(MyPE==0)  WRITE(*,*)' SUM OF NCY = ', J1
         call mfix_exit(myPE)
      ENDIF


!======================================================================
! Z-DIRECTION
!======================================================================

      IF(NO_K) RETURN

! Step 1.  Input verification
!      1.1 Shift control points arrays such that the user only needs to enter 
!          CPZ(1) and above, and CPZ(0) is automatically set to zero.

      DO N = MAX_CP,1,-1
         CPZ(N) = CPZ(N-1)
      ENDDO

      CPZ(0) = ZERO

!      1.2. Last control point must match domain length.

      NZ = 0
      DO N = 1,MAX_CP
         IF(CPZ(N)>ZERO) NZ = NZ + 1
      ENDDO

      IF(NZ>0) THEN
         IF(MyPE==0)  WRITE(*,*)' INFO: DEFINING GRID SPACING IN Z-DIRECTION... '
         IF(MyPE==0)  WRITE(*,*)' INFO: NUMBER OF CONTROL POINTS IN Z-DIRECTION = ',NZ
         IF(CPZ(NZ)/=ZLENGTH) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: LAST CONTROL POINT MUST BE EQUAL TO ZLENGTH.'
            IF(MyPE==0)  WRITE(*,*)' ZLENGTH = ',ZLENGTH
            IF(MyPE==0)  WRITE(*,*)' LAST CONTROL POINT = ',CPZ(NZ)
            call mfix_exit(myPE)
         ENDIF
      ENDIF
 
!      1.3. Check for acceptable values, and identify independent segments. If 
!           the first or last cell dimension is given, it is converted into an
!           expansion ratio.

      INDEPENDENT_SEGMENT = .TRUE.

      DO N = 1,NZ   ! For each segment

         IF(CPZ(N) <= CPZ(N-1)) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: CONTROL POINTS ALONG Z MUST BE SORTED IN ASCENDING ORDER.'
            IF(MyPE==0)  WRITE(*,*)' CPZ = ',CPZ(0:NZ)
            call mfix_exit(myPE)
         ENDIF

         IF(NCZ(N) <= 1) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: NUMBER OF CELLS MUST BE LARGER THAN 1 IN Z-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' NCZ = ',NCZ(N)
            call mfix_exit(myPE)
         ENDIF

         IF(ERZ(N) <= ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: EXPANSION RATIO MUST BE POSITIVE IN Z-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' ERZ = ',ERZ(N)
            call mfix_exit(myPE)
         ENDIF

      ENDDO

      DO N = 1,NZ   ! For each segment

         IF(FIRST_DZ(N)/=ZERO.AND.LAST_DZ(N)/=ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST AND LAST DZ ARE DEFINED, WHICH IS NOT ALLOWED IN Z-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' FIRST DZ = ',FIRST_DZ(N)
            IF(MyPE==0)  WRITE(*,*)' LAST  DZ = ',LAST_DZ(N)
            call mfix_exit(myPE)
         ELSEIF(FIRST_DZ(N)>ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' INFO: FIRST DZ DEFINED IN Z-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' FIRST DZ = ',FIRST_DZ(N)
            L = CPZ(N) - CPZ(N-1)  ! Size of the current segment
            IF(L<=FIRST_DZ(N)+TOL_F) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST DZ IS NOT SMALLER THAN SEGMENT LENGTH IN Z-SEGMENT :',N
               IF(MyPE==0)  WRITE(*,*)' FIRST DZ = ',FIRST_DZ(N)
               IF(MyPE==0)  WRITE(*,*)' SEGMENT LENGTH = ',L
               call mfix_exit(myPE)
            ENDIF
            CALL FIND_CELL_RATIO('FIRST',FIRST_DZ(N),L,NCZ(N),CELL_RATIO)
            ERZ(N) = CELL_RATIO**(NCZ(N)-1)
            IF(MyPE==0)  WRITE(*,*)' CORRESPONDING EXPANSION RATIO = ',ERZ(N)
         ELSEIF(LAST_DZ(N)>ZERO) THEN
            IF(MyPE==0)  WRITE(*,*)' INFO: LAST DZ DEFINED IN Z-SEGMENT :',N
            IF(MyPE==0)  WRITE(*,*)' LAST DZ = ',LAST_DZ(N)
            L = CPZ(N) - CPZ(N-1)  ! Size of the current segment
            IF(L<=LAST_DZ(N)+TOL_F) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: LAST DZ IS NOT SMALLER THAN SEGMENT LENGTH IN Z-SEGMENT :',N
               IF(MyPE==0)  WRITE(*,*)' LAST DZ = ',LAST_DZ(N)
               IF(MyPE==0)  WRITE(*,*)' SEGMENT LENGTH = ',L
               call mfix_exit(myPE)
            ENDIF
            CALL FIND_CELL_RATIO('LAST ',LAST_DZ(N),L,NCZ(N),CELL_RATIO)
            ERZ(N) = CELL_RATIO**(NCZ(N)-1)
            IF(MyPE==0)  WRITE(*,*)' CORRESPONDING EXPANSION RATIO = ',ERZ(N)
         ELSEIF(FIRST_DZ(N)<ZERO) THEN
            IF(N==1) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: FIRST DZ CANNOT MATCH PREVIOUS DZ FOR FIRST SEGMENT.'
               call mfix_exit(myPE)
            ELSE
               IF(MyPE==0)  WRITE(*,*)' INFO: FIRST DZ WILL ATTEMPT TO MATCH PREVIOUS DZ FOR SEGMENT :',N
               INDEPENDENT_SEGMENT(N) = .FALSE.
            ENDIF
         ELSEIF(LAST_DZ(N)<ZERO) THEN
            IF(N==NZ) THEN
               IF(MyPE==0)  WRITE(*,*)' ERROR: LAST DZ CANNOT MATCH NEXT DZ FOR LAST SEGMENT.'
               call mfix_exit(myPE)
            ELSE
               IF(MyPE==0)  WRITE(*,*)' INFO: LAST DZ WILL ATTEMPT TO MATCH NEXT DZ FOR SEGMENT :',N
               INDEPENDENT_SEGMENT(N) = .FALSE.
            ENDIF
         ENDIF

      ENDDO

! Step 3.  Computation of cell sizes.

!      3.1 First pass: Set-up all independent segments


      K1 = 0  ! First index of segment
      K2 = 0  ! Last index of segment

      DO N = 1,NZ   ! For each segment

         K2 = K1 + NCZ(N) - 1 

         IF(INDEPENDENT_SEGMENT(N)) THEN

            L = CPZ(N) - CPZ(N-1)  ! Size of the current segment

            IF(ERZ(N)/=ONE) THEN
               CELL_RATIO = ERZ(N)**(ONE/DFLOAT(NCZ(N)-1))                     ! Ratio between two consecutive cells
               DZ(K1) = L * (ONE - CELL_RATIO) / (ONE - CELL_RATIO**NCZ(N))     ! First cell size

               DO K = K1+1,K2                                                   ! All other cell sizes, geometric series
                 DZ(K) = DZ(K-1) * CELL_RATIO
               ENDDO

            ELSE
               DZ(K1:K2) = L / NCZ(N)                                           ! Uniform size if expansion ratio is unity.
            ENDIF

         ENDIF

         K1 = K2 + 1                                                            ! Prepare First index for next segment

      ENDDO

!      3.2 Second pass: Set-up all dependent segments


      K1 = 0  ! First index of segment
      K2 = 0  ! Last index of segment

      DO N = 1,NZ   ! For each segment

         K2 = K1 + NCZ(N) - 1 

         IF(.NOT.INDEPENDENT_SEGMENT(N)) THEN

            L = CPZ(N) - CPZ(N-1)  ! Size of the current segment

            IF(FIRST_DZ(N)<ZERO) THEN
               DZ(K1) = DZ(K1-1)                                                ! First cell size
               CALL FIND_CELL_RATIO('FIRST',DZ(K1),L,NCZ(N),CELL_RATIO)
               DO K = K1+1,K2                                                   ! All other cell sizes, geometric series
                 DZ(K) = DZ(K-1) * CELL_RATIO
               ENDDO
            ELSEIF(LAST_DZ(N)<ZERO) THEN
               DZ(K2) = DZ(K2+1)                                                ! Last cell size
               CALL FIND_CELL_RATIO('LAST ',DZ(K2),L,NCZ(N),CELL_RATIO)
               DO K = K2-1,K1,-1                                                ! All other cell sizes, geometric series
                 DZ(K) = DZ(K+1) / CELL_RATIO
               ENDDO
            ENDIF

         ENDIF

         K1 = K2 + 1                                                  ! Prepare First index for next segment

      ENDDO


! Step 4. Verify that the sum of cells among all segment matches the total number of cells

      IF(K1>0.AND.K1/=KMAX) THEN
         IF(MyPE==0)  WRITE(*,*)' ERROR: KMAX MUST BE EQUAL TO THE SUM OF NCZ.'
         IF(MyPE==0)  WRITE(*,*)' KMAX = ', KMAX
         IF(MyPE==0)  WRITE(*,*)' SUM OF NCZ = ', K1
         call mfix_exit(myPE)
      ENDIF



      RETURN

      END SUBROUTINE GET_DXYZ_FROM_CONTROL_POINTS





      DOUBLE PRECISION Function F(POS,ALPHAC,D_Target,L,N,TOL)
      USE constant 
      USE mpi_utility    

      IMPLICIT NONE
      DOUBLE PRECISION:: ALPHAC,D,D_Target,DU,L,TOL
      INTEGER:: N
      CHARACTER (LEN=5) :: POS

      DU = L / DFLOAT(N)    ! Cell size if uniform distribution

      IF(ALPHAC==ONE) THEN
         D = DU
      ELSE
         IF(TRIM(POS)=='FIRST') THEN
            D = L * (ONE - ALPHAC) / (ONE -ALPHAC**N)
         ELSEIF(TRIM(POS)=='LAST') THEN
            D = L * (ONE - ALPHAC) / (ONE -ALPHAC**N) * ALPHAC**(N-1)
         ELSE
            IF(MyPE==0) WRITE(*,*)' ERROR, IN FUNCTION F: POS MUST BE FIRST OR LAST.'
            call mfix_exit(myPE)
         ENDIF
      ENDIF
       

      F = D - D_Target

      RETURN

      END FUNCTION F



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: FIND_CELL_RATIO                                        C
!  Purpose: Given the interval length L, number of cells N, and the    C
!           target value of D_target, find the cell ratio alpha3       C
!           such that D(POS) matches D_Target. POS can be either       C
!           FIRST or LAST cell in the segment.                         C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE FIND_CELL_RATIO(POS,D_Target,L,N,ALPHA3)
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      
      IMPLICIT NONE
      LOGICAL :: CLIP_FLAG,CLIP_FLAG1,CLIP_FLAG2,CLIP_FLAG3,INTERSECT_FLAG,SOLUTION_FOUND

      DOUBLE PRECISION :: f1,f2,f3
      DOUBLE PRECISION :: ALPHA1,ALPHA2,ALPHA3,D_Target,L,TOL,DU
      DOUBLE PRECISION, PARAMETER :: ALPHAMAX = 1000.0D0  ! maximum  value of cell ratio
      INTEGER :: N,niter
      DOUBLE PRECISION, EXTERNAL :: F
      CHARACTER (LEN=5) :: POS


      DU = L / DFLOAT(N)                  ! Cell size if uniform distribution

      IF(DU==D_TARGET) THEN
         ALPHA3 = 1.0 
         SOLUTION_FOUND = .TRUE.
         RETURN
      ELSE

         IF(TRIM(POS)=='FIRST') THEN     ! Determine two initial guesses
            IF(D_TARGET<DU) THEN
               ALPHA1 = ONE
               ALPHA2 = ALPHAMAX
            ELSE
               ALPHA1 = ONE/ALPHAMAX
               ALPHA2 = ONE
            ENDIF
         ELSEIF(TRIM(POS)=='LAST') THEN
            IF(D_TARGET>DU) THEN
               ALPHA1 = ONE
               ALPHA2 = ALPHAMAX
            ELSE
               ALPHA1 = ONE/ALPHAMAX
               ALPHA2 = ONE
            ENDIF
         ELSE
            IF(MyPE==0) WRITE(*,*)' ERROR, IN FUNCTION F: POS MUST BE FIRST OR LAST.'
            call mfix_exit(myPE)
         ENDIF

      ENDIF


      f1 = F(POS,ALPHA1,D_Target,L,N,TOL)
      f2 = F(POS,ALPHA2,D_Target,L,N,TOL)

!======================================================================
!  The cell ratio is solution of F(alpha) = zero. The root is found by 
!  the secant method, based on two inital guesses.
!======================================================================

      niter = 1
      SOLUTION_FOUND = .FALSE.

        if(DABS(f1)<TOL_F) then         ! First guess is solution
           SOLUTION_FOUND = .TRUE.
           ALPHA3 = ALPHA1 
        elseif(DABS(f2)<TOL_F) then    ! Second guess is solution
           SOLUTION_FOUND = .TRUE.
           ALPHA3 = ALPHA2
        elseif(f1*f2 < ZERO) then       ! Solution is between two guesses
          niter = 0
          f3 = 2.0d0*TOL_F
          do while (   (abs(f3) > TOL_F)   .AND.   (niter<ITERMAX_INT)       )
           
            ALPHA3 = ALPHA1 - f1*(ALPHA2-ALPHA1)/(f2-f1)  ! secant point

            f3 = F(POS,ALPHA3,D_Target,L,N,TOL)

            if(f1*f3<0) then            ! Reduce size of interval 
              ALPHA2 = ALPHA3
              f2 = f3
            else 
              ALPHA1 = ALPHA3
              f1 = f3
            endif   
            niter = niter + 1

          end do
          if (niter < ITERMAX_INT) then
            SOLUTION_FOUND = .TRUE.
          else
             WRITE(*,*)   'Unable to find a solution'
             WRITE(*,1000)'between ALPHA1 = ', ALPHA1
             WRITE(*,1000)'   and  ALPHA2 = ', ALPHA2
             WRITE(*,1000)'Current value of ALPHA3 = ', ALPHA3
             WRITE(*,1000)'Current value of abs(f) = ', DABS(f3)
             WRITE(*,1000)'Tolerance = ', TOL_F
             WRITE(*,*)'Maximum number of iterations = ', ITERMAX_INT
             WRITE(*,*)   'Please increase the intersection tolerance, '
             WRITE(*,*)   'or the maximum number of iterations, and try again.'
             WRITE(*,*)   'MFiX will exit now.'             
             CALL MFIX_EXIT(myPE) 
             SOLUTION_FOUND = .FALSE.
          endif
        else
          WRITE(*,*)   'Unable to find a solution'
          WRITE(*,*)   'MFiX will exit now.'             
          CALL MFIX_EXIT(myPE) 
          SOLUTION_FOUND = .FALSE.
        endif


 1000 FORMAT(A,3(2X,G12.5)) 


      RETURN

      END SUBROUTINE FIND_CELL_RATIO
