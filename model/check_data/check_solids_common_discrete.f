!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_SOLIDS_COMMON_DISCRETE                            !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE


! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Invoke gas/solids coupled simulation.
      USE discretelement, only: DES_CONTINUUM_COUPLED
! Runtime Flag: Generate initial particle configuation.
      USE discretelement, only: GENER_PART_CONFIG
! Runtime Flag: Invoke MPPIC model.
      USE mfix_pic, only: MPPIC
! Runtime Flag: Store DES_*_OLD arrays.
      USE discretelement, only: DO_OLD
! Runtime Flag: Invoke TFM/DEM hybrid model.
      USE discretelement, only: DES_CONTINUUM_HYBRID
! Runtime Flag: Solve energy equations
      USE run, only: ENERGY_EQ
! Number of DEM solids phases.
      USE discretelement, only: DES_MMAX
! DEM solid phase diameters and densities.
      USE discretelement, only: DES_D_p0, DES_RO_s
! TFM solids phase diameters and densities. (DEM default)
      USE physprop, only: D_p0, RO_s0

! User specified integration method.
      USE discretelement, only: DES_INTG_METHOD
      USE discretelement, only: INTG_ADAMS_BASHFORTH
      USE discretelement, only: INTG_EULER
! User specified neighbor search method.
      USE discretelement, only: DES_NEIGHBOR_SEARCH
! User specified data out format (VTP, TecPlot)
      USE discretelement, only: DES_OUTPUT_TYPE
! Max/Min particle radi
      USE discretelement, only: MAX_RADIUS, MIN_RADIUS
! Runtime Flag: Periodic boundaries
      USE discretelement, only: DES_PERIODIC_WALLS
      USE discretelement, only: DES_PERIODIC_WALLS_X
      USE discretelement, only: DES_PERIODIC_WALLS_Y
      USE discretelement, only: DES_PERIODIC_WALLS_Z
! Flag: Solve variable solids density.
      use run, only: SOLVE_ROs
! Calculated baseline variable solids density.
      use physprop, only: BASE_ROs

! Number of ranks.
      use run, only: SOLIDS_MODEL

! Subroutine access.
      use physprop, only: MMAX

      USE run, only: MOMENTUM_X_EQ
      USE run, only: MOMENTUM_Y_EQ
      USE run, only: MOMENTUM_Z_EQ

      use run, only: RUN_TYPE
      use discretelement, only: GENER_PART_CONFIG

      USE physprop, only: CLOSE_PACKED

      USE mpi_utility


! Global Parameters:
!---------------------------------------------------------------------//
!      use param1, only: UNDEFINED_I

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Loop index.
      INTEGER :: M, lM  ! Solids phase Index

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE")


      DES_D_p0 = UNDEFINED
      DES_RO_s = UNDEFINED

      MAX_RADIUS = -UNDEFINED
      MIN_RADIUS =  UNDEFINED

      M = 0
      DO lM=1, MMAX+DES_MMAX

! The accounts for an offset between the DEM and TFM phase indices
         IF(SOLIDS_MODEL(lM) == 'TFM') CYCLE
         M = M+1

! Copy of the input keyword values into discrete solids arrays. We may be
! able to remove the DES_ specific variables moving forward.
         DES_D_p0(M) = D_p0(lM)
         DES_RO_s(M) = merge(BASE_ROs(lM), RO_s0(lM), SOLVE_ROs(lM))
! Determine the maximum particle size in the system (MAX_RADIUS), which
! in turn is used for various tasks
         MAX_RADIUS = MAX(MAX_RADIUS, 0.5d0*DES_D_P0(M))
         MIN_RADIUS = MIN(MIN_RADIUS, 0.5d0*DES_D_P0(M))
      ENDDO


! Set close_packed to true to prevent possible issues stemming from the
! pressure correction equation.  Specifically, if closed_packed is false
! then a mixture pressure correction equation is invoked and this is not
! correctly setup for DEM.  To do so would require ensuring that
! 1) the solids phase continuum quantities used in these equations are
!    correctly set based on their DEM counterparts and
! 2) the pressure correction coefficients for such solids phases are
!    also calculated (currently these calculations are turned off
!    when using DEM)
      CLOSE_PACKED((MMAX+1):DIM_M) = .TRUE.


! Turn off the 'continuum' equations for discrete solids if the user
! specified them.  We could make use of these flags.
      MOMENTUM_X_EQ((MMAX+1):DIM_M) = .FALSE.
      MOMENTUM_Y_EQ((MMAX+1):DIM_M) = .FALSE.
      MOMENTUM_Z_EQ((MMAX+1):DIM_M) = .FALSE.

! Derive periodicity from cyclic boundary flags.
      DES_PERIODIC_WALLS_X = CYCLIC_X .OR. CYCLIC_X_PD
      DES_PERIODIC_WALLS_Y = CYCLIC_Y .OR. CYCLIC_Y_PD
      DES_PERIODIC_WALLS_Z = CYCLIC_Z .OR. CYCLIC_Z_PD

      DES_PERIODIC_WALLS = (DES_PERIODIC_WALLS_X .OR.                  &
        DES_PERIODIC_WALLS_Y .OR. DES_PERIODIC_WALLS_Z)


! Overrite for restart cases.
      IF(TRIM(RUN_TYPE) .NE. 'NEW') GENER_PART_CONFIG = .FALSE.

! Check for valid neighbor search option.
      SELECT CASE(DES_NEIGHBOR_SEARCH)
      CASE (1) ! N-Square
      CASE (2)
         WRITE(ERR_MSG,2001) 2, 'QUADTREE'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      CASE (3)
         WRITE(ERR_MSG,2001) 3, 'OCTREE'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      CASE (4) ! Grid based
      CASE DEFAULT
         WRITE(ERR_MSG,2001) DES_NEIGHBOR_SEARCH,'UNKNOWN'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 2001 FORMAT('Error 2001:Invalid DES_NEIGHBOR_SEARCH method: ',I2,1X,  &
         A,/'Please correct the mfix.dat file.')

      END SELECT


! Check the output file format
      IF(DES_OUTPUT_TYPE == UNDEFINED_C) DES_OUTPUT_TYPE = 'PARAVIEW'
      SELECT CASE(trim(DES_OUTPUT_TYPE))
      CASE ('PARAVIEW')
      CASE ('TECPLOT')
      CASE DEFAULT
         WRITE(ERR_MSG,2010) trim(DES_OUTPUT_TYPE)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 2010 FORMAT('Error 2010:Invalid DES_OUTPUT_TYPE: ',A,/'Please ',       &
         'correct the mfix.dat file.')

      END SELECT


! Check for valid integration method
      SELECT CASE(trim(DES_INTG_METHOD))
      CASE ('EULER')
         INTG_EULER = .TRUE.
         INTG_ADAMS_BASHFORTH = .FALSE.
         !DES_INTG_METHOD_ENUM = 1
      CASE ('ADAMS_BASHFORTH')
         INTG_EULER = .FALSE.
         INTG_ADAMS_BASHFORTH = .TRUE.
         !DES_INTG_METHOD_ENUM = 2
      CASE DEFAULT
         WRITE(ERR_MSG,2020) trim(DES_INTG_METHOD)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 2020 FORMAT('Error 2020:Invalid DES_INGT_METHOD: ',A,/'Please ',      &
         'correct the mfix.dat file.')

      END SELECT

      DO_OLD = INTG_ADAMS_BASHFORTH .OR. MPPIC

! Set flags for energy equations
      IF(ENERGY_EQ) CALL CHECK_SOLIDS_COMMON_DISCRETE_ENERGY

! Check thermodynamic properties of discrete solids.
      CALL CHECK_SOLIDS_COMMON_DISCRETE_THERMO

! Check geometry constrains.
      CALL CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY

! Check interpolation input.
      CALL CHECK_SOLIDS_COMMON_DISCRETE_INTERP


      CALL FINL_ERR_MSG


      RETURN

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_ENERGY                      !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!  Purpose: Check input parameters for solving discrete solids phase   !
!  energy equations.  Only DEM simulations (neither hybrid nor MPPIC)  !
!  can invoke particle-particle heat transfer. Therefore, checks for   !
!  those functions are reseved for later.                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_ENERGY


! Global Variables:
!---------------------------------------------------------------------//
      use run, only: ENERGY_EQ
      use run, only: UNITS

      use discretelement, only: DES_MMAX

      use physprop, only: MMAX
      Use compar, only: NODESI, NODESJ, NODESK

      use des_thermo, only: DES_CONV_CORR
      use des_thermo, only: DES_CONV_CORR_ENUM
      use des_thermo, only: RANZ_1952

!      use des_thermo, only: DES_COND_EQ
!      use des_thermo, only: DES_COND_EQ_PFP
!      use des_thermo, only: DES_COND_EQ_PP

      Use des_thermo, only: SB_CONST
      Use des_thermo, only: DES_Em

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: UNDEFINED


! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: M
! Number of processors used. (Heat transfer is limited to serial runs!)
      INTEGER :: CHECK_MPI


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE_ENERGY")


!      IF(.NOT.ENERGY_EQ)THEN
! Reinitialize the heat transfer logicals to false.
!         DES_COND_EQ     = .FALSE.
!         DES_COND_EQ_PFP = .FALSE.
!         DES_COND_EQ_PP  = .FALSE.
!         CALL FINL_ERR_MSG
!         RETURN
!      ENDIF


! Check the number of processors. DES reactive chemistry is currently
! limited to serial runs.
      CHECK_MPI = NODESI * NODESJ * NODESK
      IF(CHECK_MPI.NE.1) THEN
         WRITE(ERR_MSG, 2000)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 2000 FORMAT('Error 2000: Currently, simulations with discrete solids',&
        ' heat transfer',/'modules are limited to serial runs. Please',&
        ' correct the mfix.dat file.')


! Gas/Solids convection:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! Verify the selected convective heat transfer coefficient model
      SELECT CASE(TRIM(DES_CONV_CORR))
! Ranz, W.E. and Marshall, W.R., "Frication and transfer coefficients
! for single particles and packed beds,"  Chemical Engineering Science,
! Vol. 48, No. 5, pp 247-253, 1952.
      CASE ('RANZ_1952')
         DES_CONV_CORR_ENUM = RANZ_1952
! If the heat transfer coefficient correlation provided by the user does
! not match one of the models outlined above, flag the error and exit.
      CASE DEFAULT
         WRITE(ERR_MSG,1001)'DES_CONV_CORR', trim(DES_CONV_CORR)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      END SELECT


! Radiation Equation:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! Verify that a emmisivity value is specified for each solids pase
      DO M = MMAX+1, MMAX+DES_MMAX
         IF(DES_Em(M) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('DES_Em',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO


! Set the value of the Stefan-Boltzman Constant based on the untis
      IF(UNITS == 'SI')THEN
         SB_CONST = 5.6704d0*(10.0d0**(-8)) ! W/((m^2).K^4)
      ELSE
         SB_CONST = 1.355282d0*(10.0d0**(-12)) ! cal/((cm^2).sec.K^4)
      ENDIF


      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_ENERGY




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON_DISCRETE_THERMO                     !
!  Author: J.Musser                                   Date: 17-Jun-10  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_THERMO

      use run, only: ANY_SPECIES_EQ
      use compar, only: NODESI, NODESJ, NODESK
      use stiff_chem, only: STIFF_CHEMISTRY

      use error_manager

      IMPLICIT NONE



!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE_THERMO")

! Check the number of processors. DES reactive chemistry is currently
! limited to serial runs.
      IF(ANY_SPECIES_EQ) THEN
         IF((NODESI*NODESJ*NODESK) /= 1) THEN
            WRITE(ERR_MSG, 9001)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

 9001 FORMAT('Error 9001: DES reactive chemistry is limited to ',      &
         'serail  runs.',/'NODESI, NODESJ, and NODESK must equal 1. ', &
         'Please correct the dat file.')

! Stiff chemistry solver is a TFM reaction model not for DES.
      IF(STIFF_CHEMISTRY) THEN
         WRITE(ERR_MSG,9003)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 9003 FORMAT('Error 9003: The stiff chemistry solver is not ',         &
      'available in DES',/'simulations. Please correct the input file.')


      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_THERMO


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY                   !
!  Author: J.Musser                                   Date: 11-DEC-13  !
!                                                                      !
!  Purpose: Check user input data                                      !
!                                                                      !
!  Comments: Geometry checks were moved here from CHECK_DES_DATA.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE geometry, only: COORDINATES
      USE geometry, only: DO_I, DO_J, DO_K
      USE geometry, only: NO_I, NO_J, NO_K
      USE geometry, only: ZLENGTH

! Flag: Use DES E-L model
      use discretelement, only: DISCRETE_ELEMENT
      USE discretelement, only: DES_CONTINUUM_COUPLED
      USE discretelement, only: MAX_RADIUS

! Flag: Use STL for DEM walls
      USE discretelement, only: USE_STL_DES
! Flag: Use Cartesian grid cut-cell implementation
      USE cutcell, only: CARTESIAN_GRID
! Flag: Use STL represenation in CG
      USE cutcell, only: USE_STL

      use param1, only: UNDEFINED_I

      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      DOUBLE PRECISION :: MIN_DEPTH

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY")


! DEM/MPPIC is restriced to CARTESIAN coordinates.
      IF(COORDINATES == 'CYLINDRICAL') THEN
         WRITE (ERR_MSG, 1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error: 1100: DES and MPPIC models only support ',        &
         'CARTESIAN coordinates.')


! Check dimension. This is redundant with check_data_03.
      IF(NO_I .OR. NO_J) THEN
         WRITE(ERR_MSG, 1200)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1200 FORMAT('Error 1200: Illegal geometry for DEM/MPPIC. 2D ',        &
         'simulations are',/'restricted to the XY plane. Please ',     &
         'correct the mfix.dat file.')


      IF(DES_CONTINUUM_COUPLED)THEN
! Check that the depth of the simulation exceeds the largest particle
! to ensure correct calculation of volume fraction. This is important
! for coupled simulations.
         MIN_DEPTH = 2.0d0*MAX_RADIUS
         IF(ZLENGTH < MIN_DEPTH)THEN
            WRITE(ERR_MSG, 1300)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

 1300 FORMAT('Error 1300: The maximum particle diameter exceeds the ', &
         'simulation',/'depth (ZLENGTH). Please correct the mfix.dat ',&
         'file.')


      IF(CARTESIAN_GRID .AND. .NOT.USE_STL) THEN
         WRITE(ERR_MSG,1400)
         CALL FLUSH_ERR_MSG(ABORT =.TRUE.)
      ENDIF

 1400 FORMAT('Error 1400: Cartesian grid and discrete models (DEM or ',&
         'PIC) only',/'support STL wall representations. Quadrics ',   &
         'and polygons are not',/'supported.')


      IF(CARTESIAN_GRID.AND.USE_STL.AND..NOT.USE_STL_DES)              &
         USE_STL_DES = .TRUE.


      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON_DISCRETE_INTERP                     !
!  Author: J.Musser                                   Date: 25-Nov-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_INTERP

! Runtime Flag: Utilize cutcell geometry.
      use cutcell, only: CARTESIAN_GRID
! Runtime Flag: Invoke gas/solids coupled simulation.
      use discretelement, only: DES_CONTINUUM_COUPLED
! Runtime Flag: Invoke MPPIC model.
      USE mfix_pic, only: MPPIC
! User input for DES interpolation scheme.
      use particle_filter, only: DES_INTERP_SCHEME
! Enumerated interpolation scheme for faster access
      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE
      use particle_filter, only: DES_INTERP_GARG
      use particle_filter, only: DES_INTERP_DPVM
      use particle_filter, only: DES_INTERP_GAUSS
! User specified filter width
      use particle_filter, only: FILTER_WIDTH
! Flag: Diffuse DES field variables.
      use particle_filter, only: DES_DIFFUSE_MEAN_FIELDS
! Flag: Interpolate continuum fields
      use particle_filter, only: DES_INTERP_MEAN_FIELDS
! Flag: Interplate variables for drag calculation.
      use particle_filter, only: DES_INTERP_ON

      use param1, only: UNDEFINED, UNDEFINED_C

      use error_manager

      IMPLICIT NONE

      DOUBLE PRECISION :: DXYZ_MIN


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE_INTERP")


! Set the interpolation ENUM value.
      SELECT CASE(trim(adjustl(DES_INTERP_SCHEME)))
      CASE ('NONE')
         DES_INTERP_SCHEME_ENUM = DES_INTERP_NONE
! Cannot use interpolation when no scheme is selected.
         IF(DES_INTERP_ON)THEN
            WRITE(ERR_MSG,2001) 'DES_INTERP_ON'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(DES_INTERP_MEAN_FIELDS)THEN
            WRITE(ERR_MSG,2001) 'DES_INTERP_MEAN_FIELDS'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF(DES_CONTINUUM_COUPLED) THEN
            IF(MPPIC) THEN
               WRITE(ERR_MSG,2002) 'MPPIC solids'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(MPPIC) THEN
               WRITE(ERR_MSG,2002) 'Cartesian grid cut-cells'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

      CASE ('GARG_2012')
         DES_INTERP_SCHEME_ENUM = DES_INTERP_GARG

      CASE ('DPVM')
         DES_INTERP_SCHEME_ENUM = DES_INTERP_DPVM

      CASE ('DPVM_GAUSS')
         DES_INTERP_SCHEME_ENUM = DES_INTERP_GAUSS

      CASE DEFAULT
         WRITE(ERR_MSG,2000) trim(adjustl(DES_INTERP_SCHEME))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      END SELECT

 2000 FORMAT('Error 2000: Invalid DES_INTERP_SCHEME: ',A,/'Please ',   &
         'correct the mfix.dat file.')

 2001 FORMAT('Error 2001: No interpolation scheme specified when ',A,/ &
         'is enabled. Please correct the mfix.dat file.')

 2002 FORMAT('Error 2002: DES simulations utilizing ',A,' require',/   &
         'interpolation (DES_INTERP_ON and DES_INTERP_MEANFIELDS). ',/ &
         'Please correct the mfix.dat file.')


      SELECT CASE(DES_INTERP_SCHEME_ENUM)

      CASE(DES_INTERP_NONE)

         IF(.NOT.DES_DIFFUSE_MEAN_FIELDS .AND.                         &
            FILTER_WIDTH /= UNDEFINED) THEN
            WRITE(ERR_MSG,2100) trim(adjustl(DES_INTERP_SCHEME))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 2100 FORMAT('Error 2100: The selected interpolation scheme (',A,') ', &
         'only',/'supports an adjustable filter size when the mean ',  &
         'fields are diffused.',/'Please correct the input file.')

      CASE(DES_INTERP_GARG)
         DES_INTERP_MEAN_FIELDS= .TRUE.

         IF(DES_DIFFUSE_MEAN_FIELDS) THEN
            WRITE(ERR_MSG,2110) trim(adjustl(DES_INTERP_SCHEME))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 2110 FORMAT('Error 2110: The selected interpolation scheme (',A,') ', &
         'does not',/'support diffusive filtering of mean field ',     &
          'quantites. Please correct',/'the input file.')

         IF(FILTER_WIDTH /= UNDEFINED) THEN
            WRITE(ERR_MSG,2111) trim(adjustl(DES_INTERP_SCHEME))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 2111 FORMAT('Error 2111: The selected interpolation scheme (',A,') ', &
         'does not',/'support an adjustable filter size. Please ',     &
         'correct',/'the input file.')


      CASE(DES_INTERP_DPVM, DES_INTERP_GAUSS)
         DES_INTERP_MEAN_FIELDS= .TRUE.

         IF(FILTER_WIDTH == UNDEFINED) THEN
            WRITE(ERR_MSG,2120) trim(adjustl(DES_INTERP_SCHEME))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 2120 FORMAT('Error 2120: The selected interpolation scheme (',A,') ', &
         'requires',/'a FILTER_WIDTH. Please correct the input file.')

      END SELECT


      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_INTERP
