      MODULE param

! Parameters describing problem size: (set from user input)
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! Number of partitions in axial direction
      INTEGER :: DIMENSION_I  ! X-Axis
      INTEGER :: DIMENSION_J  ! Y-Axis
      INTEGER :: DIMENSION_K  ! Z-Axis
! Number of computational cells
      INTEGER :: DIMENSION_3  !
      INTEGER :: DIMENSION_3L ! Local subdomain
      INTEGER :: DIMENSION_3G ! Global
      INTEGER :: DIMENSION_3P ! POST_MFIX
! Number of computational cells for higher order schemes
      INTEGER :: DIMENSION_4
! Total number of solids phases
      INTEGER :: DIMENSION_M
! Maximum number of species.
      INTEGER :: DIMENSION_N_g ! Gas
      INTEGER :: DIMENSION_N_s ! Solids
! Maximum number of user-defined scalars
      INTEGER :: DIMENSION_SCALAR
      INTEGER :: DIM_SCALAR2
      

! Parameters limiting user-specifed input.
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! Maximum number of reactions defined in data file
      INTEGER, PARAMETER :: DIMENSION_RXN = 100
! Number of user defined constants
      INTEGER, PARAMETER :: DIMENSION_C = 500
! Maximum number of items for specifying initial conditions
      INTEGER, PARAMETER :: DIMENSION_IC = 500
! Maximum number of items for specifying boundary conditions
      INTEGER, PARAMETER :: DIMENSION_BC = 500
! Maximum number of items for specifying internal surfaces.
      INTEGER, PARAMETER :: DIMENSION_IS = 500
! Maximum number of items for specifying point sources
      INTEGER, PARAMETER :: DIMENSION_PS = 500
! Maximum number of solids phases
      INTEGER, PARAMETER :: DIM_M = 10
! Maximum number of gas species
      INTEGER, PARAMETER :: DIM_N_g = 100
! Maximum number of solids species per phase.
      INTEGER, PARAMETER :: DIM_N_s = 100
! Maximum of DIM_N_g and DIM_N_s
      INTEGER, PARAMETER :: DIM_N = max(DIM_N_g, DIM_N_s)
! Maximum number of species.
      INTEGER, PARAMETER :: DIM_N_ALL = 2*DIM_N
! Maximum of the number of cells in the x direction.
      INTEGER, PARAMETER :: DIM_I = 5000
! Maximum of the number of cells in the y direction.
      INTEGER, PARAMETER :: DIM_J = 5000
! Maximum of the number of cells in the z direction.
      INTEGER, PARAMETER :: DIM_K = 5000
! Maximum number of user-defined output files
      INTEGER, PARAMETER :: DIMENSION_USR = 5
! Maximum of the number of scalars that can be read
      INTEGER, PARAMETER :: DIM_SCALAR = 100

      END MODULE param                                                                           
