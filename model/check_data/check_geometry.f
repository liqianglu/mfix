!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_GEOMETRY                                          !
!  Purpose: Check the distributed parallel namelist variables.         !
!                                                                      !
!  Author: P. Nicoletti                               Date: 14-DEC-99  !
!  Reviewer: J.Musser                                 Date: 16-Jan-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_GEOMETRY(SHIFT)


! Global Variables:
!---------------------------------------------------------------------//
! Domain partitions in various directions.
      use geometry, only: DX, XLENGTH, XMIN
      use geometry, only: DY, YLENGTH
      use geometry, only: DZ, ZLENGTH

      use geometry, only: NO_I, IMIN1, IMAX, IMAX1, IMAX3
      use geometry, only: NO_J, JMIN1, JMAX, JMAX1, JMAX3
      use geometry, only: NO_K, KMIN1, KMAX, KMAX1, KMAX3

      use bc, only: Flux_g

! Runtime flag specifying 2D simulations
!      use geometry, only: NO_K

      use geometry, only: COORDINATES, CYLINDRICAL
      use geometry, only: CYCLIC_X, CYCLIC_X_PD, CYCLIC_X_MF
      use geometry, only: CYCLIC_Y, CYCLIC_Y_PD, CYCLIC_Y_MF
      use geometry, only: CYCLIC_Z, CYCLIC_Z_PD, CYCLIC_Z_MF
!      use geometry, only: COORDINATES

      use cutcell, only: PARTIAL_CHECK_03

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ONE, ZERO, UNDEFINED_I, UNDEFINED

      use param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none


      LOGICAL, intent(IN) :: SHIFT
      LOGICAL, external :: COMPARE

! Local Variables:
!---------------------------------------------------------------------//


      
! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_GEOMETRY")

      CALL GET_DXYZ_FROM_CONTROL_POINTS

      CALL CHECK_AXIS(IMAX, IMAX3, XLENGTH, DX, 'X', 'I', NO_I, SHIFT)
      CALL CHECK_AXIS(JMAX, JMAX3, YLENGTH, DY, 'Y', 'J', NO_J, SHIFT)
      CALL CHECK_AXIS(KMAX, KMAX3, ZLENGTH, DZ, 'Z', 'K', NO_K, SHIFT)

      IF(SHIFT) CALL SHIFT_DXYZ

!  Ensure that the cell sizes across cyclic boundaries are comparable
      IF(CYCLIC_X .OR. CYCLIC_X_PD) THEN 
         IF(DX(IMIN1) /= DX(IMAX1)) THEN 
            WRITE(ERR_MSG,1100) 'DX(IMIN1)',DX(IMIN1),'DX(IMAX1)',DX(IMAX1)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDIF 

      IF(CYCLIC_Y .OR. CYCLIC_Y_PD) THEN 
         IF(DY(JMIN1) /= DY(JMAX1)) THEN 
            WRITE(ERR_MSG,1100) 'DY(JMIN1)',DY(JMIN1),'DY(JMAX1)',DY(JMAX1)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDIF 

      IF(CYCLIC_Z .OR. CYCLIC_Z_PD .OR. CYLINDRICAL) THEN 
         IF (DZ(KMIN1) /= DZ(KMAX1)) THEN 
            WRITE(ERR_MSG,1100) 'DZ(KMIN1)',DZ(KMIN1),'DZ(KMAX1)',DZ(KMAX1)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF 
      ENDIF 

 1100 FORMAT('Error 1100: Cells adjacent to cyclic boundaries must ',  &
         'be of same size:',/2X,A,' = ',G12.5,/2x,A,' = ',G12.5,/      &
         'Please correct the mfix.dat file.')
 

      CALL FINL_ERR_MSG

      RETURN  



      END SUBROUTINE CHECK_GEOMETRY