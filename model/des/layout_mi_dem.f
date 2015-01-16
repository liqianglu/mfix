!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_DEM                                              !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose: Check the data provided for the des mass inflow boundary   !
!  condition and flag errors if the data is improper.  This module is  !
!  also used to convert the proveded information into the format       !
!  necessary for the dependent subrountines to function properly.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE LAYOUT_MI_DEM(BCV, BCV_I, MAX_DIA)

      use bc, only: BC_PLANE
      USE run, only: RUN_TYPE

      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: BCV
      INTEGER, INTENT(IN) :: BCV_I      ! BC loop counter
      DOUBLE PRECISION, INTENT(IN) :: MAX_DIA

      CALL INIT_ERR_MSG("LAYOUT_MI_DEM")

! This subroutine determines the pattern that the particles will need to
! enter the system, if any. This routine only needs to be called if a
! run is new.  If a run is a RESTART_1, all of the setup information
! provided by this subroutine is will be obtained from the *_DES.RES file.
! This is done due to this routine's strong dependence on the
! RANDOM_NUMBER() subroutine.
      IF(RUN_TYPE == 'RESTART_1') THEN
         CALL SET_DEM_MI_OWNER(BCV, BCV_I)
      ELSE
         SELECT CASE (BC_PLANE(BCV))
         CASE('N','S'); CALL LAYOUT_DEM_MI_NS(BCV, BCV_I, MAX_DIA)
         CASE('E','W'); CALL LAYOUT_DEM_MI_EW(BCV, BCV_I, MAX_DIA)
         CASE('T','B'); CALL LAYOUT_DEM_MI_TB(BCV, BCV_I, MAX_DIA)
         END SELECT
      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE LAYOUT_MI_DEM

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: LAYOUT_DEM_MI_NS                                        !
!                                                                      !
!  Purpose:  This routine determines the layout of the mass inlet as   !
!  either ordered or random based upon the inlet conditions.  This     !
!  routine also verifies that the specified inlet conditions for mass  !
!  or volumetric flow rates along with inlet size (length or area) and !
!  particle inlet velocity will work.  If not an error is flagged and  !
!  the program is exited.                                              !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE LAYOUT_DEM_MI_NS(BCV, BCV_I, MAX_DIA)

      use bc, only: BC_PLANE, BC_Y_s, BC_AREA
      use bc, only: BC_X_w, BC_X_e
      use bc, only: BC_Z_b, BC_Z_t

      use des_bc, only: DEM_MI

      use compar
      use geometry
      use indices

      use funits, only: DMP_LOG

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager
      use functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! passed arguments giving index information of the boundary
      INTEGER, INTENT(IN) :: BCV
      INTEGER, INTENT(IN) :: BCV_I
! Max diameter of incoming particles at bc
      DOUBLE PRECISION, INTENT(IN) :: MAX_DIA
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters.
      INTEGER LL, LC
! Indices for mapping to fluid grid.
      INTEGER IJK, I, J, K
! Local MESH indices
      INTEGER H, W
! Temp variable for double precision
      DOUBLE PRECISION :: TMP_DP
! Temporary variable for integers
      INTEGER :: TMP_INT

      INTEGER, allocatable :: MESH_H(:)
      INTEGER, allocatable :: MESH_W(:)

      DOUBLE PRECISION, allocatable :: MESH_P(:)
      DOUBLE PRECISION, allocatable :: MESH_Q(:)

      INTEGER, allocatable :: RAND_MAP(:)
      INTEGER, allocatable :: FULL_MAP(:,:)

! max number of partitions along length of inlet
      INTEGER :: WMAX, HMAX

! the length of each side of the inlet boundary
      DOUBLE PRECISION :: PLEN, QLEN

! Number of occupied mesh cells
      INTEGER :: OCCUPANTS

      DOUBLE PRECISION :: SHIFT, WINDOW, OFFSET

      LOGICAL, EXTERNAL :: COMPARE

      LOGICAL, parameter :: setDBG = .FALSE.
      LOGICAL :: dFlag

!-----------------------------------------------

! Initialize the error manager.
      CALL INIT_ERR_MSG('LAYOUT_DEM_MI_NS')

      dFlag = (DMP_LOG .AND. setDBG)
      if(dFlag) write(*,"(2/,'Building DEM_MI: ',I3)") BCV_I

! Store the index that maps back to the user input.

      OCCUPANTS = 0

! Calculate number of partitions in first direction.
      PLEN = BC_X_e(BCV) - BC_X_w(BCV)
      WMAX = FLOOR(real(PLEN/MAX_DIA))
      allocate( MESH_W(WMAX) )
      allocate( MESH_P(0:WMAX) )

! Calculate number of partitions in the second direction.
      QLEN = merge(BC_Z_t(BCV) - BC_Z_b(BCV), MAX_DIA, DO_K)
      HMAX = FLOOR(real(QLEN/MAX_DIA))
      allocate( MESH_H(HMAX) )
      allocate( MESH_Q(0:HMAX) )

! Allocate the full map.
      allocate( FULL_MAP(WMAX, HMAX))

! Set the value of the boundary condtion offset value used in the
! placement of new particles.
      CALL CALC_CELL_INTERSECT(ZERO, BC_Y_s(BCV), DY, JMAX, J)
      SHIFT = merge(-ONE, ONE, BC_PLANE(BCV) == 'N')
      DEM_MI(BCV_I)%OFFSET = BC_Y_s(BCV) + MAX_DIA*SHIFT
      DEM_MI(BCV_I)%L = J + int(SHIFT)
      if(dFlag) write(*,"(2x,'Offset: ',3x,I4,3x,g12.5)") &
         DEM_MI(BCV_I)%L, DEM_MI(BCV_I)%OFFSET


! Dimension of grid cell for seeding particles; this may be larger than
! than the particle diameter but not smaller:
      DEM_MI(BCV_I)%WINDOW = MIN(PLEN/WMAX, QLEN/HMAX)
      WINDOW = DEM_MI(BCV_I)%WINDOW
      if(dFlag) write(*,"(2x,'Windows size: ',g12.5)") WINDOW

! Setup the first direction.
      SHIFT = HALF*(PLEN - WMAX*WINDOW)
      MESH_P(0) = BC_X_w(BCV) + SHIFT
      if(dFlag) write(*,8005) 'P', SHIFT, 'P', MESH_P(0)
      DO LC=1,WMAX
         MESH_P(LC) = MESH_P(0) + dble(LC-1)*WINDOW
         SHIFT = MESH_P(LC) + HALF*WINDOW
         CALL CALC_CELL_INTERSECT(XMIN, SHIFT, DX, IMAX, MESH_W(LC))
         IF(dFlag)WRITE(*,8006) LC, 'W', MESH_W(LC), 'P', MESH_P(LC)
      ENDDO

! Setup the second direction.
      IF(DO_K) THEN
         SHIFT = HALF*(QLEN - HMAX*WINDOW)
         MESH_Q(0) = BC_Z_b(BCV) + SHIFT
         if(dFlag) write(*,8005) 'Q',SHIFT, 'Q',MESH_Q(0)
         DO LC=1,HMAX
            MESH_Q(LC) = MESH_Q(0) + dble(LC-1)*WINDOW
            SHIFT = MESH_Q(LC) + HALF*WINDOW
            CALL CALC_CELL_INTERSECT(ZERO, SHIFT, DZ, KMAX, MESH_H(LC))
            IF(dFlag)WRITE(*,8006) LC, 'H', MESH_H(LC), 'Q', MESH_Q(LC)
         ENDDO
      ELSE
         MESH_H(1) = 1
         MESH_Q(1) = 0.0d0
      ENDIF

! Get the Jth index of the fluid cell
      CALL CALC_CELL_INTERSECT(ZERO, BC_Y_s(BCV), DY, JMAX, J)

! If the computationsl cell adjacent to the DEM_MI mesh cell is a
! fluid cell and has not been cut, store the ID of the cell owner.
      DO H=1,HMAX
      DO W=1,WMAX

         FULL_MAP(W,H) = 0

         IF(.NOT.IS_ON_myPE_owns(MESH_W(W),J,MESH_H(H))) CYCLE

         IF(DO_K) THEN

            CALL CALC_CELL_INTERSECT(XMIN, MESH_P(W), DX, IMAX, I)
            CALL CALC_CELL_INTERSECT(ZERO, MESH_Q(H), DZ, KMAX, K)
            IF(EXCLUDE_DEM_MI_CELL(I, J, K, DX(I)*DZ(K))) CYCLE

            SHIFT = MESH_P(W)+WINDOW
            CALL CALC_CELL_INTERSECT(XMIN, SHIFT, DX, IMAX, I)
            IF(EXCLUDE_DEM_MI_CELL(I, J, K, DX(I)*DZ(K))) CYCLE

            SHIFT = MESH_Q(W)+WINDOW
            CALL CALC_CELL_INTERSECT(ZERO, SHIFT, DZ, KMAX, K)
            IF(EXCLUDE_DEM_MI_CELL(I, J, K, DX(I)*DZ(K))) CYCLE

            CALL CALC_CELL_INTERSECT(XMIN, MESH_P(W), DX, IMAX, I)
            IF(EXCLUDE_DEM_MI_CELL(I, J, K, DX(I)*DZ(K))) CYCLE

         ELSE

            K = MESH_H(1)
            CALL CALC_CELL_INTERSECT(XMIN, MESH_P(W), DX, IMAX, I)
            IF(EXCLUDE_DEM_MI_CELL(I, J, K, DX(I)*DZ(K))) CYCLE

            SHIFT = MESH_P(W)+WINDOW
            CALL CALC_CELL_INTERSECT(XMIN, SHIFT, DX, IMAX, I)
            IF(EXCLUDE_DEM_MI_CELL(I, J, K, DX(I)*DZ(K))) CYCLE
         ENDIF

         I = MESH_W(W)
         K = MESH_H(H)
         IF(EXCLUDE_DEM_MI_CELL(I,J,K, DX(I)*DZ(K))) CYCLE

         OCCUPANTS = OCCUPANTS + 1
         FULL_MAP(W,H) = myPE+1
      ENDDO
      ENDDO

! Sync the full map across all ranks.
      CALL GLOBAL_ALL_SUM(OCCUPANTS)
      CALL GLOBAL_ALL_SUM(FULL_MAP)

! Throw an error and exit if there are no occupants.
      IF(OCCUPANTS == 0) THEN
         WRITE(ERR_MSG, 1100) BCV_I
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: No un-cut fluid cells adjacent to DEM_MI ',  &
         'staging area.',/'Unable to setup the discrete solids mass ', &
         'inlet for BC:',I3)

! Store the number of occupants.
      DEM_MI(BCV_I)%OCCUPANTS = OCCUPANTS

! Display the fill map if debugging
      IF(dFlag) THEN
         WRITE(*,"(2/,2x,'Displaying Fill Map:')")
         DO H=HMAX,1,-1
            WRITE(*,"(2x,'H =',I3)",advance='no')H
            DO W=1,WMAX
               IF(FULL_MAP(W,H) == 0) then
                  WRITE(*,"(' *')",advance='no')
               ELSE
                  WRITE(*,"(' .')",advance='no')
               ENDIF
            ENDDO
            WRITE(*,*)' '
         ENDDO
      ENDIF

! Construct an array of integers randomly ordered.
      if(dFLAG) write(*,"(2/,2x,'Building RAND_MAP:')")
      allocate( RAND_MAP(OCCUPANTS) )
      RAND_MAP = 0

! Only Rank 0 will randomize the layout and broadcast it to the other
! ranks. This will ensure that all ranks see the same layout.
      IF(myPE == 0) THEN
         LL = 1
         DO WHILE (RAND_MAP(OCCUPANTS) .EQ. 0)
            CALL RANDOM_NUMBER(TMP_DP)
            TMP_INT = CEILING(real(TMP_DP*dble(OCCUPANTS)))
            DO LC = 1, LL
              IF(TMP_INT .EQ. RAND_MAP(LC) )EXIT
              IF(LC .EQ. LL)THEN
                 if(dFlag) WRITE(*,"(4x,'LC:',I9,' : ',I9)") LC, TMP_INT
                 RAND_MAP(LC) = TMP_INT
                 LL = LL + 1
              ENDIF
            ENDDO
         ENDDO
      ENDIF

      CALL GLOBAL_ALL_SUM(RAND_MAP)

! Initialize the vacancy pointer.
      DEM_MI(BCV_I)%VACANCY = 1

! Allocate the mass inlet storage arrays.
      allocate( DEM_MI(BCV_I)%W(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%P(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%H(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%Q(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%OWNER(OCCUPANTS) )

      if(dFlag) write(*,8010)
! Store the MI layout in the randomized order.
      LC = 0
      DO H=1,HMAX
      DO W=1,WMAX
         IF(FULL_MAP(W,H) == 0) CYCLE
         LC = LC + 1
         LL = RAND_MAP(LC)
         DEM_MI(BCV_I)%OWNER(LL) = FULL_MAP(W,H) - 1

         DEM_MI(BCV_I)%W(LL) = MESH_W(W)
         DEM_MI(BCV_I)%H(LL) = MESH_H(H)

         DEM_MI(BCV_I)%P(LL) = MESH_P(W)
         DEM_MI(BCV_I)%Q(LL) = MESH_Q(H)

         if(dFlag) write(*,8011) DEM_MI(BCV_I)%OWNER(LL), &
            DEM_MI(BCV_I)%W(LL), DEM_MI(BCV_I)%H(LL), DEM_MI(BCV_I)%L, &
            DEM_MI(BCV_I)%P(LL), DEM_MI(BCV_I)%Q(LL), DEM_MI(BCV_I)%OFFSET

      ENDDO
      ENDDO


 8010 FORMAT(2/,2x,'Storing DEM_MI data:',/4X,'OWNER',5X,'W',5X,'H',   &
         5X,'L',7X,'P',12X,'Q',12X,'R')
 8011 FORMAT(4x,I5,3(2X,I4),3(2x,g12.5))


      if(dFlag) write(*,"(2/,2x,'Inlet area sizes:')")
      if(dFlag) write(*,9000) 'mfix.dat: ', PLEN * QLEN
      if(dFlag) write(*,9000) 'BC_AREA:  ', BC_AREA(BCV)
      if(dFlag) write(*,9000) 'DEM_MI:   ', OCCUPANTS * (WINDOW**2)
 9000 FORMAT(2x,A,g12.5)

! House keeping.
      IF( allocated(MESH_H)) deallocate(MESH_H)
      IF( allocated(MESH_W)) deallocate(MESH_W)
      IF( allocated(MESH_P)) deallocate(MESH_P)
      IF( allocated(MESH_Q)) deallocate(MESH_Q)

      IF( allocated(RAND_MAP)) deallocate(RAND_MAP)
      IF( allocated(FULL_MAP)) deallocate(FULL_MAP)

      CALL FINL_ERR_MSG

      RETURN

 8005 FORMAT(2/,2x,'Building MESH_',A1,':',/4x,'Shift:',f8.4,/4x,      &
         'MESH_',A1,'(0) = ',f8.4,/)

 8006 FORMAT(4x,'LC = ',I4,3x,A1,' =',I3,3x,A1,' =',f8.4)

      contains

!----------------------------------------------------------------------!
!  Function to exclude cells from DEM mass inlet.                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION EXCLUDE_DEM_MI_CELL(lI, lJ, lK, lAREA)
! Indicies of cell to check
      INTEGER, INTENT(IN) :: lI, lJ, lK
! Calculated face area based on grid dimensions
      DOUBLE PRECISION, INTENT(IN) :: lAREA
! Local value for IJK
      INTEGER :: IJK

      EXCLUDE_DEM_MI_CELL = .TRUE.

      IF(.NOT.IS_ON_myPE_plus2layers(lI,lJ,lK)) RETURN
      IF(DEAD_CELL_AT(lI,lJ,lK)) RETURN
      IJK = FUNIJK(lI,lJ,lK)
      IF(.NOT.FLUID_AT(IJK)) RETURN
      IF(.NOT.COMPARE(AXZ(IJK),lAREA)) RETURN

      EXCLUDE_DEM_MI_CELL = .FALSE.
      RETURN
      END FUNCTION EXCLUDE_DEM_MI_CELL

      END SUBROUTINE LAYOUT_DEM_MI_NS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: LAYOUT_DEM_MI_EW                                        !
!                                                                      !
!  Purpose:  This routine determines the layout of the mass inlet as   !
!  either ordered or random based upon the inlet conditions.  This     !
!  routine also verifies that the specified inlet conditions for mass  !
!  or volumetric flow rates along with inlet size (length or area) and !
!  particle inlet velocity will work.  If not an error is flagged and  !
!  the program is exited.                                              !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE LAYOUT_DEM_MI_EW(BCV, BCV_I, MAX_DIA)

      use bc, only: BC_PLANE, BC_X_w, BC_AREA
      use bc, only: BC_Y_s, BC_Y_n
      use bc, only: BC_Z_b, BC_Z_t

      use des_bc, only: DEM_MI

      use compar
      use geometry
      use indices

      use funits, only: DMP_LOG

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager
      use functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! passed arguments giving index information of the boundary
      INTEGER, INTENT(IN) :: BCV
      INTEGER, INTENT(IN) :: BCV_I
! Max diameter of incoming particles at bc
      DOUBLE PRECISION, INTENT(IN) :: MAX_DIA
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters.
      INTEGER LL, LC
! Indices for mapping to fluid grid.
      INTEGER IJK, I, J, K
! Local MESH indices
      INTEGER H, W
! Temp variable for double precision
      DOUBLE PRECISION :: TMP_DP
! Temporary variable for integers
      INTEGER :: TMP_INT

      INTEGER, allocatable :: MESH_H(:)
      INTEGER, allocatable :: MESH_W(:)

      DOUBLE PRECISION, allocatable :: MESH_P(:)
      DOUBLE PRECISION, allocatable :: MESH_Q(:)

      INTEGER, allocatable :: RAND_MAP(:)
      INTEGER, allocatable :: FULL_MAP(:,:)

! max number of partitions along length of inlet
      INTEGER :: WMAX, HMAX

! the length of each side of the inlet boundary
      DOUBLE PRECISION :: PLEN, QLEN

! Number of occupied mesh cells
      INTEGER :: OCCUPANTS

      DOUBLE PRECISION :: SHIFT, WINDOW, OFFSET

      LOGICAL, EXTERNAL :: COMPARE

      LOGICAL, parameter :: setDBG = .FALSE.
      LOGICAL :: dFlag

!-----------------------------------------------

! Initialize the error manager.
      CALL INIT_ERR_MSG('LAYOUT_DEM_MI_EW')

      dFlag = (DMP_LOG .AND. setDBG)
      if(dFlag) write(*,"(2/,'Building DEM_MI: ',I3)") BCV_I

      OCCUPANTS = 0

! Calculate number of partitions in first direction.
      PLEN = BC_Y_n(BCV) - BC_Y_s(BCV)
      WMAX = FLOOR(real(PLEN/MAX_DIA))
      allocate( MESH_W(WMAX) )
      allocate( MESH_P(0:WMAX) )

! Calculate number of partitions in the second direction.
      QLEN = merge(BC_Z_t(BCV) - BC_Z_b(BCV), MAX_DIA, DO_K)
      HMAX = FLOOR(real(QLEN/MAX_DIA))
      allocate( MESH_H(HMAX) )
      allocate( MESH_Q(0:HMAX) )

! Allocate the full map.
      allocate( FULL_MAP(WMAX, HMAX))

! Set the value of the boundary condtion offset value used in the
! placement of new particles.
      CALL CALC_CELL_INTERSECT(XMIN, BC_X_w(BCV), DX, IMAX, I)
      SHIFT = merge(-ONE, ONE, BC_PLANE(BCV) == 'E')
      DEM_MI(BCV_I)%OFFSET = BC_X_w(BCV) + MAX_DIA*SHIFT
      DEM_MI(BCV_I)%L = I + int(SHIFT)
      if(dFlag) write(*,"(2x,'Offset: ',3x,I4,3x,g12.5)") &
         DEM_MI(BCV_I)%L, DEM_MI(BCV_I)%OFFSET


! Dimension of grid cell for seeding particles; this may be larger than
! than the particle diameter but not smaller:
      DEM_MI(BCV_I)%WINDOW = MIN(PLEN/WMAX, QLEN/HMAX)
      WINDOW = DEM_MI(BCV_I)%WINDOW
      if(dFlag) write(*,"(2x,'Windows size: ',g12.5)") WINDOW

! Setup the first direction.
      SHIFT = HALF*(PLEN - WMAX*WINDOW)
      MESH_P(0) = BC_Y_s(BCV) + SHIFT
      if(dFlag) write(*,8005) 'P', SHIFT, 'P', MESH_P(0)
      DO LC=1,WMAX
         MESH_P(LC) = MESH_P(0) + dble(LC-1)*WINDOW
         SHIFT = MESH_P(LC) + HALF*WINDOW
         CALL CALC_CELL_INTERSECT(ZERO, SHIFT, DY, JMAX, MESH_W(LC))
         IF(dFlag)WRITE(*,8006) LC, 'W', MESH_W(LC), 'P', MESH_P(LC)
      ENDDO

! Setup the second direction.
      IF(DO_K) THEN
         SHIFT = HALF*(QLEN - HMAX*WINDOW)
         MESH_Q(0) = BC_Z_b(BCV) + SHIFT
         if(dFlag) write(*,8005) 'Q',SHIFT, 'Q',MESH_Q(0)
         DO LC=1,HMAX
            MESH_Q(LC) = MESH_Q(0) + dble(LC-1)*WINDOW
            SHIFT = MESH_Q(LC) + HALF*WINDOW
            CALL CALC_CELL_INTERSECT(ZERO, SHIFT, DZ, KMAX, MESH_H(LC))
            IF(dFlag)WRITE(*,8006) LC, 'H', MESH_H(LC), 'Q', MESH_Q(LC)
         ENDDO
      ELSE
         MESH_H(1) = 1
         MESH_Q(1) = 0.0d0
      ENDIF

! Get the Jth index of the fluid cell
      CALL CALC_CELL_INTERSECT(XMIN, BC_X_w(BCV), DX, IMAX, I)

! If the computationsl cell adjacent to the DEM_MI mesh cell is a
! fluid cell and has not been cut, store the ID of the cell owner.
      DO H=1,HMAX
      DO W=1,WMAX
         J = MESH_W(W)
         K = MESH_H(H)
         FULL_MAP(W,H) = 0
         IF(.NOT.IS_ON_myPE_owns(I,J,K)) CYCLE
         IF(DEAD_CELL_AT(I,J,K)) CYCLE
         IJK = FUNIJK(I,J,K)
         IF(.NOT.FLUID_AT(IJK)) CYCLE
         IF(.NOT.COMPARE(AYZ(IJK),DY(J)*DZ(K))) CYCLE
         OCCUPANTS = OCCUPANTS + 1
         FULL_MAP(W,H) = myPE+1
      ENDDO
      ENDDO

! Sync the full map across all ranks.
      CALL GLOBAL_ALL_SUM(OCCUPANTS)
      CALL GLOBAL_ALL_SUM(FULL_MAP)

! Throw an error and exit if there are no occupants.
      IF(OCCUPANTS == 0) THEN
         WRITE(ERR_MSG, 1100) BCV_I
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: No un-cut fluid cells adjacent to DEM_MI ',  &
         'staging area.',/'Unable to setup the discrete solids mass ', &
         'inlet for BC:',I3)

! Store the number of occupants.
      DEM_MI(BCV_I)%OCCUPANTS = OCCUPANTS

! Display the fill map if debugging
      IF(dFlag) THEN
         WRITE(*,"(2/,2x,'Displaying Fill Map:')")
         DO H=HMAX,1,-1
            WRITE(*,"(2x,'H =',I3)",advance='no')H
            DO W=1,WMAX
               IF(FULL_MAP(W,H) == 0) then
                  WRITE(*,"(' *')",advance='no')
               ELSE
                  WRITE(*,"(' .')",advance='no')
               ENDIF
            ENDDO
            WRITE(*,*)' '
         ENDDO
      ENDIF

! Construct an array of integers randomly ordered.
      if(dFLAG) write(*,"(2/,2x,'Building RAND_MAP:')")
      allocate( RAND_MAP(OCCUPANTS) )
      RAND_MAP = 0

! Only Rank 0 will randomize the layout and broadcast it to the other
! ranks. This will ensure that all ranks see the same layout.
      IF(myPE == 0) THEN
         LL = 1
         DO WHILE (RAND_MAP(OCCUPANTS) .EQ. 0)
            CALL RANDOM_NUMBER(TMP_DP)
            TMP_INT = CEILING(real(TMP_DP*dble(OCCUPANTS)))
            DO LC = 1, LL
              IF(TMP_INT .EQ. RAND_MAP(LC) )EXIT
              IF(LC .EQ. LL)THEN
                 if(dFlag) WRITE(*,"(4x,'LC:',I3,' : ',I3)") LC, TMP_INT
                 RAND_MAP(LC) = TMP_INT
                 LL = LL + 1
              ENDIF
            ENDDO
         ENDDO
      ENDIF

      CALL GLOBAL_ALL_SUM(RAND_MAP)

! Initialize the vacancy pointer.
      DEM_MI(BCV_I)%VACANCY = 1

! Allocate the mass inlet storage arrays.
      allocate( DEM_MI(BCV_I)%W(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%P(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%H(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%Q(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%OWNER(OCCUPANTS) )

      if(dFlag) write(*,8010)
! Store the MI layout in the randomized order.
      LC = 0
      DO H=1,HMAX
      DO W=1,WMAX
         IF(FULL_MAP(W,H) == 0) CYCLE
         LC = LC + 1
         LL = RAND_MAP(LC)
         DEM_MI(BCV_I)%OWNER(LL) = FULL_MAP(W,H) - 1

         DEM_MI(BCV_I)%W(LL) = MESH_W(W)
         DEM_MI(BCV_I)%H(LL) = MESH_H(H)

         DEM_MI(BCV_I)%P(LL) = MESH_P(W)
         DEM_MI(BCV_I)%Q(LL) = MESH_Q(H)

         if(dFlag) write(*,8011) DEM_MI(BCV_I)%OWNER(LL), &
            DEM_MI(BCV_I)%W(LL), DEM_MI(BCV_I)%H(LL), DEM_MI(BCV_I)%L, &
            DEM_MI(BCV_I)%P(LL), DEM_MI(BCV_I)%Q(LL), DEM_MI(BCV_I)%OFFSET

      ENDDO
      ENDDO


 8010 FORMAT(2/,2x,'Storing DEM_MI data:',/4X,'OWNER',5X,'W',5X,'H',   &
         5X,'L',7X,'P',12X,'Q',12X,'R')
 8011 FORMAT(4x,I5,3(2X,I4),3(2x,g12.5))


      if(dFlag) write(*,"(2/,2x,'Inlet area sizes:')")
      if(dFlag) write(*,9000) 'mfix.dat: ', PLEN * QLEN
      if(dFlag) write(*,9000) 'BC_AREA:  ', BC_AREA(BCV)
      if(dFlag) write(*,9000) 'DEM_MI:   ', OCCUPANTS * (WINDOW**2)
 9000 FORMAT(2x,A,g12.5)

! House keeping.
      IF( allocated(MESH_H)) deallocate(MESH_H)
      IF( allocated(MESH_W)) deallocate(MESH_W)
      IF( allocated(MESH_P)) deallocate(MESH_P)
      IF( allocated(MESH_Q)) deallocate(MESH_Q)

      IF( allocated(RAND_MAP)) deallocate(RAND_MAP)
      IF( allocated(FULL_MAP)) deallocate(FULL_MAP)

      CALL FINL_ERR_MSG

      RETURN

 8005 FORMAT(2/,2x,'Building MESH_',A1,':',/4x,'Shift:',f8.4,/4x,      &
         'MESH_',A1,'(0) = ',f8.4,/)

 8006 FORMAT(4x,'LC = ',I4,3x,A1,' =',I3,3x,A1,' =',f8.4)

      END SUBROUTINE LAYOUT_DEM_MI_EW

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: LAYOUT_DEM_MI_TB                                        !
!                                                                      !
!  Purpose:  This routine determines the layout of the mass inlet as   !
!  either ordered or random based upon the inlet conditions.  This     !
!  routine also verifies that the specified inlet conditions for mass  !
!  or volumetric flow rates along with inlet size (length or area) and !
!  particle inlet velocity will work.  If not an error is flagged and  !
!  the program is exited.                                              !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE LAYOUT_DEM_MI_TB(BCV, BCV_I, MAX_DIA)

      use bc, only: BC_PLANE, BC_Z_b, BC_AREA
      use bc, only: BC_X_w, BC_X_e
      use bc, only: BC_Y_s, BC_Y_n

      use des_bc, only: DEM_MI

      use compar
      use geometry
      use indices

      use funits, only: DMP_LOG

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager
      use functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! passed arguments giving index information of the boundary
      INTEGER, INTENT(IN) :: BCV
      INTEGER, INTENT(IN) :: BCV_I
! Max diameter of incoming particles at bc
      DOUBLE PRECISION, INTENT(IN) :: MAX_DIA
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters.
      INTEGER LL, LC
! Indices for mapping to fluid grid.
      INTEGER IJK, I, J, K
! Local MESH indices
      INTEGER H, W
! Temp variable for double precision
      DOUBLE PRECISION :: TMP_DP
! Temporary variable for integers
      INTEGER :: TMP_INT

      INTEGER, allocatable :: MESH_H(:)
      INTEGER, allocatable :: MESH_W(:)

      DOUBLE PRECISION, allocatable :: MESH_P(:)
      DOUBLE PRECISION, allocatable :: MESH_Q(:)

      INTEGER, allocatable :: RAND_MAP(:)
      INTEGER, allocatable :: FULL_MAP(:,:)

! max number of partitions along length of inlet
      INTEGER :: WMAX, HMAX

! the length of each side of the inlet boundary
      DOUBLE PRECISION :: PLEN, QLEN

! Number of occupied mesh cells
      INTEGER :: OCCUPANTS

      DOUBLE PRECISION :: SHIFT, WINDOW, OFFSET

      LOGICAL, EXTERNAL :: COMPARE

      LOGICAL, parameter :: setDBG = .FALSE.
      LOGICAL :: dFlag

!-----------------------------------------------

! Initialize the error manager.
      CALL INIT_ERR_MSG('LAYOUT_DEM_MI_TB')

      dFlag = (DMP_LOG .AND. setDBG)
      if(dFlag) write(*,"(2/,'Building DEM_MI: ',I3)") BCV_I

! Store the index that maps back to the user input.

      OCCUPANTS = 0

! Calculate number of partitions in first direction.
      PLEN = BC_X_e(BCV) - BC_X_w(BCV)
      WMAX = FLOOR(real(PLEN/MAX_DIA))
      allocate( MESH_W(WMAX) )
      allocate( MESH_P(0:WMAX) )

! Calculate number of partitions in the second direction.
      QLEN = BC_Y_n(BCV) - BC_Y_s(BCV)
      HMAX = FLOOR(real(QLEN/MAX_DIA))
      allocate( MESH_H(HMAX) )
      allocate( MESH_Q(0:HMAX) )

! Allocate the full map.
      allocate( FULL_MAP(WMAX, HMAX))

! Set the value of the boundary condtion offset value used in the
! placement of new particles.
      CALL CALC_CELL_INTERSECT(ZERO, BC_Z_b(BCV), DZ, KMAX, K)
      SHIFT = merge(-ONE, ONE, BC_PLANE(BCV) == 'T')
      DEM_MI(BCV_I)%OFFSET = BC_Z_b(BCV) + MAX_DIA*SHIFT
      DEM_MI(BCV_I)%L = K + int(SHIFT)
      if(dFlag) write(*,"(2x,'Offset: ',3x,I4,3x,g12.5)") &
         DEM_MI(BCV_I)%L, DEM_MI(BCV_I)%OFFSET


! Dimension of grid cell for seeding particles; this may be larger than
! than the particle diameter but not smaller:
      DEM_MI(BCV_I)%WINDOW = MIN(PLEN/WMAX, QLEN/HMAX)
      WINDOW = DEM_MI(BCV_I)%WINDOW
      if(dFlag) write(*,"(2x,'Windows size: ',g12.5)") WINDOW

! Setup the first direction.
      SHIFT = HALF*(PLEN - WMAX*WINDOW)
      MESH_P(0) = BC_X_w(BCV) + SHIFT
      if(dFlag) write(*,8005) 'P', SHIFT, 'P', MESH_P(0)
      DO LC=1,WMAX
         MESH_P(LC) = MESH_P(0) + dble(LC-1)*WINDOW
         SHIFT = MESH_P(LC) + HALF*WINDOW
         CALL CALC_CELL_INTERSECT(XMIN, SHIFT, DX, IMAX, MESH_W(LC))
         IF(dFlag)WRITE(*,8006) LC, 'W', MESH_W(LC), 'P', MESH_P(LC)
      ENDDO

! Setup the second direction.
      SHIFT = HALF*(QLEN - HMAX*WINDOW)
      MESH_Q(0) = BC_Y_s(BCV) + SHIFT
      if(dFlag) write(*,8005) 'Q',SHIFT, 'Q',MESH_Q(0)
      DO LC=1,HMAX
         MESH_Q(LC) = MESH_Q(0) + dble(LC-1)*WINDOW
         SHIFT = MESH_Q(LC) + HALF*WINDOW
         CALL CALC_CELL_INTERSECT(ZERO, SHIFT, DY, JMAX, MESH_H(LC))
         IF(dFlag)WRITE(*,8006) LC, 'H', MESH_H(LC), 'Q', MESH_Q(LC)
      ENDDO

! Get the Jth index of the fluid cell
      CALL CALC_CELL_INTERSECT(ZERO, BC_Z_b(BCV), DZ, KMAX, K)

! If the computationsl cell adjacent to the DEM_MI mesh cell is a
! fluid cell and has not been cut, store the ID of the cell owner.
      DO H=1,HMAX
      DO W=1,WMAX
         I = MESH_W(W)
         J = MESH_H(H)
         FULL_MAP(W,H) = 0
         IF(.NOT.IS_ON_myPE_owns(I,J,K)) CYCLE
         IF(DEAD_CELL_AT(I,J,K)) CYCLE
         IJK = FUNIJK(I,J,K)
         IF(.NOT.FLUID_AT(IJK)) CYCLE
         IF(.NOT.COMPARE(AXY(IJK),DX(I)*DY(J))) CYCLE
         OCCUPANTS = OCCUPANTS + 1
         FULL_MAP(W,H) = myPE+1
      ENDDO
      ENDDO

! Sync the full map across all ranks.
      CALL GLOBAL_ALL_SUM(OCCUPANTS)
      CALL GLOBAL_ALL_SUM(FULL_MAP)

! Throw an error and exit if there are no occupants.
      IF(OCCUPANTS == 0) THEN
         WRITE(ERR_MSG, 1100) BCV_I
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: No un-cut fluid cells adjacent to DEM_MI ',  &
         'staging area.',/'Unable to setup the discrete solids mass ', &
         'inlet for BC:',I3)

! Store the number of occupants.
      DEM_MI(BCV_I)%OCCUPANTS = OCCUPANTS

! Display the fill map if debugging
      IF(dFlag) THEN
         WRITE(*,"(2/,2x,'Displaying Fill Map:')")
         DO H=HMAX,1,-1
            WRITE(*,"(2x,'H =',I3)",advance='no')H
            DO W=1,WMAX
               IF(FULL_MAP(W,H) == 0) then
                  WRITE(*,"(' *')",advance='no')
               ELSE
                  WRITE(*,"(' .')",advance='no')
               ENDIF
            ENDDO
            WRITE(*,*)' '
         ENDDO
      ENDIF

! Construct an array of integers randomly ordered.
      if(dFLAG) write(*,"(2/,2x,'Building RAND_MAP:')")
      allocate( RAND_MAP(OCCUPANTS) )
      RAND_MAP = 0

! Only Rank 0 will randomize the layout and broadcast it to the other
! ranks. This will ensure that all ranks see the same layout.
      IF(myPE == 0) THEN
         LL = 1
         DO WHILE (RAND_MAP(OCCUPANTS) .EQ. 0)
            CALL RANDOM_NUMBER(TMP_DP)
            TMP_INT = CEILING(real(TMP_DP*dble(OCCUPANTS)))
            DO LC = 1, LL
              IF(TMP_INT .EQ. RAND_MAP(LC) )EXIT
              IF(LC .EQ. LL)THEN
                 if(dFlag) WRITE(*,"(4x,'LC:',I3,' : ',I3)") LC, TMP_INT
                 RAND_MAP(LC) = TMP_INT
                 LL = LL + 1
              ENDIF
            ENDDO
         ENDDO
      ENDIF

      CALL GLOBAL_ALL_SUM(RAND_MAP)

! Initialize the vacancy pointer.
      DEM_MI(BCV_I)%VACANCY = 1

! Allocate the mass inlet storage arrays.
      allocate( DEM_MI(BCV_I)%W(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%P(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%H(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%Q(OCCUPANTS) )
      allocate( DEM_MI(BCV_I)%OWNER(OCCUPANTS) )

      if(dFlag) write(*,8010)
! Store the MI layout in the randomized order.
      LC = 0
      DO H=1,HMAX
      DO W=1,WMAX
         IF(FULL_MAP(W,H) == 0) CYCLE
         LC = LC + 1
         LL = RAND_MAP(LC)
         DEM_MI(BCV_I)%OWNER(LL) = FULL_MAP(W,H) - 1

         DEM_MI(BCV_I)%W(LL) = MESH_W(W)
         DEM_MI(BCV_I)%H(LL) = MESH_H(H)

         DEM_MI(BCV_I)%P(LL) = MESH_P(W)
         DEM_MI(BCV_I)%Q(LL) = MESH_Q(H)

         if(dFlag) write(*,8011) DEM_MI(BCV_I)%OWNER(LL), &
            DEM_MI(BCV_I)%W(LL), DEM_MI(BCV_I)%H(LL), DEM_MI(BCV_I)%L, &
            DEM_MI(BCV_I)%P(LL), DEM_MI(BCV_I)%Q(LL), DEM_MI(BCV_I)%OFFSET

      ENDDO
      ENDDO


 8010 FORMAT(2/,2x,'Storing DEM_MI data:',/4X,'OWNER',5X,'W',5X,'H',   &
         5X,'L',7X,'P',12X,'Q',12X,'R')
 8011 FORMAT(4x,I5,3(2X,I4),3(2x,g12.5))


      if(dFlag) write(*,"(2/,2x,'Inlet area sizes:')")
      if(dFlag) write(*,9000) 'mfix.dat: ', PLEN * QLEN
      if(dFlag) write(*,9000) 'BC_AREA:  ', BC_AREA(BCV)
      if(dFlag) write(*,9000) 'DEM_MI:   ', OCCUPANTS * (WINDOW**2)
 9000 FORMAT(2x,A,g12.5)

! House keeping.
      IF( allocated(MESH_H)) deallocate(MESH_H)
      IF( allocated(MESH_W)) deallocate(MESH_W)
      IF( allocated(MESH_P)) deallocate(MESH_P)
      IF( allocated(MESH_Q)) deallocate(MESH_Q)

      IF( allocated(RAND_MAP)) deallocate(RAND_MAP)
      IF( allocated(FULL_MAP)) deallocate(FULL_MAP)

      CALL FINL_ERR_MSG

      RETURN

 8005 FORMAT(2/,2x,'Building MESH_',A1,':',/4x,'Shift:',f8.4,/4x,      &
         'MESH_',A1,'(0) = ',f8.4,/)

 8006 FORMAT(4x,'LC = ',I4,3x,A1,' =',I3,3x,A1,' =',f8.4)

      END SUBROUTINE LAYOUT_DEM_MI_TB


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_DEM_MI_OWNER                                        !
!  Author: J.Musser                                   Date: 08-Aug-14  !
!                                                                      !
!  Purpose: Set the owner of the background mesh used for seeding new  !
!  particles into the domain. The mesh is stored and read from the RES !
!  file, however, the owers need to be set per-run in case the DMP     !
!  partitions changed.                                                 !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_DEM_MI_OWNER(BCV, BCV_I)

      use bc, only: BC_PLANE
      use des_bc, only: DEM_MI

      use compar
      use geometry
      use indices

      use funits, only: DMP_LOG

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager
      use functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! passed arguments giving index information of the boundary
      INTEGER, INTENT(IN) :: BCV
      INTEGER, INTENT(IN) :: BCV_I

! Generic loop counter
      INTEGER :: LC1, OCCUPANTS
      INTEGER :: IJK, I, J, K

      OCCUPANTS = DEM_MI(BCV_I)%OCCUPANTS
      allocate(DEM_MI(BCV_I)%OWNER(OCCUPANTS))

      DEM_MI(BCV_I)%OWNER = 0

      SELECT CASE (BC_PLANE(BCV))
      CASE('N','S')

         J = DEM_MI(BCV_I)%L
         DO LC1=1,OCCUPANTS
            I = DEM_MI(BCV_I)%W(LC1)
            K = DEM_MI(BCV_I)%H(LC1)
            IF(IS_ON_myPE_owns(I,J,K)) &
               DEM_MI(BCV_I)%OWNER(LC1) = myPE
         ENDDO

      CASE('E','W')

         I = DEM_MI(BCV_I)%L
         DO LC1=1,OCCUPANTS
            J = DEM_MI(BCV_I)%W(LC1)
            K = DEM_MI(BCV_I)%H(LC1)
            IF(IS_ON_myPE_owns(I,J,K)) &
               DEM_MI(BCV_I)%OWNER(LC1) = myPE
         ENDDO

      CASE('T','B')

         K = DEM_MI(BCV_I)%L
         DO LC1=1,OCCUPANTS
            I = DEM_MI(BCV_I)%W(LC1)
            J = DEM_MI(BCV_I)%H(LC1)
            IF(IS_ON_myPE_owns(I,J,K)) &
               DEM_MI(BCV_I)%OWNER(LC1) = myPE
         ENDDO

      END SELECT

      CALL GLOBAL_ALL_SUM(DEM_MI(BCV_I)%OWNER(:))

      RETURN
      END SUBROUTINE SET_DEM_MI_OWNER
