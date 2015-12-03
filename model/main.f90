! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MFIX                                                    !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: The main module in the MFIX program                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!
!> \mainpage Multiphase Flow with Interphase eXchanges
!!
!! MFIX is a general-purpose computer code developed at the National
!! Energy Technology Laboratory, NETL, for describing the hydrodynamics,
!! heat transfer, and chemical reactions in fluid-solid systems.
!!
!! It has been used for describing bubbling and circulating fluidized
!! beds and spouted beds. MFiX calculations give transient data on the
!! three-dimensional distribution of pressure, velocity, temperature,
!! and species mass fractions. MFiX code is based on a generally
!! accepted set of multiphase flow equations. The code is used as a
!! "test-stand" for testing and developing multiphase flow constitutive
!!  equations.
!!
!! \section Notice
!! Neither the United States Government nor any agency thereof, nor any
!! of their employees, makes any warranty, expressed or implied, or
!! assumes any legal liability or responsibility for the accuracy,
!! completeness, or usefulness of any information, apparatus, product,
!! or process disclosed or represents that its use would not infringe
!! privately owned rights.
!!
!! * MFIX is provided without any user support for applications in the
!!   user's immediate organization. It should not be redistributed in
!!   whole or in part.
!!
!! * The use of MFIX is to be acknowledged in any published paper based
!!   on computations using this software by citing the MFIX theory
!!   manual. Some of the submodels are being developed by researchers
!!   outside of NETL. The use of such submodels is to be acknowledged
!!   by citing the appropriate papers of the developers of the submodels.
!!
!! * The authors would appreciate receiving any reports of bugs or other
!!   difficulties with the software, enhancements to the software, and
!!   accounts of practical applications of this software.
!!
!! \section Disclaimer
!! This report was prepared as an account of work sponsored by an agency
!! of the United States Government. Neither the United States Government
!! nor any agency thereof, nor any of their employees, makes any
!! warranty, express or implied, or assumes any legal liability or
!! responsibility for the accuracy, completeness, or usefulness of any
!! information, apparatus, product, or process disclosed, or represents
!! that its use would not infringe privately owned rights. Reference
!! herein to any specific commercial product, process, or service by
!! trade name, trademark, manufacturer, or otherwise does not
!! necessarily constitute or imply its endorsement, recommendation, or
!! favoring by the United States Government or any agency thereof. The
!! views and opinions of authors expressed herein do not necessarily
!! state or reflect those of the United States Government or any
!! agency thereof.

MODULE MAIN

  CHARACTER*512 :: out_buffer
  CHARACTER*512 :: in_buffer

  LOGICAL :: in_semaphore(1) = .false.
  LOGICAL :: out_semaphore(1) = .false.

  CONTAINS

SUBROUTINE MFIX()
!f2py threadsafe

  !-----------------------------------------------
  ! Modules
  !-----------------------------------------------
  USE MFIX_netcdf
  USE cdist
  USE compar
  USE cutcell
  USE dashboard
  USE discretelement
  USE fldvar
  USE functions
  USE funits
  USE machine
  USE mfix_pic
  USE mpi_utility
  USE output
  USE parallel_mpi
  USE param
  USE param1
  USE read_input
  USE qmom_kinetic_equation
  USE quadric
  USE run
  USE time_cpu

  USE vtk, only : WRITE_VTK_FILES

  use error_manager

  IMPLICIT NONE
  !-----------------------------------------------
  ! Local variables
  !-----------------------------------------------
  ! Final value of CPU time.
  DOUBLE PRECISION :: CPU1
  ! time used for computations.
  DOUBLE PRECISION :: CPUTIME_USED, WALLTIME_USED
  ! CPU time unit.
  CHARACTER(LEN=4) :: TUNIT
  ! Save TIME in input file for RESTART_2
  DOUBLE PRECISION :: TIME_SAVE
  ! Temporary storage for DT
  DOUBLE PRECISION :: DT_tmp
  ! loop counter
  INTEGER :: L
  ! DISTIO variable for specifying the mfix version
  CHARACTER(LEN=512) :: version
  ! environment variable
  !$      CHARACTER(LEN=512) :: omp_num_threads
  !$      INTEGER :: length
  !$      INTEGER :: status

  !$      INTEGER num_threads, threads_specified, omp_id
  !$      INTEGER omp_get_num_threads
  !$      INTEGER omp_get_thread_num

  ! ! C Function
  !       INTERFACE
  !          SUBROUTINE INIT_CMD_SOCKET(port) BIND ( C )
  !            use, INTRINSIC :: iso_c_binding
  !            CHARACTER(KIND=C_CHAR), INTENT(IN) :: port(*)
  !          END SUBROUTINE INIT_CMD_SOCKET
  !          SUBROUTINE INIT_LOG_SOCKET(port) BIND ( C )
  !            use, INTRINSIC :: iso_c_binding
  !            CHARACTER(KIND=C_CHAR), INTENT(IN) :: port(*)
  !          END SUBROUTINE INIT_LOG_SOCKET
  !       END INTERFACE

  !-----------------------------------------------

  ! DISTIO
  ! If you change the value below in this subroutine, you must also
  ! change it in write_res0.f and the value should also be consistent
  ! with the check in read_res0
  version = 'RES = 01.6'

  bDoing_postmfix = .false.

  ! Invoke MPI initialization routines and get rank info.
  CALL PARALLEL_INIT
  CALL GEN_LOG_BASENAME

  ! we want only PE_IO to write out common error messages
  DMP_LOG = (myPE == PE_IO)

  ! set the version.release of the software
  ID_VERSION = '2015-2'

  ! set automatic restart flag to false
  !      AUTOMATIC_RESTART = .FALSE.
  !      ITER_RESTART      = 1

  !       ! create libmsockets server
  ! #ifdef socket
  !       CALL INIT_CMD_SOCKET("7777"//CHAR(0))
  !       CALL INIT_LOG_SOCKET("8888"//CHAR(0))
  ! #endif

  ! specify the number of processors to be used
  !$        call get_environment_variable("OMP_NUM_THREADS",omp_num_threads,length,status, .true.)
  !$      if (status.eq.0 .and. length.ne.0) then
  !$        read(omp_num_threads,*) threads_specified
  !$      else
  !$        WRITE(*,'(A,$)') 'Enter the number of threads to be used for SMP: '
  !$        READ(*,*) threads_specified
  !$      endif

  !$      call omp_set_num_threads(threads_specified)

  ! Find the number of processors used
  !$omp  parallel
  !$      num_threads = omp_get_num_threads()
  !$      omp_id = omp_get_thread_num()
  !$      if(omp_id.eq.0) Write(*,*)' Number of threads used for SMP = ',  num_threads
  !$omp  end parallel


  ! Set machine dependent constants
  CALL MACHINE_CONS

  ! Get the date and time. They give the unique run_id in binary output
  ! files
  CALL GET_RUN_ID

  ! AEOLUS: stop trigger mechanism to terminate MFIX normally before batch
  ! queue terminates. timestep at the beginning of execution
  CALL CPU_TIME (CPU00)
  WALL0 = WALL_TIME()

  ! Read input data, check data, do computations for IC and BC locations
  ! and flows, and set geometry parameters such as X, X_E, DToDX, etc.
  CALL GET_DATA

  ! Write the initial part of the standard output file
  CALL WRITE_OUT0
  IF(.NOT.CARTESIAN_GRID)  CALL WRITE_FLAGS

  ! Write the initial part of the special output file(s)
  CALL WRITE_USR0

  !$    CALL START_LOG
  !$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' '
  !$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' Number of processors used = ', threads_specified
  !$    IF(DMP_LOG)WRITE (UNIT_LOG, *) ' '
  !$    CALL END_LOG

  !  setup for PC quickwin application
  CALL PC_QUICKWIN


  CALL INIT_ERR_MSG('MFIX')


101 CONTINUE


  ! if not netcdf writes asked for ... globally turn off netcdf
  if (MFIX_usingNETCDF()) then
     bGlobalNetcdf = .false.
     do L = 1,20
        if (bWrite_netcdf(L)) bGlobalNetcdf = .true.
     enddo
  endif


  IF(AUTOMATIC_RESTART) THEN
     RUN_TYPE = 'RESTART_1'
     AUTOMATIC_RESTART = .FALSE.
     ITER_RESTART = ITER_RESTART + 1
     CALL CHECK_INITIAL_CONDITIONS
     CALL CHECK_BOUNDARY_CONDITIONS
     CALL CHECK_INTERNAL_SURFACES
     CALL CHECK_POINT_SOURCES
     CALL CHECK_CHEMICAL_RXNS
     CALL SET_FLAGS
     CALL SET_CONSTPROP
  ENDIF

  DT_TMP = DT
  SELECT CASE (TRIM(RUN_TYPE))

  CASE ('NEW')
     ! Write the initial part of the restart files
     CALL WRITE_RES0
     DO L = 1, N_SPX
        CALL WRITE_SPX0 (L, 0)
     ENDDO

  CASE ('RESTART_1')
     ! Read the time-dependent part of the restart file
     CALL READ_RES1
     WRITE(ERR_MSG, 1010) TIME, NSTEP
     CALL FLUSH_ERR_MSG()

  CASE ('RESTART_2')
     TIME_SAVE = TIME
     ! DISTIO
     if (myPE .ne. PE_IO .and. bDist_IO .and. bStart_with_one_res) then
        write (unit_res,rec=1) version
        write (unit_res,rec=2) 4
        write (unit_res,rec=3) 4
     endif

     CALL READ_RES1
     TIME = TIME_SAVE

     WRITE(ERR_MSG, 1010) TIME, NSTEP
     CALL FLUSH_ERR_MSG()

     CALL WRITE_RES0

     ! Writing the RES1 and SPX1 can only be done here when re-indexing is turned off
     ! This will be done after the cell re-indexing is done later in this file.
     ! This allows restarting independently of the re-indexing setting between
     ! the previous and current run.
     IF(.NOT.RE_INDEXING) THEN
        CALL WRITE_RES1
        DO L = 1, N_SPX
           CALL WRITE_SPX0 (L, 0)
           CALL WRITE_SPX1 (L, 0)
        END DO
        call write_netcdf(0,0,time)
     ENDIF

  CASE DEFAULT
     CALL START_LOG
     IF(DMP_LOG)WRITE (UNIT_LOG, *) &
          ' MFIX: Do not know how to process'
     IF(DMP_LOG)WRITE (UNIT_LOG, *) ' RUN_TYPE in data file'
     CALL END_LOG
     call mfix_exit(myPE)

  END SELECT

  call MPI_Barrier(MPI_COMM_WORLD,mpierr)

  IF (DT_TMP /= UNDEFINED) THEN
     DT = MAX(DT_MIN,MIN(DT_MAX,DT))
  ELSE
     DT = DT_TMP
  ENDIF

  ! Set arrays for computing indices. A secondary call is made
  ! after cut cell-preprocessing to update array indices.
  IF(CARTESIAN_GRID) THEN
     CALL SET_INCREMENTS
     CALL SET_INCREMENTS3
  ENDIF


  !      IF(.NOT.RE_INDEXING) CALL WRITE_IJK_VALUES


  ! Set the flags for wall surfaces impermeable and identify flow
  ! boundaries using FLAG_E, FLAG_N, and FLAG_T
  CALL SET_FLAGS1

  !  Update flags for Cartesian_GRID.
  IF(CARTESIAN_GRID) CALL CHECK_BC_FLAGS

  ! Calculate cell volumes and face areas
  IF(.NOT.CARTESIAN_GRID)  THEN
     CALL SET_GEOMETRY1
     !       ELSE
     !         CALL SET_GEOMETRY
  ENDIF

  ! Find corner cells and set their face areas to zero
  IF(.NOT.CARTESIAN_GRID)  THEN
     CALL GET_CORNER_CELLS()
  ELSE
     IF (SET_CORNER_CELLS)  CALL GET_CORNER_CELLS ()
  ENDIF

  ! Set constant physical properties
  CALL SET_CONSTPROP

  ! Set initial conditions
  CALL SET_IC

  ! Set point sources.
  CALL SET_PS

  ! Set boundary conditions
  CALL ZERO_NORM_VEL
  CALL SET_BC0
  !      IF(DISCRETE_ELEMENT) CALL MAKE_ARRAYS_DES

  ! JFD: cartesian grid implementation
  IF(CARTESIAN_GRID) CALL CG_SET_BC0
  !      IF(DEM_SOLIDS) CALL SET_BC_DEM

  ! Set gas mixture molecular weight
  CALL SET_MW_MIX_G

  ! Set the pressure field for a fluidized bed
  IF (RUN_TYPE == 'NEW') CALL SET_FLUIDBED_P

  ! Initialize densities.
  IF (RUN_TYPE == 'NEW') CALL SET_RO_G
  IF (RUN_TYPE == 'NEW') CALL SET_RO_S

  ! Initialize time dependent boundary conditions
  CALL SET_BC1

  ! Check the field variable data and report errors.
  IF(.NOT.CARTESIAN_GRID)  CALL CHECK_DATA_20


  !=======================================================================
  ! JFD: START MODIFICATION FOR RE-INDEXING CELLS
  !=======================================================================
  IF(CARTESIAN_GRID.AND.RE_INDEXING) THEN

     IF(myPE == PE_IO) THEN
        WRITE(*,"(72('='))")
        WRITE(*,*)' RE-INDEXING CELLS FOR CARTESIAN GRID...'
     ENDIF
     CALL RE_INDEX_ARRAYS


     !      IF(myPE == PE_IO)print*,'Calling REPORT_BEST_IJK_SIZE:'
     !       CALL REPORT_BEST_IJK_SIZE
     CALL REPORT_BEST_PROCESSOR_SIZE
     !      IF(myPE == PE_IO)print*,'Exiting MFIX after REPORT_BEST_IJK_SIZE.'


     IF(myPE == PE_IO) WRITE(*,"(72('='))")

     ! In case of a RESTART_2, write the RES1 and SPX1 files here
     ! This was commented out earlier in this file.
     IF(RUN_TYPE == 'RESTART_2') THEN
        CALL WRITE_RES1
        DO L = 1, N_SPX
           CALL WRITE_SPX0 (L, 0)
           CALL WRITE_SPX1 (L, 0)
        END DO
        call write_netcdf(0,0,time)
     ENDIF
  ENDIF

  !=======================================================================
  ! JFD: END MODIFICATION FOR RE-INDEXING CELLS
  !=======================================================================

  ! Setup VTK data for regular (no cut cells) grid
  IF(.NOT.CARTESIAN_GRID.AND.WRITE_VTK_FILES) CALL SETUP_VTK_NO_CUTCELL

  IF(DISCRETE_ELEMENT) CALL MAKE_ARRAYS_DES
  IF(QMOMK) CALL QMOMK_MAKE_ARRAYS

  ! Set the inflow/outflow BCs for DEM solids
  IF(DEM_SOLIDS) CALL SET_BC_DEM
  ! Set the inflow/outflow BC for PIC solids
  IF(PIC_SOLIDS) CALL SET_BC_PIC

  ! Set the inital properties of each particle.
  IF(DEM_SOLIDS) CALL SET_IC_DEM

  ! AEOLUS: debug prints
  if (DBGPRN_LAYOUT .or. bdist_io) then
     !     write (*,*) myPE , ' E.4 ... version = ' , version(1:33)
     call debug_write_layout()
     call write_parallel_info()
  endif

  ! Initializations for CPU time calculations in iterate
  CPUOS = 0.
  CALL CPU_TIME (CPU1)
  CPU_NLOG = CPU1
  TIME_NLOG = TIME - DT

  ! Get the initial value of CPU time
  CALL CPU_TIME (CPU0)

  ! Find the solution of the equations from TIME to TSTOP at
  ! intervals of DT
  CALL TIME_MARCH
  IF(AUTO_RESTART.AND.AUTOMATIC_RESTART&
       .AND.ITER_RESTART.LE.10) GOTO 101

  ! Call user-defined subroutine after time-loop.
  IF (CALL_USR) CALL USR3

  ! Get the final value of CPU time.  The difference gives the
  ! CPU time used for the computations.
  CALL CPU_TIME (CPU1)

  ! Compute the CPU time and write it out in the .OUT file.
  CPUTIME_USED = CPU1 - CPU0 - CPU_IO
  WALLTIME_USED = WALL_TIME() - WALL0
  CALL WRITE_OUT3 (CPUTIME_USED, WALLTIME_USED, CPU_IO)

  ! JFD: cartesian grid implementation
  IF(WRITE_DASHBOARD) THEN
     IF(DT>=DT_MIN) THEN
        RUN_STATUS = 'Complete.'
     ELSE
        RUN_STATUS = 'DT < DT_MIN.  Recovery not possible!'
     ENDIF
     CALL GET_TUNIT(CPUTIME_USED,TUNIT)
     CALL UPDATE_DASHBOARD(0,CPUTIME_USED,TUNIT)
  ENDIF
  IF(CARTESIAN_GRID)  CALL CLOSE_CUT_CELL_FILES

  ! Finalize and terminate MPI
  call parallel_fin

  CALL FINL_ERR_MSG

  STOP

1000 FORMAT(/1X,'MFIX ',A,' Simulation:'/)

1010 FORMAT('Message 1010: Read in data from .RES file for TIME = ',&
       G12.5,/'Time step number (NSTEP) =',I7)

END SUBROUTINE MFIX

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GEN_LOG_BASENAME                                        !
!  Author: Aytekin Gel                                Date: 19-SEP-03  !
!                                                                      !
!  Purpose: Generate the file base for DMP logs.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
SUBROUTINE GEN_LOG_BASENAME

  use compar, only: myPE
  use compar, only: fbname

  implicit none

  ! Variables for generating file basename with processor id
  INTEGER :: i1, i10, i100, i1000, i10000

  ! PAR_I/O Generate file basename for LOG files
  i10000 = int(myPE/10000)
  i1000  = int((myPE-i10000*10000)/1000)
  i100   = int((myPE-i10000*10000-i1000*1000)/100)
  i10    = int((myPE-i10000*10000-i1000*1000-i100*100)/10)
  i1     = int((myPE-i10000*10000-i1000*1000-i100*100-i10*10)/1)

  i10000 = i10000 + 48
  i1000  = i1000  + 48
  i100   = i100   + 48
  i10    = i10    + 48
  i1     = i1     + 48

  fbname=char(i10000)//char(i1000)//char(i100)//char(i10)//char(i1)

  RETURN
END SUBROUTINE GEN_LOG_BASENAME




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: debug_write()                                          C
!  Purpose: Write out full geometry index setup information for the
!  case
!                                                                      C
!  Author: Aytekin Gel                                Date: 19-SEP-03  C
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

SUBROUTINE debug_write_layout()

  !-----------------------------------------------
  ! Modules
  !-----------------------------------------------
  USE param
  USE param1
  USE parallel
  USE geometry
  USE compar
  USE mpi_utility
  USE sendrecv
  USE sendrecv3
  USE indices
  USE leqsol
  USE cdist
  USE funits
  USE run
  USE time_cpu
  USE functions
  IMPLICIT NONE
  !-----------------------------------------------
  ! Local Variables
  !-----------------------------------------------
  ! phase index
  INTEGER :: M
  ! indices
  INTEGER :: i, j, k, ijk, ijk_GL, ijk_PROC, ijk_IO
  !
  integer :: indxA, indxA_gl, indxB, indxB_gl, indxC, indxC_gl
  integer :: indxD, indxD_gl, indxE, indxE_gl, indxF, indxF_gl
  integer :: indxG, indxG_gl, indxH, indxH_gl
  !
  logical :: amgdbg = .TRUE.

  character(LEN=80) :: fname

  !DISTIO
  !      fname = "layout_xxxx.txt"
  !      write (fname(8:11),'(i4.4)') myPE
  fname = "layout_xxxxx.txt"
  write (fname(8:12),'(i5.5)') myPE
  open (unit=11,file=fname,status='unknown')

  write (11,*) ' ********************************************'
  write (11,*) ' ********************************************'
  write (11,*) ' ********************************************'
  write (11,*) ' ********************************************'
  write (11,*) ' '
  write (11,*) ' '
  write (11,*) ' myPE =           ' , myPE
  write (11,*) ' '
  write (11,*) ' '


  IF (AMGDBG .OR. bDist_IO) THEN
     write(11,"('BLK1: Running from istart3,iend3 .AND. jstart3, jend3 .AND. kstart3, kend3')")
     write(11,"(' (   i ,    j,     k) =>    ijk      ijk_GL     ijk_PROC    ijk_IO')")
     write(11,"(' ====================      =====     =======    ========    ======')")
     DO k = kstart3, kend3
        DO i = istart3,iend3
           DO j = jstart3, jend3
              ijk = FUNIJK(i,j,k)
              ijk_GL = FUNIJK_GL(i,j,k)
              ijk_PROC = FUNIJK_PROC(i,j,k,myPE)
              ijk_IO = FUNIJK_IO(i,j,k)
              write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',4(I8,' , '))") &
                   i,j,k,ijk,ijk_GL,ijk_PROC,ijk_IO
           ENDDO
        ENDDO
     ENDDO

     write(11,"(/,/,'BLK2: Print out Bottom, South, West, East, North, Top neighbors')")
     write(11,"(' (   i ,    j,     k) =>    ijk    ijk_GL    B_of    S_of    W_of    E_of    N_of    T_of')")
     write(11,"(' ====================      =====   =======  ======  ======  ======  ======  ======  ======')")
     DO k = kstart3, kend3
        DO i = istart3,iend3
           DO j = jstart3, jend3
              ijk = FUNIJK(i,j,k)
              ijk_GL = FUNIJK_GL(i,j,k)
              write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',2(I7,' , '),6(I7,2X))") &
                   i,j,k,ijk,ijk_GL,bottom_of(ijk),south_of(ijk),west_of(ijk),&
                   east_of(ijk),north_of(ijk),top_of(ijk)
           ENDDO
        ENDDO
     ENDDO

     write(11,"(/,/,'BLK3: Print out km, jm, im, ip, jp, kp neighbors')")
     write(11,"(' (   i ,    j,     k) =>    ijk    ijk_GL    km_of   jm_of   im_of   ip_of   jp_of   kp_of')")
     write(11,"(' ====================      =====   =======  ======  ======  ======  ======  ======  ======')")
     DO k = kstart3, kend3
        DO i = istart3,iend3
           DO j = jstart3, jend3
              ijk = FUNIJK(i,j,k)
              ijk_GL = FUNIJK_GL(i,j,k)
              write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',2(I7,' , '),6(I7,2X))") &
                   i,j,k,ijk,ijk_GL,km_of(ijk),jm_of(ijk),im_of(ijk),&
                   ip_of(ijk),jp_of(ijk),kp_of(ijk)
           ENDDO
        ENDDO
     ENDDO

     write(11,"(/,'BLK4a: Active Fluid Cells:FLUID_AT(ijk)=.T.',/,&
          &           ' (   i ,    j,     k) =>    ijk  [   x ,     ,     z]')")
     write(11,"(' ====================      =====  ====================')")
     DO ijk = ijkstart3, ijkend3
        I = I_OF(IJK)
        J = J_OF(IJK)
        K = K_OF(IJK)

        !         IF (FLOW_AT_E(IJK)) THEN
        IF (FLUID_AT(IJK)) THEN
           !          write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',I8)") I,J,K,ijk
           write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',I8,' [',E12.5,',',E12.5,' ]')") I,J,K,ijk,X(i),Z(k)
        ENDIF
     ENDDO

     write(11,"(/,'BLK4b: Cells that are (.NOT.WALL_AT(IJK)) = .T.',/,&
          &           ' (   i ,    j,     k) =>    ijk  [   x ,     ,     z]')")
     write(11,"(' ====================      =====  ====================')")
     DO ijk = ijkstart3, ijkend3
        I = I_OF(IJK)
        J = J_OF(IJK)
        K = K_OF(IJK)

        IF (.NOT.WALL_AT(IJK)) THEN
           !          write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',I8)") I,J,K,ijk
           write(11,"(' (',I4,' , ',I4,' , ',I4,') => ',I8,' [',E12.5,',',E12.5,' ]')") I,J,K,ijk,X(i),Z(k)
        ENDIF
     ENDDO

     DO k = kstart3, kend3
        DO i = istart3,iend3
           DO j = jstart3, jend3
              ijk = FUNIJK(i,j,k)
              ijk_GL = FUNIJK_GL(i,j,k)

              if (i == istart2 .AND. j == jstart2) then
                 indxA = ijk
                 indxA_gl = ijk_GL
              endif
              if (i == istart1 .AND. j == jstart1) then
                 indxE = ijk
                 indxE_gl = ijk_GL
              endif
              if (i == istart2 .AND. j == jend2) then
                 indxB = ijk
                 indxB_gl = ijk_GL
              endif
              if (i == istart1 .AND. j == jend1) then
                 indxF = ijk
                 indxF_gl = ijk_GL
              endif
              if (i == iend1 .AND. j == jstart1) then
                 indxH = ijk
                 indxH_gl = ijk_GL
              endif
              if (i == iend2 .AND. j == jstart2) then
                 indxD = ijk
                 indxD_gl = ijk_GL
              endif
              if (i == iend1 .AND. j == jend1) then
                 indxG = ijk
                 indxG_gl = ijk_GL
              endif
              if (i == iend2 .AND. j == jend2) then
                 indxC = ijk
                 indxC_gl = ijk_GL
              endif
           ENDDO
        ENDDO
        write(11,"('BLK5:')")
        write(11,"(57('='))")
        write(11,"('k= ',I5,/,57('='))") k
        write(11,"('B= ',I5,' (',I7,')',20X,'C= ',I5,' (',I7,')',/)") indxB, indxB_gl, &
             indxC, indxC_gl
        !        write(UNIT_LOG,"(' \',34X,'/')")
        !        write(UNIT_LOG,"(2X,'\',32X,'/')")
        write(11,"(3X,'F= ',I5,' (',I7,')',12X,'G= ',I5,' (',I7,')')") indxF, indxF_gl, &
             indxG, indxG_gl
        write(11,"(4(9X,'|',29X,'|',/))")
        write(11,"(3X,'E= ',I5,' (',I7,')',12X,'H= ',I5,' (',I7,')',/)") indxE, indxE_gl, &
             indxH, indxH_gl
        !        write(UNIT_LOG,"(2X,'/',32X,'\')")
        !        write(UNIT_LOG,"('/',34X,'\')")
        write(11,"('A= ',I5,' (',I7,')',20X,'D= ',I5,' (',I7,')',/,/)") indxA, indxA_gl, &
             indxD, indxD_gl

        !        write(UNIT_LOG,"(' (',I4,' , ',I4,' , ',I4,') => ',2(I7,' , '),6(I7,2X))") &
        !                                         i,j,k,ijk,ijk_GL,bottom_of(ijk),south_of(ijk),west_of(ijk),&
        !                                        east_of(ijk),north_of(ijk),top_of(ijk)

     ENDDO

     !      write(UNIT_LOG,"(/,' (   i ,    j,     k) =>    ijk (Active Fluid)')")
     !      write(UNIT_LOG,"(' ====================      =====')")
     !       DO ijk = ijkstart3, ijkend3
     !         I = I_OF(IJK)
     !         J = J_OF(IJK)
     !         K = K_OF(IJK)

     !         IF (FLOW_AT_E(IJK)) THEN
     !         IF (FLUID_AT(IJK)) THEN
     !           write(UNIT_LOG,"(' (',I4,' , ',I4,' , ',I4,') => ',I8)") I,J,K,ijk
     !         ENDIF
     !      END DO


  endif   ! end if(amgdbg .or. bdist_io)

  M = 0
  !      CALL WRITE_AB_M (A_M, B_M, IJKMAX2, M, IER)

  IF (AMGDBG .OR. bDist_IO) THEN
     write(11,"(/,/,'BLK6: ========= ORIGINAL MFIX VARIABLES ===========')")
     write(11,"('PE ',I5,': imin1  = ',I6,3X,'imax1= ',I6,/,'PE ',I5,': jmin1  = ',I6,3X,'jmax1= ',I6)") &
          myPE,imin1,imax1,myPE,jmin1,jmax1
     write(11,"('PE ',I5,': kmin1  = ',I6,3X,'kmax1= ',I6)") myPE,kmin1,kmax1
     write(11,"('-----')")
     write(11,"('PE ',I5,': imin2  = ',I6,3X,'imax2= ',I6,/,'PE ',I5,': jmin2  = ',I6,3X,'jmax2= ',I6)") &
          myPE,imin2,imax2,myPE,jmin2,jmax2
     write(11,"('PE ',I5,': kmin2  = ',I6,3X,'kmax2= ',I6)") myPE,kmin2,kmax2
     write(11,"('----- Below xxx3 set is DMP extension ------------')")
     write(11,"('PE ',I5,': imin3  = ',I6,3X,'imax3= ',I6,/,'PE ',I5,': jmin3  = ',I6,3X,'jmax3= ',I6)") &
          myPE,imin3,imax3,myPE,jmin3,jmax3
     write(11,"('PE ',I5,': kmin3  = ',I6,3X,'kmax3= ',I6)") myPE,kmin3,kmax3
     write(11,"('----- End of Below xxx3 set is DMP extension -----')")
     !      write(11,"('PE ',I5,': ijkmax2= ',I6)") myPE,ijkmax2
     write(11,"('PE ',I5,': ijmax2 = ',I6)") myPE,ijmax2
     write(11,"('PE ',I5,': ijkmin1= ',I6,' ijkmax1= ',I12)") myPE,ijkmin1, ijkmax1
     write(11,"('PE ',I5,':          ',6X,' ijkmax2= ',I12)") myPE,ijkmax2
     write(11,"('PE ',I5,':          ',6X,' ijkmax3= ',I12)") myPE,ijkmax3
     write(11,"('PE ',I5,': ijkmin4= ',I6,' ijkmax4= ',I12)") myPE,ijkmin4, ijkmax4


     write(11,"(/,/,' ========= DMP EXTENSION VARIABLES ===========')")
     !      write(UNIT_LOG,"('PE ',I5,': ijksize  = ',I6)") myPE,ijksize
     write(11,"('PE ',I5,': ijksize3 = ',I6,3X,'ijksize3_all = ',I6)") myPE,ijksize3,ijksize3_all(myPE)
     write(11,"('PE ',I5,': ijksize4 = ',I6,3X,'ijksize4_all = ',I6)") myPE,ijksize4,ijksize4_all(myPE)
     write(11,"('PE ',I5,': ijkstart3  = ',I6,3X,'ijkend3  = ',I6)") myPE,ijkstart3, ijkend3
     write(11,"('PE ',I5,': ijkstart3_all = ',I6,3X,'ijkstart4_all = ',I6)") myPE,ijkstart3_all(myPE),ijkstart4_all(myPE)
     write(11,"('PE ',I5,': istart_all = ',I6,3X,'iend_all = ',I6,/,'PE ',I5,': jstart_all = ',I6,3X,'jend_all = ',I6)") &
          myPE,istart_all(myPE),iend_all(myPE),myPE,jstart_all(myPE),jend_all(myPE)
     write(11,"('PE ',I5,': kstart_all = ',I6,3X,'kend_all = ',I6,/,'----------------------')") &
          myPE,kstart_all(myPE),kend_all(myPE)

     write(11,"('PE ',I5,': istart1_all= ',I6,3X,'iend1_all= ',I6,/,'PE ',I5,': jstart1_all= ',I6,3X,'jend3_all= ',I6)") &
          myPE,istart1_all(myPE),iend1_all(myPE),myPE,jstart1_all(myPE),jend1_all(myPE)
     write(11,"('PE ',I5,': kstart1_all= ',I6,3X,'kend1_all= ',I6,/,'----------------------')") &
          myPE,kstart1_all(myPE),kend1_all(myPE)

     write(11,"('PE ',I5,': istart2_all= ',I6,3X,'iend2_all= ',I6,/,'PE ',I5,': jstart2_all= ',I6,3X,'jend3_all= ',I6)") &
          myPE,istart2_all(myPE),iend2_all(myPE),myPE,jstart2_all(myPE),jend2_all(myPE)
     write(11,"('PE ',I5,': kstart2_all= ',I6,3X,'kend2_all= ',I6,/,'----------------------')") &
          myPE,kstart2_all(myPE),kend2_all(myPE)

     write(11,"('PE ',I5,': istart3_all= ',I6,3X,'iend3_all= ',I6,/,'PE ',I5,': jstart3_all= ',I6,3X,'jend3_all= ',I6)") &
          myPE,istart3_all(myPE),iend3_all(myPE),myPE,jstart3_all(myPE),jend3_all(myPE)
     write(11,"('PE ',I5,': kstart3_all= ',I6,3X,'kend3_all= ',I6,/,'----------------------')") &
          myPE,kstart3_all(myPE),kend3_all(myPE)

     write(11,"('PE ',I5,': istart1= ',I6,3X,'iend1= ',I6,/,'PE ',I5,': jstart1= ',I6,3X,'jend1= ',I6)") &
          myPE,istart1,iend1,myPE,jstart1,jend1
     write(11,"('PE ',I5,': kstart1= ',I6,3X,'kend1= ',I6,/,'----------------------')") &
          myPE,kstart1,kend1
     write(11,"('PE ',I5,': istart2= ',I6,3X,'iend2= ',I6,/,'PE ',I5,': jstart2= ',I6,3X,'jend2= ',I6)") &
          myPE,istart2,iend2,myPE,jstart2,jend2
     write(11,"('PE ',I5,': kstart2= ',I6,3X,'kend2= ',I6,/,'----------------------')") &
          myPE,kstart2,kend2
     write(11,"('PE ',I5,': istart3= ',I6,3X,'iend3= ',I6,/,'PE ',I5,': jstart3= ',I6,3X,'jend3= ',I6)") &
          myPE,istart3,iend3,myPE,jstart3,jend3
     write(11,"('PE ',I5,': kstart3= ',I6,3X,'kend3= ',I6,/,'----------------------')") &
          myPE,kstart3,kend3

  ENDIF   ! end if(amgdbg .or. bdist_io)

  close(unit=11)


  RETURN
END SUBROUTINE DEBUG_WRITE_LAYOUT



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

SUBROUTINE write_parallel_info()

  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE param
  USE param1
  USE parallel
  USE geometry
  USE compar
  USE mpi_utility
  USE sendrecv
  USE sendrecv3
  USE indices
  USE leqsol
  USE funits
  USE run
  USE time_cpu
  USE functions
  IMPLICIT NONE
  !-----------------------------------------------
  ! Dummy arguments
  !-----------------------------------------------
  ! Local Variables
  !-----------------------------------------------
  ! phase index
  INTEGER :: M
  ! indices
  INTEGER :: i, j, k, ijk, ijk_GL, ijk_PROC, ijk_IO
  !
  character(LEN=80) :: fname
  !-----------------------------------------------

  !DISTIO
  !      fname = "p_info_xxxx.txt"
  !      write (fname(8:11),'(i4.4)') myPE
  fname = "p_info_xxxxx.txt"
  write (fname(8:12),'(i5.5)') myPE
  open (unit=11,file=fname,status='unknown')

  write (11,*) myPe , ' = myPE'

  write (11,*) myPE , istart3,iend3
  write (11,*) myPE , jstart3,jend3
  write (11,*) myPE , kstart3,kend3

  write(11,"('BLK1: Running from istart3,iend3 .AND. jstart3, jend3 .AND. kstart3, kend3')")
  write(11,"(' (   i ,    j,     k)       ijk      ijk_GL     ijk_PROC    ijk_IO')")
  write(11,"(' ====================      =====     =======    ========    ======')")
  DO k = kstart3, kend3
     DO i = istart3,iend3
        DO j = jstart3, jend3
           ijk = FUNIJK(i,j,k)
           ijk_GL = FUNIJK_GL(i,j,k)
           ijk_PROC = FUNIJK_PROC(i,j,k,myPE)
           ijk_IO = FUNIJK_IO(i,j,k)
           write(11,"('  ',I4,'   ',I4,'   ',I4,'     ',4(I8,'   '))" ) &
                i,j,k,ijk,ijk_GL,ijk_PROC,ijk_IO
        ENDDO
     ENDDO
  ENDDO

  M = 0
  !      CALL WRITE_AB_M (A_M, B_M, IJKMAX2, M, IER)

  close(unit=11)

  RETURN
END SUBROUTINE write_parallel_info

END MODULE MAIN
