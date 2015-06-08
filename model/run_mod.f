!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: run                                                    C
!  Purpose: Common block containing run control data                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE run

! Modules
!---------------------------------------------------------------------//
      Use param, only: dim_M, dim_eqs
!---------------------------------------------------------------------//


! Main filename to be used for output files  Name must
! still be legal after extensions are added to it.
      CHARACTER(LEN=60) :: RUN_NAME

! Brief description of the problem.
      CHARACTER(LEN=60) :: DESCRIPTION

! Units for data input and output: CGS.
      CHARACTER(LEN=16) :: UNITS

! Type of run: NEW, RESTART
      CHARACTER(LEN=16) :: RUN_TYPE

! Variable which triggers automatic restart
      LOGICAL :: AUTOMATIC_RESTART

! counter to keep track of how many auto_retart were performed
      INTEGER :: ITER_RESTART

! version.release of software
      CHARACTER(LEN=10) :: ID_VERSION

! Start-time of the run.
      DOUBLE PRECISION :: TIME

! Stop-time of the run.
      DOUBLE PRECISION :: TSTOP

! Time step.
      DOUBLE PRECISION :: DT

! 1./Time step.
      DOUBLE PRECISION :: oDT

! Number of times steps completed.
      INTEGER :: NSTEP

! Declare a new variable to use on CN with RESTART cases
! Number of time steps when restart file was read
      INTEGER :: NSTEPRST

! Discretization scheme for different equations
      INTEGER :: DISCRETIZE(DIM_EQS)

! Use Chi scheme for discretizing certain equation sets
!  (species mass fractions)
      LOGICAL :: Chi_scheme

! If .TRUE. solve X momentum equations
      LOGICAL :: MOMENTUM_X_EQ(0:DIM_M)

! If .TRUE. solve Y momentum equations
      LOGICAL :: MOMENTUM_Y_EQ(0:DIM_M)

! If .TRUE. solve Z momentum equations
      LOGICAL :: MOMENTUM_Z_EQ(0:DIM_M)

! IF .TRUE. use Jackson form momentum equations
      LOGICAL :: JACKSON
! IF .TRUE. use Ishii form momentum equations
      LOGICAL :: ISHII

! If .TRUE. use Model-B momentum equations
      LOGICAL :: Model_B

! If .TRUE. include added (virtual) mass in momentum eq.
      LOGICAL :: Added_Mass

! phase number where added mass is applied.
      INTEGER :: M_AM

! If .TRUE. solve K_Epsilon turbulence eq.
      LOGICAL :: K_Epsilon

! If .TRUE. solve energy equations
      LOGICAL :: ENERGY_EQ

! If .TRUE. use the deferred correction method
      LOGICAL :: DEF_COR

! If .TRUE. use the fourth order interpolation
      LOGICAL :: FPFOI

! If .TRUE. activate 2nd order accurate time implementation
      LOGICAL :: CN_ON

! If .TRUE. solve granular energy equations
      LOGICAL :: GRANULAR_ENERGY

! If .TRUE. solve species balance equations
      LOGICAL :: SPECIES_EQ(0:DIM_M)

! If .TRUE. one of the species equations is being solved
      LOGICAL :: ANY_SPECIES_EQ

! If .TRUE. call user-defined subroutines
      LOGICAL :: CALL_USR

! If .TRUE. solve population balance  equations
      LOGICAL :: Call_DQMOM

! Drag model options (see drag_gs for full details)
! default is syam_obrien (may enforce a corrected Umf by defining
! drag_c1 and drag_d1 accordingly)
      CHARACTER(64) :: DRAG_TYPE
      INTEGER :: DRAG_TYPE_ENUM
      INTEGER,PARAMETER :: SYAM_OBRIEN=0
      INTEGER,PARAMETER :: GIDASPOW=1
      INTEGER,PARAMETER :: GIDASPOW_PCF=2
      INTEGER,PARAMETER :: GIDASPOW_BLEND=3
      INTEGER,PARAMETER :: GIDASPOW_BLEND_PCF=4
      INTEGER,PARAMETER :: WEN_YU=5
      INTEGER,PARAMETER :: WEN_YU_PCF=6
      INTEGER,PARAMETER :: KOCH_HILL=7
      INTEGER,PARAMETER :: KOCH_HILL_PCF=8
      INTEGER,PARAMETER :: BVK=9
      INTEGER,PARAMETER :: HYS=10
      INTEGER,PARAMETER :: USER_DRAG=11

! filtered/subgrid corrections to the drag coefficient & granular
! stress terms including granular viscosity and solids pressure
! current options are 'igci' and 'milioli'
      CHARACTER(64) :: SUBGRID_TYPE

      INTEGER :: SUBGRID_TYPE_ENUM
      INTEGER,PARAMETER :: UNDEFINED_SUBGRID_TYPE=0
      INTEGER,PARAMETER :: IGCI=1
      INTEGER,PARAMETER :: MILIOLI=2

! If .TRUE. incorporate the wall effects upon the calculation of the
! subgrid solids viscosity, solids pressure, and gas-solids drag
      LOGICAL :: SUBGRID_Wall
! the ratio of the FilterSize to the GridSize
      DOUBLE PRECISION :: filter_size_ratio

! Single particle drag correlation
      CHARACTER(64) :: CD_FUNCTION

! Parameter used to calculate lubrication interactions between
! different particles in HYS drag model
      DOUBLE PRECISION :: LAM_HYS

! Kinetic theory model options (see calc_mu_s for details)
! for m > 1 : IA_nonep, GHD, LUN_1984
! for m = 1 : LUN_1984, simonin, ahmadi, or
!             GD_99 for granular flow or GTSH for gas-solids flow
      CHARACTER(64) :: KT_TYPE
      INTEGER :: KT_TYPE_ENUM
      INTEGER,PARAMETER :: LUN_1984=0
      INTEGER,PARAMETER :: SIMONIN_1996=1
      INTEGER,PARAMETER :: AHMADI_1995=2
      INTEGER,PARAMETER :: GD_1999=3
      INTEGER,PARAMETER :: GTSH_2012=4
      INTEGER,PARAMETER :: IA_2005=5
      INTEGER,PARAMETER :: GHD_2007=6

! If .TRUE. use Simonin model (k_epsilon must also be true)
      LOGICAL :: SIMONIN

! If .TRUE. use Ahmadi model (k_epsilon must also be true)
      LOGICAL :: AHMADI

! If .TRUE. calculate frictional stress terms
      LOGICAL :: FRICTION
! Form of friction model:
!             If 0: use S:S
!             If 1: use the form of Savage to compute S:S
!             If 2: use combination of both for frictional stress terms
      INTEGER :: SAVAGE

! If .TRUE. use Scheffer frictional stress (default set to .TRUE.)
      LOGICAL :: SCHAEFFER

! If .TRUE. use blending frictional/kinetic stresses
! (default set to .FALSE. do not blend)
      LOGICAL :: BLENDING_STRESS
      LOGICAL :: TANH_BLEND ! default set to true
      LOGICAL :: SIGM_BLEND ! default set to false

! If .TRUE. use Jenkins small friction BC
      LOGICAL :: JENKINS
! If .TRUE. use revised phip for JJ BC
      LOGICAL :: BC_JJ_M
! If .TRUE. output PHIP to JJ_PHIP.dat
      LOGICAL :: PHIP_OUT_JJ
! to write specularity
      INTEGER :: PHIP_OUT_ITER

! If .TRUE. treat system as if shearing
      LOGICAL :: SHEAR
! Shear Vel
      DOUBLE PRECISION :: V_sh

! Radial distribution function options (see g_0 for details)
! for m > 1 options are lebowitz, modified_lebowitz,
! mansoori, modified_mansoori.  default = lebowitz
! for m = 1 then carnahan and starling rdf used
      CHARACTER(64) :: RDF_TYPE
      INTEGER :: RDF_TYPE_ENUM
      INTEGER, PARAMETER :: LEBOWITZ=0
      INTEGER, PARAMETER :: MODIFIED_LEBOWITZ=1
      INTEGER, PARAMETER :: MANSOORI=2
      INTEGER, PARAMETER :: MODIFIED_MANSOORI=3
      INTEGER, PARAMETER :: CARNAHAN_STARLING=4

! If .TRUE. use Yu and Standish correlation to compute ep_star
      LOGICAL :: YU_STANDISH

! If .TRUE. use Fedors and Landel correlation to compute ep_star
      LOGICAL :: FEDORS_LANDEL

! STOP Trigger mechanism to terminate MFIX normally before batch
! queue terminates flag variable to check for end of batch queue when
! set to TRUE check performed at the beginning of each time step and
! termination of mfix triggered after saving all files if condition
! is met
      LOGICAL :: CHK_BATCHQ_END
! variable to store the total wall clock duration of the batch queue
! session wall clock time specified in seconds
! for jaguarcnl@NCCS max wall clock limit is 2.5 hr limit up to 512
! processors
      DOUBLE PRECISION :: BATCH_WALLCLOCK
! variable to set a buffer time before the batch queue session ends to
! make sure once MFIX is triggered to shutdown, there is sufficient
! time to save files, make copies to HPSS storage before batch queue
! time runs out. Current logic in MFIX checks for:
!    if CPU_TIME > (BATCH_WALLCLOCK - TERM_BUFFER) then
!    save all .RES .SP files and trigger shutdown
      DOUBLE PRECISION :: TERM_BUFFER

! If .TRUE. code will automatically restart for DT < DT_MIN
      LOGICAL :: AUTO_RESTART

! parameters for dynamically adjusting time step
! +1 -> increase dt; -1 decrease dt
      INTEGER :: DT_dir

! number of steps since last adjustment
      INTEGER :: STEPS_tot

! number of iterations since last adjustment
      INTEGER :: NIT_tot

! iterations per second for last dt
      DOUBLE PRECISION :: NITos

! Maximum Time step.
      DOUBLE PRECISION :: DT_MAX

! Minimum Time step.
      DOUBLE PRECISION :: DT_MIN

! Time step adjustment factor (<1.0)
      DOUBLE PRECISION :: DT_FAC

! The previous time step used in iterate (before it is
! changed by adjust_dt)
      DOUBLE PRECISION :: DT_prev

! in case iterations converged and DT modified, use old dt
! to advance time in time_march.
      LOGICAL :: use_DT_prev

! Slope limiter parameter (0 < C _FAC <= 1.0)
      DOUBLE PRECISION :: C_FAC

! If .TRUE. reduce time step when residuals do not decrease
      LOGICAL :: DETECT_STALL

! String which controls reduction of global sums for residual
! calculations
      LOGICAL :: DEBUG_RESID

! Generate log files when negative gas density is detected.
      LOGICAL :: REPORT_NEG_DENSITY

       common /run_dp/ time      !for Linux

! Flag to invoke the variable solids density model.
!      LOGICAL          SOLID_RO_V


! Flags indicating variable solids density.
      LOGICAL :: SOLVE_ROs(DIM_M), ANY_SOLVE_ROs

! Specifies the type of solids: TFM, DEM, MPPIC
      CHARACTER(len=3), DIMENSION(DIM_M) :: SOLIDS_MODEL

! Flags for various solids phase models.
      LOGICAL :: TFM_SOLIDS
      LOGICAL :: DEM_SOLIDS
      LOGICAL :: PIC_SOLIDS
! The number of the various solids phases.
      INTEGER :: TFM_COUNT = 0
      INTEGER :: DEM_COUNT = 0
      INTEGER :: PIC_COUNT = 0

    CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function name: ADJUST_DT(IER, NIT)                                  C
!  Purpose: Automatically adjust time step                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: FEB-10-97  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:
!                                                                      C
!  Local variables: NONE                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      LOGICAL FUNCTION ADJUST_DT (IER, NIT)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE machine, only: start_log, end_log
      USE mpi_utility
      USE output
      USE param
      USE param1

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER IER, NIT
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: STEPS_MIN = 5
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      DOUBLE PRECISION :: NITOS_NEW
      CHARACTER(LEN=70), DIMENSION(1) :: LINE
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
!//  only do this routine if the root processor
!//SP
!     if (myPE.ne.PE_IO) return

      ADJUST_DT = .FALSE.                     !No need to iterate again
      IF (DT==UNDEFINED .OR. DT<ZERO) RETURN

      use_DT_prev = .FALSE.  ! when false use current dt to advance time.


!     Initialize
      IF (IER == 100) THEN
         DT_DIR = -1
         STEPS_TOT = 0
         NITOS = 0.
         NIT_TOT = 0

! Iterate successfully converged.
!---------------------------------------------------------------------//
      ELSE IF(IER == 0) THEN

!AE TIME 041601: Set back the timestep to original size which was
! halved previously for 2nd order accurate time implementation.
         IF((CN_ON.AND.NSTEP>1 .AND. RUN_TYPE == 'NEW') .OR. &
            (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) &
            DT = 2.D0*DT

         IF(STEPS_TOT >= STEPS_MIN) THEN
            NITOS_NEW = DBLE(NIT_TOT)/(STEPS_TOT*DT)
            IF (NITOS_NEW > NITOS) DT_DIR = DT_DIR*(-1)
            STEPS_TOT = 0
            NITOS = NITOS_NEW
            NIT_TOT = 0
            IF (DT_DIR > 0) THEN
               DT = MIN(DT_MAX,DT/DT_FAC)
            ELSE
               DT = DT*DT_FAC
            ENDIF
! in case IER = 0 and DT is modified, use DT_prev to advance time
           use_DT_prev = .TRUE.
            IF (FULL_LOG.and.myPE.eq.PE_IO) &
               WRITE (*, *) 'DT = ', DT, '  NIT/s = ', NINT(NITOS)
         ELSE
            STEPS_TOT = STEPS_TOT + 1
            NIT_TOT = NIT_TOT + NIT
         ENDIF
! No need to iterate again
         ADJUST_DT = .FALSE.

!AE TIME 041601: Cut the timestep into half for 2nd order accurate
! time implementation.
         IF((CN_ON.AND.NSTEP>1 .AND. RUN_TYPE == 'NEW') .OR. &
            (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) &
            DT = 0.5D0*DT

! Iterate failed to converge.
!---------------------------------------------------------------------//
      ELSE
!AE TIME 041601: Set back the timestep to original size which was
! halved previously for 2nd order accurate time implementation
         IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. &
            (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) &
            DT = 2.D0*DT

! The step size has decreased to the minimum.
         IF (DT < DT_MIN) THEN
            IF(AUTO_RESTART) THEN
               LINE(1) = 'DT < DT_MIN.  Recovery not possible! Trying Automatic Restart'
              ELSE
                 LINE(1) = 'DT < DT_MIN.  Recovery not possible!'
              ENDIF
            IF (FULL_LOG.and.myPE.eq.PE_IO) WRITE (*, *) LINE(1)
            CALL WRITE_ERROR ('ADJUST_DT', LINE, 1)
!           CALL MFIX_EXIT(myPE)

! The time step cannot be reduced the step size is fixed.
         ELSE IF (DT_FAC >= ONE) THEN
            LINE(1) = 'DT_FAC >= 1.  Recovery not possible!'
            IF (FULL_LOG.and.myPE.eq.PE_IO) WRITE (*, *) LINE(1)
            CALL WRITE_ERROR ('ADJUST_DT', LINE, 1)
            CALL MFIX_EXIT(myPE)
         ELSE

! Reduce the step size.
            DT_DIR = -1
            STEPS_TOT = 0
            NITOS = 0.
            NIT_TOT = 0
            DT = DT*DT_FAC

            IF (FULL_LOG) THEN
               IF(myPE.eq.PE_IO) then
                  WRITE (*, '(12X,A,G12.5,9X,A)') &
                     ' Dt=', DT, ' Recovered            :-)'
               ENDIF
               CALL START_LOG
               IF(DMP_LOG)WRITE (UNIT_LOG, '(12X,A,G12.5,9X,A)') &
                  ' Dt=', DT, ' Recovered            :-)'
               CALL END_LOG
            ENDIF
            CALL RESET_NEW
! Iterate again with new dt
            ADJUST_DT = .TRUE.


!AE TIME 041601: Cut the timestep into half for 2nd order accurate
! time implementation.
            IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. &
               (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP>=(NSTEPRST+1)))&
                 DT = 0.5D0*DT
         ENDIF

! Get courant numbers and max p_star for the diverged iteration.
!        call get_stats(IER)

      ENDIF

!AIKEDEBUG Get the min. Courant number displayed.
!      call get_stats(IER)

      ODT = ONE/DT

      call bcast (dt,root)
      call bcast (odt,root)
!
      RETURN
      END FUNCTION ADJUST_DT

      END MODULE RUN
