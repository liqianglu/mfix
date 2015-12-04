! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: TIME_MARCH                                              !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Controlling module for time marching and finding the       !
!           solution of equations from TIME to TSTOP at intervals of   !
!           DT, updating the b.c.'s, and creating output.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE TIME_MARCH

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE MFIX_netcdf
      USE bc
      USE main
      USE cdist
      USE check
      USE coeff
      USE compar
      USE cont
      USE cutcell
      USE dashboard
      USE discretelement
      USE drag
      USE fldvar
      USE funits
      USE geometry
      USE indices
      USE leqsol, only: SOLVER_STATISTICS, REPORT_SOLVER_STATS
      USE output, only: RES_DT
      USE param
      USE param1
      USE pgcor
      USE physprop
      USE pscor
      USE qmom_kinetic_equation
      USE run
      USE scalars
      USE stiff_chem, only : STIFF_CHEMISTRY, STIFF_CHEM_SOLVER
      USE tau_g
      USE tau_s
      USE time_cpu
      USE toleranc
      USE utilities, ONLY: MAX_VEL_INLET
      USE visc_g
      USE visc_s
      USE vshear
      USE vtk
      USE vtp
      use adjust_dt
      use interactive, only: INTERACT, INIT_INTERACTIVE_MODE
      use output, only: RES_DT, NLOG

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Flag to indicate one pass through iterate for steady
! state conditions.
      LOGICAL :: FINISH

! Loop indices
      INTEGER :: M
! Error index
      INTEGER :: IER
! Number of iterations
      INTEGER :: NIT, NIT_TOTAL
! used for activating check_data_30
      INTEGER :: NCHECK, DNCHECK
! dummy logical variable for initializing adjust_dt
      LOGICAL :: dummy

! Flag to save results and cleanly exit.
      LOGICAL :: EXIT_SIGNAL = .FALSE.

! C Function
      INTERFACE
         SUBROUTINE CHECK_SOCKETS() BIND ( C )
           use, INTRINSIC :: iso_c_binding
         END SUBROUTINE CHECK_SOCKETS
      END INTERFACE

!-----------------------------------------------
! External functions
!-----------------------------------------------
! use function vavg_v_g to catch NaN's
      DOUBLE PRECISION, EXTERNAL :: VAVG_U_G, VAVG_V_G, VAVG_W_G, X_vavg

!-----------------------------------------------

      IF(AUTOMATIC_RESTART) RETURN

      FINISH  = .FALSE.
      NCHECK  = NSTEP
      DNCHECK = 1
      CPU_IO  = ZERO
      NIT_TOTAL = 0

      CALL INIT_OUTPUT_VARS

      IF(INTERACTIVE_MODE) CALL INIT_INTERACTIVE_MODE

! Parse residual strings
      CALL PARSE_RESID_STRING ()

! Call user-defined subroutine to set constants, check data, etc.
      IF (CALL_USR) CALL USR0

      CALL RRATES_INIT()

! Calculate all the coefficients once before entering the time loop
      CALL INIT_COEFF(IER)

      DO M=1, MMAX
         F_gs(1,M) = ZERO
      ENDDO

! Remove undefined values at wall cells for scalars
      CALL UNDEF_2_0 (ROP_G)
      DO M = 1, MMAX
         CALL UNDEF_2_0 (ROP_S(1,M))
      ENDDO

! Initialize d's and e's to zero
      DO M = 0, MMAX
         D_E(1,M) = ZERO
         D_N(1,M) = ZERO
         D_T(1,M) = ZERO
      ENDDO
      E_E(:) = ZERO
      E_N(:) = ZERO
      E_T(:) = ZERO

! calculate shear velocities if periodic shear BCs are used
      IF(SHEAR) CALL CAL_D(V_sh)

! Initialize check_mass_balance.  This routine is not active by default.
! Specify a reporting interval (hard-wired in the routine) to activate
! the routine.
      Call check_mass_balance (0)

! sof modification: now it's only needed to do this once before time-loop
! Mark the phase whose continuity will be solved and used to correct
! void/volume fraction in calc_vol_fr (see subroutine for details)
      CALL MARK_PHASE_4_COR (PHASE_4_P_G, PHASE_4_P_S, DO_CONT, MCP,&
          DO_P_S, SWITCH_4_P_G, SWITCH_4_P_S)

! The TIME loop begins here.............................................
 100  CONTINUE

      ! wait for request to finish
      do while(1.eq.request_pending)
         mfix_waiting = 1
         ! without a write, this loop gets optimized away
         write(*,fmt='(a)',advance='no') ''
      enddo
      mfix_waiting = 0

#ifdef socket
      IF(INTERACTIVE_MODE) THEN
         write(*,*) "INTERACTIVE_MODE: ",INTERACTIVE_MODE
         CALL CHECK_SOCKETS()
         CALL INTERACT(EXIT_SIGNAL, ABORT_SIGNAL)
         IF(ABORT_SIGNAL) RETURN
      ENDIF
#endif

! Terminate MFIX normally before batch queue terminates.
      IF (CHK_BATCHQ_END) CALL CHECK_BATCH_QUEUE_END(EXIT_SIGNAL)

      IF (CALL_USR) CALL USR1

! Remove solids from cells containing very small quantities of solids
      IF(.NOT.(DISCRETE_ELEMENT .OR. QMOMK) .OR. &
         DES_CONTINUUM_HYBRID) THEN
         IF(KT_TYPE_ENUM == GHD_2007) THEN
            CALL ADJUST_EPS_GHD
         ELSE
            CALL ADJUST_EPS
         ENDIF
      ENDIF

! sof modification: uncomment code below and modify MARK_PHASE_4_COR to
! use previous MFIX algorithm. Nov 22 2010.
! Mark the phase whose continuity will be solved and used to correct
! void/volume fraction in calc_vol_fr (see subroutine for details)
!      CALL MARK_PHASE_4_COR (PHASE_4_P_G, PHASE_4_P_S, DO_CONT, MCP,&
!           DO_P_S, SWITCH_4_P_G, SWITCH_4_P_S, IER)

! Set wall boundary conditions and transient flow b.c.'s
      CALL SET_BC1

! Uncoupled discrete element simulations
      IF(DISCRETE_ELEMENT .AND. .NOT.DES_CONTINUUM_COUPLED) THEN
         IF (DEM_SOLIDS) CALL DES_TIME_MARCH
         IF (PIC_SOLIDS) CALL PIC_TIME_MARCH
         RETURN
      ENDIF

      CALL OUTPUT_MANAGER(EXIT_SIGNAL, FINISH)

      IF (DT == UNDEFINED) THEN
         IF (FINISH) THEN
            RETURN
         ELSE
            FINISH = .TRUE.
         ENDIF

! Mechanism to terminate MFIX normally before batch queue terminates.
      ELSEIF (TIME + 0.1d0*DT >= TSTOP .OR. EXIT_SIGNAL) THEN
         IF(SOLVER_STATISTICS) &
            CALL REPORT_SOLVER_STATS(NIT_TOTAL, NSTEP)
         RETURN
      ENDIF

! Update previous-time-step values of field variables
      CALL UPDATE_OLD

! Calculate coefficients
      CALL CALC_COEFF_ALL (0, IER)

! Calculate the stress tensor trace and cross terms for all phases.
      CALL CALC_TRD_AND_TAU()

! Calculate additional solid phase momentum source terms
! that arise from kinetic theory constitutive relations
      IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
         CALL CALC_KTMOMSOURCE_U_S ()
         CALL CALC_KTMOMSOURCE_V_S ()
         CALL CALC_KTMOMSOURCE_W_S ()
      ENDIF

! Check rates and sums of mass fractions every NLOG time steps
      IF (NSTEP == NCHECK) THEN
         IF (DNCHECK < 256) DNCHECK = DNCHECK*2
         NCHECK = NCHECK + DNCHECK
! Upate the reaction rates for checking
         CALL CALC_RRATE(IER)
         CALL CHECK_DATA_30
      ENDIF

! Double the timestep for 2nd order accurate time implementation
      IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. &
          (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
         DT = 0.5d0*DT
         ODT = ODT * 2.0d0
      ENDIF

! Check for maximum velocity at inlet to avoid convergence problems
      MAX_INLET_VEL = 100.0d0*MAX_VEL_INLET()
! if no inlet velocity is specified, use an upper limit defined in
! toleranc_mod.f
      IF(MAX_INLET_VEL <= SMALL_NUMBER) THEN
         MAX_INLET_VEL = MAX_ALLOWED_VEL
         IF (UNITS == 'SI') MAX_INLET_VEL = 1D-2 * MAX_ALLOWED_VEL
      ENDIF
! Scale the value using a user defined scale factor
      MAX_INLET_VEL = MAX_INLET_VEL * MAX_INLET_VEL_FAC

! Advance the solution in time by iteratively solving the equations
 150  CALL ITERATE (IER, NIT)

      IF(AUTOMATIC_RESTART) RETURN

! Just to Check for NaN's, Uncomment the following lines and also lines
! of code in  VAVG_U_G, VAVG_V_G, VAVG_W_G to use.
!      X_vavg = VAVG_U_G ()
!      X_vavg = VAVG_V_G ()
!      X_vavg = VAVG_W_G ()
!      IF(AUTOMATIC_RESTART) RETURN

      DO WHILE (ADJUSTDT(IER,NIT))
         CALL ITERATE (IER, NIT)
      ENDDO

! A single to interupt the time step was sent for interactive mode.
      IF(INTERUPT) GOTO 100

      IF(DT < DT_MIN) THEN
         AUTO_RESTART = .FALSE.
         IF(AUTO_RESTART) THEN
            IF(TIME <= RES_DT) THEN

 1000 FORMAT('Automatic restart not possible as Total Time < RES_DT')

               WRITE(ERR_MSG,1000)
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
               AUTOMATIC_RESTART = .TRUE.
            ENDIF
            RETURN

         ELSE

 1100 FORMAT('DT < DT_MIN.  Recovery not possible!')

            WRITE(ERR_MSG,1100)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

            IF(WRITE_DASHBOARD) THEN
               RUN_STATUS = 'DT < DT_MIN.  Recovery not possible!'
               CALL UPDATE_DASHBOARD(NIT,0.0d0,'    ')
            ENDIF
            CALL MFIX_EXIT(MyPE)
         ENDIF
      ENDIF


! Stiff Chemistry Solver.
      IF(STIFF_CHEMISTRY) THEN
         CALL STIFF_CHEM_SOLVER(DT, IER)
         IF(IER /= 0) THEN
            dummy = ADJUSTDT(IER, NIT)
            GOTO 150
         ENDIF
      ENDIF

! Check over mass and elemental balances.  This routine is not active by default.
! Edit the routine and specify a reporting interval to activate it.
      CALL CHECK_MASS_BALANCE (1)

! Other solids model implementations
      IF (DEM_SOLIDS) CALL DES_TIME_MARCH
      IF (PIC_SOLIDS) CALL PIC_TIME_MARCH
      IF (QMOMK) CALL QMOMK_TIME_MARCH
      IF (CALL_DQMOM) CALL USR_DQMOM

! Advance the time step and continue
      IF((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. &
         (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
! Double the timestep for 2nd order accurate time implementation
         DT = 2.d0*DT
         ODT = ODT * 0.5d0
! Perform the explicit extrapolation for CN implementation
         CALL CN_EXTRAPOL
      ENDIF

      IF (DT /= UNDEFINED) THEN
         IF(USE_DT_PREV) THEN
            TIME = TIME + DT_PREV
         ELSE
            TIME = TIME + DT
         ENDIF
         USE_DT_PREV = .FALSE.
         NSTEP = NSTEP + 1
      ENDIF

      NIT_TOTAL = NIT_TOTAL+NIT

! write (*,"('Compute the Courant number')")
! call get_stats(IER)

      FLUSH (6)
! The TIME loop ends here....................................................
      GOTO 100

      IF(SOLVER_STATISTICS) CALL REPORT_SOLVER_STATS(NIT_TOTAL, NSTEP)

      END SUBROUTINE TIME_MARCH
