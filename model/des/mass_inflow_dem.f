!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_MASS_INLET                                          !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose:  This routine fills in the necessary information for new   !
!  particles entereing the system.                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MASS_INFLOW_DEM

      use discretelement

      use bc
      use des_bc

      implicit none

      INTEGER :: IP, LS, M, NP, IJK, LC
      INTEGER :: BCV, BCV_I
      INTEGER :: IERR
      LOGICAL :: CHECK_FOR_ERRORS, OWNS

! I/J/K index of fluid cell containing the new particle.
      INTEGER :: IJKP(3)

      DOUBLE PRECISION :: DIST, POS(3)

      CHECK_FOR_ERRORS = .FALSE.

      IERR = 0
      DO BCV_I = 1, DEM_BCMI
         BCV = DEM_BCMI_MAP(BCV_I)


         DO LC=DEM_BCMI_IJKSTART(BCV_I), DEM_BCMI_IJKEND(BCV_I)
           IJK = DEM_BCMI_IJK(LC)

            DO LS= 1,PINC(IJK)
               NP = PIC(IJK)%p(LS)
               IF(PEA(NP,3)) CYCLE
               SELECT CASE (BC_PLANE(BCV))
               CASE('N'); DIST = DES_POS_NEW(2,NP) - YN(BC_J_s(BCV))
               CASE('S'); DIST = YN(BC_J_s(BCV)-1) - DES_POS_NEW(2,NP)
               CASE('E'); DIST = DES_POS_NEW(1,NP) - XE(BC_I_w(BCV))
               CASE('W'); DIST = XE(BC_I_w(BCV)-1) - DES_POS_NEW(1,NP)
               CASE('T'); DIST = DES_POS_NEW(3,NP) - ZT(BC_K_b(BCV))
               CASE('B'); DIST = ZT(BC_K_b(BCV)-1) - DES_POS_NEW(3,NP)
               END SELECT
! The particle is still inside the domain
               IF(DIST > DES_RADIUS(NP)) PEA(NP,2) = .FALSE.
            ENDDO
         ENDDO


! Check if any particles need seeded this time step.
         IF(DEM_MI_TIME(BCV_I) > S_TIME) CYCLE

         LS = 1

! Loop over the particles being injected
         PLoop: DO IP = 1, PI_COUNT(BCV_I)

! Increment the global max particle ID (all ranks).
            iMAX_GLOBAL_ID = iMAX_GLOBAL_ID + 1

! Determine the location and phase of the incoming particle.
            CALL SEED_NEW_PARTICLE(BCV, BCV_I, M, POS, IJKP, OWNS)

! Only the rank receiving the new particle needs to continue
            IF(.NOT.OWNS) CYCLE PLoop

! Increment the number of particle on the processor by one. If the max
! number of particles is exceeded, set the error flag and cycle.
            PIP = PIP + 1
            IF(PIP  >= MAX_PIP) IERR = 1
            IF(IERR /= 0) CYCLE PLoop

! Find the first free space in the particle existance array.
            NP_LP: DO NP = LS, MAX_PIP
               IF(.NOT.PEA(NP,1)) THEN
                  LS = NP
                  EXIT NP_LP
               ENDIF
            ENDDO NP_LP

! Set the particle's global ID.
            iGLOBAL_ID(NP) = iMAX_GLOBAL_ID

! Set the properties of the new particle.
            CALL SET_NEW_PARTICLE_PROPS(BCV, M, NP, POS, IJKP)

         ENDDO PLoop

! Update the time for seeding the next particle.
         DEM_MI_TIME(BCV_I) = S_TIME + PI_FACTOR(BCV_I)*DTSOLID
! Set the flag for error checking.
         CHECK_FOR_ERRORS = .TRUE.
      ENDDO


     IF(CHECK_FOR_ERRORS) THEN
     ENDIF

 1000 FORMAT(/1X,70('*')//,' From: DES_MASS_INLET -',/&
         ' Message: Maximum number of particles in the system MAX_PIS',&
         /10X,' has been exceeded. Increase the value in mfix.dat',/&
         1X,70('*')/)


      RETURN
      END SUBROUTINE MASS_INFLOW_DEM



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SEED_NEW_PARTICLE                                       !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Purpose:  This routine uses the classification information to place !
!  a new particle in the proper location.                              !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SEED_NEW_PARTICLE(lBCV, lBCV_I, lM, lPOS, lIJKP, lOWNS)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE bc
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE param1
      USE physprop
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! The associated bc no.
      INTEGER, INTENT(IN) :: lBCV, lBCV_I
! Phase of incoming particle.
      INTEGER, INTENT(OUT) :: lM
! Position of incoming particle.
      DOUBLE PRECISION, INTENT(OUT) :: lPOS(3)
! I/J/K index of fluid cell containing the new particle.
      INTEGER, INTENT(OUT) :: lIJKP(3)
! Logical indicating that the current rank is the owner
      LOGICAL, INTENT(OUT) :: lOWNS

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! the associated bc no.
!      INTEGER :: BCV
! a random number between 0-1
      DOUBLE PRECISION RAND1, RAND2
! Array index of vacant position
      INTEGER :: VACANCY
! Number of array position.
      INTEGER :: OCCUPANTS

      INTEGER :: RAND_I

      DOUBLE PRECISION :: WINDOW


!      IF(PARTICLE_PLCMNT(lBCV_I) == 'ORDR')THEN
         VACANCY = DEM_MI(lBCV_I)%VACANCY
         OCCUPANTS = DEM_MI(lBCV_I)%OCCUPANTS
         DEM_MI(lBCV_I)%VACANCY = MOD(VACANCY,OCCUPANTS) + 1
!      ELSE
!         call bcast(lpar_rad)
!         call bcast(lpar_pos)
!         call bcast(m)
!      ENDIF

      lOWNS = (myPE == DEM_MI(lBCV_I)%OWNER(VACANCY))

! Only the rank that owns the new particle needs to continue.
      IF(.NOT.lOWNS) RETURN

! Obtain a random numbers from (0,1]
      CALL RANDOM_NUMBER(RAND1)
      CALL RANDOM_NUMBER(RAND2)

! Assign the phase of the incoming particle.
      IF(DEM_MI(lBCV_I)%POLYDISPERSE) THEN
         RAND_I = ceiling(dble(NUMFRAC_LIMIT)*RAND1)
         lM = DEM_BC_POLY_LAYOUT(lBCV_I, RAND_I)
      ELSE
         lM = DEM_BC_POLY_LAYOUT(lBCV_I,1)
      ENDIF

      WINDOW = DEM_MI(lBCV_I)%WINDOW
      RAND1 = HALF*DES_D_p0(lM) + (WINDOW - DES_D_p0(lM))*RAND1
      RAND2 = HALF*DES_D_p0(lM) + (WINDOW - DES_D_p0(lM))*RAND2


! Set the physical location and I/J/K location of the particle.
      SELECT CASE(BC_PLANE(lBCV))
      CASE('N','S')
         lPOS(1) = DEM_MI(lBCV_I)%P(VACANCY) + RAND1
         lPOS(3) = DEM_MI(lBCV_I)%Q(VACANCY) + RAND2
         lPOS(2) = DEM_MI(lBCV_I)%OFFSET

         lIJKP(1) = DEM_MI(lBCV_I)%W(VACANCY)
         lIJKP(3) = DEM_MI(lBCV_I)%H(VACANCY)
         lIJKP(2) = DEM_MI(lBCV_I)%L

      CASE('E','W')
         lPOS(2) = DEM_MI(lBCV_I)%P(VACANCY) + RAND1
         lPOS(3) = DEM_MI(lBCV_I)%Q(VACANCY) + RAND2
         lPOS(1) = DEM_MI(lBCV_I)%OFFSET

         lIJKP(2) = DEM_MI(lBCV_I)%W(VACANCY)
         lIJKP(3) = DEM_MI(lBCV_I)%H(VACANCY)
         lIJKP(1) = DEM_MI(lBCV_I)%L

      CASE('T','B')
         lPOS(1) = DEM_MI(lBCV_I)%P(VACANCY) + RAND1
         lPOS(2) = DEM_MI(lBCV_I)%Q(VACANCY) + RAND2
         lPOS(3) = DEM_MI(lBCV_I)%OFFSET

         lIJKP(1) = DEM_MI(lBCV_I)%W(VACANCY)
         lIJKP(2) = DEM_MI(lBCV_I)%H(VACANCY)
         lIJKP(3) = DEM_MI(lBCV_I)%L
      END SELECT

! Easier and cleaner to clear out the third component at the end.
      IF(NO_K) lPOS(3) = ZERO

      RETURN
      END SUBROUTINE SEED_NEW_PARTICLE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_NEW_PARTICLE_PROPS                                  !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Purpose:  Set the properties of the new particle.                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_NEW_PARTICLE_PROPS(lBCV, lM, lNP, lPOS, lIJKP)

      USE compar
      USE bc
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE param1
      USE physprop
      USE des_thermo
      use des_rxns

      use run, only: ENERGY_EQ
      use run, only: ANY_SPECIES_EQ
      use constant, only: PI

      use indices
      use functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! The associated bc no.
      INTEGER, INTENT(IN) :: lBCV
! Phase of incoming particle.
      INTEGER, INTENT(IN) :: lM
! Index of new particle
      INTEGER, INTENT(IN) :: lNP
! Position of incoming particle.
      DOUBLE PRECISION, INTENT(IN) :: lPOS(3)
! I/J/K index of fluid cell containing the new particle.
      INTEGER, INTENT(IN) :: lIJKP(3)

      INTEGER :: IJK

! Global phase index
      INTEGER :: BC_M

! Shift the phase index by SMAX to refernece global variables.
      BC_M = lM + SMAX

! Set the PEA Flags:
      PEA(lNP,1:2) = .TRUE.  ! The particle exists and is entering
      PEA(lNP,3:4) = .FALSE. ! It is not exiting nor a ghost particle

! Set the initial position values based on mass inlet class
      DES_POS_NEW(:,lNP) = lPOS(:)
      DES_VEL_NEW(1,lNP) = BC_U_s(lBCV,BC_M)
      DES_VEL_NEW(2,lNP) = BC_V_s(lBCV,BC_M)
      DES_VEL_NEW(3,lNP) = BC_W_s(lBCV,BC_M)

! Set the initial velocity values
      IF (DO_OLD) THEN
         DES_POS_OLD(:,lNP) = lPOS(:)
         DES_VEL_OLD(:,lNP) = DES_VEL_NEW(:,lNP)
         OMEGA_OLD(:,lNP) = 0
      ENDIF

! Set the initial angular velocity values
      OMEGA_NEW(:,lNP) = 0

! Set the particle radius value
      DES_RADIUS(lNP) = HALF * DES_D_P0(lM)

! Set the particle density value
      RO_Sol(lNP) = DES_RO_S(lM)

! Store the I/J/K indices of the particle.
      PIJK(lNP,1:3) = lIJKP(:)
      PIJK(lNP,4) = FUNIJK(lIJKP(1), lIJKP(2), lIJKP(3))

! Set the particle mass phase
      PIJK(lNP,5) = lM

! Calculate the new particle's Volume, Mass, OMOI
      PVOL(lNP) = (4.0d0/3.0d0) * PI * DES_RADIUS(lNP)**3
      PMASS(lNP) = PVOL(lNP) * RO_Sol(lNP)
      OMOI(lNP) = 5.d0 / (2.d0 * PMASS(lNP) * DES_RADIUS(lNP)**2)

! If solving the energy equations, set the temperature
      IF(ANY_SPECIES_EQ .OR. ENERGY_EQ ) THEN
         DES_T_s_NEW(lNP) = BC_T_s(lBCV,BC_M)
         DES_T_s_OLD(lNP) = DES_T_s_NEW(lNP)
      ENDIF

! Set species mass fractions
      IF((ENERGY_EQ .AND. C_PS0(BC_M)/=UNDEFINED) .OR. ANY_SPECIES_EQ)&
         DES_X_s(lNP,1:NMAX(BC_M)) = BC_X_s(lBCV,BC_M,1:NMAX(BC_M))

! Calculate time dependent physical properties
      CALL DES_PHYSICAL_PROP(lNP, .FALSE.)


      RETURN
      END SUBROUTINE SET_NEW_PARTICLE_PROPS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine:  DES_NEW_PARTICLE_TEST                                  !
!                                                                      !
!  Purpose:  This routine checks if a new particle placed using the    !
!  random inlet was placed in contact with an existing particle.  If   !
!  so a flag is set indicating contact, and the new particle is        !
!  repositioned within the inlet domain.                               !
!                                                                      !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Purpose: This routine has to be modified for parallel version       !
!           the parameter now accepts the lpar_rad and lpar_pos and    !
!           tests if it touches any particles                          !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_NEW_PARTICLE_TEST(BCV_I,ppar_rad,ppar_pos,TOUCHING)

      USE compar
      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param1
      USE physprop
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! index of boundary condition
      INTEGER, INTENT(IN) :: BCV_I
      DOUBLE PRECISION, INTENT(IN) :: ppar_pos(DIMN)
      DOUBLE PRECISION, INTENT(IN) :: ppar_rad
      LOGICAL, INTENT(INOUT) :: TOUCHING
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! particle number id of a potential overlapping/contacting particle
      INTEGER NP2
! total number of particles in current ijk cell and loop counter
      INTEGER NPG, LL
! i, j, k indices along boundary used for loop counters
      INTEGER I, J, K, IJK
! for parallel processing
      integer listart,liend,ljstart,ljend,lkstart,lkend

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT
      DOUBLE PRECISION  DISTVEC(DIMN), DIST, R_LM
!-----------------------------------------------

      TOUCHING = .FALSE.

! For parallel processing the arrays has to be limited
!      select case (des_mi_class(bcv_i))
!      case ('XW','XE', 'YZw','YZe')
!         listart = gs_array(bcv_i,1)
!         liend = gs_array(bcv_i,2)
!         ljstart = max(gs_array(bcv_i,3),jstart)
!         ljend = min(gs_array(bcv_i,4),jend)
!         lkstart = max(gs_array(bcv_i,5),jstart)
!         lkend = min(gs_array(bcv_i,6),jend)
!      case ('YN','YS', 'XZn','XZs')
!         listart = max(gs_array(bcv_i,1),istart)
!         liend = min(gs_array(bcv_i,2),iend)
!         ljstart = gs_array(bcv_i,3)
!         ljend = gs_array(bcv_i,4)
!         lkstart = max(gs_array(bcv_i,5),jstart)
!         lkend = min(gs_array(bcv_i,6),jend)
!      case ('ZT','ZB', 'XYt','XYb')
!         listart = max(gs_array(bcv_i,1),istart)
!         liend = min(gs_array(bcv_i,2),iend)
!         ljstart = max(gs_array(bcv_i,3),jstart)
!         ljend = min(gs_array(bcv_i,4),jend)
!         lkstart = gs_array(bcv_i,5)
!         lkend = gs_array(bcv_i,6)
!      end select

      DO k = lkstart,lkend
      DO j = ljstart,ljend
      DO i = listart,liend
!      DO K = GS_ARRAY(BCV_I,5), GS_ARRAY(BCV_I,6)
!         DO J = GS_ARRAY(BCV_I,3), GS_ARRAY(BCV_I,4)
!           DO I =  GS_ARRAY(BCV_I,1), GS_ARRAY(BCV_I,2)
             IJK = FUNIJK(I,J,K)
             IF(ASSOCIATED(PIC(IJK)%P)) THEN
               NPG =  SIZE(PIC(IJK)%P)
               DO LL = 1, NPG
                  NP2 = PIC(IJK)%P(LL)
                  DISTVEC(:) = ppar_pos(:) - DES_POS_NEW(:,NP2)
                  DIST = SQRT(DES_DOTPRDCT(DISTVEC,DISTVEC))
                  R_LM = ppar_rad + DES_RADIUS(NP2)
                  IF(DIST .LE. R_LM) TOUCHING = .TRUE.
               ENDDO
             ENDIF
           ENDDO
         ENDDO
       ENDDO

      RETURN
      END SUBROUTINE DES_NEW_PARTICLE_TEST
