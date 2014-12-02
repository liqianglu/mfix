!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS

      use particle_filter, only: DES_INTERP_MEAN_FIELDS
      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE
      use particle_filter, only: DES_INTERP_GARG


      IF(DES_INTERP_MEAN_FIELDS) THEN
         SELECT CASE(DES_INTERP_SCHEME_ENUM)
         CASE(DES_INTERP_NONE) ; CALL COMP_MEAN_FIELDS_ZERO_ORDER
         CASE(DES_INTERP_GARG) ; CALL COMP_MEAN_FIELDS_INTERP0
         CASE DEFAULT; CALL COMP_MEAN_FIELDS_INTERP1
         END SELECT
      ELSE
         CALL COMP_MEAN_FIELDS_ZERO_ORDER
      ENDIF

!      CALL DIFFUSE_MEAN_FIELDS

      CALL CALC_EPG_DES


      END SUBROUTINE COMP_MEAN_FIELDS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE indices
      USE compar
      USE parallel
      USE sendrecv
      USE discretelement
      use desgrid
      use desmpi
      USE mfix_pic
      USE functions

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Loop counters: partciles, filter cells, phases
      INTEGER NP, LC, M
! Fluid cell index
      INTEGER IJK
! Total Mth solids phase volume in IJK
      DOUBLE PRECISION :: SOLVOLINC(DIMENSION_3,DES_MMAX)
! One divided by the total solids volume.
      DOUBLE PRECISION :: OoSOLVOL
! PVOL times statistical weight
      DOUBLE PRECISION :: VOL_WT
! Loop bound for filter
      INTEGER :: LP_BND

      SOLVOLINC(:,:) = ZERO

      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO

! Calculate the gas phae forces acting on each particle.
      DO NP=1,MAX_PIP
         IF(.NOT.PEA(NP,1)) CYCLE

         VOL_WT = PVOL(NP)
         IF(MPPIC) VOL_WT = VOL_WT*DES_STAT_WT(NP)
! Fluid cell containing the particle
         IJK = PIJK(NP,4)
! Particle phase for data binning.
         M = PIJK(NP,5)
! Accumulate total solids volume (by phase)
         SOLVOLINC(IJK,M) = SOLVOLINC(IJK,M) + VOL_WT
! Accumulate total solids momenum-ish (by phase)
         DES_U_S(IJK,M) = DES_U_S(IJK,M) +                             &
            DES_VEL_NEW(1,NP)*VOL_WT
         DES_V_S(IJK,M) = DES_V_S(IJK,M) +                             &
            DES_VEL_NEW(2,NP)*VOL_WT
         IF(DO_K) DES_W_S(IJK,M) = DES_W_S(IJK,M) +                    &
            DES_VEL_NEW(3,NP)*VOL_WT
      ENDDO

! Calculate the cell average solids velocity, the bulk density,
! and the void fraction.
!----------------------------------------------------------------//
!$omp parallel do if(ijkend3 .ge. 2000) default(shared)        &
!$omp private(IJK,M,OoSOLVOL)
      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE

! calculating the cell average solids velocity for each solids phase
         DO M = 1, DES_MMAX
            IF(SOLVOLINC(IJK,M).GT.ZERO) THEN
               OoSOLVOL = ONE/SOLVOLINC(IJK,M)
               DES_U_s(IJK,M) = DES_U_s(IJK,M)*OoSOLVOL
               DES_V_s(IJK,M) = DES_V_s(IJK,M)*OoSOLVOL
               IF(DO_K) DES_W_s(IJK,M) = DES_W_s(IJK,M)*OoSOLVOL
            ENDIF

! calculating the bulk density of solids phase m based on the total
! number of particles having their center in the cell
            DES_ROP_S(IJK,M) = DES_RO_S(M)*SOLVOLINC(IJK,M)/VOL(IJK)

         ENDDO   ! end loop over M=1,DES_MMAX

      ENDDO     ! end loop over IJK=ijkstart3,ijkend3
!$omp end parallel do

! Halo exchange of solids volume fraction data.
      CALL SEND_RECV(DES_ROP_S,2)

      RETURN
      END SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER
