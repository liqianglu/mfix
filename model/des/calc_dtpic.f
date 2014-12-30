!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Module name: MAKE_ARRAYS_DES                                        !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: DES - allocating DES arrays
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DTPIC

      USE param1
      USE funits
      USE run
      USE compar
      USE discretelement
      USE cutcell
      use desmpi
      use mpi_utility
      USE geometry
      USE des_rxns
      USE des_thermo
      USE des_stl_functions

      use mfix_pic, only: CFL_PIC, DTPIC_CFL, DTPIC_TAUP
      use mfix_pic, only: DTPIC_MAX

      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: I, J, K, L, IJK, PC, SM_CELL


! MPPIC related quantities
      DOUBLE PRECISION :: DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------

      DTPIC_CFL = LARGE_NUMBER

      PC = 1
      DO L = 1, MAX_PIP
      IF(PC.GT.PIP) EXIT
         IF(.NOT.PEA(L,1)) CYCLE
         PC = PC+1
         IF(PEA(L,4)) CYCLE

         DTPIC_TMPX = (CFL_PIC*DX(PIJK(L,1)))/&
            (ABS(DES_VEL_NEW(1,L))+SMALL_NUMBER)
         DTPIC_TMPY = (CFL_PIC*DY(PIJK(L,2)))/&
            (ABS(DES_VEL_NEW(2,L))+SMALL_NUMBER)
         DTPIC_TMPZ = LARGE_NUMBER
         IF(DO_K) DTPIC_TMPZ = (CFL_PIC*DZ(PIJK(L,3)))/&
            (ABS(DES_VEL_NEW(3,L))+SMALL_NUMBER)

         DTPIC_CFL = MIN(DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ)
      ENDDO

      CALL global_all_max(DTPIC_CFL)

      DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)
      DTSOLID = DTPIC_MAX

      WRITE(ERR_MSG,2000) dtpic_cfl, dtpic_taup, DTSOLID
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 2000 FORMAT('DTPIC BASED ON CFL AND TAUP:', 2x, 2(2x,g11.4),          &
         /'DTSOLID set to ', g11.4)

      RETURN
      END SUBROUTINE CALC_DTPIC
