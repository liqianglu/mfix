!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: Mass_OUTFLOW_PIC                                        !
!  Author: R. Garg                                   Date: 23-Jun-14   !
!                                                                      !
!  Purpose:  Routine to delete out of domain parcels for PIC           !
!  implementation                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE MASS_OUTFLOW_PIC

      USE error_manager
      USE mpi_utility
      use bc
      use discretelement
      use functions
      use pic_bc

      implicit none

      INTEGER :: IJK
      INTEGER :: LC, LP, NP
      INTEGER :: BCV, BCV_I

      DOUBLE PRECISION :: DIST


      DO BCV_I = 1, PIC_BCMO

         BCV = PIC_BCMO_MAP(BCV_I)

         DO LC=PIC_BCMO_IJKSTART(BCV_I), PIC_BCMO_IJKEND(BCV_I)
            IJK = PIC_BCMO_IJK(LC)
            DO LP= 1,PINC(IJK)

               NP = PIC(IJK)%p(LP)
               IF(IS_NONEXISTENT(NP)) CYCLE

               SELECT CASE (BC_PLANE(BCV))
               CASE('S'); DIST = YN(BC_J_s(BCV)-1) - DES_POS_NEW(2,NP)
               CASE('N'); DIST = DES_POS_NEW(2,NP) - YN(BC_J_s(BCV))
               CASE('W'); DIST = XE(BC_I_w(BCV)-1) - DES_POS_NEW(1,NP)
               CASE('E'); DIST = DES_POS_NEW(1,NP) - XE(BC_I_w(BCV))
               CASE('B'); DIST = ZT(BC_K_b(BCV)-1) - DES_POS_NEW(3,NP)
               CASE('T'); DIST = DES_POS_NEW(3,NP) - ZT(BC_K_b(BCV))
               END SELECT

               IF(DIST < DES_RADIUS(NP)) CALL DELETE_PARCEL(NP)

            ENDDO
         ENDDO
      ENDDO


      RETURN
      END SUBROUTINE MASS_OUTFLOW_PIC




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DELETE_PARCEL                                           !
!  Author: R. Garg                                    Date: 23-Jun-14  !
!                                                                      !
!  Purpose:  Routine to delete parcel                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DELETE_PARCEL(NP)

      USE compar
      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param1
      USE physprop
      USE mfix_pic
      USE functions

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NP

      CALL SET_NONEXISTENT(NP)

      DES_POS_OLD(:,NP) = ZERO
      DES_POS_NEW(:,NP) = ZERO
      DES_VEL_OLD(:,NP) = ZERO
      DES_VEL_NEW(:,NP) = ZERO
      DES_RADIUS(NP) = ZERO
      PMASS(NP) = ZERO
      PVOL(NP) = ZERO
      RO_Sol(NP) = ZERO
      OMOI(NP) = ZERO

      DES_STAT_WT(NP) = ZERO

      FC(:,NP) = ZERO

      PIP = PIP - 1

      RETURN
      END SUBROUTINE DELETE_PARCEL
