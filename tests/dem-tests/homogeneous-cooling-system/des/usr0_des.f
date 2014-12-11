!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS0_DES                                               !
!                                                                      !
!  Purpose: This routine is called before the discrete phase time loop !
!  and is user-definable. The user may insert code in this routine or  !
!  call appropriate user defined subroutines.                          !
!                                                                      !
!  This routien is not called from a loop, hence all indicies are      !
!  undefined.                                                          !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR0_DES

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use run
      Use usr

      IMPLICIT NONE

      CALL WRITE_GRAN_TEMP

      RETURN
      END SUBROUTINE USR0_DES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS0_DES                                               !
!                                                                      !
!  Purpose: This routine is called before the discrete phase time loop !
!  and is user-definable. The user may insert code in this routine or  !
!  call appropriate user defined subroutines.                          !
!                                                                      !
!  This routien is not called from a loop, hence all indicies are      !
!  undefined.                                                          !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_GRAN_TEMP

      use discretelement, only: OVERLAP_MAX

      use discretelement, only: DES_MMAX
      use discretelement, only: GLOBAL_GRAN_TEMP
      use discretelement, only: DES_VEL_AVG
      use discretelement, only: des_ke,des_rote
      use discretelement, only: S_TIME
      use geometry, only: NO_K

      Use run, only: RUN_NAME


      use desmpi
      use mpi_utility
      use parallel
      use sendrecv
      
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! file name
      CHARACTER*64 :: FNAME
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: GT_UNIT = 2030

! Variables used for calculation in output file
      DOUBLE PRECISION :: GRAN_TEMP, AVG_VEL
! Working dimension
      INTEGER :: wDIMN

! Set the working dimension for the problem.
      wDIMN = merge(2, 3, NO_K)

      CALL DES_GRANULAR_TEMPERATURE


      if(mype.eq.pe_io)then

      GRAN_TEMP = (1.0/wDIMN)*SUM( GLOBAL_GRAN_TEMP(1:wDIMN) )
      AVG_VEL = (1.0/wDIMN)*SUM( DES_VEL_AVG(1:wDIMN) )
   

      FNAME = 'POST_GT.dat'
      INQUIRE(FILE=FNAME,EXIST=F_EXISTS)

      IF(.NOT.F_EXISTS) THEN
         OPEN(UNIT=GT_UNIT,FILE=FNAME,STATUS='NEW')
         WRITE(GT_UNIT,"(7(2X,A))") 'SOLIDS-TIME', 'GRAN-ENERGY', &
            'KNTC-ENERGY','ROT-ENERGY', 'VG-VELOCTY','MAX-OVERLAP','MAX-NEIGH'
      ELSE
         OPEN(UNIT=GT_UNIT,FILE=FNAME,&
            POSITION="APPEND",STATUS='OLD')
      ENDIF
      endif

!  SOLIDS-TIME  GRAN-ENERGY  KNTC-ENERGY  AVG-VELOCTY  MAX-OVERLAP  MAX-NEIGH
!xxFFFFFFFFFFFxxFFFFFFFFFFFxxFFFFFFFFFFFxxFFFFFFFFFFFxxFFFFFFFFFFFxxiiiiiiiii

      WRITE(GT_UNIT,"(6(2x,g11.5),2X,I5)") S_TIME, GRAN_TEMP, DES_KE,  &
          DES_ROTE,AVG_VEL, OVERLAP_MAX, NEIGH_MAX

      CLOSE(GT_UNIT)

      RETURN
      END SUBROUTINE WRITE_GRAN_TEMP


