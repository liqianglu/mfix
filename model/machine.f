!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: MACHINE_CONS                                           C
!  Purpose: set the machine constants    ( SGI ONLY )                  C
!                                                                      C
!  Author: P. Nicoletti                               Date: 28-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: OPEN_N1, NWORDS_DP, NWORDS_R, N_WORDS_I         C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE MACHINE_CONS
!
!
      Use machine
      IMPLICIT NONE
!
      OPEN_N1   = 512
      NWORDS_DP =  64
      NWORDS_R  = 128
      NWORDS_I  = 128
      JUST_FLUSH = .TRUE.
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_RUN_ID                                             C
!  Purpose: get the run id for this run                                C
!                                                                      C
!  Author: P. Nicoletti                               Date: 16-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: add ndoe name                                              C
!  Author: P.Nicoletti                                Date: 07-FEB-92  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR, ID_MINUTE   C
!                      ID_SECOND, ID_NODE                              C
!                                                                      C
!  Local variables: TIME_ARRAY                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_RUN_ID
!
      USE param
      USE run
      IMPLICIT NONE
!
!             temporary array to hold time data
      INTEGER DAT(8)
      CHARACTER*10 DATE, TIM, ZONE

      CALL DATE_AND_TIME(DATE, TIM, ZONE, DAT)
      ID_YEAR   = DAT(1)
      ID_MONTH  = DAT(2)
      ID_DAY    = DAT(3)
      ID_HOUR   = DAT(5)
      ID_MINUTE = DAT(6)
      ID_SECOND = DAT(7)

!     For SGI only
!      CALL GETHOSTNAME(ID_NODE,64)
!     For Linux with Portland Group compilers
      call hostnm(ID_NODE)
!
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function name: WALL_TIME (CPU)                                      C
!  Purpose: returns wall time since start of the run                   C
!                                                                      C
!  Author: P. Nicoletti                               Date: 10-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: TA, XT                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION WALL_TIME()
!
      IMPLICIT NONE

      INTEGER, SAVE :: COUNT_OLD=0, WRAP=0
!
! local variables
!                       clock cycle
      INTEGER           COUNT

!                       number of cycles per second
      INTEGER           COUNT_RATE

!                       max number of cycles, after which count is reset to 0
      INTEGER           COUNT_MAX

      CALL SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)
      IF(COUNT_OLD .GT. COUNT) THEN
        WRAP = WRAP + 1
      ENDIF
      COUNT_OLD = COUNT

      WALL_TIME      = DBLE(COUNT)/DBLE(COUNT_RATE) &
                     + DBLE(WRAP) * DBLE(COUNT_MAX)/DBLE(COUNT_RATE)
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: START_LOG                                              C
!  Purpose: does nothing ... for VAX/VMS compatibility (SGI ONLY)      C
!                                                                      C
!  Author: P. Nicoletti                               Date: 28-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE START_LOG
      IMPLICIT NONE
      RETURN
      END
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: END_LOG                                                C
!  Purpose: flushes the log file                                       C
!                                                                      C
!  Author: P. Nicoletti                               Date: 28-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE END_LOG
      USE funits
      IMPLICIT NONE
      IF(DMP_LOG)CALL FLUSH (UNIT_LOG)
      RETURN
      END
!
      subroutine slumber
      return
      end
!
      subroutine pc_quickwin
      return
      end
!


      subroutine ran
      return
      end

!
      subroutine flush_bin(iunit)
      implicit none
      integer :: iunit
      call flush(iunit)
      return
      end
      subroutine flush_res(iunit)
      implicit none
      integer :: iunit
      call flush(iunit)
      return
      end
