!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: READ_NAMELIST(POST)                                    C
!  Purpose: Read in the NAMELIST variables                             C
!                                                                      C
!  Author: P. Nicoletti                               Date: 25-NOV-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:  Modifications for generating 3 character filebasename,    C
!            for each PE, i.e., XXX000.LOG, XXX001.LOG                 C
!            Each PE reads it's own mfix.dat, e.g.                     C
!                for PE 0 : mfix000.dat, for PE 1 : mfix001.dat etc.   C
!            Added checks to inquire whether mfixXXX.dat exists or not C           
!            Modifications for preliminary checks for NO_K = TRUE AND  C
!            number of processors (numPEs) > 1 situation               C
!            (see comment line : !// 300 0809)                         C
!                                                                      C
!  Author:   Aeolus Res. Inc.                         Date: 05-AUG-99  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: RUN_NAME, DESCRIPTION, UNITS, RUN_TYPE, TIME,   C
!                      TSTOP, DT, RES_DT, SPX_DT , OUT_DT, USR_DT,     C
!                      USR_TYPE, USR_VAR, USR_FORMAT, USR_EXT, USR_X_w C
!                      USR_X_e , USR_Y_s , USR_Y_n , USR_Z_b , USR_Z_t,C
!                      NLOG, COORDINATES, IMAX, DX, XLENGTH, JMAX, DY, C
!                      YLENGTH, KMAX, DZ, ZLENGTH, MMAX, D_p, RO_s,    C
!                      EP_star, MU_g0, MW_AVG, IC_X_w, IC_X_e, IC_Y_s, C
!                      IC_Y_n, IC_Z_b, IC_Z_t, IC_I_w, IC_I_e, IC_J_s, C
!                      IC_J_n, IC_K_b, IC_K_t, IC_EP_g, IC_P_g,        C
!                      IC_ROP_s, IC_T_g, IC_T_s,  IC_U_g,     C
!                      IC_U_s, IC_V_g, IC_V_s, IC_W_g, IC_W_s, BC_X_w, C
!                      BC_X_e, BC_Y_s, BC_Y_n, BC_Z_b, BC_Z_t, BC_I_w, C
!                      BC_I_e, BC_J_s, BC_J_n, BC_K_b, BC_K_t, BC_EP_g,C
!                      BC_P_g, BC_RO_g, BC_ROP_g, BC_ROP_s, BC_T_g,    C
!                      BC_T_s,  BC_U_g, BC_U_s,BC_V_g, BC_V_s,C
!                      BC_W_g, BC_W_s, BC_TYPE, BC_VOLFLOW_g,          C
!                      BC_VOLFLOW_s, BC_MASSFLOW_g, BC_MASSFLOW_s,     C
!                      BC_DT_0, BC_Jet_g0, BC_DT_h, BC_Jet_gh, BC_DT_l,C
!                      BC_Jet_gl, RO_g0                                C
!                                                                      C
!  Local variables: LINE_STRING, MAXCOL, LINE_TOO_BIG, COMMENT_INDEX   C
!                   SEEK_COMMENT                                       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE READ_NAMELIST(POST) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run
      USE output
      USE physprop
      USE geometry
      USE ic
      USE is
      USE bc
      USE fldvar
      USE constant
      USE indices
      USE toleranc 
      USE funits 
      USE scales 
      USE ur_facs 
      USE leqsol 
      USE residual
      USE rxns
      USE scalars
      USE compar      !// 001 Include MPI header file
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!               This routine is called from: 0 -  mfix; 1 - post_mfix
      INTEGER POST 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!               LINE_STRING(1:MAXCOL) has valid input data
      INTEGER, PARAMETER :: MAXCOL = 80 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!               holds one line in the input file
      CHARACTER LINE_STRING*132
!
!               integer function ... if > 0 ... data past column
!               MAXCOL in input file
      INTEGER   LINE_TOO_BIG 
!
!               Length of noncomment string
      INTEGER   LINE_LEN
!
!               integer function which returns COMMENT_INDEX
      INTEGER   SEEK_COMMENT 
!
!                      Coefficient of restitution (old symbol)
      DOUBLE PRECISION e

!
!                      Indicate whether currently reading rxns or rate
      LOGICAL          RXN_FLAG
!
!                      Indicate whether to do a namelist read on the line
      LOGICAL          READ_FLAG

!
!                      blank line function
      LOGICAL          BLANK_LINE
      INTEGER   L
!// 500 0805 define local variables, i1, i10, i100 for generating filebasename
      INTEGER :: i1, i10, i100
      LOGICAL :: present

!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!      INTEGER , EXTERNAL :: LINE_TOO_BIG, SEEK_COMMENT 
!      LOGICAL , EXTERNAL :: BLANK_LINE 
!-----------------------------------------------
!

      INCLUDE 'usrnlst.inc' 
      INCLUDE 'namelist.inc'
!
!
      E = UNDEFINED 
!
      RXN_FLAG = .FALSE. 
!
      READ_FLAG = .TRUE. 
!
      NO_OF_RXNS = 0 

!// 500 0805 Generate the decimals for file basename for each PE (max 999 PEs!)
      i100 = int(myPE/100)
      i10  = int((myPE-i100*100)/10)
      i1   = int((myPE-i100*100-i10*10)/1)

      i100 = i100 + 48
      i10  = i10  + 48
      i1   = i1   + 48
!// 500 0805 generate filebasename for each active PE
      fbname=char(i100)//char(i10)//char(i1)

!//AIKEPARDBG
!      if(idbgpar.ge.1)&                     !//AIKEPARDBG
!       write(*,FMT="('File basename for PE ',I3,' is ',A3)") myPE,fbname  !//AIKEPARDBG

!// 500 0805 ensure that mfixXXX.dat exists
      inquire(file='mfix.dat',exist=present)
      if(.not.present) then
        write(*,"('(PE ',I3,'): input data file, ',A11,' is missing: run aborted')") &
            myPE,'mfix.dat'
        call mfix_exit(myPE) !// 990 0807 Abort all PEs, not only the current one
      endif
!
!
! OPEN MFIX ASCII INPUT FILE, AND LOOP THRU IT
!
      OPEN(UNIT=UNIT_DAT, FILE='mfix.dat', STATUS='OLD', ERR=910) !// 500 modified filename
!
  100 CONTINUE 
      READ (UNIT_DAT, 1100, END=500) LINE_STRING 
!
      LINE_LEN = SEEK_COMMENT(LINE_STRING,LEN(LINE_STRING)) - 1 
      CALL REMOVE_COMMENT (LINE_STRING, LINE_LEN + 1, LEN(LINE_STRING)) 
!
      IF (LINE_LEN <= 0) GO TO 100               !comment line 
      IF (BLANK_LINE(LINE_STRING)) GO TO 100     !blank line 
!
      IF (LINE_TOO_BIG(LINE_STRING,LINE_LEN,MAXCOL) > 0) THEN 
         WRITE (*, 1300) LINE_STRING 
!//  MPI_Exit
         CALL MFIX_EXIT(myPE) 
      ENDIF 
!
!
      CALL MAKE_UPPER_CASE (LINE_STRING, LINE_LEN) 
!
!
!  Complete arithmetic operations and expand line
!
!
      CALL PARSE_LINE (LINE_STRING, LINE_LEN, RXN_FLAG, READ_FLAG) 
!
!
! Write the current line to a scratch file
! and read the scratch file in NAMELIST format
!
      IF (READ_FLAG) THEN 
         OPEN(UNIT=UNIT_TMP, STATUS='SCRATCH', ERR=900) 
         WRITE (UNIT_TMP, 1000) 
         WRITE (UNIT_TMP, 1100) LINE_STRING(1:LINE_LEN) 
         WRITE (UNIT_TMP, 1200) 
         REWIND (UNIT=UNIT_TMP) 
!       READ (UNIT_TMP,NML=INPUT_DATA,ERR=930,END=930)  ! Use this for FPS
         READ (UNIT_TMP, NML=INPUT_DATA, ERR=400, END=930) 
         GO TO 410 
  400    CONTINUE 
         IF (POST == 1) GO TO 410 
         REWIND (UNIT=UNIT_TMP)                  !If called from POST_MFIX ignore the error 
         WRITE (UNIT_TMP, 1010) 
         WRITE (UNIT_TMP, 1100) LINE_STRING(1:LINE_LEN) 
         WRITE (UNIT_TMP, 1200) 
         REWIND (UNIT=UNIT_TMP) 
         READ (UNIT_TMP, NML=USR_INPUT_DATA, ERR=930, END=930) 
  410    CONTINUE 
         CLOSE(UNIT=UNIT_TMP) 
      ENDIF 
!
      GO TO 100 
!
  500 CONTINUE 

!//AIKEPARDBGSTOP 0922
!      write(*,"('(PE ',I2,'): SPECIES_EQ : ',12L2)") myPE,SPECIES_EQ !//AIKEPARDBG
!      call mfix_exit(myPE) !//AIKEPARDBG	 

      CLOSE(UNIT=UNIT_DAT) 
!
      IF (E /= UNDEFINED) C_E = E 
!
      DO L = 1, NO_OF_RXNS 
!
!  Do all the rxns have a rate defined?
!
         IF (.NOT.GOT_RATE(L)) THEN 
!
            WRITE (*, 1610) myPE,RXN_NAME(L) !//PAR_I/O added processor id for output
!            STOP  
             call mfix_exit(myPE) !// 990 0807 Abort all PEs, not only the current one
!
         ENDIF 
!
!  Do all the rxns have a stoichiometry defined?
!
         IF (.NOT.GOT_RXN(L)) THEN 
!
            WRITE (*, 1620) myPE,RXN_NAME(L) !//PAR_I/O added processor id for output
!            STOP  
             call mfix_exit(myPE) !// 990 0807 Abort all PEs, not only the current one
!
         ENDIF 
      END DO 


!//AIKEPARDBG debug output to make sure data read from mfix.dat is same on all PEs
!      open(unit=UNIT_OUT,file='test'//fbname//'.out',status='UNKNOWN')  !//AIKEPARDBG
!      CALL WRITE_OUT0                                  !//AIKEPARDBG
!      call MPI_Barrier(MPI_COMM_WORLD,mpierr)			!//AIKEPARDBG
!      write(*,"('(PE ',I3,'): reached end of read_namelist')") myPE	!//AIKEPARDBG
!      call mfix_exit(myPE)								!//AIKEPARDBG

      RETURN  
!
! HERE IF AN ERROR OCCURED OPENNING/READING THE FILES
!
  900 CONTINUE 
      WRITE (*, 1400) 
      CALL MFIX_EXIT(myPE) 
  910 CONTINUE 
      WRITE (*, 1500) 
      CALL MFIX_EXIT(myPE) 
  930 CONTINUE 
      WRITE (*, 1600) LINE_STRING(1:LINE_LEN) 
      CALL MFIX_EXIT(myPE) 
!
 1000 FORMAT(1X,'$INPUT_DATA') 
 1010 FORMAT(1X,'$USR_INPUT_DATA') 
 1100 FORMAT(A) 
 1200 FORMAT(1X,'$END') 
 1300 FORMAT(/1X,70('*')//' From: READ_NAMELIST',/' Message: ',&
         'The following line in mfix.dat has data past column 80',/1X,A80,/1X,&
         70('*')/) 
 1400 FORMAT(/1X,70('*')//' From: READ_NAMELIST',/' Message: ',&
         'Unable to open scratch file',/1X,70('*')/) 
 1500 FORMAT(/1X,70('*')//' From: READ_NAMELIST',/' Message: ',&
         'Unable to open mfix.dat file',/1X,70('*')/) 
 1600 FORMAT(/1X,70('*')//' From: READ_NAMELIST',/' Message: ',&
         'Error while reading the following line in mfix.dat',/1X,A80,/1X,&
         'Possible causes are 1. incorrect format, 2. unknown name,',/1X,&
         '3. the item is dimensioned too small (see PARAM.INC file).',/1X,70(&
         '*')/) 
 1610 FORMAT(/1X,70('*')//'(PE ',I3,'): From: READ_NAMELIST',/' Message: ',&
         'No rxn rate defined for rxn: ',A,/1X,70('*')/) 
 1620 FORMAT(/1X,70('*')//'(PE ',I3,'): From: READ_NAMELIST',/' Message: ',&
         'No stoichiometry defined for rxn: ',A,/1X,70('*')/) 
!
      END SUBROUTINE READ_NAMELIST 
!
!
      logical function blank_line (line) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character line*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l 
!-----------------------------------------------
!
!
      blank_line = .FALSE. 
      do l = 1, len(line) 
       if (line(l:l)/=' ' .and. line(l:l)/='    ') return
      end do 
!
      blank_line = .TRUE. 
      return  
!
      end function blank_line 
