!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: read_database(Ier)                                     C
!  Purpose: read thermochemical database                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-OCT-05  C
!                                                                      C
!  Modification 1: J. Musser                          Date: 02-May-11  C
!  Purpose: Provided support for DEM access to database.               C
!                                                                      C
!  Modification 2: J. Musser                          Date: 02-Oct-12  C
!  Purpose: Calls to READ_DATABASE were moved to CHECK_DATA_04/05      C
!  duing input data integrity checks for the gas and solids phases.    C
!  Rather than looping through all species for each phase, the model   C
!  (TFM/DEM), phase index, species index, species name, and molecular  C
!  weight are passed as dummy arguments so that only infomration for   C
!  referenced species (lName) is obtained.                             C                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE READ_DATABASE(MODEL, lM, lN, lName, lMW)

      USE param 
      USE param1 
      USE physprop
      USE constant
      USE compar
      USE rxns
      USE funits 
      USE discretelement
      USE des_rxns

      IMPLICIT NONE

! Indicates the model. TFM, DEM
      CHARACTER(len=*), INTENT(IN) :: MODEL
! Phase and species indices
      INTEGER, INTENT(IN) :: lM, lN
! Species name from data file.
      CHARACTER(len=*), INTENT(IN) :: lName
! Species molecular weight from the data file (if any)
      DOUBLE PRECISION, INTENT(INOUT) :: lMW

! Molecular weight read from database.
      DOUBLE PRECISION dbMW

! Error message returned from Read_Therm and sent to calling routine
      INTEGER IER
! File unit of databases searched
      INTEGER FUNIT
! Loop counter for checking all three locations for data
      INTEGER FILE
! Input/Output error status ( 0 is no error)
      INTEGER IOS

! Identifies if an error was detected.
      LOGICAL ErrorFlag

! File name of Burcat and Ruscic database
      CHARACTER(len=10) :: THERM = 'BURCAT.THR'
! Full path to model directory
      INCLUDE 'mfix_directory_path.inc'
! Full path to Burcat and Ruscic database
      CHARACTER(len=147) FILENAME

! External function. Integrates the temperature-dependent specific 
! heat from zero to Tref.
      DOUBLE PRECISION, EXTERNAL :: calc_ICpoR
      DOUBLE PRECISION, EXTERNAL :: DES_calc_ICpoR
! Tcom +/- SMALL_NUMBER: This is done so that the specific heat poly-
! nomail can be evaluate at Tcom with the high and low coefficients.
      DOUBLE PRECISION :: xTc

! Various integrations of the specific heat polynomials:
      DOUBLE PRECISION :: ICpoR_TrL  ! 0.0 --> Tref using Alow
      DOUBLE PRECISION :: ICpoR_TcL  ! 0.0 --> Tcom using Alow
      DOUBLE PRECISION :: ICpoR_TcH  ! 0.0 --> Tcom using Ahigh

      LOGICAL :: testCp = .FALSE.

! Initialize the file unit to be used.
      FUNIT = UNIT_DAT  ! .dat file unit
! Read data from mfix.dat or from BURCAT.THR in run directory or
! mfix directory.
      FILE = 0
      DO
         FILE = FILE + 1
! Check for thermochemical data in the mfix.dat file.
         IF(FILE == 1) then
           OPEN(UNIT=FUNIT, FILE='mfix.dat', STATUS='OLD', IOSTAT= IOS)
	          IF(IOS /= 0) CYCLE       ! Cycle on opening error
           IF(myPE == PE_IO) THEN
              WRITE(*,1000)'mfix.dat'  ! Identify read file. (screen)
              WRITE(UNIT_LOG,1000)'mfix.dat' ! (log file)
           ENDIF
! Read thermochemical data from the BURCAT.THR database in the local
! run directory.
	        ELSEIF(FILE == 2) THEN
            OPEN(UNIT=FUNIT,FILE=TRIM(THERM), STATUS='OLD', IOSTAT= IOS)
	           IF(IOS /= 0) CYCLE         ! Cycle on opening error
            IF(myPE == PE_IO) THEN
               WRITE(*,1000) TRIM(THERM)  ! Identify read file.
               WRITE(UNIT_LOG,1000) TRIM(THERM)  ! (log file)
            ENDIF
! Read thermochemical data from the BURCAT.THR database in the model
! directory (model/thermochemical/BURCAT.THR).
     	   ELSEIF(file == 3) then
            FILENAME = trim(MFIX_PATH)//'/thermochemical/'//TRIM(THERM)
            OPEN(UNIT=FUNIT,FILE=TRIM(FILENAME), STATUS='OLD',IOSTAT= IOS)
	           IF(IOS /= 0) CYCLE            ! Cycle on opening error
            IF(myPE == PE_IO) THEN
               WRITE(*,1000) ('/thermochemical/'//TRIM(THERM)) ! (screen)
               WRITE(UNIT_LOG,1000) ('/thermochemical/'//TRIM(THERM))  ! (log file)
            ENDIF
        	ELSE
            EXIT 
        	ENDIF

         REWIND(UNIT=funit)
 
! Initialize the error flag
         IER = 0

         IF(MODEL == "TFM") THEN
            CALL READ_THERM(FUNIT, lName, Thigh(lM,lN), Tlow(lM,lN),&
               Tcom(lM,lN), dbMW, Ahigh(:,lM,lN), Alow(:,lM,lN),    &
               HfrefoR(lM,lN), IER)

            IF(IER == 0) THEN
! If the user did not supply a value for the gas phase molecular weight
! in the mfix.dat file, use the value from the database.
               IF(lMW == UNDEFINED) lMW = dbMW
! There are a number of species with Tlow as 300, for which the
! following calculation will produce an error because T_ref = 298.  So
! slightly extend validity of the correaltion.
               IF(ABS(Tlow(lM,lN)-T_ref)<=2.0D0 .AND. &
                  Tlow(lM,lN) > T_ref) Tlow(lM,lN) = T_ref

! Initialize the reference integrals.
               ICpoR_l(lM,lN) = ZERO
               ICpoR_h(lM,lN) = ZERO

! Calculate the integral of specific heat from zero to Tref using the
! Alow coefficients.
               ICpoR_TrL = calc_ICpoR(T_ref, lM, lN, IER)
! Calculate the integral of specific heat from zero to Tcom using the 
! Alow coefficients.
               xTc = Tcom(lM,lN)-SMALL_NUMBER
               ICpoR_TcL = calc_ICpoR(xTc, lM, lN, IER)
! Calculate the integral of specific heat from zero to Tcom using the 
! Ahigh coefficients.
               xTc = Tcom(lM,lN)+SMALL_NUMBER
               ICpoR_TcH = calc_ICpoR(xTc, lM, lN, IER)
! Store the integrals in global variables.
               ICpoR_l(lM,lN) = ICpoR_TrL
               ICpoR_h(lM,lN) = ICpoR_TcL - ICpoR_TcH - ICpoR_TrL
            ENDIF
! Get DEM solids phase species data.
         ELSEIF(MODEL == "DEM") THEN
! Read the database.
            CALL Read_Therm(FUNIT, lName, DES_Thigh(lM,lN),            &
               DES_Tlow(lM,lN), DES_Tcom(lM,lN), dbMW,                 &
               DES_Ahigh(1,lM,lN), DES_Alow(1,lM,lN),                  &
               DES_HfrefoR(lM,lN), IER)
            IF(IER == 0)THEN
! If the user did not supply a value for the solids phase molecular
! weight of species in the mfix.dat file, use the value from the database.
               IF(lMW == UNDEFINED) lMW = dbMW
! There are a number of species with Tlow as 300, for which the
! following calculation will produce an error because T_ref = 298.  So
! slightly extend validity of the correaltion.
               IF( abs(DES_Tlow(lM,lN)-T_ref) <= 2.0D0 .AND.           &
                  DES_Tlow(lM,lN) > T_ref )                            &
                  DES_Tlow(lM,lN) = T_ref

! Initialize the reference integrals.
               DES_ICpoR_l(lM,lN) = ZERO
               DES_ICpoR_h(lM,lN) = ZERO

! Calculate the integral of specific heat from zero to Tref using the
! Alow coefficients.
               ICpoR_TrL = DES_calc_ICpoR(T_ref, lM, lN, IER)
! Calculate the integral of specific heat from zero to Tcom using the 
! Alow coefficients.
               xTc = DES_Tcom(lM,lN)-SMALL_NUMBER
               ICpoR_TcL = DES_calc_ICpoR(xTc, lM, lN, IER)
! Calculate the integral of specific heat from zero to Tcom using the 
! Ahigh coefficients.
               xTc = DES_Tcom(lM,lN)+SMALL_NUMBER
               ICpoR_TcH = DES_calc_ICpoR(xTc, lM, lN, IER)

! Store the integrals in global variables.
               DES_ICpoR_l(lM,lN) = ICpoR_TrL
               DES_ICpoR_h(lM,lN) = ICpoR_TcL - ICpoR_TcH - ICpoR_TrL
            ENDIF
         ELSE
! No other models have been set to use the thermochemical database.
            IF(myPE == PE_IO) THEN
               WRITE(*,1020) TRIM(ADJUSTL(MODEL))
               WRITE(UNIT_LOG,1020) TRIM(ADJUSTL(MODEL))
            ENDIF
            CALL MFIX_EXIT(myPE)
         ENDIF

         ErrorFlag = .TRUE.
         IF(IER == 0) THEN
            IF(myPE == PE_IO) THEN
               WRITE(*,1001) 'Found!'
               WRITE(UNIT_LOG,1001) 'Found!'
            ENDIF
            ErrorFlag = .FALSE.

            if(testCP) CALL writeCp(MODEL, lM, lN, lName, lMW)
            EXIT
        ELSEIF(myPE == PE_IO) THEN
            WRITE(*,1001) 'Not found.'
            WRITE(UNIT_LOG,1001) 'Not found.'
        ENDIF

      ENDDO

! Error message control.
!-----------------------------------------------------------------------
! Write remaining error message if needed.
      IF(ErrorFlag) THEN
         IF(myPE .EQ. PE_IO) THEN
            WRITE(*,1010) trim(lName)
            WRITE(*,1011) TRIM(FILENAME)
            WRITE(UNIT_LOG,1010) trim(lName)
            WRITE(UNIT_LOG,1011) TRIM(FILENAME)
         ENDIF
         CALL MFiX_Exit(myPE)
      ENDIF

      RETURN  


! Messages
!-----------------------------------------------------------------------
 1000 FORMAT(5X,'Checking ',A,$)
 1001 FORMAT(1X,' :: ',A)

! Error Flags
!-----------------------------------------------------------------------
 1010 FORMAT(//1X,70('*')/, ' From: READ_DATABASE',/, ' Message:',     &
         ' Species "',A,'" was not matched to any entry in the',/      &
         ' thermochemical databases.')
 1011 FORMAT(/' SUGGESTION: Search the database for the exact species',&
         ' name. The',/' species names are case sensitive and should', &
         ' match the names in',/' BURCAT.THR exactly excluding',       &
         ' trailing blanks and tabs. Also verify',/' that the data',   &
         ' section in the mfix.dat file (if any) is below a line',/    &
         ' that starts with THERMO DATA.',//' Database location:',     &
         /1X,A,/1X,70('*')/)
 1012 FORMAT(1X,70('*')/)

 1020 FORMAT(//1X,70('*')/' From: READ_DATABASE',/                     &
         ' Message: This routine was entered with an invalid model',   &
         ' argument.',/' Acceptable values are TFM and DEM.',/         &
         ' Passed Argument: MODEL = ',A,/1X,70('*')//)

      END SUBROUTINE READ_DATABASE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: READ_DATABASE0(IER)                                    C
!  Purpose: Provides legacy support for rrates files.                  C
!                                                                      C
!  Author: J. Musser                                  Date: 02-Oct-12  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE READ_DATABASE0(Ier)

      USE param 
      USE param1 
      USE physprop
      USE constant
      USE compar
      USE rxns
      USE funits 
      USE discretelement
      USE des_rxns

      IMPLICIT NONE

! Error message returned from Read_Therm and sent to calling routine
      INTEGER IER
! Loop indicies for mass phase and species
      INTEGER M, N
! Loop counter for continuum and discrete species
      INTEGER Nsp, DES_Nsp

! Return to the calling routine if the database has already been called.
      IF(database_read)RETURN

! Set the flag identifying that the database has been read.
      database_read = .TRUE.

! Initialize counters
      Nsp = 0
      DES_Nsp = 0

! Read species data for the gas phase.
!-----------------------------------------------------------------------
      DO N = 1, NMAX(0)
         Nsp = Nsp + 1
! If a species name was not specified in mfix.dat, flag error and exit.
	        IF(SPECIES_NAME(Nsp) == UNDEFINED_C) THEN
            WRITE(*,1010) N         ! screen
            IF(DMP_LOG) WRITE(UNIT_LOG,1010) N  ! log file
	           CALL MFIX_EXIT(mypE)
	        ENDIF
! Read the database.
         CALL READ_DATABASE('TFM', 0, N, SPECIES_NAME(Nsp), MW_g(N))
	     ENDDO

! Read species data for the continuum solids phases.
!-----------------------------------------------------------------------
! Skip reading the database for the continuum solids phase if the 
! simulation is only employing discrete solids.
      IF(.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID)THEN
       	 DO M = 1, MMAX
            DO N = 1, NMAX(M)
	              Nsp = Nsp + 1
! If a species name was not specified in mfix.dat, flag error and exit.
	              IF(SPECIES_NAME(Nsp) == UNDEFINED_C)THEN
                  WRITE(*,1011)'continuum', M, N ! screen
                  IF(DMP_LOG) WRITE(UNIT_LOG,1011)'continuum', M, N
	                 CALL MFIX_EXIT(mypE)
	              ENDIF
               CALL READ_DATABASE('TFM', M, N, &
                  SPECIES_NAME(Nsp), MW_s(M,N))
	           ENDDO   ! N=1, NMAX(M)
	        ENDDO   ! M=1, MMAX
      ENDIF

      RETURN  

! Error Messages
!-----------------------------------------------------------------------
 1010 FORMAT(/1X,70('*')/, ' From: READ_DATABASE0',/, ' Message: ',    &
         'No SPECIES_NAME provided for gas phase species ',I3,'.',/' ',&
         'Check mfix.dat.',/1X,70('*')/)
 1011 FORMAT(/1X,70('*')/, ' From: READ_DATABASE0',/, ' Message: ',    &
         'No SPECIES_NAME provided for ',A,' solids phase ',I2,', ',/  &
         ' species ',I3,'.',/' Check mfix.dat.',/1X,70('*')/)

      END SUBROUTINE READ_DATABASE0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: read_database(Ier)                                     C
!  Purpose: read thermochemical database                               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-OCT-05  C
!                                                                      C
!  Modification 1: J. Musser                          Date: 02-May-11  C
!  Purpose: Provided support for DEM access to database.               C
!                                                                      C
!  Modification 2: J. Musser                          Date: 02-Oct-12  C
!  Purpose: Calls to READ_DATABASE were moved to CHECK_DATA_04/05      C
!  duing input data integrity checks for the gas and solids phases.    C
!  Rather than looping through all species for each phase, the model   C
!  (TFM/DEM), phase index, species index, species name, and molecular  C
!  weight are passed as dummy arguments so that only infomration for   C
!  referenced species (lName) is obtained.                             C                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE writeCp(MODEL, lM, lN, lName, lMW)


      USE physprop
      USE compar

! Universal gas constant in cal/mol.K
      use constant, only: RGAS => GAS_CONST_cal

      IMPLICIT NONE

! Indicates the model. TFM, DEM
      CHARACTER(len=*), INTENT(IN) :: MODEL

! Phase and species indices
      INTEGER, INTENT(IN) :: lM, lN
! Species name from data file.
      CHARACTER(len=*), INTENT(IN) :: lName
! Species molecular weight from the data file (if any)
      DOUBLE PRECISION, INTENT(in) :: lMW

      INTEGER :: IER1, IER2, lc

      DOUBLE PRECISION :: T
      DOUBLE PRECISION :: lCP, lICP

      DOUBLE PRECISION, EXTERNAL :: calc_CpoR
      DOUBLE PRECISION, EXTERNAL :: calc_ICpoR


      IF(MODEL == "DEM") return


      write(*,"(2/,3x,'Specific Heat report for ',A)")trim(lName)

      write(*,"(/12x,'Low',9x,'High')")
      write(*,"(6x,'T',3x,g11.5,2x,g11.5)") Tlow(lM,lN), Thigh(lM,lN)
      DO lc=1,5
         write(*,"(4x,'A(',I1,')',2(2x,g11.5))") lc, &
            Alow(lc,lM,lN), Ahigh(lc,lM,lN)
      ENDDO
      write(*,"('')")
      write(*,"(5x,'Tcom: ',g11.5)")Tcom(lM,lN)
      write(*,"('')")

      write(*,"(5x,'Temperature',8x,'Cp',11x,'ICp')")

      T = 300
      DO WHILE(T <= Thigh(lM,lN))

         IER1 = 0
         IER2 = 0

         write(*,"(7x,g11.5,$)") T
         lCP  = calc_CpoR(T, lM, lN, IER1) * RGAS / lMW
         lICP = calc_ICpoR(T, lM, lN, IER2) * RGAS / lMW
         write(*,"(2(3x,g11.5),$)")lCP, lICP

         IF(IER1 /= 0) write(*,"(3x,'Cp Error!',$)")
         IF(IER2 /= 0) write(*,"(3x,'ICp Error!',$)")
         write(*,"('')")

         T = T + 100.0
      ENDDO

      write(*,"('')")
      END SUBROUTINE writeCp
