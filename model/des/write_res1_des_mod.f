      MODULE WRITE_RES1_DES

      use desmpi
      use compar, only: myPE
      use compar, only: PE_IO

      use cdist, only: bDist_IO

      use mpi_comm_des, only: DES_GATHER
      use mpi_comm_des, only: DESMPI_GATHERV

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: INIT_WRITE_RES_DES
      PUBLIC :: FINL_WRITE_RES_DES

      PUBLIC :: WRITE_RES_DES
      PUBLIC :: WRITE_RES_pARRAY
      PUBLIC :: WRITE_RES_cARRAY


! Write scalar and data WITHOUT MPI data collection.
      INTERFACE WRITE_RES_DES
         MODULE PROCEDURE WRITE_RES_DES_0I, WRITE_RES_DES_1I
         MODULE PROCEDURE WRITE_RES_DES_0D, WRITE_RES_DES_1D
         MODULE PROCEDURE WRITE_RES_DES_0L, WRITE_RES_DES_1L
      END INTERFACE

! Write particle array data.
      INTERFACE WRITE_RES_pARRAY
         MODULE PROCEDURE WRITE_RES_pARRAY_1I
         MODULE PROCEDURE WRITE_RES_pARRAY_1D
         MODULE PROCEDURE WRITE_RES_pARRAY_1L
      END INTERFACE

! Write neighbor/collision array data.
      INTERFACE WRITE_RES_cARRAY
         MODULE PROCEDURE WRITE_RES_cARRAY_1I
         MODULE PROCEDURE WRITE_RES_cARRAY_1D
         MODULE PROCEDURE WRITE_RES_cARRAY_1L
      END INTERFACE


      INTEGER, PARAMETER :: RDES_UNIT = 901

! Send/Recv parameters for Particle arrays:
      INTEGER :: pROOTCNT, pPROCCNT
      INTEGER :: pSEND
      INTEGER, allocatable :: pGATHER(:)
      INTEGER, allocatable :: pDISPLS(:)

! Send/Recv parameters for Collision/Neighbor arrays:
      INTEGER :: cROOTCNT, cPROCCNT
      INTEGER :: cSEND
      INTEGER, allocatable :: cGATHER(:)
      INTEGER, allocatable :: cDISPLS(:)


      CONTAINS

!``````````````````````````````````````````````````````````````````````!
! Subroutine: OPEN_RES_DES                                             !
!                                                                      !
! Purpose: Construct the file name and open the DES RES file.          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE OPEN_RES_DES(BASE)

      use machine, only: OPEN_N1

      CHARACTER(len=*), INTENT(IN)  :: BASE
      CHARACTER(len=32) :: lFNAME

      IF(bDIST_IO) THEN
         WRITE(lFNAME,'(A,I4.4,A)') BASE//'_DES_',myPE,'.RES'
         OPEN(UNIT=RDES_UNIT, FILE=lFNAME, FORM='UNFORMATTED',         &
            STATUS='UNKNOWN', ACCESS='DIRECT', RECL=OPEN_N1)

      ELSEIF(myPE == PE_IO) THEN
         WRITE(lFNAME,'(A,A)') BASE//'_DES.RES'
         OPEN(UNIT=RDES_UNIT, FILE=lFNAME, FORM='UNFORMATTED',         &
            STATUS='UNKNOWN', ACCESS='DIRECT', RECL=OPEN_N1)
      ENDIF

      END SUBROUTINE OPEN_RES_DES


!``````````````````````````````````````````````````````````````````````!
! Subroutine: INIT_WRITE_RES_DES                                       !
!                                                                      !
! Purpose: Construct the file name and open the DES RES file.          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE INIT_WRITE_RES_DES(BASE, lVERSION, lNEXT_REC)


      use mpi_utility, only: GLOBAL_SUM

      use desmpi, only: iProcBuf, dProcBuf
      use desmpi, only: iRootBuf, dRootBuf

      use desmpi, only: iGath_SendCnt

      use desmpi, only: iDisPls

      use discretelement, only: PIP, iGHOST_CNT

      use machine, only: OPEN_N1

      CHARACTER(len=*), INTENT(IN)  :: BASE
      DOUBLE PRECISION, INTENT(IN) :: lVERSION
      INTEGER, INTENT(OUT) :: lNEXT_REC

! Number of real particles on local rank
      INTEGER :: lParCnt, lPROC
! Total number of real particles.
      INTEGER :: lGHOST_CNT

      INTEGER :: lGatherCnts(0:NUMPEs-1)


      CALL OPEN_RES_DES(BASE)


      allocate(pGATHER(0:numPEs-1))
      allocate(pDISPLS(0:numPEs-1))

      allocate(cGATHER(0:numPEs-1))
      allocate(cDISPLS(0:numPEs-1))


      IF(bDIST_IO) THEN

         pROOTCNT = PIP
         pPROCCNT = pROOTCNT

         lGHOST_CNT = iGHOST_CNT

         cROOTCNT = PAIR_NUM
         cPROCCNT = cROOTCNT
      ELSE

! Setup data for particle array data collection:
         pROOTCNT = 10
         pPROCCNT = PIP - iGHOST_CNT

! Rank 0 gets the total number of gloabl particles.
         CALL GLOBAL_SUM(pPROCCNT, pROOTCNT)

! Serial IO does not store ghost particle data.
         lGHOST_CNT = 0

! Construct an array for the Root process that states the number of
! (real) particles on each process.
         pSEND = pPROCCNT

         lGatherCnts = 0
         lGatherCnts(myPE) = pPROCCNT

         CALL GLOBAL_SUM(lGatherCnts, pGATHER)

! Calculate the displacements for each process in the global array.
         pDISPLS(0) = 0
         DO lPROC = 1,NUMPES-1
            pDISPLS(lPROC) = pDISPLS(lPROC-1) + pGATHER(lPROC-1)
         ENDDO


! Setup data for neighbor/collision arrays
         cROOTCNT = 10
         cPROCCNT = PAIR_NUM

! Rank 0 gets the total number of gloabl particles.
         CALL GLOBAL_SUM(cPROCCNT, cROOTCNT)

! Construct an array for the Root process that states the number of
! (real) particles on each process.
         cSEND = cPROCCNT

         lGatherCnts = 0
         lGatherCnts(myPE) = cPROCCNT

         CALL GLOBAL_SUM(lGatherCnts, cGATHER)

! Calculate the displacements for each process in the global array.
         cDISPLS(0) = 0
         DO lPROC = 1,NUMPES-1
            cDISPLS(lPROC) = cDISPLS(lPROC-1) + cGATHER(lPROC-1)
         ENDDO

      ENDIF


! Write out the initial data.
      lNEXT_REC = 1
      CALL WRITE_RES_DES(lNEXT_REC, lVERSION)    ! RES file version
      CALL WRITE_RES_DES(lNEXT_REC, pROOTCNT)    ! Number of Particles
      CALL WRITE_RES_DES(lNEXT_REC, lGHOST_CNT)  ! Number of Ghosts
      CALL WRITE_RES_DES(lNEXT_REC, cROOTCNT)    ! Number of neighbors


      RETURN
      END SUBROUTINE INIT_WRITE_RES_DES


!``````````````````````````````````````````````````````````````````````!
! Subroutine: CLOSE_RES_DES                                            !
!                                                                      !
! Purpose: Close the DES RES file.                                     !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE FINL_WRITE_RES_DES

      IF(bDIST_IO .OR. myPE == PE_IO) close(RDES_UNIT)

      IF(allocated(dPROCBUF)) deallocate(dPROCBUF)
      IF(allocated(dROOTBUF)) deallocate(dROOTBUF)
      IF(allocated(iPROCBUF)) deallocate(iPROCBUF)
      IF(allocated(iROOTBUF)) deallocate(iROOTBUF)

      if(allocated(pGATHER)) deallocate(pGATHER)
      if(allocated(pDISPLS)) deallocate(pDISPLS)

      if(allocated(cGATHER)) deallocate(cGATHER)
      if(allocated(cDISPLS)) deallocate(cDISPLS)

      RETURN
      END SUBROUTINE FINL_WRITE_RES_DES



!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_0I                                         !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_0I(lNEXT_REC, INPUT_I)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(IN) :: INPUT_I

      IF(bDIST_IO .OR. myPE == PE_IO) &
         WRITE(RDES_UNIT, REC=lNEXT_REC) INPUT_I

      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE WRITE_RES_DES_0I


!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_1I                                         !
!                                                                      !
! Purpose: Write an array of integers to RES file. Note that data is   !
! not collected and hsould be on rank 0.                               !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_1I(lNEXT_REC, INPUT_I)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(IN) :: INPUT_I(:)

      INTEGER :: lSIZE

      lSIZE = size(INPUT_I)

      IF(bDIST_IO .OR. myPE == PE_IO) &
         CALL OUT_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)

      RETURN
      END SUBROUTINE WRITE_RES_DES_1I



!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_0D                                         !
!                                                                      !
! Purpose: Write scalar double percision values to RES file.           !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_0D(lNEXT_REC, INPUT_D)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(IN) :: INPUT_D

      IF(bDIST_IO .OR. myPE == PE_IO) &
         WRITE(RDES_UNIT, REC=lNEXT_REC) INPUT_D

      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE WRITE_RES_DES_0D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_1D                                         !
!                                                                      !
! Purpose: Write an array of double percision values to RES file. Note !
! that data is not collected and should be on rank 0.                  !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_1D(lNEXT_REC, INPUT_D)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(IN) :: INPUT_D(:)

      INTEGER :: lSIZE

      lSIZE = size(INPUT_D)

      IF(bDIST_IO .OR. myPE == PE_IO) &
         CALL OUT_BIN_512(RDES_UNIT, INPUT_D, lSIZE, lNEXT_REC)

      RETURN
      END SUBROUTINE WRITE_RES_DES_1D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_0L                                         !
!                                                                      !
! Purpose: Write scalar logical values to RES file.                    !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_0L(lNEXT_REC, INPUT_L)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(IN) :: INPUT_L

      INTEGER :: INPUT_I

      INPUT_I = merge(1,0,INPUT_L)

      IF(bDIST_IO .OR. myPE == PE_IO) &
         WRITE(RDES_UNIT, REC=lNEXT_REC) INPUT_I

      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE WRITE_RES_DES_0L


!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_DES_1L                                         !
!                                                                      !
! Purpose: Write an array of integers to RES file. Note that data is   !
! not collected and hsould be on rank 0.                               !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_DES_1L(lNEXT_REC, INPUT_L)

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(IN) :: INPUT_L(:)

      INTEGER, ALLOCATABLE :: INPUT_I(:)

      INTEGER :: lSIZE, LC1

      lSIZE = size(INPUT_L)
      ALLOCATE(INPUT_I(lSIZE))

      DO LC1=1, lSIZE
         INPUT_I(LC1) = merge(1,0,INPUT_L(LC1))
      ENDDO

      IF(bDIST_IO .OR. myPE == PE_IO) &
         CALL OUT_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)

      IF(allocated(INPUT_I)) deallocate(INPUT_I)

      RETURN
      END SUBROUTINE WRITE_RES_DES_1L



!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_PARRAY_1I                                      !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_PARRAY_1I(lNEXT_REC, INPUT_I, pLOC2GLB)

      use desmpi, only: iProcBuf
      use discretelement, only: PEA

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(IN) :: INPUT_I(:)
      LOGICAL, INTENT(IN), OPTIONAL :: pLOC2GLB

      LOGICAL :: lLOC2GLB
! Loop counters
      INTEGER :: LC1, LC2

      lLOC2GLB = .FALSE.
      IF(present(pLOC2GLB)) lLOC2GLB = pLOC2GLB

      allocate(iPROCBUF(pPROCCNT))
      allocate(iROOTBUF(pROOTCNT))

      iDISPLS = pDISPLS
      iGath_SendCnt = pSEND
      iGatherCnts   = pGATHER

      IF(bDIST_IO) THEN
         LC1 = 1

         IF(lLOC2GLB) THEN
            DO LC2 = 1, MAX_PIP
               IF(LC1 > PIP) EXIT
               IF(.NOT. PEA(LC1,1)) CYCLE
               iProcBuf(LC1) = iGLOBAL_ID(INPUT_I(LC2))
               LC1 = LC1 + 1
            ENDDO
         ELSE
            DO LC2 = 1, MAX_PIP
               IF(LC1 > PIP) EXIT
               IF(.NOT. PEA(LC1,1)) CYCLE
               iProcBuf(LC1) = INPUT_I(LC2)
               LC1 = LC1 + 1
            ENDDO
         ENDIF
         CALL OUT_BIN_512i(RDES_UNIT, iProcBuf, pROOTCNT, lNEXT_REC)

      ELSE
         CALL DES_GATHER(INPUT_I, lLOC2GLB)
         IF(myPE == PE_IO) &
            CALL OUT_BIN_512i(RDES_UNIT,iROOTBUF, pROOTCNT, lNEXT_REC)
      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE WRITE_RES_PARRAY_1I



!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_PARRAY_1D                                      !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_PARRAY_1D(lNEXT_REC, INPUT_D)

      use desmpi, only: iProcBuf
      use discretelement, only: PEA

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(IN) :: INPUT_D(:)

! Loop counters
      INTEGER :: LC1, LC2


      allocate(dPROCBUF(pPROCCNT))
      allocate(dROOTBUF(pROOTCNT))

      iDISPLS = pDISPLS
      iGath_SendCnt = pSEND
      iGatherCnts   = pGATHER

      IF(bDIST_IO) THEN
         LC1 = 1
         DO LC2 = 1, MAX_PIP
            IF(LC1 > PIP) EXIT
            IF(.NOT. PEA(LC1,1)) CYCLE
            dProcBuf(LC1) = INPUT_D(LC2)
            LC1 = LC1 + 1
         ENDDO
         CALL OUT_BIN_512(RDES_UNIT, dProcBuf, pROOTCNT, lNEXT_REC)
      ELSE
         CALL DES_GATHER(INPUT_D)
         IF(myPE == PE_IO) &
            CALL OUT_BIN_512(RDES_UNIT,dRootBuf, pROOTCNT, lNEXT_REC)
      ENDIF

      deallocate(dPROCBUF)
      deallocate(dROOTBUF)

      RETURN
      END SUBROUTINE WRITE_RES_PARRAY_1D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_PARRAY_1D                                      !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_PARRAY_1L(lNEXT_REC, INPUT_L)

      use desmpi, only: iProcBuf
      use discretelement, only: PEA

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(IN) :: INPUT_L(:)

! Loop counters
      INTEGER :: LC1, LC2

      allocate(iPROCBUF(pPROCCNT))
      allocate(iROOTBUF(pROOTCNT))

      iDISPLS = pDISPLS
      iGath_SendCnt = pSEND
      iGatherCnts   = pGATHER

      IF(bDIST_IO) THEN
         LC1 = 1
         DO LC2 = 1, MAX_PIP
            IF(LC1 > PIP) EXIT
            IF(.NOT. PEA(LC1,1)) CYCLE
            iProcBuf(LC1) = merge(1,0,INPUT_L(LC2))
            LC1 = LC1 + 1
         ENDDO
         CALL OUT_BIN_512i(RDES_UNIT, iProcBuf, pROOTCNT, lNEXT_REC)
      ELSE
         CALL DES_GATHER(INPUT_L)
         IF(myPE == PE_IO) &
            CALL OUT_BIN_512i(RDES_UNIT,iRootBuf, pROOTCNT, lNEXT_REC)
      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE WRITE_RES_PARRAY_1L

!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_cARRAY_1I                                      !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_cARRAY_1I(lNEXT_REC, INPUT_I, pLOC2GLB)

      use desmpi, only: iProcBuf
      use discretelement, only: PEA

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(IN) :: INPUT_I(:)
      LOGICAL, INTENT(IN), OPTIONAL :: pLOC2GLB

      LOGICAL :: lLOC2GLB
! Loop counters
      INTEGER :: LC1

      lLOC2GLB = .FALSE.
      IF(present(pLOC2GLB)) lLOC2GLB = pLOC2GLB

      allocate(iPROCBUF(cPROCCNT))
      allocate(iROOTBUF(cROOTCNT))

      iDISPLS = cDISPLS
      iGath_SendCnt = cSEND
      iGatherCnts   = cGATHER

      IF(lLOC2GLB) THEN
         DO LC1 = 1, cPROCCNT
            iProcBuf(LC1) = iGLOBAL_ID(INPUT_I(LC1))
         ENDDO
      ELSE
         iProcBuf(1:cPROCCNT) = INPUT_I(1:cPROCCNT)
      ENDIF

      IF(bDIST_IO) THEN
         CALL OUT_BIN_512i(RDES_UNIT, iProcBuf, cPROCCNT, lNEXT_REC)

      ELSE
         CALL DESMPI_GATHERV(pTYPE=1)
         IF(myPE == PE_IO) &
            CALL OUT_BIN_512i(RDES_UNIT,iROOTBUF, cROOTCNT, lNEXT_REC)
      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE WRITE_RES_cARRAY_1I

!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_cARRAY_1D                                      !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_cARRAY_1D(lNEXT_REC, INPUT_D)

      use desmpi, only: dPROCBUF ! Local process buffer
      use desmpi, only: dROOTBUF ! Root process buffer

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(IN) :: INPUT_D(:)

! Loop counters
      INTEGER :: LC1

      allocate(dPROCBUF(cPROCCNT))
      allocate(dROOTBUF(cROOTCNT))

      iDISPLS = cDISPLS
      iGath_SendCnt = cSEND
      iGatherCnts   = cGATHER

      dProcBuf(1:cPROCCNT) = INPUT_D(1:cPROCCNT)

      IF(bDIST_IO) THEN
         CALL OUT_BIN_512(RDES_UNIT, dProcBuf, cPROCCNT, lNEXT_REC)

      ELSE
         CALL DESMPI_GATHERV(pTYPE=2)
         IF(myPE == PE_IO) &
            CALL OUT_BIN_512(RDES_UNIT, dROOTBUF, cROOTCNT, lNEXT_REC)
      ENDIF

      deallocate(dPROCBUF)
      deallocate(dROOTBUF)

      RETURN
      END SUBROUTINE WRITE_RES_cARRAY_1D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: WRITE_RES_cARRAY_1L                                      !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE WRITE_RES_cARRAY_1L(lNEXT_REC, INPUT_L)

      use desmpi, only: iProcBuf
      use discretelement, only: PEA

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(IN) :: INPUT_L(:)

! Loop counters
      INTEGER :: LC1

      allocate(iPROCBUF(cPROCCNT))
      allocate(iROOTBUF(cROOTCNT))

      iDISPLS = cDISPLS
      iGath_SendCnt = cSEND
      iGatherCnts   = cGATHER

      DO LC1 = 1, cPROCCNT
         iProcBuf(LC1) = merge(1,0,INPUT_L(LC1))
      ENDDO

      IF(bDIST_IO) THEN
         CALL OUT_BIN_512i(RDES_UNIT, iProcBuf, cPROCCNT, lNEXT_REC)

      ELSE
         CALL DESMPI_GATHERV(pTYPE=1)
         IF(myPE == PE_IO) &
            CALL OUT_BIN_512i(RDES_UNIT,iROOTBUF, cROOTCNT, lNEXT_REC)
      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE WRITE_RES_cARRAY_1L

      END MODULE WRITE_RES1_DES
