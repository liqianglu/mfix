!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
!
!  module name: des_write_restart
!  purpose: writing des data for restart
!
!  Author : Pradeep G
!  Purpose : Reads either single restart file or multiple restart files
!            (based on bdist_io) flag
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^c
      subroutine WRITE_RES0_DES

      use param1
      use compar
      use discretelement
      use run
      use des_bc
      use des_rxns
      use des_thermo
      use mfix_pic, only: MPPIC
      use mfix_pic, only: DES_STAT_WT

      use error_manager
      use write_res1_des
      use mpi_utility

      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      INTEGER :: LC1
      INTEGER :: lNEXT_REC
      INTEGER :: lDIMN

      DOUBLE PRECISION :: VERSION


! Set the version of the DES RES file.
      VERSION = 1.0d0
! Set the output dimension.
      lDIMN = merge(2,3,NO_K)

      CALL INIT_WRITE_RES_DES(trim(RUN_NAME), VERSION, lNEXT_REC)

      CALL WRITE_RES_DES(lNEXT_REC, VTP_FINDEX)
      CALL WRITE_RES_DES(lNEXT_REC, TECPLOT_FINDEX)
      CALL WRITE_RES_DES(lNEXT_REC, DTSOLID)

      DO LC1 = 1, lDIMN
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_POS_NEW(LC1,:))
      ENDDO

      CALL WRITE_RES_pARRAY(lNEXT_REC, iGLOBAL_ID)

      DO LC1 = 2, 4
         CALL WRITE_RES_pARRAY(lNEXT_REC, PEA(:,LC1))
      ENDDO

      DO LC1 = 1, lDIMN
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_VEL_NEW(LC1,:))
      ENDDO

      DO LC1 = 1, merge(1,3,NO_K)
         CALL WRITE_RES_pARRAY(lNEXT_REC, OMEGA_NEW(LC1,:))
      ENDDO

      CALL WRITE_RES_pARRAY(lNEXT_REC, DES_RADIUS)
      CALL WRITE_RES_pARRAY(lNEXT_REC, RO_SOL)

      IF(MPPIC) &
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_STAT_WT)

      IF(ENERGY_EQ) &
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_T_s_NEW)

      IF(ANY_SPECIES_EQ) THEN
         DO LC1=1, DIMENSION_N_S
            CALL WRITE_RES_pARRAY(lNEXT_REC, DES_X_s(:,LC1))
         ENDDO
      ENDIF

      DO LC1=1,2
         CALL WRITE_RES_cARRAY(lNEXT_REC, COLLISIONS(LC1,:), pLOC2GLB=.TRUE.)
      ENDDO

      CALL WRITE_RES_cARRAY(lNEXT_REC, PV_COLL(:))
      DO LC1=1, lDIMN
         CALL WRITE_RES_cARRAY(lNEXT_REC,PFN_COLL(LC1,:))
         CALL WRITE_RES_cARRAY(lNEXT_REC,PFT_COLL(LC1,:))
      ENDDO

      CALL WRITE_RES_DES(lNEXT_REC, DEM_BCMI)
      DO LC1=1, DEM_BCMI
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI_TIME(LC1))
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%VACANCY)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%OCCUPANTS)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%WINDOW)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%OFFSET)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%L)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%W(:))
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%H(:))
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%P(:))
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%Q(:))
      ENDDO

      CALL FINL_WRITE_RES_DES

! Notify that output was written.
      WRITE(ERR_MSG, 1000) trim(iVal(S_TIME))
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1000 FORMAT('DES restart data written at time = ',A)

      RETURN
      END SUBROUTINE WRITE_RES0_DES
