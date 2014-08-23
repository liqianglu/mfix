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

      use write_res1_des

      use mpi_utility

      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      INTEGER :: LC1, LC2
      INTEGER :: lNEXT_REC
      INTEGER :: lDIMN

      DOUBLE PRECISION :: VERSION 

!-----------------------------------------------
      VERSION = 1.7d0

      lDIMN = merge(2,3,NO_K)

      CALL INIT_WRITE_RES_DES(trim(RUN_NAME), lNEXT_REC)

      CALL WRITE_RES_DES(lNEXT_REC, VERSION)
      CALL WRITE_RES_DES(lNEXT_REC, VTP_FINDEX)
      CALL WRITE_RES_DES(lNEXT_REC, TECPLOT_FINDEX)
      CALL WRITE_RES_DES(lNEXT_REC, DTSOLID)

      DO LC1 = 1, lDIMN
         CALL WRITE_RES_DES(lNEXT_REC, DES_POS_NEW(LC1,:))
      ENDDO

      CALL WRITE_RES_DES(lNEXT_REC, iGLOBAL_ID)

      DO LC1 = 2, 4
         CALL WRITE_RES_DES(lNEXT_REC, PEA(:,LC1))
      ENDDO

      DO LC1 = 1, lDIMN
         CALL WRITE_RES_DES(lNEXT_REC, DES_VEL_NEW(LC1,:))
      ENDDO

      DO LC1 = 1, merge(1,3,NO_K)
         CALL WRITE_RES_DES(lNEXT_REC, OMEGA_NEW(LC1,:))
      ENDDO

      CALL WRITE_RES_DES(lNEXT_REC, DES_RADIUS)
      CALL WRITE_RES_DES(lNEXT_REC, RO_SOL)

      IF(MPPIC) &
         CALL WRITE_RES_DES(lNEXT_REC, DES_STAT_WT)

      IF(ENERGY_EQ) &
         CALL WRITE_RES_DES(lNEXT_REC, DES_T_s_NEW)

      IF(ANY_SPECIES_EQ) THEN
         DO LC1=1, DIMENSION_N_S 
            CALL WRITE_RES_DES(lNEXT_REC, DES_X_s(:,LC1))
         ENDDO
      ENDIF

      CALL WRITE_RES_DES(lNEXT_REC, NEIGHBOURS(:,1))
      CALL WRITE_RES_DES(lNEXT_REC, PN(1,:))
      CALL WRITE_RES_DES(lNEXT_REC, PN_WALL(1,:))

      DO LC1=2, MAXNEIGHBORS
         CALL WRITE_RES_DES(lNEXT_REC, NEIGHBOURS(:,LC1), pLOC2GLB=.TRUE.)
         CALL WRITE_RES_DES(lNEXT_REC, PN(LC1,:), pLOC2GLB=.TRUE.)
         CALL WRITE_RES_DES(lNEXT_REC, PV(LC1,:))
      ENDDO

      DO LC1=1, 6
         CALL WRITE_RES_DES(lNEXT_REC, PN_WALL(LC1,:), pLOC2GLB=.TRUE.)
         CALL WRITE_RES_DES(lNEXT_REC, PV_WALL(LC1,:))

         DO LC2=1, lDIMN
            CALL WRITE_RES_DES(lNEXT_REC, PFT(:,LC1,LC2))
            CALL WRITE_RES_DES(lNEXT_REC, PFT_WALL(:,LC1,LC2))
         ENDDO
      ENDDO

      CALL WRITE_RES_DES(lNEXT_REC, COLLISION_NUM)
      DO LC1=1, COLLISION_NUM
         CALL WRITE_RES_DES(lNEXT_REC, PV_COLL(LC1))
         DO LC2=1, lDIMN
            CALL WRITE_RES_DES(lNEXT_REC, PFN_COLL(LC2,LC1))
            CALL WRITE_RES_DES(lNEXT_REC, PFT_COLL(LC2,LC1))
         ENDDO
      ENDDO


      CALL WRITE_RES_DES(lNEXT_REC, DEM_BCMI)
      DO LC1=1, DEM_BCMI
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI_TIME(LC1))
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%VACANCY)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%OCCUPANTS)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%WINDOW)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%OFFSET)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%L)
         CALL WRITE_RES_DES_NPA(lNEXT_REC, DEM_MI(LC1)%W(:))
         CALL WRITE_RES_DES_NPA(lNEXT_REC, DEM_MI(LC1)%H(:))
         CALL WRITE_RES_DES_NPA(lNEXT_REC, DEM_MI(LC1)%P(:))
         CALL WRITE_RES_DES_NPA(lNEXT_REC, DEM_MI(LC1)%Q(:))
      ENDDO

      CALL FINL_WRITE_RES_DES

      RETURN
      END SUBROUTINE WRITE_RES0_DES