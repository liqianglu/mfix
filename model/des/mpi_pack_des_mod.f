!----------------------------------------------------------------------!
!  Module: MPI_PACK_DES                                                !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
!  Purpose: Contains routines for packing real and ghost particles     !
!     into the MPI send buffers.                                       !
!----------------------------------------------------------------------!
      MODULE MPI_PACK_DES

      PRIVATE
      PUBLIC :: DESMPI_PACK_PARCROSS, DESMPI_PACK_GHOSTPAR

      interface pack_dbuf
        module procedure pack_db0
        module procedure pack_db1
        module procedure pack_i0
        module procedure pack_i1
        module procedure pack_l0
      end interface pack_dbuf

      CONTAINS

!----------------------------------------------------------------------!
!  Subroutine: DESMPI_PACK_GHOSTPAR                                    !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
! Purpose: Packs ghost particle in the send buffer.                    !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_PACK_GHOSTPAR(PFACE)

! Global Variables:
!---------------------------------------------------------------------//
! Size of ghost particle packets
      use desmpi, only: iGhostPacketSize
! Flag indicating that the ghost particle was updated
      use discretelement, only: iGHOST_UPDATED
! The MPI send buffer
      use desmpi, only: dSENDBUF
! Buffer offset
      use desmpi, only: iBUFOFFSET
! Runtime flag for solving the energy equations
      use run, only: ENERGY_EQ
! The neighbor processor's rank
      use desmpi, only: iNEIGHPROC
! DES grid cell containing each particle: current/previous
      use discretelement, only: DG_PIJKPRV
! The global ID for each particle
      use discretelement, only: iGLOBAL_ID
! Particle positions: current/previous
      use discretelement, only: DES_POS_NEW
! Particle tangential velocities: current/previous
      use discretelement, only: DES_VEL_NEW
! Particle rotational velocities: current/previous
      use discretelement, only: OMEGA_NEW
! Particle tempertures. current/previous
      use des_thermo, only: DES_T_s_NEW
! Particle radius, volume
      use discretelement, only: DES_RADIUS
! Number of cells used in interpolation
      use particle_filter, only: FILTER_SIZE
! Cells and weights for interpolation
      use particle_filter, only: FILTER_CELL, FILTER_WEIGHT
! Map to fluid grid cells and solids phase (I,J,K,IJK,M)
      use discretelement, only: PIJK
! Number of particles on the process (max particle array size)
      use discretelement, only: MAX_PIP
! User-defined variables for each particle.
      use discretelement, only: DES_USR_VAR, DES_USR_VAR_SIZE
! Function to convert DES grid IJK to new proc value.
      use desgrid, only: dg_ijkconv
! Size of the send buffer
      use desmpi, only: isendcnt
! Offset for particles with cycle BCs (otherwise zero)
      use desmpi, only: dcycl_offset
! Map of particles to DES grid
      use derived_types, only: DG_PIC
! Cell number of ghost particles
      use desmpi, only: iSENDINDICES

! Global Constants:
!---------------------------------------------------------------------//
      use constant, only: PI
! Dimension of particle spatial arrays.
      use discretelement, only: DIMN
      use functions, only: is_exiting
      use functions, only: is_ghost, is_entering_ghost, is_exiting_ghost

      IMPLICIT NONE

! Dummy arguments:
!---------------------------------------------------------------------//
! Processor boundary being packed (Top/Bottom/North/South/East/West)
      INTEGER, INTENT(IN) :: PFACE

! Local variables
!---------------------------------------------------------------------//
      integer :: lijk,lindx,ltot_ind,lpicloc,lpar_cnt,lcurpar
      integer :: lbuf
!......................................................................!

      lpar_cnt = 0
      ltot_ind = isendindices(1,pface)
      do lindx = 2,ltot_ind+1
         lijk = isendindices(lindx,pface)
         do lpicloc =1,dg_pic(lijk)%isize
            lbuf = lpar_cnt*iGhostPacketSize + ibufoffset
            lcurpar = dg_pic(lijk)%p(lpicloc)

! Do not send particle data for a ghost particle whose owner has not yet
! updated the particle's data on this processor.
            if((is_ghost(lcurpar) .or. &
                is_entering_ghost(lcurpar) .or. &
                is_exiting_ghost(lcurpar)) .and. &
                .not.ighost_updated(lcurpar)) cycle

! 1) Global ID
            call pack_dbuf(lbuf,iglobal_id(lcurpar),pface)
! 2) DES grid IJK
            call pack_dbuf(lbuf,dg_ijkconv(lijk,pface,                 &
               ineighproc(pface)),pface)
! 3) DES grid IJK - previous
            call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(lcurpar),pface,  &
               ineighproc(pface)),pface)
! 4) Radius
            call pack_dbuf(lbuf,des_radius(lcurpar),pface)
! 5) Phase index
            call pack_dbuf(lbuf,pijk(lcurpar,5),pface)
! 6) Position
            call pack_dbuf(lbuf,des_pos_new(lcurpar,:)+                &
               dcycl_offset(pface,:),pface)
! 7) Translational Velocity
            call pack_dbuf(lbuf,des_vel_new(lcurpar,:),pface)
! 8) Rotational Velocity
            call pack_dbuf(lbuf,omega_new(lcurpar,:),pface)
! 9) Exiting particle flag
            call pack_dbuf(lbuf,merge(1,0,is_exiting(lcurpar).or.is_exiting_ghost(lcurpar)),pface)
! 10) Temperature
            IF(ENERGY_EQ) &
               call pack_dbuf(lbuf,des_t_s_new(lcurpar),pface)
! 11) User Variable
            IF(DES_USR_VAR_SIZE > 0) &
               call pack_dbuf(lbuf,des_usr_var(:,lcurpar),pface)
! 12) Interpolation weights
            IF(FILTER_SIZE > 0) THEN
               call pack_dbuf(lbuf,filter_cell(:,lcurpar),pface)
               call pack_dbuf(lbuf,filter_weight(:,lcurpar),pface)
            ENDIF

            lpar_cnt = lpar_cnt + 1
         end do
      end do
      dsendbuf(1+mod(pface,2))%facebuf(1)=lpar_cnt
      isendcnt(pface) = lpar_cnt*iGhostPacketSize + ibufoffset

      end subroutine desmpi_pack_ghostpar

!----------------------------------------------------------------------!
!  Subroutine: DESMPI_PACK_PARCROSS                                    !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
! Purpose: Packs real particle in the send buffer.                     !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_PACK_PARCROSS(PFACE)

! Global Variables:
!---------------------------------------------------------------------//
! The MPI send buffer
      use desmpi, only: dSENDBUF
! Buffer offset
      use desmpi, only: iBUFOFFSET
! Runtime flag for solving the energy equations
      use run, only: ENERGY_EQ
! Runtime flag for solving species equations
      use run, only: ANY_SPECIES_EQ
! Runtime flag for MPPIC solids
      use mfix_pic, only: MPPIC
! DES grid cell containing each particle: current/previous
      use discretelement, only: DG_PIJKPRV
! The neighbor processor's rank
      use desmpi, only: iNEIGHPROC
! The statistical weight of each particle.
      use mfix_pic, only: DES_STAT_WT
! The global ID for each particle
      use discretelement, only: iGLOBAL_ID
! Particle positions: current/previous
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
! Particle tangential velocities: current/previous
      use discretelement, only: DES_VEL_NEW, DES_VEL_OLD
! Particle rotational velocities: current/previous
      use discretelement, only: OMEGA_NEW, OMEGA_OLD
! Particle orientation
      use discretelement, only: PARTICLE_ORIENTATION,ORIENTATION
! Particle radius, volume, density, mass
      use discretelement, only: DES_RADIUS, PVOL, RO_SOL, PMASS
! Previous value for particle acceleration (tangential/rotational)
      use discretelement, only: DES_ACC_OLD, ROT_ACC_OLD
! Particle species composition
      use des_rxns, only: DES_X_s
! Particle tempertures. current/previous
      use des_thermo, only: DES_T_s_NEW, DES_T_s_OLD
! Force arrays acting on the particle
      use discretelement, only: FC, TOW
! One of the moment of inertia
      use discretelement, only: OMOI
! Map to fluid grid cells and solids phase (I,J,K,IJK,M)
      use discretelement, only: PIJK
! Flag to send/recv old (previous) values
      use discretelement, only: DO_OLD
! Number of particles on the process (max particle array size)
      use discretelement, only: PIP, MAX_PIP
! Number of ghost particles on the current process
      use discretelement, only: iGHOST_CNT
! User-defined variables for each particle.
      use discretelement, only: DES_USR_VAR, DES_USR_VAR_SIZE
! Particle pair (neighborhood) arrays:
      use discretelement, only: NEIGHBORS, NEIGHBOR_INDEX, NEIGH_NUM
! Pair collision history information
      use discretelement, only: PFT_NEIGHBOR
! Dimension of particle spatial arrays.
      use discretelement, only: DIMN
! Flag indicating the the fluid-particle drag is explicitly coupled.
      use discretelement, only: DES_EXPLICITLY_COUPLED
! Explicit particle drag force
      use discretelement, only: DRAG_FC
! Cells and weights for interpolation
      use particle_filter, only: FILTER_WEIGHT

      use desgrid, only: dg_ijkconv, icycoffset
      use desmpi, only: dcycl_offset, isendcnt
      use desmpi, only: irecvindices

      use desmpi, only: iParticlePacketSize
      use desmpi, only: iPairPacketSize
      use derived_types, only: dg_pic

      use functions

      implicit none

! Dummy arguments:
!---------------------------------------------------------------------//
! Processor boundary being packed (Top/Bottom/North/South/East/West)
      INTEGER, INTENT(IN) :: PFACE

! Local variables
!---------------------------------------------------------------------//
      integer :: li, lj, lk
      integer :: ltot_ind,lindx,cc
      integer :: lneigh,lijk,&
                 lpicloc,lparcnt,lcurpar
      integer :: lbuf,num_neighborlists_to_send

      logical, allocatable, dimension(:) :: going_to_send

! Location in the buffer where the number of pair data is specified.
      integer :: num_neighborlists_send_buf_loc
!......................................................................!

! pack the particle crossing the boundary
      ltot_ind = irecvindices(1,pface)
      lparcnt = 0

      allocate(going_to_send(max_pip))
      going_to_send(:) = .false.

      do lindx = 2,ltot_ind+1
         lijk = irecvindices(lindx,pface)
         do lpicloc = 1,dg_pic(lijk)%isize
            lcurpar = dg_pic(lijk)%p(lpicloc)

! if ghost particle then cycle
            if(is_ghost(lcurpar) .or. &
               is_entering_ghost(lcurpar) .or. &
               is_exiting_ghost(lcurpar)) cycle

            going_to_send(lcurpar) = .true.
            lbuf = lparcnt*iParticlePacketSize + ibufoffset

! 1) Global ID
            call pack_dbuf(lbuf,iglobal_id(lcurpar),pface)
! 2) DES Grid IJK
            call pack_dbuf(lbuf,dg_ijkconv(lijk,pface,                 &
               ineighproc(pface)),pface)
! 3) DES grid IJK - previous
            call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(lcurpar),pface,  &
               ineighproc(pface)),pface)
! 4) Radius
            call pack_dbuf(lbuf,des_radius(lcurpar),pface)
! 5) Fluid cell I index with cycle offset
            li = pijk(lcurpar,1) + icycoffset(pface,1)
            call pack_dbuf(lbuf,li,pface)
! 6) Fluid cell J index with cycle offset
            lj = pijk(lcurpar,2) + icycoffset(pface,2)
            call pack_dbuf(lbuf,lj,pface)
! 7) Fluid cell K index with cycle offset
            lk = pijk(lcurpar,3) + icycoffset(pface,3)
            call pack_dbuf(lbuf,lk,pface)
! 8) Fluid cell IJK on destination process
            call pack_dbuf(lbuf,funijk_proc(li,lj,lk,                  &
               ineighproc(pface)),pface)
! 9) Particle solids phase index
            call pack_dbuf(lbuf,pijk(lcurpar,5),pface)
! 10) Entering particle flag.
            call pack_dbuf(lbuf, is_entering(lcurpar).or.is_entering_ghost(lcurpar), pface)
! 11) Exiting particle flag.
            call pack_dbuf(lbuf, is_exiting(lcurpar).or.is_exiting_ghost(lcurpar), pface)
! 12) Density
            call pack_dbuf(lbuf,ro_sol(lcurpar),pface)
! 13) Volume
            call pack_dbuf(lbuf,pvol(lcurpar),pface)
! 14) Mass
            call pack_dbuf(lbuf,pmass(lcurpar),pface)
! 15) 1/Moment of Inertia
            call pack_dbuf(lbuf,omoi(lcurpar),pface)
! 16) Position with cyclic shift
            call pack_dbuf(lbuf,des_pos_new(lcurpar,:) +               &
               dcycl_offset(pface,:),pface)
! 17) Translational velocity
            call pack_dbuf(lbuf,des_vel_new(lcurpar,:),pface)
! 18) Rotational velocity
            call pack_dbuf(lbuf,omega_new(lcurpar,:),pface)
! 19) Accumulated translational forces
            call pack_dbuf(lbuf,fc(lcurpar,:),pface)
! 20) Accumulated torque forces
            call pack_dbuf(lbuf,tow(lcurpar,:),pface)
! 21) Temperature
            IF(ENERGY_EQ) &
               call pack_dbuf(lbuf,des_t_s_new(lcurpar),pface)
! 22) Species composition
            IF(ANY_SPECIES_EQ) &
               call pack_dbuf(lbuf,des_x_s(lcurpar,:),pface)
! 23) Explicit drag force
            IF(DES_EXPLICITLY_COUPLED) &
               call pack_dbuf(lbuf, drag_fc(lcurpar,:),pface)
! 24) User defined variable
            IF(DES_USR_VAR_SIZE > 0) &
               call pack_dbuf(lbuf, des_usr_var(:,lcurpar),pface)
! 25) Particle orientation
            IF(PARTICLE_ORIENTATION) &
               call pack_dbuf(lbuf,orientation(:,lcurpar),pface)

! -- Higher order integration variables
            IF (DO_OLD) THEN
! 26) Position (previous)
               call pack_dbuf(lbuf,des_pos_old(lcurpar,:) +            &
                  dcycl_offset(pface,:),pface)
! 27) Translational velocity (previous)
               call pack_dbuf(lbuf,des_vel_old(lcurpar,:),pface)
! 28) Rotational velocity (previous)
               call pack_dbuf(lbuf,omega_old(lcurpar,:),pface)
! 29) Translational acceleration (previous)
               call pack_dbuf(lbuf,des_acc_old(lcurpar,:),pface)
! 30) Rotational acceleration (previous)
               call pack_dbuf(lbuf,rot_acc_old(lcurpar,:),pface)
! 31) Temperature (previous)
               IF(ENERGY_EQ) &
                  call pack_dbuf(lbuf,des_t_s_old(lcurpar),pface)
            ENDIF

! PIC particles are removed and the number of particles on the processor
! is decremented.
            IF (MPPIC) THEN
! 32) Statistical weight
               call pack_dbuf(lbuf,des_stat_wt(lcurpar),pface)
               call set_nonexistent(lcurpar)
               pip = pip - 1

! DEM particles are converted to ghost particles. This action does not
! change the number of particles.
            ELSE
               if (is_entering(lcurpar)) then
                  call set_entering_ghost(lcurpar)
               elseif (is_exiting(lcurpar)) then
                  call set_exiting_ghost(lcurpar)
               else
                  call set_ghost(lcurpar)
               endif
               ighost_cnt = ighost_cnt + 1
            END IF

! Clear out the force array.
            fc(lcurpar,:) = 0.
            lparcnt = lparcnt + 1
         end do
      end do

! Calculate the location in buffer where the number of pair data is
! stored and skip specifying the entry. After all the pair data is
! packed, then this value is set.
      lbuf = lparcnt*iParticlePacketSize + ibufoffset
      num_neighborlists_send_buf_loc = lbuf
      lbuf = lbuf+1

       num_neighborlists_to_send = 0
       lcurpar = 1
       do cc = 1, NEIGH_NUM
          IF (0 .eq. NEIGHBORS(cc)) EXIT

          IF (cc.eq.NEIGHBOR_INDEX(lcurpar)) THEN
             lcurpar = lcurpar + 1
          ENDIF

! Only packup pairing data for particles being transfered.
          if (.not. going_to_send(lcurpar)) cycle

! Do not send pairing data if the pair no longer exists or if the
! particle is exiting as it may be locatable during unpacking.
          lneigh = neighbors(lcurpar)
          if(is_nonexistent(lneigh)) cycle
          if(is_exiting(lneigh)) cycle

! 34) Global ID of particle being packed.
          call pack_dbuf(lbuf,iglobal_id(lcurpar),pface)
! 35) DES grid IJK of cell receiving the particle.
          call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(lcurpar),pface,    &
               ineighproc(pface)),pface)
! 36) Global ID of neighbor particle.
          call pack_dbuf(lbuf,iglobal_id(lneigh),pface)
! 37) DES grid IJK of cell containing the neighbor particle.
          call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(lneigh),pface,     &
               ineighproc(pface)),pface)
! 38) Tangential collision history.
          call pack_dbuf(lbuf,PFT_NEIGHBOR(:,CC),pface) 
! Increment the number of pairs being sent.
          num_neighborlists_to_send = num_neighborlists_to_send + 1
       enddo

! Store the number of pair datasets being sent. This information is
! stored before the pairing data so the receiving process knows the
! amount of data to 'unpack.'
      lbuf = num_neighborlists_send_buf_loc
! 33) Number of pair datasets.
      call pack_dbuf(lbuf,num_neighborlists_to_send,pface)

      dsendbuf(1+mod(pface,2))%facebuf(1) = lparcnt
      isendcnt(pface) = lparcnt*iParticlePacketSize +                  &
         num_neighborlists_to_send*iPairPacketSize + ibufoffset + 3

      deallocate(going_to_send)

      RETURN
      END SUBROUTINE DESMPI_PACK_PARCROSS

!----------------------------------------------------------------------!
! PACK SUBROUTINE FOR SINGLE REAL VARIABLES                            !
!----------------------------------------------------------------------!
      subroutine pack_db0(lbuf,idata,pface)
      use desmpi, only: dSENDBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(in) :: idata

      dsendbuf(1+mod(pface,2))%facebuf(lbuf) = idata
      lbuf = lbuf + 1

      return
      end subroutine pack_db0

!----------------------------------------------------------------------!
! Pack subroutine for real arrays                                      !
!----------------------------------------------------------------------!
      subroutine pack_db1(lbuf,idata,pface)
      use desmpi, only: dSENDBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(in) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      dsendbuf(1+mod(pface,2))%facebuf(lbuf:lbuf+lsize-1) = idata
      lbuf = lbuf + lsize

      return
      end subroutine pack_db1

!----------------------------------------------------------------------!
! Pack subroutine for single integer variables                         !
!----------------------------------------------------------------------!
      subroutine pack_i0(lbuf,idata,pface)
      use desmpi, only: dSENDBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(in) :: idata

      dsendbuf(1+mod(pface,2))%facebuf(lbuf) = idata
      lbuf = lbuf + 1

      return
      end subroutine pack_i0

!----------------------------------------------------------------------!
! Pack subroutine for integer arrays                                   !
!----------------------------------------------------------------------!
      subroutine pack_i1(lbuf,idata,pface)
      use desmpi, only: dSENDBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(in) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      dsendbuf(1+mod(pface,2))%facebuf(lbuf:lbuf+lsize-1) = idata
      lbuf = lbuf + lsize

      return
      end subroutine pack_i1

!----------------------------------------------------------------------!
! Pack subroutine for logical scalars                                  !
!----------------------------------------------------------------------!
      subroutine pack_l0(lbuf,ldata,pface)
      use desmpi, only: dSENDBUF

      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      logical, intent(in) :: ldata

      dsendbuf(1+mod(pface,2))%facebuf(lbuf) = merge(1.0, 0.0, ldata)
      lbuf = lbuf + 1

      return
      end subroutine pack_l0

      end module mpi_pack_des
