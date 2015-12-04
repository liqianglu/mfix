! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!   Module name: DISCRETELEMENT                                        C
!   Purpose: DES mod file                                              C
!            Common Block containing DEM conditions                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE DERIVED_TYPES

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE multi_sweep_and_prune, only: multisap_t, boxhandlelist_t
      IMPLICIT NONE
!-----------------------------------------------
      ! the global multisap
      type(multisap_t) multisap

      type(boxhandlelist_t), DIMENSION(:),  ALLOCATABLE :: boxhandle         !(PARTICLES)

! Dynamic information related to computational (eulerian) fluid grid
!----------------------------------------------------------------->>>
! Dynamic variable. for each ijk computational fluid cell store the
! total number of particles and the id's of the particles in that cell
      TYPE iap1
         INTEGER, DIMENSION(:), POINTER:: p
      END TYPE iap1

      ! in order to facilitate the parallel processing the PIC is defined
      ! as single array IJK
      TYPE(iap1), DIMENSION(:), ALLOCATABLE:: pic  ! (DIMENSION_3)

! particle in cell related variable
      type iap2
         integer :: isize
         integer, dimension(:), pointer:: p
      end type iap2

      type(iap2), dimension(:),allocatable:: dg_pic

      END MODULE DERIVED_TYPES
