!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Write_Ab_m_var(A_m, b_m, var, IER)                  C 
!  Purpose: Write the sparse matrix coefficients and the               C
!           source vector.                                             C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_AB_M_VAR(A_M, B_M, VAR, IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE matrix 

      USE geometry      
      USE compar      
      USE mpi_utility  
      USE indices       
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Local index
      INTEGER          L, IJK1
!
!                      Phase index
      INTEGER          M
!
!                      cell index
      INTEGER          IJK
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3)
!
!                      Source vector
      DOUBLE PRECISION b_m(DIMENSION_3)

!                      Source vector
      DOUBLE PRECISION var(DIMENSION_3)

      double precision, allocatable :: array1(:) , array2(:)   !//
      double precision, allocatable :: am(:,:)                !//
!
!-----------------------------------------------
!
      integer i, j, k
      include 'function.inc'
      
      if (myPE == PE_IO) then
         allocate (array1(ijkmax3))    
         allocate (array2(ijkmax3))    
         allocate (am(ijkmax3,-3:3))    
      else
         allocate (array1(1))          
         allocate (array2(1))          
         allocate (am(1,-3:3))          
      end if

      if (myPE == PE_IO) then
         CALL START_LOG 
         WRITE (UNIT_LOG,*) ' Note : write_am_m is VERY inefficient '
         WRITE (*,*) ' Note : write_am_m is VERY inefficient '
         WRITE (UNIT_LOG,*) '  '
         WRITE (UNIT_LOG,*) ' A_m and B_m arrays below are in the '
         WRITE (UNIT_LOG,*) ' mfix INTERNAL order'
         WRITE (UNIT_LOG,*) ' '
         WRITE (UNIT_LOG, '(A,A)') &
           '  IJK  I  J  K   b         s         w         p         e       ', &
           '  n         t        Source     Variable' 
      end if


      do L = -3,3

      call gather(a_m(:,L),array1,root)

      DO K = Kmin2, Kmax2
      DO I = Imin2, Imax2
      DO J = Jmin2, Jmax2

      IJK = FUNIJK_GL(IMAP_C(I),JMAP_C(J),KMAP_C(K))
!     IJK = FUNIJK_GL(I,J,K)

      if (myPE == PE_IO) am(ijk,l) = array1(ijk)


      END DO 
      END DO 
      END DO 

      end do

      call gather(var(:),array1,root) 
      call gather(b_m(:),array2,root) 

      DO K = Kmin2, Kmax2
      DO I = Imin2, Imax2
      DO J = Jmin2, Jmax2

!     IJK = FUNIJK_GL(I,J,K)
      IJK = FUNIJK_GL(IMAP_C(I),JMAP_C(J),KMAP_C(K))

	if (myPE == PE_IO) WRITE (UNIT_LOG, '(I5, 3(I3), 9(1X,G9.2))') FUNIJK_IO(I,J,K), I, J, K,&
                                    (AM(ijk,L),L=-3,3), array2(IJK), array1(IJK)

      END DO
      END DO
      END DO

      if (myPE == PE_IO) CALL END_LOG 


      deallocate (array1)    !//
      deallocate (array2)    !//

      RETURN  
      END SUBROUTINE WRITE_AB_M_VAR 
      
!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 020 New local variables for parallelization: array1,array2,i,j,k
!// 400 Added mpi_utility module and other global reduction (gather) calls
