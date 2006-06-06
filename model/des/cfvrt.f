!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFVRT(Vtan, VRl, TANGNT)                               C
!  Purpose: DES - Calculate the tangential component of                C
!           relative velocity                                          C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFVRT(Vtan, VRl, TANGNT)

      USE param1      
      USE discretelement
      IMPLICIT NONE

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

      INTEGER K
      DOUBLE PRECISION Vtan, VRl(DIMN), TANGNT(DIMN), Vno, NORM(DIMN)

!-----------------------------------------------------------------------

      Vtan = DES_DOTPRDCT(VRl,TANGNT)
      
      RETURN
      END SUBROUTINE CFVRT


