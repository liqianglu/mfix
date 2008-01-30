!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GAS_DRAG(A_M, B_M, VXF_GS, IER, UV, VV, WV)            C
!  Purpose: DES - Accounting for the equal and opposite drag force     C
!           on gas due to particles by introducing the drag            C      
!           as a source term. Face centered                            C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!     
      SUBROUTINE GAS_DRAG(A_M, B_M, VXF_GS, IER, UV, VV, WV)
!-----------------------------------------------
!     M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_g
      USE rxns
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is
      USE tau_g
      USE bc
      USE compar    
      USE sendrecv 
      USE discretelement
      USE drag 
      
      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     D u m m y   A r g u m e n t s
!-----------------------------------------------
!     
!     
!     Error index
      INTEGER          IER
!     
!     Indices
      INTEGER          IJK, NC
!     
!     Phase index
      INTEGER          M, UV, VV, WV
!     Averaging Factor 
      DOUBLE PRECISION :: AVG_FACTOR
!     Grid indices
      INTEGER          I,J,K
!     
      DOUBLE PRECISION USFCM, VSFCM, WSFCM, VCELL
      DOUBLE PRECISION A_M(DIMENSION_3, -3:3, 0:DIMENSION_M)
      DOUBLE PRECISION B_M(DIMENSION_3, 0:DIMENSION_M)
      DOUBLE PRECISION VXF_GS(DIMENSION_3, DIMENSION_M)  
      Logical DES_ON_TMP
      DOUBLE PRECISION tmp_A, tmp_B
 !     
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      AVG_FACTOR = 0.25D0*(DIMN-2) + 0.5D0*(3-DIMN)
	DES_ON_TMP = DES_INTERP_ON
      IF(UV.EQ.1) THEN
         DO M = 1, MMAX
            DO IJK = IJKSTART3, IJKEND3
               IF(FLUID_AT(IJK)) THEN
                  if(DES_ON_TMP) then 
                     I = I_OF(IJK)
                     J = J_OF(IJK)
                     K = K_OF(IJK)
                     
                     tmp_A =  - AVG_FACTOR*(DRAG_AM(I&
                          &,J,K,M) + DRAG_AM(I,J-1,K,M))
                     tmp_B =  - AVG_FACTOR*(DRAG_BM(I&
                          &,J,K,1,M) + DRAG_BM(I,J-1,K,1,M))
                         
                     IF(DIMN.eq.3) then 
                        tmp_A = tmp_A - AVG_FACTOR*(DRAG_AM(I&
                             &,J,K-1,M) + DRAG_AM(I,J-1,K-1,M))
                        tmp_B = tmp_B - AVG_FACTOR*(DRAG_BM(I&
                             &,J,K-1,1,M) + DRAG_BM(I,J-1,K-1,1,M))
                     end IF
                     
                     VCELL = VOL(IJK)
                     
                     tmp_A = tmp_A*VCELL
                     tmp_B = tmp_B*VCELL
                     
                     else 
                        USFCM = AVG_X(DES_U_S(IJK,M),DES_U_S(EAST_OF(IJK),M),I_OF(IJK))
                        tmp_A =  - VXF_GS(IJK,M)
                        tmp_B =  - VXF_GS(IJK,M)*USFCM
                     end if
                     
                     A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
                     B_M(IJK,0) = B_M(IJK,0) + tmp_B
               END IF
            END DO
         END DO


      ELSE IF(VV.EQ.1) THEN

         DO M = 1, MMAX
            DO IJK = IJKSTART3, IJKEND3
               IF(FLUID_AT(IJK)) THEN
                  I = I_OF(IJK)
                  J = J_OF(IJK)
                  K = K_OF(IJK)
                  
                  if(DES_ON_TMP) then 
                     
                 tmp_A =  - AVG_FACTOR*(DRAG_AM(I&
                      &,J,K,M) + DRAG_AM(I-1,J,K,M))
                 tmp_B =  - AVG_FACTOR*(DRAG_BM(I&
                      &,J,K,2,M) + DRAG_BM(I-1,J,K,2,M))
                 
                 IF(DIMN.eq.3) then 
                    tmp_A = tmp_A - AVG_FACTOR*(DRAG_AM(I&
                         &,J,K-1,M) + DRAG_AM(I-1,J,K-1,M))
                    tmp_B = tmp_B - AVG_FACTOR*(DRAG_BM(I&
                         &,J,K-1,2,M) + DRAG_BM(I-1,J,K-1,2,M))
                 end IF
                 
                 VCELL = VOL(IJK)

                                      
                 tmp_A = tmp_A*VCELL
                 tmp_B = tmp_B*VCELL
                 
                 else 
                    I = I_OF(IJK)
                    J = J_OF(IJK)
                    K = K_OF(IJK)
                    
                    VCELL = VOL(IJK)
                    VSFCM = AVG_Y(DES_V_S(IJK,M),DES_V_S(NORTH_OF(IJK),M),J_OF(IJK))
                    tmp_A =  - VXF_GS(IJK,M)
                    tmp_B =  - VXF_GS(IJK,M)*VSFCM
                    
                 end if
                 A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
                 B_M(IJK,0) = B_M(IJK,0) + tmp_B
                 
               END IF
             END DO
          END DO
90                  FORMAT(3(1x,i2),10(2x,E12.5))
          close(900)

      ELSE IF(WV.EQ.1) THEN
         DO M = 1, MMAX
            DO IJK = IJKSTART3, IJKEND3
               IF(FLUID_AT(IJK)) THEN
                  I = I_OF(IJK)
                  J = J_OF(IJK)
                  K = K_OF(IJK)
                  
                  if(DES_ON_TMP) then 
                     
                     tmp_A =  - 0.25d0*(DRAG_AM(I&
                          &,J,K,M) + DRAG_AM(I-1,J,K,M)+DRAG_AM(I,J-1,K&
                          &,M)+DRAG_AM(I-1,J-1,K,M))
                     tmp_B =  - 0.25d0*(DRAG_BM(I&
                          &,J,K,3,M) + DRAG_BM(I-1,J,K,3,M)+ DRAG_BM(I,J&
                          &-1,K,3,M)+ DRAG_BM(I-1,J-1,K,3,M))
                     
                     VCELL = VOL(IJK)
                                                           
                     tmp_A = tmp_A*VCELL
                     tmp_B = tmp_B*VCELL

                  else 
                     WSFCM = AVG_Z(DES_W_S(IJK,M),DES_W_S(TOP_OF(IJK),M),K_OF(IJK))
                     tmp_A = - VXF_GS(IJK,M)
                     tmp_B = - VXF_GS(IJK,M)*WSFCM
                     
                  end if
                  A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
                  B_M(IJK,0) = B_M(IJK,0) + tmp_B
                  
               END IF
            END DO
         END DO
         
      END IF
      !stop
      RETURN
      END SUBROUTINE GAS_DRAG
