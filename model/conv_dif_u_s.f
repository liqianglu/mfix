!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CONV_DIF_U_s(A_m, B_m, IER)                            C
!  Purpose: Determine convection diffusion terms for U_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative;              C
!  See source_u_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!  Reviewer:                                          Date:            C
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
      SUBROUTINE CONV_DIF_U_S(A_M, B_M, IER) 
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
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Solids phase index 
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
 
!-----------------------------------------------
!
!
      DO M = 1, MMAX 
        IF (MOMENTUM_X_EQ(M)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!	  IF DEFERRED CORRECTION IS USED TO SOLVE U_S
!
	  IF(DEF_COR)THEN
	    CALL STORE_A_U_S0 (A_M(1,-3,M), M, IER)
	    IF (DISCRETIZE(3) > 1) CALL STORE_A_U_SDC (A_M(1,-3,M), M, B_M, IER)
	  ELSE
!
!	  NO DEFERRED CORRECTION IS TO BE USED TO SOLVE FOR U_S
!  
            IF (DISCRETIZE(3) == 0) THEN         ! 0 & 1 => FOUP 
               CALL STORE_A_U_S0 (A_M(1,-3,M), M, IER) 
            ELSE 
               CALL STORE_A_U_S1 (A_M(1,-3,M), M, IER) 
            ENDIF 
!
          ENDIF
	  	  
	  CALL DIF_U_IS (MU_S(1,M), A_M, B_M, M, IER)
	ENDIF 
      END DO 
      
      RETURN  
      END SUBROUTINE CONV_DIF_U_S 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_U_s0(A_U_s, M, IER)                            C
!  Purpose: Determine convection diffusion terms for U_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; FOUP         C
!  See source_u_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-APR-96  C
!  Reviewer:                                          Date:            C
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
      SUBROUTINE STORE_A_U_S0(A_U_S, M, IER) 
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
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc 
      USE fldvar
      USE output
      USE compar 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
! 
!                      Indices 
      INTEGER          I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,& 
                       IJKNE, IJKP, IJKT, IJKTE 
      INTEGER          IMJK, IM, IJKW 
      INTEGER          IJMK, JM, IPJMK, IJKS, IJKSE 
      INTEGER          IJKM, KM, IPJKM, IJKB, IJKBE 
! 
!                      Solids phase 
      INTEGER          M 
! 
!                      Face velocity 
      DOUBLE PRECISION V_f 
! 
!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_U_s 
      DOUBLE PRECISION A_U_s(DIMENSION_3, -3:3, M:M) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!
!$omp      parallel do 	&
!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,	&
!$omp&                    IJKNE, IJKP, IJKT, IJKTE,  V_f, D_f,	&
!$omp&                    IMJK, IM, IJKW,	&
!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,	&
!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_E(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKE 
            ELSE 
               IJKC = IJK 
            ENDIF 
            IP = IP1(I) 
            IJKN = NORTH_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           East face (i+1, j, k)
            V_F = AVG_X_E(U_S(IJK,M),U_S(IPJK,M),IP) 
            D_F = MU_S(IJKE,M)*ODX(IP)*AYZ_U(IJK) 
            IF (V_F >= ZERO) THEN 
               A_U_S(IJK,E,M) = D_F 
               A_U_S(IPJK,W,M) = D_F + AVG_X(ROP_S(IJKC,M),ROP_S(IJKE,M),I)*V_F&
                  *AYZ_U(IJK) 
            ELSE 
               A_U_S(IJK,E,M) = D_F - AVG_X(ROP_S(IJKE,M),ROP_S(EAST_OF(IJKE),M&
                  ),IP)*V_F*AYZ_U(IJK) 
               A_U_S(IPJK,W,M) = D_F 
            ENDIF 
!
!           North face (i+1/2, j+1/2, k)
            V_F = AVG_X(V_S(IJK,M),V_S(IPJK,M),I) 
            D_F = AVG_X_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(&
               IJKE,M),MU_S(IJKNE,M),J),I)*ODY_N(J)*AXZ_U(IJK) 
            IF (V_F >= ZERO) THEN 
               A_U_S(IJK,N,M) = D_F 
               A_U_S(IJPK,S,M) = D_F + AVG_X(ROP_S(IJKC,M),ROP_S(IJKE,M),I)*V_F&
                  *AXZ_U(IJK) 
            ELSE 
               A_U_S(IJK,N,M) = D_F - AVG_X(ROP_S(IJKN,M),ROP_S(IJKNE,M),I)*V_F&
                  *AXZ_U(IJK) 
               A_U_S(IJPK,S,M) = D_F 
            ENDIF 
!
!           Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 
               V_F = AVG_X(W_S(IJK,M),W_S(IPJK,M),I) 
               D_F = AVG_X_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),AVG_Z_H(MU_S(&
                  IJKE,M),MU_S(IJKTE,M),K),I)*OX_E(I)*ODZ_T(K)*AXY_U(IJK) 
               IF (V_F >= ZERO) THEN 
                  A_U_S(IJK,T,M) = D_F 
                  A_U_S(IJKP,B,M) = D_F + AVG_X(ROP_S(IJKC,M),ROP_S(IJKE,M),I)*&
                     V_F*AXY_U(IJK) 
               ELSE 
                  A_U_S(IJK,T,M) = D_F - AVG_X(ROP_S(IJKT,M),ROP_S(IJKTE,M),I)*&
                     V_F*AXY_U(IJK) 
                  A_U_S(IJKP,B,M) = D_F 
               ENDIF 
            ENDIF 
!
!           West face (i, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
!
               V_F = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I) 
               D_F = MU_S(IJKC,M)*ODX(I)*AYZ_U(IMJK) 
               IF (V_F >= ZERO) THEN 
                  A_U_S(IJK,W,M) = D_F + AVG_X(ROP_S(IJKW,M),ROP_S(IJKC,M),IM)*&
                     V_F*AYZ_U(IMJK) 
               ELSE 
                  A_U_S(IJK,W,M) = D_F 
               ENDIF 
            ENDIF 
!
!           South face (i+1/2, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IJMK)) THEN 
               JM = JM1(J) 
               IPJMK = IP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKSE = EAST_OF(IJKS) 
               V_F = AVG_X(V_S(IJMK,M),V_S(IPJMK,M),I) 
               D_F = AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),AVG_Y_H(MU_S&
                  (IJKSE,M),MU_S(IJKE,M),JM),I)*ODY_N(JM)*AXZ_U(IJMK) 
               IF (V_F >= ZERO) THEN 
                  A_U_S(IJK,S,M) = D_F + AVG_X(ROP_S(IJKS,M),ROP_S(IJKSE,M),I)*&
                     V_F*AXZ_U(IJMK) 
               ELSE 
                  A_U_S(IJK,S,M) = D_F 
               ENDIF 
            ENDIF 
!
!           Bottom face (i+1/2, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLOW_AT_E(IJKM)) THEN 
                  KM = KM1(K) 
                  IPJKM = IP_OF(IJKM) 
                  IJKB = BOTTOM_OF(IJK) 
                  IJKBE = EAST_OF(IJKB) 
                  V_F = AVG_X(W_S(IJKM,M),W_S(IPJKM,M),I) 
                  D_F = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),AVG_Z_H(&
                     MU_S(IJKBE,M),MU_S(IJKE,M),KM),I)*OX_E(I)*ODZ_T(KM)*AXY_U(&
                     IJKM) 
                  IF (V_F >= ZERO) THEN 
                     A_U_S(IJK,B,M) = D_F + AVG_X(ROP_S(IJKB,M),ROP_S(IJKBE,M),&
                        I)*V_F*AXY_U(IJKM) 
                  ELSE 
                     A_U_S(IJK,B,M) = D_F 
                  ENDIF 
               ENDIF 
            ENDIF 
!
         ENDIF 
      END DO 
      
      RETURN  
      END SUBROUTINE STORE_A_U_S0 

!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_U_sdc(A_U_s, M, B_M, IER)                           C
!  Purpose:TO USE DEFERRED CORRECTION METHOD TO SOLVE THE U-MOMENTUM   C
!  EQUATION. THIS METHOD COMBINES FIRST ORDER UPWIND AND A USER        C
!  SPECIFIED HIGH ORDER METHOD                                         C
!                                                                      C
!  Author: C. GUENTHER                                 Date: 8-APR-99  C
!  Reviewer:                                          Date:            C
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
      SUBROUTINE STORE_A_U_SDC(A_U_S, M, B_M, IER) 
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
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc 
      USE fldvar
      USE output
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar  
      USE sendrecv
      USE sendrecv3
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
! 
!                      Indices 
      INTEGER          I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,& 
                       IJKNE, IJKP, IJKT, IJKTE 
      INTEGER          IMJK, IM, IJKW 
      INTEGER          IJMK, JM, IPJMK, IJKS, IJKSE 
      INTEGER          IJKM, KM, IPJKM, IJKB, IJKBE 
      INTEGER          IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER          IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
! 
!                      Solids phase 
      INTEGER          M 
! 
! loezos
	INTEGER  incr
! loezos

!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
!
!                      Septadiagonal matrix A_U_s 
      DOUBLE PRECISION A_U_s(DIMENSION_3, -3:3, M:M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
!
!	DEFERRED CORRCTION CONTRIBUTION FORM HIGH ORDER METHOD
	DOUBLE PRECISION MOM_HO
!
!	LOW ORDER APPROXIMATION 
	DOUBLE PRECISION MOM_LO
!
!	CONVECTION FACTOR AT THE FACE
	DOUBLE PRECISION CONV_FAC
!
!	DEFERRED CORRECTION CONTRIBUTIONS FROM EACH FACE
	DOUBLE PRECISION  EAST_DC
	DOUBLE PRECISION  WEST_DC
	DOUBLE PRECISION  NORTH_DC
	DOUBLE PRECISION  SOUTH_DC
        DOUBLE PRECISION  TOP_DC
        DOUBLE PRECISION  BOTTOM_DC
!
!
!
!---------------------------------------------------------------
!	EXTERNAL FUNCTIONS
!---------------------------------------------------------------
	DOUBLE PRECISION , EXTERNAL :: FPFOI_OF
!---------------------------------------------------------------
! 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'function3.inc'


      call lock_tmp4_array

      
      call lock_tmp_array
      call lock_xsi_array      
!
!  Calculate convection factors
!
!
! Send recv the third ghost layer
      IF ( FPFOI ) THEN
         Do IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = U_S(IJK,M)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

!
!$omp parallel do private(IJK,I,IP,IPJK,IJKE)
      DO IJK = ijkstart3, ijkend3
!
         I = I_OF(IJK) 
         IP = IP1(I) 
         IPJK = IP_OF(IJK) 
         IJKE = EAST_OF(IJK) 
!
!
!           East face (i+1, j, k)
         U(IJK) = AVG_X_E(U_S(IJK,M),U_S(IPJK,M),IP) 
!
!
!           North face (i+1/2, j+1/2, k)
         V(IJK) = AVG_X(V_S(IJK,M),V_S(IPJK,M),I) 
!
!
!           Top face (i+1/2, j, k+1/2)
         IF (DO_K) WW(IJK) = AVG_X(W_S(IJK,M),W_S(IPJK,M),I) 
      END DO 

! loezos
	incr=1		
! loezos

      CALL CALC_XSI (DISCRETIZE(3), U_S(1,M), U, V, WW, XSI_E, XSI_N,&
	XSI_T,incr) 
!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!$omp      parallel do 	&
!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,	&
!$omp&                    IJKNE, IJKP, IJKT, IJKTE,  D_f,	&
!$omp&                    IMJK, IM, IJKW,	&
!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,	&
!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE, &
!$omp&              MOM_HO, MOM_LO, CONV_FAC,EAST_DC,WEST_DC,NORTH_DC,&
!$omp&              SOUTH_DC, TOP_DC,BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3
!
         IF (FLOW_AT_E(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKE = EAST_OF(IJK) 
!
!           Third Ghost layer information
            IPPP  = IP_OF(IP_OF(IPJK))
            IPPP4 = funijk3(I_OF(IPPP), J_OF(IPPP), K_OF(IPPP))
            IMMM  = IM_OF(IM_OF(IMJK))
            IMMM4 = funijk3(I_OF(IMMM), J_OF(IMMM), K_OF(IMMM))
            JPPP  = JP_OF(JP_OF(IJPK))
            JPPP4 = funijk3(I_OF(JPPP), J_OF(JPPP), K_OF(JPPP))
            JMMM  = JM_OF(JM_OF(IJMK))
            JMMM4 = funijk3(I_OF(JMMM), J_OF(JMMM), K_OF(JMMM))
            KPPP  = KP_OF(KP_OF(IJKP))
            KPPP4 = funijk3(K_OF(IPPP), J_OF(KPPP), K_OF(KPPP))
            KMMM  = KM_OF(KM_OF(IJKM))
            KMMM4 = funijk3(I_OF(KMMM), J_OF(KMMM), K_OF(KMMM))
!
!
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKE 
            ELSE 
               IJKC = IJK 
            ENDIF 
            IP = IP1(I) 
            IJKN = NORTH_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE East face (i+1, j, k)
!           
		IF(U(IJK) >= ZERO)THEN
	           CONV_FAC = AVG_X(ROP_S(IJK,M),ROP_S(IJKE,M),I)&
		               *U(IJK)*AYZ_U(IJK) 
		   MOM_LO = U_S(IJK,M)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(U_S(IPJK,M), U_S(IJK,M), & 
                            U_S(IMJK,M), U_S(IM_OF(IMJK),M))
                     ELSE
                     ENDIF
		ELSE
		   CONV_FAC = AVG_X(ROP_S(IJKE,M),ROP_S(EAST_OF(IJKE),M),IP)&
                              *U(IJK)*AYZ_U(IJK)
		   MOM_LO = U_S(IPJK,M)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IPJK,M), & 
                            U_S(IP_OF(IPJK),M), TMP4(IPPP4))
                     ELSE
                     ENDIF
		ENDIF
                     IF (.NOT. FPFOI ) THEN
		      MOM_HO = XSI_E(IJK)*U_S(IPJK,M)+ &
                                (1.0-XSI_E(IJK))*U_S(IJK,M)
                     ELSE
                     ENDIF
		EAST_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE North face (i+1/2, j+1/2, k)
!            
		IF(V(IJK) >= ZERO)THEN
	            CONV_FAC = AVG_X(ROP_S(IJKC,M),ROP_S(IJKE,M),I)&
		              *V(IJK)*AXZ_U(IJK) 
		    MOM_LO = U_S(IJK,M)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(U_S(IJPK,M), U_S(IJK,M), & 
                            U_S(IJMK,M), U_S(JM_OF(IJMK),M))
                     ELSE
                     ENDIF
		ELSE
		    CONV_FAC = AVG_X(ROP_S(IJKN,M),ROP_S(IJKNE,M),I)&
		              *V(IJK)*AXZ_U(IJK)
		    MOM_LO = U_S(IJPK,M)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJPK,M), & 
                            U_S(JP_OF(IJPK),M), TMP4(JPPP4))
                     ELSE
                     ENDIF
		ENDIF
                     IF (.NOT. FPFOI ) THEN
		      MOM_HO = XSI_N(IJK)*U_S(IJPK,M)+ &
                              (1.0-XSI_N(IJK))*U_S(IJK,M)
                     ELSE
                     ENDIF
		NORTH_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Top face (i+1/2, j, k+1/2)
!
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 
	       IF(WW(IJK) >= ZERO)THEN
	          CONV_FAC = AVG_X(ROP_S(IJKC,M),ROP_S(IJKE,M),I)&
		             *WW(IJK)*AXY_U(IJK) 
		  MOM_LO = U_S(IJK,M)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(U_S(IJKP,M), U_S(IJK,M), & 
                            U_S(IJKM,M), U_S(KM_OF(IJKM),M))
                     ELSE
                     ENDIF
	       ELSE
		 CONV_FAC = AVG_X(ROP_S(IJKT,M),ROP_S(IJKTE,M),I)&
		           *WW(IJK)*AXY_U(IJK)
		 MOM_LO = U_S(IJKP,M)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJKP,M), & 
                            U_S(KP_OF(IJKP),M), TMP4(KPPP4))
                     ELSE
                     ENDIF
	       ENDIF
                     IF (.NOT. FPFOI ) THEN
	              MOM_HO = XSI_T(IJK)*U_S(IJKP,M)+ &
                                 (1.0-XSI_T(IJK))*U_S(IJK,M)
                     ELSE
                     ENDIF
	       TOP_DC = CONV_FAC*(MOM_LO-MOM_HO)
	    ELSE
	       TOP_DC = ZERO
	    
            ENDIF
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE West face (i, j, k)
!
            IMJK = IM_OF(IJK) 
            IM = IM1(I) 
            IJKW = WEST_OF(IJK) 
	    IF(U(IMJK) >= ZERO)THEN
	      CONV_FAC = AVG_X(ROP_S(IJKW,M),ROP_S(IJK,M),IM)&
	                 *U(IMJK)*AYZ_U(IMJK) 
	      MOM_LO = U_S(IMJK,M)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IMJK,M), & 
                            U_S(IM_OF(IMJK),M), TMP4(IMMM4))
                     ELSE
                     ENDIF	  
	    ELSE
	      CONV_FAC=AVG_X(ROP_S(IJK,M),ROP_S(IJKE,M),I)*U(IMJK)*AYZ_U(IMJK)
	      MOM_LO = U_S(IJK,M)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(U_S(IMJK,M), U_S(IJK,M), & 
                            U_S(IPJK,M), U_S(IP_OF(IPJK),M))
                     ELSE
                     ENDIF	
	    ENDIF
                     IF (.NOT. FPFOI ) THEN
	              MOM_HO = XSI_E(IMJK)*U_S(IJK,M)+ &
                               (1.0-XSI_E(IMJK))*U_S(IMJK,M)
                     ELSE
                     ENDIF
	    WEST_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE South face (i+1/2, j-1/2, k)
!
            IJMK = JM_OF(IJK) 
            JM = JM1(J) 
            IPJMK = IP_OF(IJMK) 
            IJKS = SOUTH_OF(IJK) 
            IJKSE = EAST_OF(IJKS) 
	    IF(V(IJMK) >= ZERO)THEN
	       CONV_FAC = AVG_X(ROP_S(IJKS,M),ROP_S(IJKSE,M),I)&
	                  *V(IJMK)*AXZ_U(IJMK) 
	       MOM_LO = U_S(IJMK,M)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJMK,M), & 
                            U_S(JM_OF(IJMK),M), TMP4(JMMM4))
                     ELSE
                     ENDIF
	    ELSE
	       CONV_FAC = AVG_X(ROP_S(IJK,M),ROP_S(IJKE,M),I)&
		           *V(IJMK)*AXZ_U(IJMK)
	       MOM_LO = U_S(IJK,M)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(U_S(IJMK,M), U_S(IJK,M), & 
                            U_S(IJPK,M), U_S(JP_OF(IJPK),M))
                     ELSE
                     ENDIF
	    ENDIF
                     IF (.NOT. FPFOI ) THEN
	              MOM_HO = XSI_N(IJMK)*U_S(IJK,M)+ &
                               (1.0-XSI_N(IJMK))*U_S(IJMK,M)
                     ELSE
                     ENDIF
	    SOUTH_DC = CONV_FAC*(MOM_LO-MOM_HO)
!
!           DEFERRED CORRECTION CONTRIBUTION AT THE Bottom face (i+1/2, j, k-1/2)
!
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               KM = KM1(K) 
               IPJKM = IP_OF(IJKM) 
               IJKB = BOTTOM_OF(IJK) 
               IJKBE = EAST_OF(IJKB) 
	       IF(WW(IJK) >= ZERO)THEN
	            CONV_FAC = AVG_X(ROP_S(IJKB,M),ROP_S(IJKBE,M),I)&
		               *WW(IJK)*AXY_U(IJKM) 
		    MOM_LO = U_S(IJKM,M)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJKM,M), & 
                            U_S(KM_OF(IJKM),M), TMP4(KMMM4))
                     ELSE
                     ENDIF
	       ELSE
		    CONV_FAC = AVG_X(ROP_S(IJK,M),ROP_S(IJKE,M),I)&
		               *WW(IJK)*AXY_U(IJKM)
		    MOM_LO = U_S(IJK,M)
                     IF ( FPFOI ) THEN
                      MOM_HO = FPFOI_OF(U_S(IJKM,M), U_S(IJK,M), & 
                            U_S(IJKP,M), U_S(KP_OF(IJKP),M))
                     ELSE
                     ENDIF
	       ENDIF
                     IF (.NOT. FPFOI ) THEN
	              MOM_HO = XSI_T(IJKM)*U_S(IJK,M)+ &
                               (1.0-XSI_T(IJKM))*U_S(IJKM,M)
                     ELSE
                     ENDIF
	       BOTTOM_DC = CONV_FAC*(MOM_LO-MOM_HO)
            ELSE
	       BOTTOM_DC = ZERO
            ENDIF
!
!		CONTRIBUTION DUE TO DEFERRED CORRECTION
!
		B_M(IJK,M) = B_M(IJK,M)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC&
				+BOTTOM_DC-TOP_DC
!
         ENDIF 
      END DO 
     
      call unlock_tmp4_array
      call unlock_tmp_array
      call unlock_xsi_array
      
      RETURN  
      END SUBROUTINE STORE_A_U_SDC 


!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_A_U_s1(A_U_s, M, IER)
!  Purpose: Determine convection diffusion terms for U_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive. The C
!  center coefficient and the source vector are negative; Higher order C
!  See source_u_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAR-97  C
!  Reviewer:                                          Date:            C
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
      SUBROUTINE STORE_A_U_S1(A_U_S, M, IER) 
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
      USE parallel 
      USE matrix 
      USE geometry
      USE indices
      USE run
      USE physprop
      USE visc_s
      USE toleranc 
      USE fldvar
      USE output
      USE vshear
      Use xsi_array
      Use tmp_array,  U => Array1, V => Array2, WW => Array3
      USE compar  
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
! 
!                      Indices 
      INTEGER          I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,& 
                       IJKNE, IJKP, IJKT, IJKTE 
      INTEGER          IMJK, IM, IJKW 
      INTEGER          IJMK, JM, IPJMK, IJKS, IJKSE 
      INTEGER          IJKM, KM, IPJKM, IJKB, IJKBE 
! 
!                      Solids phase 
      INTEGER          M 
! 
! loezos                     
      INTEGER incr    
!loezos

!                      Diffusion parameter 
      DOUBLE PRECISION D_f 
! 
!                      Septadiagonal matrix A_U_s 
      DOUBLE PRECISION A_U_s(DIMENSION_3, -3:3, M:M) 
! 
!                      Convection weighting factors 
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),& 
!                       XSI_t(DIMENSION_3) 
!      DOUBLE PRECISION U(DIMENSION_3),& 
!                       V(DIMENSION_3), WW(DIMENSION_3) 
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'


      call lock_tmp_array
      call lock_xsi_array

!
!  Calculate convection factors
!
!$omp parallel do private(IJK,I,IP,IPJK,IJKE)
      DO IJK = ijkstart3, ijkend3
!
         I = I_OF(IJK) 
         IP = IP1(I) 
         IPJK = IP_OF(IJK) 
         IJKE = EAST_OF(IJK) 


!
!
!           East face (i+1, j, k)
         U(IJK) = AVG_X_E(U_S(IJK,M),U_S(IPJK,M),IP) 
!
!
!           North face (i+1/2, j+1/2, k)
         V(IJK) = AVG_X(V_S(IJK,M),V_S(IPJK,M),I) 
!
!
!           Top face (i+1/2, j, k+1/2)
         IF (DO_K) WW(IJK) = AVG_X(W_S(IJK,M),W_S(IPJK,M),I) 
      END DO 

! loezos
	incr=1		
! loezos

      CALL CALC_XSI (DISCRETIZE(3), U_S(1,M), U, V, WW, XSI_E, XSI_N, XSI_T,&
	incr) 

! loezos      
! update to true velocity
      IF (SHEAR) THEN
!$omp  parallel do private(IJK)
	 DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN  
	   V(IJK)=V(IJK)+VSHE(IJK)		
          END IF
        END DO
      END IF
! loezos

!
!
!  Calculate convection-diffusion fluxes through each of the faces
!
!$omp      parallel do 	&
!$omp&     private(I,  J, K, IP, IJK, IJKC, IPJK, IJPK, IJKE, IJKN,	&
!$omp&                    IJKNE, IJKP, IJKT, IJKTE,   D_f,	&
!$omp&                    IMJK, IM, IJKW,	&
!$omp&                    IJMK, JM, IPJMK, IJKS, IJKSE,	&
!$omp&                    IJKM, KM, IPJKM, IJKB, IJKBE)
      DO IJK = ijkstart3, ijkend3 
!
         IF (FLOW_AT_E(IJK)) THEN 
            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK) 
            IPJK = IP_OF(IJK) 
            IJPK = JP_OF(IJK) 
            IJKE = EAST_OF(IJK) 
            IF (WALL_AT(IJK)) THEN 
               IJKC = IJKE 
            ELSE 
               IJKC = IJK 
            ENDIF 
            IP = IP1(I) 
            IJKN = NORTH_OF(IJK) 
            IJKNE = EAST_OF(IJKN) 
!
!           East face (i+1, j, k)
            D_F = MU_S(IJKE,M)*ODX(IP)*AYZ_U(IJK) 
!
            A_U_S(IJK,E,M) = D_F - XSI_E(IJK)*AVG_X(ROP_S(IJKE,M),ROP_S(EAST_OF&
               (IJKE),M),IP)*U(IJK)*AYZ_U(IJK) 
!
            A_U_S(IPJK,W,M) = D_F + (ONE - XSI_E(IJK))*AVG_X(ROP_S(IJKC,M),&
               ROP_S(IJKE,M),I)*U(IJK)*AYZ_U(IJK) 
!
!
!           North face (i+1/2, j+1/2, k)
            D_F = AVG_X_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),AVG_Y_H(MU_S(&
               IJKE,M),MU_S(IJKNE,M),J),I)*ODY_N(J)*AXZ_U(IJK) 
!
            A_U_S(IJK,N,M) = D_F - XSI_N(IJK)*AVG_X(ROP_S(IJKN,M),ROP_S(IJKNE,M&
               ),I)*V(IJK)*AXZ_U(IJK) 
!
            A_U_S(IJPK,S,M) = D_F + (ONE - XSI_N(IJK))*AVG_X(ROP_S(IJKC,M),&
               ROP_S(IJKE,M),I)*V(IJK)*AXZ_U(IJK) 
!
!
!           Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN 
               IJKP = KP_OF(IJK) 
               IJKT = TOP_OF(IJK) 
               IJKTE = EAST_OF(IJKT) 
!
               D_F = AVG_X_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),AVG_Z_H(MU_S(&
                  IJKE,M),MU_S(IJKTE,M),K),I)*OX_E(I)*ODZ_T(K)*AXY_U(IJK) 
!
               A_U_S(IJK,T,M) = D_F - XSI_T(IJK)*AVG_X(ROP_S(IJKT,M),ROP_S(&
                  IJKTE,M),I)*WW(IJK)*AXY_U(IJK) 
!
               A_U_S(IJKP,B,M) = D_F + (ONE - XSI_T(IJK))*AVG_X(ROP_S(IJKC,M),&
                  ROP_S(IJKE,M),I)*WW(IJK)*AXY_U(IJK) 
            ENDIF 
!
!           West face (i, j, k)
            IMJK = IM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IMJK)) THEN 
               IM = IM1(I) 
               IJKW = WEST_OF(IJK) 
!
               D_F = MU_S(IJKC,M)*ODX(I)*AYZ_U(IMJK) 
!
               A_U_S(IJK,W,M) = D_F + (ONE - XSI_E(IMJK))*AVG_X(ROP_S(IJKW,M),&
                  ROP_S(IJKC,M),IM)*U(IMJK)*AYZ_U(IMJK) 
            ENDIF 
!
!           South face (i+1/2, j-1/2, k)
            IJMK = JM_OF(IJK) 
            IF (.NOT.FLOW_AT_E(IJMK)) THEN 
               JM = JM1(J) 
               IPJMK = IP_OF(IJMK) 
               IJKS = SOUTH_OF(IJK) 
               IJKSE = EAST_OF(IJKS) 
!
               D_F = AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),AVG_Y_H(MU_S&
                  (IJKSE,M),MU_S(IJKE,M),JM),I)*ODY_N(JM)*AXZ_U(IJMK) 
!
               A_U_S(IJK,S,M) = D_F + (ONE - XSI_N(IJMK))*AVG_X(ROP_S(IJKS,M),&
                  ROP_S(IJKSE,M),I)*V(IJMK)*AXZ_U(IJMK) 
            ENDIF 
!
!           Bottom face (i+1/2, j, k-1/2)
            IF (DO_K) THEN 
               IJKM = KM_OF(IJK) 
               IF (.NOT.FLOW_AT_E(IJKM)) THEN 
                  KM = KM1(K) 
                  IPJKM = IP_OF(IJKM) 
                  IJKB = BOTTOM_OF(IJK) 
                  IJKBE = EAST_OF(IJKB) 
!
                  D_F = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),AVG_Z_H(&
                     MU_S(IJKBE,M),MU_S(IJKE,M),KM),I)*OX_E(I)*ODZ_T(KM)*AXY_U(&
                     IJKM) 
!
                  A_U_S(IJK,B,M) = D_F + (ONE - XSI_T(IJKM))*AVG_X(ROP_S(IJKB,M&
                     ),ROP_S(IJKBE,M),I)*WW(IJKM)*AXY_U(IJKM) 
               ENDIF 
            ENDIF 
!
         ENDIF 
      END DO 
      
      call unlock_tmp_array
      call unlock_xsi_array

      
      RETURN  
      END SUBROUTINE STORE_A_U_S1 
!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3

