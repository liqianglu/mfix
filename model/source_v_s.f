!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_V_s(A_m, B_m, IER)                              C
!  Purpose: Determine source terms for V_s momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.    The center      C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-JUN-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Allow for partial-slip boundary conditions proposed by     C
!           by Johnson & Jackson (1987) if the Granular Temperature    C
!           equation is used.                                          C
!  Author: K. Agrawal, Princeton University           Date: 24-JAN-98  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SOURCE_V_S(A_M, B_M, IER) 
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
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is
      USE tau_s
      USE bc
      USE vshear
      USE compar 
      USE sendrecv 
      use kintheory
      USE ghdtheory
      USE drag
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE cutcell
      USE quadric
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! 
! 
!                      Error index 
      INTEGER          IER 
! 
!                      Indices 
      INTEGER          I, J, K, IJK,IMJK, IJMK, IJKM, IJKN
! 
!                      Phase index 
      INTEGER          M, MM, L
      DOUBLE PRECISION   SUM_EPS_CP 
! 
!                      Internal surface 
      INTEGER          ISV 
! 
!                      Pressure at north cell 
      DOUBLE PRECISION PgN 
! 
!                      Average volume fraction 
      DOUBLE PRECISION EPSA, EPStmp, epse, epsw, epsn, epss, &
                       epst, epsb, epsMix, epsMixN
! 
!                      Average density 
      DOUBLE PRECISION ROPSA 
! 
!                      Average density difference 
      DOUBLE PRECISION dro1, dro2, droa 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      Source terms (Surface) 
      DOUBLE PRECISION Sdp, Sdps 
! 
!                      Source terms (Volumetric) 
      DOUBLE PRECISION V0, Vmt, Vbf, Vmttmp 
!
!                      Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION Ghd_drag, avgRop

!			Source terms for HYS drag relation
      DOUBLE PRECISION HYS_drag, avgDrag
!
! loezos
      DOUBLE PRECISION VSH_n,VSH_s,VSH_e,VSH_w,VSH_p,Source_conv
      DOUBLE PRECISION SRT
! loezos
 
!                      error message 
      CHARACTER*80     LINE(2) 
!
!                      FOR CALL_DI and CALL_ISAT = .true.
      DOUBLE PRECISION SUM_R_S_temp(DIMENSION_3, DIMENSION_M)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      INTEGER ::          IM,JM,IP,JP,KM,KP
      INTEGER ::          IPJK,IJPK,IJKP,IJKC,IJKE,IJKNE,IJKW,IJKWN,IMJPK
      INTEGER ::          IJKT,IJKTN,IJKB,IJKBN
      DOUBLE PRECISION :: Vn,Vs,Ve,Vw, Vt,Vb
      DOUBLE PRECISION :: B_NOC
      DOUBLE PRECISION :: MU_S_E,MU_S_W,MU_S_N,MU_S_S,MU_S_T,MU_S_B,MU_S_CUT
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!-----------------------------------------------


      DO M = 1, MMAX 
        IF(TRIM(KT_TYPE) /= 'GHD' .OR. (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN

          IF (MOMENTUM_Y_EQ(M)) THEN 

! CHEM & ISAT (nan xie)
! Set the source terms zero
            IF (CALL_DI .or. CALL_ISAT) THEN
               SUM_R_S_temp = SUM_R_S
               SUM_R_S = ZERO
            ENDIF

!$omp  parallel do private( I, J, K, IJK, IJKN, ISV, Sdp, Sdps, V0, Vmt, &
!$omp&  PGN,DRO1,DRO2,DROA, Vbf, ROPSA, EPSA, EPStmp, VSH_n,VSH_s,VSH_e,&
!$omp&  VSH_w,VSH_p,Source_conv, SRT,SUM_EPS_CP,MM) &
!$omp&  schedule(static)
            DO IJK = ijkstart3, ijkend3
                I = I_OF(IJK) 
                J = J_OF(IJK) 
                K = K_OF(IJK)
                IMJK = IM_OF(IJK)
                IJMK = JM_OF(IJK)
                IJKM = KM_OF(IJK) 
                IJKN = NORTH_OF(IJK)
                IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                  EPStmp = ZERO     
                  epsMix = ZERO
                  epsMixN= ZERO  
                  DO L = 1, SMAX
                    EPStmp = EPStmp + AVG_Y(EP_S(IJK,L),EP_S(IJKN,L),J) 
                    epsMix  = epsMix  + EP_S(IJK,L) ! epsMix, epsMixN to be used for modelB
                    epsMixN = epsMixN + EP_S(IJKN,L)
                  ENDDO                        
                  EPSA = EPStmp
                ELSE                  
                  EPSA = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J) 
                ENDIF 
                IF (IP_AT_N(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
                ELSEIF (SIP_AT_N(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  ISV = IS_ID_AT_N(IJK) 
                  B_M(IJK,M) = -IS_VEL_S(ISV,M) 

! dilute flow
                ELSEIF (EPSA <= DIL_EP_S) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
                  IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                      EPSw = ZERO
                      EPSe = ZERO
                      EPSn = ZERO
                      EPSs = ZERO
                      EPSt = ZERO
                      EPSb = ZERO
                      DO L = 1, SMAX
                        EPSw = EPSw + EP_S(WEST_OF(IJK),L)
                        EPSe = EPSe + EP_S(EAST_OF(IJK),L)
                        EPSn = EPSn + EP_S(NORTH_OF(IJK),L)
                        EPSs = EPSs + EP_S(SOUTH_OF(IJK),L)
                        IF(.NOT. NO_K) THEN
                          EPSt = EPSt + EP_S(TOP_OF(IJK),L)
                          EPSb = EPSb + EP_S(BOTTOM_OF(IJK),L)
                        ENDIF
                      ENDDO
                  ELSE
                      EPSw = EP_S(WEST_OF(IJK),M)
                      EPSe = EP_S(EAST_OF(IJK),M)
                      EPSn = EP_S(NORTH_OF(IJK),M)
                      EPSs = EP_S(SOUTH_OF(IJK),M)
                      IF(.NOT. NO_K) THEN
                        EPSt = EP_S(TOP_OF(IJK),M)
                        EPSb = EP_S(BOTTOM_OF(IJK),M)
                      ENDIF
                  ENDIF                  
! using the average boundary cell values to compute V_s (sof, Aug 23 2005)
                  IF (EPSw > DIL_EP_S .AND. .NOT.IS_AT_E(IMJK)) A_M(IJK,W,M) = ONE 
                  IF (EPSe > DIL_EP_S .AND. .NOT.IS_AT_E(IJK)) A_M(IJK,E,M) = ONE 
                  IF (EPSs > DIL_EP_S .AND. .NOT.IS_AT_N(IJMK)) A_M(IJK,S,M) = ONE 
                  IF (EPSn > DIL_EP_S .AND. .NOT.IS_AT_N(IJK)) A_M(IJK,N,M) = ONE
                  IF(.NOT. NO_K) THEN
                    IF (EPSb > DIL_EP_S .AND. .NOT.IS_AT_T(IJKM)) A_M(IJK,B,M) = ONE 
                    IF (EPSt > DIL_EP_S .AND. .NOT.IS_AT_T(IJK)) A_M(IJK,T,M) = ONE 
                  ENDIF
               
                  IF((A_M(IJK,W,M)+A_M(IJK,E,M)+A_M(IJK,S,M)+A_M(IJK,N,M)+ &
                    A_M(IJK,B,M)+A_M(IJK,T,M)) == ZERO) THEN
                    B_M(IJK,M) = -V_S(IJK,M)           
                  ELSE
                    A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+ &
                                     A_M(IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M))
                  ENDIF
! Normal case
                ELSE 

! Surface forces

! Pressure term
                  PGN = P_G(IJKN) 
                  IF (CYCLIC_Y_PD) THEN 
                    IF (CYCLIC_AT_N(IJK)) PGN = P_G(IJKN) - DELP_Y 
                  ENDIF 

                  IF (MODEL_B) THEN 
                    SDP = ZERO 
                  ELSE 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                     IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
                        SDP = -P_SCALE*EPSA*(PGN - P_G(IJK))*AXZ(IJK) 
                     ELSE
                        SDP = -P_SCALE*EPSA*(PGN * A_VPG_N(IJK)  - P_G(IJK) * A_VPG_S(IJK) )
                     ENDIF
! Original terms
!                     SDP = -P_SCALE*EPGA*(PGN - P_G(IJK))*AXZ(IJK) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

                  ENDIF 

                  IF (CLOSE_PACKED(M)) THEN
                    IF(SMAX > 1 .AND. TRIM(KT_TYPE) /= 'GHD') THEN
                      SUM_EPS_CP=0.0 
                      DO MM=1,SMAX
                        IF (CLOSE_PACKED(MM))&
                          SUM_EPS_CP=SUM_EPS_CP+AVG_Y(EP_S(IJK,MM),EP_S(IJKN,MM),J)
                      ENDDO
                       SUM_EPS_CP = Max(SUM_EPS_CP, small_number)
                       SDPS = - ((P_S(IJKN,M)-P_S(IJK,M))+(EPSA/SUM_EPS_CP)* &
                           (P_STAR(IJKN)-P_STAR(IJK)))*AXZ(IJK) 
                    ELSE
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                     IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
                        SDPS = - ((P_S(IJKN,M)-P_S(IJK,M))+(P_STAR(IJKN)-P_STAR(IJK)))*AXZ(IJK)
                     ELSE
                        SDPS = - ((P_S(IJKN,M)* A_VPG_N(IJK)-P_S(IJK,M)* A_VPG_S(IJK))+(P_STAR(IJKN)* A_VPG_N(IJK)-P_STAR(IJK)* A_VPG_S(IJK)))
                     ENDIF
! Original terms
!                     SDPS = - ((P_S(IJKN,M)-P_S(IJK,M))+(P_STAR(IJKN)-P_STAR(IJK)))*AXZ(IJK)
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                    ENDIF
                  ELSE 
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                     IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
                        SDPS = -(P_S(IJKN,M)-P_S(IJK,M))*AXZ(IJK) 
                     ELSE
                        SDPS = -(P_S(IJKN,M) * A_VPG_N(IJK)-P_S(IJK,M) * A_VPG_S(IJK))
                     ENDIF
! Original terms
!                     SDPS = -(P_S(IJKN,M)-P_S(IJK,M))*AXZ(IJK) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

                  ENDIF 

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN
! Volumetric forces
                     ROPSA = AVG_Y(ROP_S(IJK,M),ROP_S(IJKN,M),J)
! Previous time step
                     V0 = AVG_Y(ROP_SO(IJK,M),ROP_SO(IJKN,M),J)*ODT 
                  ELSE
! Volumetric forces
                     ROPSA =  (VOL(IJK)*ROP_S(IJK,M) + VOL(IJKN)*ROP_S(IJKN,M))/(VOL(IJK) + VOL(IJKN))
! Previous time step
                     V0 = (VOL(IJK)*ROP_SO(IJK,M) + VOL(IJKN)*ROP_SO(IJKN,M))*ODT/(VOL(IJK) + VOL(IJKN))
                  ENDIF
! Original terms
! Volumetric forces
!                  ROPSA = AVG_Y(ROP_S(IJK,M),ROP_S(IJKN,M),J)
! Previous time step
!                  V0 = AVG_Y(ROP_SO(IJK,M),ROP_SO(IJKN,M),J)*ODT 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

! Interphase mass transfer
                  IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                    VMTtmp = ZERO
                    DO L = 1,SMAX
                      VMTtmp = VMTtmp + AVG_Y(SUM_R_S(IJK,L),SUM_R_S(IJKN,L),J) 
                    ENDDO
                    VMT = VMTtmp
                  ELSE
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
                  IF(.NOT.CUT_V_TREATMENT_AT(IJK)) THEN   
                     VMT = AVG_Y(SUM_R_S(IJK,M),SUM_R_S(IJKN,M),J) 
                  ELSE
                     VMT = (VOL(IJK)*SUM_R_S(IJK,M) + VOL(IJKN)*SUM_R_S(IJKN,M))/(VOL(IJK) + VOL(IJKN))
                  ENDIF
! Original terms
!                  VMT = AVG_Y(SUM_R_S(IJK,M),SUM_R_S(IJKN,M),J) 
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================

                  ENDIF

! Body force
                  IF (MODEL_B) THEN 
                    IF (TRIM(KT_TYPE) /= 'GHD') THEN
                      DRO1 = (RO_S(M)-RO_G(IJK))*EP_S(IJK,M) 
                      DRO2 = (RO_S(M)-RO_G(IJKN))*EP_S(IJKN,M) 
                      DROA = AVG_Y(DRO1,DRO2,J) 
                      VBF = DROA*BFY_S(IJK,M) 
                    ELSE ! GHD and M = MMAX
                      DRO1 = ROP_S(IJK,M)  - RO_G(IJK) *epsMix
                      DRO2 = ROP_S(IJKN,M) - RO_G(IJKN)*epsMixN 
                      DROA = AVG_Y(DRO1,DRO2,J) 
                      VBF = DROA*BFY_S(IJK,M) 
                    ENDIF
                  ELSE 
                    VBF = ROPSA*BFY_S(IJK,M) 
                  ENDIF 

! Additional force for GHD from darg force sum(beta_ig * Joi/rhop_i)
                  Ghd_drag = ZERO
                  IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                    DO L = 1,SMAX
                      avgRop = AVG_Y(ROP_S(IJK,L),ROP_S(IJKN,L),J)
                      if(avgRop > ZERO) Ghd_drag = Ghd_drag -&
                           AVG_Y(F_GS(IJK,L),F_GS(IJKN,L),J) * JoiY(IJK,L) / avgRop
                    ENDDO
                  ENDIF
! end of modifications for GHD theory

! Additional force for HYS drag force, do not use with mixture GHD theory
                  HYS_drag = ZERO
                  avgDrag = ZERO
                  IF (TRIM(DRAG_TYPE) .EQ. 'HYS' .AND. TRIM(KT_TYPE) /= 'GHD') THEN
                     DO L = 1,MMAX
                        IF (L /= M) THEN
                           avgDrag = AVG_Y(beta_ij(IJK,M,L),beta_ij(IJKN,M,L),J)
                           HYS_drag = HYS_drag - avgDrag * (V_g(ijk) - V_s(IJK,L))
                        ENDIF
                     ENDDO
                  ENDIF
! end of modifications for HYS drag

! loezos, Source terms from convective mom. flux
                  IF (SHEAR) THEN
                    SRT=(2d0*V_sh/XLENGTH)        
                    VSH_p=VSH(IJK)
                    VSH_n=VSH_p
                    VSH_s=VSH_p
                    VSH_e=VSH(IJK)+SRT*1d0/oDX_E(I)
                    VSH_w=VSH(IJK)-SRT*1d0/oDX_E(IM1(I))
                    Source_conv=A_M(IJK,N,m)*VSH_n+A_M(IJK,S,m)*VSH_s&
                      +A_M(IJK,W,m)*VSH_w+A_M(IJK,E,m)*VSH_e&
                      -(A_M(IJK,N,m)+A_M(IJK,S,m)+A_M(IJK,W,m)+A_M(IJK,E,m))&
                      *VSH_p
                  ELSE 
                    Source_conv=0d0
                  END IF

! Collect the terms
                  A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+A_M(&
                     IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M)+(V0+ZMAX(VMT))*VOL_V(&
                     IJK)) 

                  IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN 
                    B_M(IJK,M) = -(SDP + KTMOM_V_S(IJK,M) + SDPS + TAU_V_S(IJK,M)&
                        +Source_conv+((V0+ZMAX((-VMT)))&
                        *V_SO(IJK,M)+VBF+HYS_drag)*VOL_V(IJK))+B_M(IJK,M)

                  ELSE
                    B_M(IJK,M) = -(SDP + SDPS + TAU_V_S(IJK,M)&
                       +Source_conv+((V0+ZMAX((-VMT)))&
                        *V_SO(IJK,M)+VBF+Ghd_drag+HYS_drag)*VOL_V(IJK))+B_M(IJK,M) 
                  ENDIF
                ENDIF   ! end if sip or ip or dilute flow branch
            ENDDO

!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(CARTESIAN_GRID) CALL CG_SOURCE_V_S(A_M, B_M, M, IER)

            CALL SOURCE_V_S_BC (A_M, B_M, M, IER) 

            IF(CARTESIAN_GRID) CALL CG_SOURCE_V_S_BC(A_M, B_M, M, IER)
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================


! CHEM & ISAT (nan xie)
            IF (CALL_DI .or. CALL_ISAT) THEN
              SUM_R_S = SUM_R_S_temp
            ENDIF
    
          ENDIF  
        ENDIF ! for GHD Theory
      ENDDO 
   
      RETURN  
      END SUBROUTINE SOURCE_V_S 
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_V_s_BC(A_m, B_m, M, IER)                        C
!  Purpose: Determine source terms for V_s momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.    The center      C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this                 C
!  stage.                                                              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-JUN-96   C
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
      SUBROUTINE SOURCE_V_S_BC(A_M, B_M, M, IER) 
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
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns 
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is 
      USE tau_s 
      USE bc
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
!                      Error index 
      INTEGER          IER 
! 
!                      Boundary condition 
      INTEGER          L 
! 
!                      Indices 
      INTEGER          I,  J, K, JM, I1, I2, J1, J2, K1, K2, IJK,& 
                       IM, KM, IJKS, IJMK, IJPK 
! 
!                      Solids phase 
      INTEGER          M 
! 
!                      Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
! 
!                      Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! 
!                      coefficient for granular bc 
      DOUBLE PRECISION hw 
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'
!
!
!  Set the default boundary conditions
!
      IF (DO_K) THEN 
         K1 = 1 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE		 
               IJK = FUNIJK(I1,J1,K1) 
               IF (NS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = -ONE 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSE IF (FS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ONE 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            END DO 
         END DO 
         K1 = KMAX2 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3 
   	       IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE	 
               IJK = FUNIJK(I1,J1,K1) 
               IF (NS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = -ONE 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSE IF (FS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ONE 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            END DO 
         END DO 
      ENDIF 
!
      I1 = 1
      DO K1 = kmin3, kmax3 
         DO J1 = jmin3, jmax3 
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE       
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = -ONE 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSE IF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ONE 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      I1 = IMAX2 
      DO K1 = kmin3, kmax3 
         DO J1 = jmin3, jmax3 
   	    IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE      
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = -ONE 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSE IF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ONE 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         END DO 
      END DO 
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 
            IF (BC_TYPE(L) == 'NO_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
!
               IF (BC_JJ_PS(L) == 0) THEN 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     		     			
                           IJK = FUNIJK(I,J,K) 
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ZERO 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                           IF (FLUID_AT(EAST_OF(IJK))) THEN 
                              A_M(IJK,E,M) = -ONE 
                           ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN 
                              A_M(IJK,W,M) = -ONE 
                           ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                              A_M(IJK,T,M) = -ONE 
                           ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                              A_M(IJK,B,M) = -ONE 
                           ENDIF 
                        END DO 
                     END DO 
                  END DO 
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_V_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!
               ENDIF 
            ELSE IF (BC_TYPE(L) == 'FREE_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
!
               IF (BC_JJ_PS(L) == 0) THEN 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     		     			
			
                           IJK = FUNIJK(I,J,K) 
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ZERO 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                           IF (FLUID_AT(EAST_OF(IJK))) THEN 
                              A_M(IJK,E,M) = ONE 
                           ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN 
                              A_M(IJK,W,M) = ONE 
                           ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                              A_M(IJK,T,M) = ONE 
                           ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                              A_M(IJK,B,M) = ONE 
                           ENDIF 
                        END DO 
                     END DO 
                  END DO 
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_V_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!
               ENDIF 
            ELSE IF (BC_TYPE(L) == 'PAR_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
!
               IF (BC_JJ_PS(L) == 0) THEN 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     		     			
			
                           IJK = FUNIJK(I,J,K) 
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                           IM = IM1(I) 
                           KM = KM1(K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ZERO 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                           IF (FLUID_AT(EAST_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,E,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_VW_S(L,M) 
                              ELSE 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODX_E(I)) 
                                 A_M(IJK,E,M) = -(HALF*BC_HW_S(L,M)-ODX_E(I)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_VW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,W,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_VW_S(L,M) 
                              ELSE 
                                 A_M(IJK,W,M) = -(HALF*BC_HW_S(L,M)-ODX_E(IM)) 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODX_E(IM)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_VW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,T,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_VW_S(L,M) 
                              ELSE 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODZ_T(K)*OX&
                                    (I)) 
                                 A_M(IJK,T,M) = -(HALF*BC_HW_S(L,M)-ODZ_T(K)*OX&
                                    (I)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_VW_S(L,M) 
                              ENDIF 
                           ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,B,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_VW_S(L,M) 
                              ELSE 
                                 A_M(IJK,B,M) = -(HALF*BC_HW_S(L,M)-ODZ_T(KM)*&
                                    OX(I)) 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODZ_T(KM)*&
                                    OX(I)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_VW_S(L,M) 
                              ENDIF 
                           ENDIF 
                        END DO 
                     END DO 
                  END DO 
               ELSE                              !Johnson and Jackson partial slip 
!
                  CALL JJ_BC_V_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
!
               ENDIF 
            ELSE IF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'S') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     		     			
			
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ZERO 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ONE 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
               ENDIF 
            ELSE IF (BC_TYPE(L) == 'OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'S') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     		     			
			
                           IJK = FUNIJK(I,J,K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ZERO 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ONE 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
!
                           IJMK = JM_OF(IJK) 
                           A_M(IJMK,E,M) = ZERO 
                           A_M(IJMK,W,M) = ZERO 
                           A_M(IJMK,N,M) = ZERO 
                           A_M(IJMK,S,M) = ONE 
                           A_M(IJMK,T,M) = ZERO 
                           A_M(IJMK,B,M) = ZERO 
                           A_M(IJMK,0,M) = -ONE 
                           B_M(IJMK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
               ELSE IF (BC_PLANE(L) == 'N') THEN 
                  I1 = BC_I_W(L) 
                  I2 = BC_I_E(L) 
                  J1 = BC_J_S(L) 
                  J2 = BC_J_N(L) 
                  K1 = BC_K_B(L) 
                  K2 = BC_K_T(L) 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     		     			
			
                           IJK = FUNIJK(I,J,K) 
!
                           IJPK = JP_OF(IJK) 
                           A_M(IJPK,E,M) = ZERO 
                           A_M(IJPK,W,M) = ZERO 
                           A_M(IJPK,N,M) = ONE 
                           A_M(IJPK,S,M) = ZERO 
                           A_M(IJPK,T,M) = ZERO 
                           A_M(IJPK,B,M) = ZERO 
                           A_M(IJPK,0,M) = -ONE 
                           B_M(IJPK,M) = ZERO 
                        END DO 
                     END DO 
                  END DO 
               ENDIF 
            ELSE 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L) 
               DO K = K1, K2 
                  DO J = J1, J2 
                     DO I = I1, I2 
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     		     			
		     
                        IJK = FUNIJK(I,J,K) 
                        A_M(IJK,E,M) = ZERO 
                        A_M(IJK,W,M) = ZERO 
                        A_M(IJK,N,M) = ZERO 
                        A_M(IJK,S,M) = ZERO 
                        A_M(IJK,T,M) = ZERO 
                        A_M(IJK,B,M) = ZERO 
                        A_M(IJK,0,M) = -ONE 
                        B_M(IJK,M) = -V_S(IJK,M) 
                        IF (BC_PLANE(L) == 'S') THEN 
                           IJKS = SOUTH_OF(IJK) 
                           A_M(IJKS,E,M) = ZERO 
                           A_M(IJKS,W,M) = ZERO 
                           A_M(IJKS,N,M) = ZERO 
                           A_M(IJKS,S,M) = ZERO 
                           A_M(IJKS,T,M) = ZERO 
                           A_M(IJKS,B,M) = ZERO 
                           A_M(IJKS,0,M) = -ONE 
                           B_M(IJKS,M) = -V_S(IJKS,M) 
                        ENDIF 
                     END DO 
                  END DO 
               END DO 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE SOURCE_V_S_BC 
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: JJ_BC_V_s(I1, I2, J1, J2, K1, K2, L, M, A_m, b_m)      C
!  Purpose: Implement Johnson and Jackson boundary condition           C
!                                                                      C
!  Author: K. Agrawal, A. Srivastava,                 Date: 14-APR-98  C
!          Princeton University                                        C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!  Modified: S. Benyahia, Fluent Inc.                 Date: 02-FEB-05  C
!      Added the argument L to calc_grbdry                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE JJ_BC_V_S(I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
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
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s 
      USE rxns 
      USE run
      USE toleranc 
      USE geometry
      USE indices
      USE is 
      USE tau_s 
      USE bc
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
!                      Boundary condition
      INTEGER          L
!
!                      Indices
      INTEGER          I,  J, K, I1, I2, J1, J2, K1, K2, IJK, &
                      IM, KM, IJPK
!
!                      Solids phase
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      coefficient for granular bc
      DOUBLE PRECISION hw, gw, cw
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!
      DO K = K1, K2 
         DO J = J1, J2 
            DO I = I1, I2 
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE		     
	    
               IJK = FUNIJK(I,J,K) 
               IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
               IM = IM1(I) 
               KM = KM1(K) 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
               IF (FLUID_AT(EAST_OF(IJK))) THEN 
	          IJPK = JP_OF(EAST_OF(IJK))
                  IF (WALL_AT(IJPK)) CYCLE
                  IF (EP_S(EAST_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,E,M) = ONE 
                  ELSE 
! start anuj 04/20
		     IF (FRICTION .AND. (ONE-EP_G(EAST_OF(IJK)))>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, EAST_OF(IJK), 'E', 'V', L, M&
                           , GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, EAST_OF(IJK), 'E', 'V', M, L,&
			      HW)
                           GW = 1D0 
                           CW = HW*BC_VW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_VW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,E,M) = -(HALF*HW - ODX_E(I)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODX_E(I)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
!
               ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN 
	          IJPK = JP_OF(WEST_OF(IJK))
                  IF (WALL_AT(IJPK)) CYCLE
                  IF (EP_S(WEST_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,W,M) = ONE 
                  ELSE 
                     IF (FRICTION .AND. (ONE-EP_G(WEST_OF(IJK)))>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, WEST_OF(IJK), 'W', 'V', L, M&
                           , GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, WEST_OF(IJK), 'W', 'V', M, L,&
                              HW) 
                           GW = 1D0 
                           CW = HW*BC_VW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_VW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,W,M) = -(HALF*HW - ODX_E(IM)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODX_E(IM)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
!
               ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN 
	          IJPK = JP_OF(TOP_OF(IJK))
                  IF (WALL_AT(IJPK)) CYCLE
                  IF (EP_S(TOP_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,T,M) = ONE 
                  ELSE 
                     IF (FRICTION .AND. (ONE-EP_G(TOP_OF(IJK)))>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, TOP_OF(IJK), 'T', 'V', L, M&
                           , GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, TOP_OF(IJK), 'T', 'V', M, L, HW) 
                           GW = 1D0 
                           CW = HW*BC_VW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_VW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,T,M) = -(HALF*HW - ODZ_T(K)*OX(I)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODZ_T(K)*OX(I)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
!
               ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
	          IJPK = JP_OF(BOTTOM_OF(IJK))
                  IF (WALL_AT(IJPK)) CYCLE
                  IF (EP_S(BOTTOM_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,B,M) = ONE 
                  ELSE 
                     IF (FRICTION .AND. (ONE-EP_G(BOTTOM_OF(IJK)))>EPS_F_MIN) THEN 
                        CALL CALC_U_FRICTION (IJK, BOTTOM_OF(IJK), 'B', 'V', L&
                           , M, GW, HW, CW) 
                     ELSE 
                        IF (BC_JJ_PS(L) == 1) THEN 
                           CALL CALC_GRBDRY (IJK, BOTTOM_OF(IJK), 'B', 'V', M, &
                              L, HW) 
                           GW = 1D0 
                           CW = HW*BC_VW_S(L,M) 
                        ELSE IF (BC_JJ_PS(L) == 2) THEN 
                           GW = 0D0 
                           HW = 1D0 
                           CW = BC_VW_S(L,M) 
                        ELSE 
                           GW = 1D0 
                           CW = 0D0 
                           HW = 0D0 
                        ENDIF 
                     ENDIF 
                     A_M(IJK,B,M) = -(HALF*HW - ODZ_T(KM)*OX(I)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODZ_T(KM)*OX(I)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
               ENDIF 
            END DO 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE JJ_BC_V_S 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,kmax2->kmin3,kmax3      
!// 360 Check if i,j,k resides on current processor

