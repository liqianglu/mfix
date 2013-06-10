!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_U_s                                             C
!  Purpose: Determine source terms for U_s momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.  The center        C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this stage           C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Allow for partial-slip boundary conditions proposed by     C
!           by Johnson & Jackson (1987) if the Granular Temperature    C
!           equation is used.                                          C
!  Author: K. Agrawal, Princeton University           Date: 24-JAN-98  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_U_S(A_M, B_M, IER) 

!-----------------------------------------------
! Modules
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
      USE compar    
      USE sendrecv  
      use kintheory
      USE ghdtheory
      USE drag
      USE cutcell
      USE quadric
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m 
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 
 
! Vector b_m 
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M) 
! Error index 
      INTEGER          IER 
!-----------------------------------------------
! Local Variables 
!-----------------------------------------------
! Indices 
      INTEGER :: I, IJK,IMJK, IJMK, IJKE, IJKM, IPJK, IPJKM, &
                 J, K, IJPK, IPJMK, IJKP
! Phase index 
      INTEGER :: M, MM, L 
! Internal surface number 
      INTEGER :: ISV 
! Pressure at east cell 
      DOUBLE PRECISION :: PgE 
! average volume fraction 
      DOUBLE PRECISION :: EPSA, EPStmp, epse, epsw, epsn, epss, &
                          epst, epsb, epsMix, epsMixE
      DOUBLE PRECISION :: SUM_EPS_CP
! Average density 
      DOUBLE PRECISION :: ROPSA 
! Average density difference 
      DOUBLE PRECISION :: dro1, dro2, droa 
! Average quantities 
      DOUBLE PRECISION :: wse, EPMUSA 
! Source terms (Surface) 
      DOUBLE PRECISION :: Sdp, Sdps 
! Source terms (Volumetric) 
      DOUBLE PRECISION :: V0, Vmt, Vbf, Vcf, Vtza, Vmttmp
! Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION :: Ghd_drag, avgRop
! Source terms for HYS drag relation
      DOUBLE PRECISION :: HYS_drag, avgDrag
! virtual (added) mass
      DOUBLE PRECISION :: F_vir, ROP_MA, Uge, Ugw, Vgw, Vge, Ugn,&
                          Ugs, Wgb, Wgt, Wge, Ugb, Ugt
      DOUBLE PRECISION :: F_vir_tmp1, F_vir_tmp2	!Handan Liu added July 19 2012						  
! error message 
      CHARACTER*80     LINE(2) 
!-----------------------------------------------      
! Include statement functions
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
          IF (MOMENTUM_X_EQ(M)) THEN 

!!!$omp  parallel do private( IJK, IJKE, ISV, Sdp, Sdps, V0, Vmt, Vbf, &
!!!$omp&  I,PGE,DRO1,DRO2,DROA, IJKM,IPJK,IPJKM,  WSE,VCF,EPMUSA,VTZA, &
!!!$omp&  EPSA, EPStmp, ROPSA, LINE,SUM_EPS_CP,MM) &
!!!$omp&  schedule(static)
!Handan Liu added on July 19 2012, then modified March 1 2013:
!===================================================================<<< Handan Liu
!$omp  parallel do default(shared)												&
!$omp  private( I, J, K, IJK, IJKE, IMJK, IPJK, IJMK, IJPK, IPJMK, IJKM,		&
!$omp			IPJKM,IJKP, ISV, epsMix, epsMixE, EPStmp, EPSA, EPSw, EPSe,		&
!$omp			EPSn, EPSs, EPSt, EPSb, PGE, Sdp, SUM_EPS_CP, Sdps, MM, 		&
!$omp			ROPSA,V0, ROP_MA, L, Vgw, Vge, 	Uge, Ugw, Ugs, Ugn, Wgt, 		&
!$omp			Wgb, Ugt, Ugb, F_vir,  VMT, VMTtmp,	F_vir_tmp1, F_vir_tmp2,		&
!$omp 			DRO1,DRO2,DROA, Vbf, avgRop, Ghd_drag,HYS_drag, avgDrag, 		&
!$omp			WSE, VCF, EPMUSA, VTZA)
            DO IJK = ijkstart3, ijkend3 
                I = I_OF(IJK) 
                J = J_OF(IJK)
                K = K_OF(IJK)
                IJKE = EAST_OF(IJK)
                IMJK = IM_OF(IJK) 
                IJMK = JM_OF(IJK)
                IJKM = KM_OF(IJK) 
                IPJK = IP_OF(IJK) 
                IJPK = JP_OF(IJK) 
                IPJMK = IP_OF(IJMK)
                IPJKM = IP_OF(IJKM)
                IJKP = KP_OF(IJK)

                IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                  EPStmp = ZERO     
                  epsMix = ZERO
                  epsMixE= ZERO          
                  DO L = 1, SMAX
                    EPStmp = EPStmp + AVG_X(EP_S(IJK,L),EP_S(IJKE,L),I) 
                    epsMix  = epsMix  + EP_S(IJK,L) ! epsMix, epsMixE to be used for modelB
                    epsMixE = epsMixE + EP_S(IJKE,L)
                    IF(IP_AT_E(IJK)) THEN
                       U_S(IJK,L) = ZERO
                    ELSEIF(SIP_AT_E(IJK)) THEN 
                       ISV = IS_ID_AT_E(IJK)
                       U_S(IJK,L) = IS_VEL_S(ISV,L)
                    ENDIF
                  ENDDO                        
                  EPSA = EPStmp
                ELSE                  
                  EPSA = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I) 
                ENDIF

! Impermeable internal surface                
                IF (IP_AT_E(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO
! Semi-permeable internal surface                  
                ELSEIF (SIP_AT_E(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  ISV = IS_ID_AT_E(IJK) 
                  B_M(IJK,M) = -IS_VEL_S(ISV,M) 
! Dilute flow
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
! using the average boundary cell values to compute U_s (sof, Aug 23 2005)
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
                    B_M(IJK,M) = -U_S(IJK,M)          
                  ELSE
                    A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+A_M(IJK,N,M)+ &
                                     A_M(IJK,S,M)+A_M(IJK,T,M)+A_M(IJK,B,M))
                  ENDIF
! Cartesian grid implementation
               ELSEIF (BLOCKED_U_CELL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 

! Normal case
               ELSE 


! Surface forces:
! Pressure terms
                  PGE = P_G(IJKE) 
                  IF (CYCLIC_X_PD) THEN 
                     IF (CYCLIC_AT_E(IJK)) PGE = P_G(IJKE) - DELP_X 
                  ENDIF 
                  IF (MODEL_B) THEN 
                     SDP = ZERO 
                  ELSE 
                     IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
                        SDP = -P_SCALE*EPSA*(PGE - P_G(IJK))*AYZ(IJK) 
                     ELSE
                        SDP = -P_SCALE*EPSA*(PGE * A_UPG_E(IJK) - P_G(IJK) * A_UPG_W(IJK) )
                     ENDIF
                  ENDIF 

                  IF (CLOSE_PACKED(M)) THEN 
                     IF(SMAX > 1 .AND. TRIM(KT_TYPE) /= 'GHD') THEN
                        SUM_EPS_CP=0.0 
                        DO MM=1,SMAX
                          IF (CLOSE_PACKED(MM))&
                            SUM_EPS_CP=SUM_EPS_CP+AVG_X(EP_S(IJK,MM),EP_S(IJKE,MM),I)
                        ENDDO
                        SUM_EPS_CP = Max(SUM_EPS_CP, small_number)
                        SDPS = -( (P_S(IJKE,M)-P_S(IJK,M))+(EPSA/SUM_EPS_CP)*&
                           (P_STAR(IJKE)-P_STAR(IJK)) )*AYZ(IJK) 
                     ELSE
                        IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
                           SDPS =-((P_S(IJKE,M)-P_S(IJK,M))+(P_STAR(IJKE)-P_STAR(IJK)))*AYZ(IJK)  
                        ELSE
                           SDPS =-((P_S(IJKE,M)* A_UPG_E(IJK)-P_S(IJK,M)* A_UPG_W(IJK)) & 
                              +(P_STAR(IJKE)* A_UPG_E(IJK)-P_STAR(IJK)* A_UPG_W(IJK)))
                        ENDIF
                     ENDIF
                  ELSE 
                     IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
                        SDPS = -(P_S(IJKE,M)-P_S(IJK,M))*AYZ(IJK) 
                     ELSE
                        SDPS = - (P_S(IJKE,M) * A_UPG_E(IJK) - P_S(IJK,M) * A_UPG_W(IJK))
                     ENDIF
                  ENDIF 


                  IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
! Volumetric forces
                     ROPSA = HALF * (VOL(IJK)*ROP_S(IJK,M) + &
                        VOL(IPJK)*ROP_S(IJKE,M))/VOL_U(IJK) 
! Previous time step
                     V0 = HALF * (VOL(IJK)*ROP_SO(IJK,M) +&
                        VOL(IPJK)*ROP_SO(IJKE,M))*ODT/VOL_U(IJK)  
! Added mass implicit transient term {Cv eps rop_g dU/dt}
                     IF(Added_Mass .AND. M==M_AM) THEN
                        ROP_MA = AVG_X(ROP_g(IJK)*EP_s(IJK,M),ROP_g(IJKE)*EP_s(IJKE,M),I)
                        V0 = V0 + Cv * ROP_MA * ODT
                     ENDIF
                  ELSE
! Volumetric forces
                     ROPSA =  (VOL(IJK)*ROP_S(IJK,M) + &
                        VOL(IPJK)*ROP_S(IJKE,M))/(VOL(IJK) + VOL(IPJK))
! Previous time step
                     V0 = (VOL(IJK)*ROP_SO(IJK,M) + &
                        VOL(IPJK)*ROP_SO(IJKE,M))*ODT/&
                        (VOL(IJK) + VOL(IPJK))
! Added mass implicit transient term {Cv eps rop_g dU/dt}
                     IF(Added_Mass .AND. M==M_AM) THEN
                       ROP_MA = (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M) +&
                          VOL(IPJK)*ROP_g(IJKE)*EP_s(IJKE,M))/&
                          (VOL(IJK) + VOL(IPJK))
                       V0 = V0 + Cv * ROP_MA * ODT
                     ENDIF
                  ENDIF

! VIRTUAL MASS SECTION (explicit terms)
! adding transient term dvg/dt - dVs/dt to virtual mass term    
                  F_vir = ZERO
!Handan Liu changed on July 19 2012, then modified on March 1 2013:
!===================================================================<< Handan Liu				  
                  IF(Added_Mass .AND. M==M_AM .AND.&
                    (.NOT.CUT_U_TREATMENT_AT(IJK))) THEN      
                     !F_vir = ( (U_G(IJK) - U_GO(IJK)) )*ODT*VOL_U(IJK)
					 F_vir_tmp1 = ( (U_G(IJK) - U_GO(IJK)) )*ODT*VOL_U(IJK)

! defining gas-particles velocity at momentum cell faces (or scalar cell center)    
                     Ugw = AVG_X_E(U_G(IMJK),U_G(IJK),I)
                     Uge = AVG_X_E(U_G(IJK),U_G(IPJK),IP1(I))  
                     Vgw = AVG_Y_N(V_G(IJMK),V_G(IJK))  
                     Vge = AVG_Y_N(V_G(IPJMK),V_G(IPJK)) 
                     Ugs = AVG_Y(U_G(IJMK),U_G(IJK),JM1(J))
                     Ugn = AVG_Y(U_G(IJK),U_G(IJPK),J)
                     IF(DO_K) THEN
                        Wgb = AVG_Z_T(W_g(IJKM),W_g(IJK))  
                        Wgt = AVG_Z_T(W_g(IPJKM),W_g(IPJK)) 
                        Wge = AVG_X(Wgb,Wgt,I)
                        Ugb = AVG_Z(U_g(IJKM),U_g(IJK),KM1(K))
                        Ugt = AVG_Z(U_g(IJK),U_g(IJKP),K)
                        !F_vir = F_vir + Wge*OX_E(I) * (Ugt - Ugb) *AXY(IJK)
                        !IF (CYLINDRICAL) F_vir = F_vir - &
                        !   Wge**2*OX_E(I) ! centrifugal force
                        F_vir_tmp2 = F_vir_tmp1 + Wge*OX_E(I) * (Ugt - Ugb) *AXY(IJK)
                        IF (CYLINDRICAL) F_vir_tmp2 = F_vir_tmp1 - &
                           Wge**2*OX_E(I)						   
                     ENDIF
! adding convective terms (U dU/dx + V dU/dy + W dU/dz) to virtual mass
                     !F_vir = F_vir + U_g(IJK)*(Uge - Ugw) *AYZ(IJK) + &
                     !   AVG_X(Vgw,Vge,I) * (Ugn - Ugs) *AXZ(IJK)
                     !F_vir = F_vir * Cv * ROP_MA
                     F_vir_tmp2 = F_vir_tmp1 + U_g(IJK)*(Uge - Ugw) *AYZ(IJK) + &
                        AVG_X(Vgw,Vge,I) * (Ugn - Ugs) *AXZ(IJK)
                     F_vir = F_vir_tmp2 * Cv * ROP_MA					 
                  ENDIF
!===================================================================>> Handan Liu

! Interphase mass transfer
                  IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                     VMTtmp = ZERO
                     DO L = 1,SMAX
                        VMTtmp = VMTtmp + HALF*(VOL(IJK)*SUM_R_S(IJK,L) + &
                           VOL(IPJK)*SUM_R_S(IJKE,L))/VOL_U(IJK)
                     ENDDO
                     VMT = VMTtmp
                  ELSE
                     IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
                        VMT = HALF * (VOL(IJK)*SUM_R_S(IJK,M) + &
                           VOL(IPJK)*SUM_R_S(IJKE,M))/VOL_U(IJK)  
                     ELSE
                        VMT = (VOL(IJK)*SUM_R_S(IJK,M) + &
                           VOL(IPJK)*SUM_R_S(IJKE,M))/(VOL(IJK) + VOL(IPJK))
                     ENDIF
                  ENDIF

! Body force
                  IF (MODEL_B) THEN 
                     IF (TRIM(KT_TYPE) /= 'GHD') THEN
!QX                     
                       DRO1 = (RO_SV(IJK,M)-RO_G(IJK))*EP_S(IJK,M) 
                       DRO2 = (RO_SV(IJK,M)-RO_G(IJKE))*EP_S(IJKE,M) 
                       DROA = AVG_X(DRO1,DRO2,I) 
                       VBF = DROA*BFX_S(IJK,M) 
                     ELSE ! GHD and M = MMAX
                       DRO1 = ROP_S(IJK,M)  - RO_G(IJK) *epsMix
                       DRO2 = ROP_S(IJKE,M) - RO_G(IJKE)*epsMixE
                       DROA = AVG_X(DRO1,DRO2,I) 
                       VBF = DROA*BFX_S(IJK,M) 
                     ENDIF
                  ELSE ! model A
                     VBF = ROPSA*BFX_S(IJK,M) 
                  ENDIF 

! Additional force for GHD from darg force sum(beta_ig * Joi/rhop_i)
                  Ghd_drag = ZERO
                  IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                    DO L = 1,SMAX
                      avgRop = AVG_X(ROP_S(IJK,L),ROP_S(IJKE,L),I)
                      if(avgRop > ZERO) Ghd_drag = Ghd_drag - &
                           AVG_X(F_GS(IJK,L),F_GS(IJKE,L),I) * JoiX(IJK,L) / avgRop
                    ENDDO
                  ENDIF

! Additional force for HYS drag force, do not use with mixture GHD theory
                  HYS_drag = ZERO
                  IF (TRIM(DRAG_TYPE) .EQ. 'HYS' .AND. TRIM(KT_TYPE) /= 'GHD') THEN
                     DO L = 1,MMAX
                        IF (L /= M) THEN
                           avgDrag = AVG_X(beta_ij(IJK,M,L),beta_ij(IJKE,M,L),I)
                           HYS_drag = HYS_drag - avgDrag * (U_g(ijk) - U_s(IJK,L))
                        ENDIF
                     ENDDO
                  ENDIF


! Special terms for cylindrical coordinates
                  IF (CYLINDRICAL) THEN 
! centrifugal force
                     WSE = AVG_X(HALF*(W_S(IJK,M)+W_S(IJKM,M)),&
                                 HALF*(W_S(IPJK,M)+W_S(IPJKM,M)),I) 
                     VCF = ROPSA*WSE**2*OX_E(I) 
                     IF(Added_Mass .AND. M==M_AM) &
                        VCF = VCF + Cv*ROP_MA*WSE**2*OX_E(I) ! virtual mass contribution
! -(2mu/x)*(u/x) part of Tau_zz/X
                     EPMUSA = AVG_X(MU_S(IJK,M),MU_S(IJKE,M),I) 
                     VTZA = 2.d0*EPMUSA*OX_E(I)*OX_E(I) 
                  ELSE 
                     VCF = ZERO 
                     VTZA = ZERO 
                  ENDIF 

! Collect the terms
                  A_M(IJK,0,M) = -(A_M(IJK,E,M)+A_M(IJK,W,M)+&
                     A_M(IJK,N,M)+A_M(IJK,S,M)+A_M(IJK,T,M)+&
                     A_M(IJK,B,M)+(V0+ZMAX(VMT)+VTZA)*VOL_U(IJK)) 

                    B_M(IJK,M) = B_M(IJK,M) - (SDP + SDPS + &
                       TAU_U_S(IJK,M) + F_vir + &
                       ( (V0+ZMAX((-VMT)))*U_SO(IJK,M)+&
                       VBF + VCF + HYS_drag)*VOL_U(IJK) )
                     
                  IF (TRIM(KT_TYPE) .EQ. 'IA_NONEP') THEN 
                     B_M(IJK,M) = B_M(IJK,M) - KTMOM_U_S(IJK,M)
                  ELSEIF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
                     B_M(IJK,M) = B_M(IJK,M) - Ghd_drag*VOL_U(IJK)
                  ENDIF

                ENDIF   ! end branching on cell type (sip/ip/dilute/block/else branches)
            ENDDO   ! end do loop over ijk
!$omp end parallel do	

! modifications for cartesian grid implementation
            IF(CARTESIAN_GRID) CALL CG_SOURCE_U_S (A_M, B_M, M, IER)
! modifications for bc            
            CALL SOURCE_U_S_BC (A_M, B_M, M, IER) 
            IF(CARTESIAN_GRID) CALL CG_SOURCE_U_S_BC (A_M, B_M, M, IER)

          ENDIF   ! end if (momentum_x_eq)  
        ENDIF   ! end if for ghd theory
      ENDDO   ! end do loop over mmax
      
      RETURN  
      END SUBROUTINE SOURCE_U_S 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_U_s_BC                                          C
!  Purpose: Determine source terms for U_g momentum eq. The terms      C
!  appear in the center coefficient and RHS vector.  The center        C
!  coefficient and source vector are negative.  The off-diagonal       C
!  coefficients are positive.                                          C
!  The drag terms are excluded from the source at this stage.          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-MAY-96  C
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

      SUBROUTINE SOURCE_U_S_BC(A_M, B_M, M, IER) 

!-----------------------------------------------
! Modules
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
! Dummy Arguments
!-----------------------------------------------    
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Solids phase index
      INTEGER, INTENT(IN) :: M 
! Error index 
      INTEGER, INTENT(INOUT) :: IER 
!-----------------------------------------------
! Local Variables      
!-----------------------------------------------
! Boundary condition 
      INTEGER ::  L 
! Indices 
      INTEGER :: I, J, K, I1, I2, J1, J2, K1, K2, IJK,& 
                 IM, JM, KM, IJKW, IMJK, IPJK, IP 
!----------------------------------------------- 
! Include statements functions
!-----------------------------------------------
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
!-----------------------------------------------

! Setting the default boundary conditions
      IF (DO_K) THEN 
         K1 = 1 
         DO J1 = jmin3,jmax3 
            DO I1 = imin3, imax3  
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IJK = FUNIJK(I1,J1,K1) 
               IF (NS_WALL_AT(IJK)) THEN 
! Setting the wall velocity to zero
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = -ONE 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ELSEIF (FS_WALL_AT(IJK)) THEN 
! Setting the wall velocity equal to the adjacent fluid velocity
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ONE 
                  A_M(IJK,B,M) = ZERO 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            ENDDO 
         ENDDO 
 
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
               ELSEIF (FS_WALL_AT(IJK)) THEN 
                  A_M(IJK,E,M) = ZERO 
                  A_M(IJK,W,M) = ZERO 
                  A_M(IJK,N,M) = ZERO 
                  A_M(IJK,S,M) = ZERO 
                  A_M(IJK,T,M) = ZERO 
                  A_M(IJK,B,M) = ONE 
                  A_M(IJK,0,M) = -ONE 
                  B_M(IJK,M) = ZERO 
               ENDIF 
            ENDDO 
         ENDDO 
      ENDIF 

      J1 = 1 
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = -ONE 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSEIF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ONE 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         ENDDO 
      ENDDO 

      J1 = JMAX2 
      DO K1 = kmin3, kmax3 
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IJK = FUNIJK(I1,J1,K1) 
            IF (NS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = -ONE 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ELSEIF (FS_WALL_AT(IJK)) THEN 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ONE 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 
            ENDIF 
         ENDDO 
      ENDDO 
! End setting the default boundary conditions
  

! Setting user specified boundary conditions
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 

! Setting wall boundary conditions
            IF (BC_TYPE(L) == 'NO_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L)
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
                           IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                              A_M(IJK,N,M) = -ONE 
                           ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                              A_M(IJK,S,M) = -ONE 
                           ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN 
                              A_M(IJK,T,M) = -ONE 
                           ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                              A_M(IJK,B,M) = -ONE 
                           ENDIF 
                        ENDDO 
                     ENDDO 
                  ENDDO   
               ELSE   ! Johnson and Jackson partial slip 
                  CALL JJ_BC_U_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
               ENDIF 

            ELSEIF (BC_TYPE(L) == 'FREE_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L)
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
                           IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                              A_M(IJK,N,M) = ONE 
                           ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                              A_M(IJK,S,M) = ONE 
                           ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN 
                              A_M(IJK,T,M) = ONE 
                           ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                              A_M(IJK,B,M) = ONE 
                           ENDIF 
                        ENDDO 
                     ENDDO 
                  ENDDO
               ELSE   ! Johnson and Jackson partial slip 
                  CALL JJ_BC_U_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
               ENDIF 

            ELSEIF (BC_TYPE(L) == 'PAR_SLIP_WALL') THEN 
               I1 = BC_I_W(L) 
               I2 = BC_I_E(L) 
               J1 = BC_J_S(L) 
               J2 = BC_J_N(L) 
               K1 = BC_K_B(L) 
               K2 = BC_K_T(L)
               IF (BC_JJ_PS(L) == 0) THEN 
                  DO K = K1, K2 
                     DO J = J1, J2 
                        DO I = I1, I2 
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IJK = FUNIJK(I,J,K) 
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                           JM = JM1(J) 
                           KM = KM1(K) 
                           A_M(IJK,E,M) = ZERO 
                           A_M(IJK,W,M) = ZERO 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                           IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,N,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODY_N(J)) 
                                 A_M(IJK,N,M) = -(HALF*BC_HW_S(L,M)-ODY_N(J)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,S,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,S,M) = -(HALF*BC_HW_S(L,M)-ODY_N(JM)) 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODY_N(JM)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,T,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODZ_T(K)*&
                                    OX_E(I)) 
                                 A_M(IJK,T,M) = -(HALF*BC_HW_S(L,M)-ODZ_T(K)*&
                                    OX_E(I)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN 
                                 A_M(IJK,B,M) = -HALF 
                                 A_M(IJK,0,M) = -HALF 
                                 B_M(IJK,M) = -BC_UW_S(L,M) 
                              ELSE 
                                 A_M(IJK,B,M) = -(HALF*BC_HW_S(L,M)-ODZ_T(KM)*&
                                    OX_E(I)) 
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODZ_T(KM)*&
                                    OX_E(I)) 
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_UW_S(L,M) 
                              ENDIF 
                           ENDIF 
                        ENDDO 
                     ENDDO 
                  ENDDO
               ELSE   ! Johnson and Jackson partial slip 
                  CALL JJ_BC_U_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 
               ENDIF 

! Setting flow boundary conditions
            ELSEIF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'W') THEN 
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
                           A_M(IJK,W,M) = ONE 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                        ENDDO 
                     ENDDO 
                  ENDDO  
               ENDIF

            ELSEIF (BC_TYPE(L) == 'OUTFLOW') THEN 
               IF (BC_PLANE(L) == 'W') THEN 
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
                           A_M(IJK,W,M) = ONE 
                           A_M(IJK,N,M) = ZERO 
                           A_M(IJK,S,M) = ZERO 
                           A_M(IJK,T,M) = ZERO 
                           A_M(IJK,B,M) = ZERO 
                           A_M(IJK,0,M) = -ONE 
                           B_M(IJK,M) = ZERO 
                           IM = IM1(I) 
                           IMJK = IM_OF(IJK) 
                           A_M(IMJK,E,M) = ZERO 
                           A_M(IMJK,W,M) = X_E(IM)/X_E(IM1(IM)) 
                           A_M(IMJK,N,M) = ZERO 
                           A_M(IMJK,S,M) = ZERO 
                           A_M(IMJK,T,M) = ZERO 
                           A_M(IMJK,B,M) = ZERO 
                           A_M(IMJK,0,M) = -ONE 
                           B_M(IMJK,M) = ZERO 
                        ENDDO 
                     ENDDO 
                  ENDDO   
               ELSEIF (BC_PLANE(L) == 'E') THEN 
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
                           IP = IP1(I) 
                           IPJK = IP_OF(IJK) 
                           A_M(IPJK,E,M) = X_E(IP)/X_E(I) 
                           A_M(IPJK,W,M) = ZERO 
                           A_M(IPJK,N,M) = ZERO 
                           A_M(IPJK,S,M) = ZERO 
                           A_M(IPJK,T,M) = ZERO 
                           A_M(IPJK,B,M) = ZERO 
                           A_M(IPJK,0,M) = -ONE 
                           B_M(IPJK,M) = ZERO 
                        ENDDO 
                     ENDDO 
                  ENDDO  
               ENDIF

! Setting bc that are not ns, fs, psw, p_inflow, p_outflow, or outflow
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
                        B_M(IJK,M) = -U_S(IJK,M) 
                        IF (BC_PLANE(L) == 'W') THEN 
                           IJKW = WEST_OF(IJK) 
                           A_M(IJKW,E,M) = ZERO 
                           A_M(IJKW,W,M) = ZERO 
                           A_M(IJKW,N,M) = ZERO 
                           A_M(IJKW,S,M) = ZERO 
                           A_M(IJKW,T,M) = ZERO 
                           A_M(IJKW,B,M) = ZERO 
                           A_M(IJKW,0,M) = -ONE 
                           B_M(IJKW,M) = -U_S(IJKW,M) 
                        ENDIF 
                     ENDDO 
                  ENDDO 
               ENDDO 
       
            ENDIF   ! end if (bc_type)
         ENDIF   ! end if (bc_defined) 
      ENDDO   ! end L do loop over dimension_bc

      RETURN  
      END SUBROUTINE SOURCE_U_S_BC 


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: JJ_BC_U_s                                              C
!  Purpose: Implement Johnson and Jackson boundary condition           C
!                                                                      C
!  Author: K. Agrawal & A. Srivastava,                Date: 14-APR-98  C
!          Princeton University                                        C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE JJ_BC_U_S(I1, I2, J1, J2, K1, K2, L, M, A_M, B_M) 

!-----------------------------------------------
! Modules
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
! Dummy Arguments
!-----------------------------------------------    
! Boundary condition
      INTEGER, INTENT(IN) :: L
! Indices
      INTEGER, INTENT(IN) :: I1, I2, J1, J2, K1, K2
! Solids phase index
      INTEGER, INTENT(IN) ::  M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local Variables      
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, JM, KM, IPJK
! coefficients for granular bc
      DOUBLE PRECISION :: hw, gw, cw
!----------------------------------------------- 
! Include statements functions
!-----------------------------------------------
      INCLUDE 'ep_s1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

      DO K = K1, K2 
         DO J = J1, J2 
            DO I = I1, I2 
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE

               IJK = FUNIJK(I,J,K) 
               IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells 
               JM = JM1(J) 
               KM = KM1(K) 
               A_M(IJK,E,M) = ZERO 
               A_M(IJK,W,M) = ZERO 
               A_M(IJK,N,M) = ZERO 
               A_M(IJK,S,M) = ZERO 
               A_M(IJK,T,M) = ZERO 
               A_M(IJK,B,M) = ZERO 
               A_M(IJK,0,M) = -ONE 
               B_M(IJK,M) = ZERO 

               IF (FLUID_AT(NORTH_OF(IJK))) THEN 
                  IPJK = IP_OF(NORTH_OF(IJK))
                  IF (WALL_AT(IPJK)) CYCLE
                  IF (EP_S(NORTH_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,N,M) = ONE 
                  ELSE 
                     IF (BC_JJ_PS(L) == 1) THEN 
! Setting the wall velocity based on Johnson and Jackson B.C which are
! also modified to include frictional effects if requested
                        CALL CALC_GRBDRY (IJK, NORTH_OF(IJK), 'N', 'U',&
                           M, L, GW, HW, CW) 
                     ELSEIF (BC_JJ_PS(L) == 2) THEN 
! Setting the wall velocity equal to the user specified value 
! (if wall velocity = 0 then it is equivalent to no slip wall) 
                        GW = 0D0 
                        HW = 1D0 
                        CW = BC_UW_S(L,M) 
                     ELSE 
! Setting the wall velocity equal to the adjacent fluid cell 
! velocity (i.e. zero flux/free slip)
                        GW = 1D0 
                        CW = 0D0 
                        HW = 0D0 
                     ENDIF 
                     A_M(IJK,N,M) = -(HALF*HW - ODY_N(J)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODY_N(J)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 

               ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                  IPJK = IP_OF(SOUTH_OF(IJK))
                  IF (WALL_AT(IPJK)) CYCLE
                  IF (EP_S(SOUTH_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,S,M) = ONE 
                  ELSE 
                     IF (BC_JJ_PS(L) == 1) THEN 
                        CALL CALC_GRBDRY (IJK, SOUTH_OF(IJK), 'S', 'U',&
                           M, L, GW, HW, CW) 
                     ELSEIF (BC_JJ_PS(L) == 2) THEN 
                        GW = 0D0 
                        HW = 1D0 
                        CW = BC_UW_S(L,M) 
                     ELSE 
                        GW = 1D0 
                        CW = 0D0 
                        HW = 0D0 
                     ENDIF 
                     A_M(IJK,S,M) = -(HALF*HW - ODY_N(JM)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODY_N(JM)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 

               ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN 
                  IPJK = IP_OF(TOP_OF(IJK))
                  IF (WALL_AT(IPJK)) CYCLE
                  IF (EP_S(TOP_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,T,M) = ONE 
                  ELSE 
                     IF (BC_JJ_PS(L) == 1) THEN 
                        CALL CALC_GRBDRY (IJK, TOP_OF(IJK), 'T', 'U',&
                           M, L, GW, HW, CW) 
                     ELSEIF (BC_JJ_PS(L) == 2) THEN 
                        GW = 0D0 
                        HW = 1D0 
                        CW = BC_UW_S(L,M) 
                     ELSE 
                        GW = 1D0 
                        CW = 0D0 
                        HW = 0D0 
                     ENDIF 
                     A_M(IJK,T,M) = -(HALF*HW - ODZ_T(K)*OX_E(I)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODZ_T(K)*OX_E(I)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 

               ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN 
                  IPJK = IP_OF(BOTTOM_OF(IJK))
                  IF (WALL_AT(IPJK)) CYCLE
                  IF (EP_S(BOTTOM_OF(IJK),M) <= DIL_EP_S) THEN 
                     A_M(IJK,B,M) = ONE 
                  ELSE 
                     IF (BC_JJ_PS(L) == 1) THEN 
                        CALL CALC_GRBDRY (IJK, BOTTOM_OF(IJK), 'B', 'U',&
                           M,L, GW, HW, CW) 
                     ELSEIF (BC_JJ_PS(L) == 2) THEN 
                        GW = 0D0 
                        HW = 1D0 
                        CW = BC_UW_S(L,M) 
                     ELSE 
                        GW = 1D0 
                        CW = 0D0 
                        HW = 0D0 
                     ENDIF 
                     A_M(IJK,B,M) = -(HALF*HW - ODZ_T(KM)*OX_E(I)*GW) 
                     A_M(IJK,0,M) = -(HALF*HW + ODZ_T(KM)*OX_E(I)*GW) 
                     B_M(IJK,M) = -CW 
                  ENDIF 
               ENDIF

            ENDDO 
         ENDDO 
      ENDDO 
     
      RETURN  
      END SUBROUTINE JJ_BC_U_S 

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_U_S                                        C
!  Purpose: Adds point sources to the solids U-Momentum equations.     C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_U_S(A_M, B_M, IER) 

      use compar    
      use constant
      use fldvar
      use geometry
      use indices
      use physprop
      use ps
      use run

! To be removed upon complete integration of point source routines.
      use bc
      use usr

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m 
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) 

! Vector b_m 
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M) 

! Error index 
      INTEGER, INTENT(INOUT) :: IER 

!----------------------------------------------- 
! Local Variables
!----------------------------------------------- 

! Indices 
      INTEGER :: IJK, I, J, K
      INTEGER :: IJKE, IPJK
      INTEGER :: BCV, M

! terms of bm expression
      DOUBLE PRECISION :: pSource

      include 'function.inc'

      do m=1 , mmax

! Calculate the mass going into each IJK cell. This is done for each 
! call in case the point source is time dependent.
         BC_LP: do BCV = 50, DIMENSION_BC
            if(POINT_SOURCES(BCV) == 0) cycle BC_LP
            if(abs(BC_U_s(BCV,M)) < small_number) cycle BC_LP

            do k = BC_K_B(BCV), BC_K_T(BCV)
            do j = BC_J_S(BCV), BC_J_N(BCV)
            do i = BC_I_W(BCV), BC_I_E(BCV)

               ijk = funijk(i,j,k)
               if(BC_U_s(BCV,M) < 0.0d0) ijk = WEST_OF(ijk)

               if(fluid_at(ijk)) then

                  if(A_M(IJK,0,M) == -ONE .AND.                        &
                    B_M(IJK,M) == -U_s(IJK,M)) then
                    B_M(IJK,M) = -BC_U_s(BCV,M) * PS_VEL_MAG_S(BCV,M)
                  else
                     pSource =  BC_MASSFLOW_S(BCV,M) *  &
                        (VOL(IJK)/PS_VOLUME(BCV))

                     A_M(IJK,0,M) = A_M(IJK,0,M) - pSource 
                     B_M(IJK,M) = B_M(IJK,M) - pSource *               &
                        BC_U_s(BCV,M) * PS_VEL_MAG_S(BCV,M)
                  endif
               endif

            enddo
            enddo
            enddo

         enddo BC_LP

      enddo ! do M=1, MMAX


      RETURN
      END SUBROUTINE POINT_SOURCE_U_S
