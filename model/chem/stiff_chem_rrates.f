!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STIFF_CHEM_RRATES(IER)                                  C
!                                                                      C
!  Purpose: Calculate reaction rates for various reactions in cell ijk C
!           using information from the data file                       C
!                                                                      C
!  Author: J. Musser                                  Date: 12-Feb-13  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STIFF_CHEM_RRATES(lNEQ, lTime, Y, YDOT)

! Global Variables:
!---------------------------------------------------------------------//
      use fldvar,     only : EP_g, RO_g, T_g, X_g     ! Gas Phase
      use fldvar,     only : ROP_S, T_s, X_s, D_p     ! Solids Phases

      use rxns,       only : REACTION, NO_OF_RXNS
      use run,        only : SPECIES_EQ, UNITS
      use physprop,   only : RO_s, C_pg, C_ps, MMAX, NMAX, RO_sv
      use stiff_chem, only : CALL_GROW, mapODEtoMFIX

      use sendrecv

      implicit none

! Passed Variables: Dummy argument format required by ODEPACK.
!---------------------------------------------------------------------//
! (1) Number of ODEs to be solve
! (2) Fluid cell index
      INTEGER, intent(in) :: lNEQ(2)
! Independent variable (not used) 
      INTEGER, intent(in) :: lTime
! Array of dependent variable initial values.
      DOUBLE PRECISION, dimension(lNEQ(1)), intent(in)  :: Y
! Rate of change of dependent variables.
      DOUBLE PRECISION, dimension(lNEQ(1)), intent(out) :: YDOT

! Local Variables:
!---------------------------------------------------------------------//
      INTEGER IJK  ! Fluid Cell Index
      INTEGER H    ! Reaction loop counter
      INTEGER L, M ! Global Phase index loop counters
      INTEGER N    ! Global species index
      INTEGER lN   ! Local reaction speices index/loop counter
      INTEGER mXfr ! Global phase index for mass transfer

! User-defined reaction rates returned from USR_RATES
      DOUBLE PRECISION RATES(NO_OF_RXNS)
! Rate of formation (+) or consumption (-) of a species (GENERIC).
      DOUBLE PRECISION lRate
! Rate of formation (+) or consumption (-) of each gas phase species.
! > Accumulative over all reactions.
      DOUBLE PRECISION lRg(DIMENSION_N_g)
! Rate of formation (+) or consumption (-) of each gas phase species.
! > Single reaction.
      DOUBLE PRECISION rRg(DIMENSION_N_g)
! Net rate of change of gas phase material.
      DOUBLE PRECISION sumlRg
! Heat of reaction assigned to gas phase for a reaction.
      DOUBLE PRECISION rHORg
! Cummulative heat of reaction assigned to gas phase.
      DOUBLE PRECISION lHORg
! Rate of formation (+) or consumption (-) of each solids phase species.
! > Accumulative over all reactions.
      DOUBLE PRECISION lRs(DIMENSION_M, DIMENSION_N_s)
! Rate of formation (+) or consumption (-) of each solids phase species.
! > Single reaction.
      DOUBLE PRECISION rRs(DIMENSION_M, DIMENSION_N_s)
! Net rate of change of solids phase material.
      DOUBLE PRECISION sumlRs(DIMENSION_M)
! Heat of reaction assigned to the solids phase for a reaction.
      DOUBLE PRECISION rHORs(DIMENSION_M)
! Cummulative heat of reaction assinged to the solids phase.
      DOUBLE PRECISION lHORs(DIMENSION_M)
! Rate of interphase enthalpy transfer due to mass transfer.
      DOUBLE PRECISION RxH(0:DIMENSION_M, 0:DIMENSION_M)

! Loop index for ODEs
      INTEGER :: Node
! Error Flag. Used to supress checks in mapODEtoMFIX.
      INTEGER :: iErr

! Reaction limiters.
      DOUBLE PRECISION, parameter :: speciesLimiter = 1.0d-6

! External Function for comparing two numbers.
      LOGICAL, external :: COMPARE
      DOUBLE PRECISION, external :: CALC_H0, CALC_H

! UDF for Reaction Rates:
      external USR_RATES

!-----------------------------------------------
      include 'ep_s1.inc'
      include 'function.inc'
      include 'ep_s2.inc'

! Initialize variables:
      IJK   = lNEQ(2)
      RATES = ZERO
      YDOT  = ZERO

      lRg = ZERO; sumlRg = ZERO; lHORg = ZERO
      lRs = ZERO; sumlRs = ZERO; lHORs = ZERO

! Map the current ODE independent variables to MFIX variables.
      CALL mapODEtoMFIX(IJK, Y)

! Calculate user defined reaction rates.
      CALL USR_RATES(IJK, RATES)

! Loop over reactions.
      RXN_LP: DO H = 1, NO_OF_RXNS

! Skip empty reactions
         IF(Reaction(H)%nSpecies == 0) CYCLE RXN_LP

! Initialize local loop arrays
         rRg = ZERO; rHORg = ZERO
         rRs = ZERO; rHORs = ZERO
         RxH = ZERO

! Calculate the rate of formation/consumption for each species.
!---------------------------------------------------------------------//
         DO lN = 1, Reaction(H)%nSpecies
! Global phase index.
            M = Reaction(H)%Species(lN)%pMap
! Global species index.
            N = Reaction(H)%Species(lN)%sMap
! Index for interphase mass transfer. For a gas/solid reaction, the 
! index is stored with the gas phase. For solid/solid mass transfer
! the index is stored with the source phase.
            mXfr = Reaction(H)%Species(lN)%mXfr
            lRate = RATES(H) * Reaction(H)%Species(lN)%MWxStoich
! Gas Phase:
            IF(M == 0) THEN
! Consumption of gas phase species.
               IF(lRate < ZERO) THEN
! Check that there is a sufficient amount of reactant.
                  IF(X_g(IJK,N) > speciesLimiter) THEN
                     rRg(N) = rRg(N) + lRate
! Enthalpy transfer associated with mass transfer. (gas/solid)
                     IF(M /= mXfr) RxH(M,mXfr) =  RxH(M,mXfr) + &
                        lRate * CALC_H0(T_G(IJK),0,N)
                  ELSE
! There is an insignificant amount of reactant. Skip this reaction.
                     CYCLE RXN_LP
                  ENDIF
               ELSE
! Formation of gas phase species.
                  rRg(N) = rRg(N) + lRate
! Enthalpy transfer associated with mass transfer. (gas/solid)
                  IF(M /= mXfr) RxH(M,mXfr) =  RxH(M,mXfr) + &
                     lRate * CALC_H0(T_s(IJK,mXfr),0,N)
               ENDIF
! Solids Phase M:
            ELSE
! Consumption of solids phase species.
               IF(lRate < ZERO) THEN
                  IF(X_s(IJK,M,N) > speciesLimiter) THEN
                     rRs(M,N) = rRs(M,N) + lRate
! Enthalpy transfer associated with mass transfer. (solid/solid) This
! is only calculated from the source (reactant) material.
                     IF(M /= mXfr) THEN
                        IF(M < mXfr) THEN
                           RxH(M,mXfr) =  RxH(M,mXfr) +                &
                              lRate * CALC_H(IJK,M,N) *                &
                              Reaction(H)%Species(lN)%xXfr
                        ELSE
                           RxH(mXfr,M) =  RxH(mXfr,M) -                &
                             lRate * CALC_H(IJK,M,N) *                 &
                             Reaction(H)%Species(lN)%xXfr
                        ENDIF
                     ENDIF
                  ELSE
! There is an insignificant amount of reactant. Skip this reaction.
                     CYCLE RXN_LP
                  ENDIF
               ELSE
! Formation of solids phase species.
                  rRs(M,N) = rRs(M,N) + lRate
               ENDIF
            ENDIF
         ENDDO ! Loop of species

! Copy the single reaction rates of formation/consumption to the total
! (accumulative) rates of formation/consumption arrays. This ensures
! that any reactions without sufficient reactants are not included.
         lRg = lRg + rRg
         lRs = lRs + rRs

! Calculate and store the heat of reaction.
!---------------------------------------------------------------------//
! Automated heat of reaction calculations
         IF(Reaction(H)%Calc_DH) THEN
! Loop over reaction species.
            DO lN = 1, Reaction(H)%nSpecies
! Global phase index.
               M = Reaction(H)%Species(lN)%pMap
! Global species index.
               N = Reaction(H)%Species(lN)%sMap
! Rate of formation/consumption for speices N
               lRate = RATES(H) * Reaction(H)%Species(lN)%MWxStoich
! Gas phase enthalpy chnage from energy equation derivation.
               IF(M == 0) THEN
                  rHORg = rHORg + lRate * CALC_H0(T_g(IJK),0,N)
! Solid phase enthalpy change from energy equation derivation.
               ELSE
                  rHORs(M) = rHORs(M) + lRate * CALC_H0(T_s(IJK,M),M,N)
               ENDIF
            ENDDO

! Complete the skew-symettric for enthalpy transfer with mass transfer
            DO M=1, MMAX
               DO L=0, M-1
                  RxH(M,L) = - RxH(L,M)
               ENDDO
            ENDDO
! Apply enthalpy transfer associated with mass transfer to get the
! complete heat of reaction of heat phse for Reaction H.
            DO L=0, MMAX
               DO M = 0, MMAX
                  IF(L == M) CYCLE
                  IF(L == 0) THEN
                     rHORg = rHORg - RxH(L,M)
                  ELSE
                     rHORs(L) = rHORs(L) - RxH(L,M)
                  ENDIF
               ENDDO
            ENDDO

! Convert the heat of reaction to the appropriate units (if SI), and 
! store in the global array.
            IF(UNITS == 'SI') THEN
               lHORg = lHORg + 4.183925d3*rHORg
               DO M=1,MMAX
                  lHORs = lHORs + 4.183925d3*rHORs(M)
               ENDDO
            ELSE
               lHORg = lHORg + rHORg
               DO M=1,MMAX
                  lHORs(M) = lHORs(M) + rHORs(M)
               ENDDO
            ENDIF
         ELSE
! User-defined heat of reaction.
            lHORg = Reaction(H)%HoR(0) * RATES(H)
            DO M=1, MMAX
               lHORs = Reaction(H)%HoR(M) * RATES(H)
            ENDDO
         ENDIF

      ENDDO RXN_LP ! Loop over reactions.

      sumlRg = sum(lRg)
      DO M=1, MMAX
         sumlRs(M) = sum(lRs(M,:))
      ENDDO



!<<<--------------------------  ODE Setup  ------------------------->>>!

! Gas Phase Variables:
!---------------------------------------------------------------------//
! Density:  RO_g
      YDOT(1) = (sumlRg + RO_g(IJK)*sum(sumlRs(:)/RO_s(:)))/EP_G(IJK)
! Temperature: T_g
      YDOT(2) = -lHORg/(EP_G(IJK)*RO_g(IJK)*C_pg(IJK))
! Species mass fractions: X_g
      DO N=1,NMAX(0)
         Node = 2 + N
         YDOT(Node) = (lRg(N) - X_g(IJK,N)*sumlRg)/(EP_G(IJK)*RO_g(IJK))
      ENDDO

      Node = 3+NMAX(0)

! Solids Phase Variables:
!---------------------------------------------------------------------//
      DO M = 1, MMAX
         IF(SPECIES_EQ(M)) THEN
            IF(COMPARE(EP_S(IJK,1),ZERO)) THEN
! Solids volume fraction.
               YDOT(Node) = ZERO
               Node = Node + 1
! Solids temperature.
               YDOT(Node) = ZERO
               Node = Node + 1
! Solids phase species mass fractions.
               DO N=1,NMAX(M)
                  YDOT(Node) = ZERO
                  Node = Node + 1
               ENDDO
! Solids phase diameter.
               IF(CALL_GROW) THEN
                  YDOT(Node) = ZERO
                  Node = Node + 1
               ENDIF
            ELSE
! Solids volume fraction: EP_s
               YDOT(Node) = sumlRs(M)/RO_S(M)
               Node = Node + 1
! Temperature: T_s
               YDOT(Node) = -lHORs(M)/(EP_S(IJK,M)*RO_s(M)*C_ps(IJK,M))
               Node = Node + 1
! Species Mass Fraction: X_s
               DO N=1, NMAX(M)
                  YDOT(Node) = (lRs(M,N) - X_s(IJK,M,N)*sumlRs(M)) /   &
                     (EP_S(IJK,M)*RO_s(M))
                  Node = Node + 1
               ENDDO
! Particle Diameter: Dp
               IF(CALL_GROW) THEN
                  YDOT(Node) = (ONE/(3.0d0*EP_S(IJK,M)*RO_S(M))) *     &
                     D_P(IJK,M)*sumlRs(M)
                  Node = Node + 1
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE STIFF_CHEM_RRATES