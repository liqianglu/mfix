!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_CHECK_PARTICLE                                     !
!                                                                      !
!  Purpose:  This routine is used to check if a new particle has fully 
!  entered the domain.  If so, the flag classifying the particle as new
!  is removed, allowing the particle to respond to contact forces from 
!  walls and other particles.                                          
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_CHECK_PARTICLE

      USE compar
      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param1
      USE physprop

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER I,J,K,IJK       ! Unused indices at this time
      INTEGER PC              ! Loop counter
      INTEGER NP              ! Particle Index
      INTEGER NPT             ! Temp Particle Index
! Tmp holder of particle position
      DOUBLE PRECISION XPOS, YPOS, ZPOS

! Loop indices used for clearing forces associated with exiting particles
      INTEGER NLIMNP, NLIM, NEIGHNP, NLIM_NEIGHNP
      
! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG
!-----------------------------------------------
      INCLUDE 'function.inc'


      DES_LOC_DEBUG = .FALSE.

      PC = 1
      DO NP = 1, MAX_PIS
         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(NP,1)) CYCLE

         XPOS = DES_POS_NEW(NP,1)
         YPOS = DES_POS_NEW(NP,2)
         IF (DIMN .EQ. 3) ZPOS = DES_POS_NEW(NP,3)
! check the status of a new injected particle. if the new particle
! fully entered the domain change its 'new' status PEA(NP,2) to false
         IF(PEA(NP,2))THEN
            IF(DIMN == 2 .AND. &
               (XPOS - DES_RADIUS(NP)) .GT. 0 .AND.&
               (XPOS + DES_RADIUS(NP)) .LT. XLENGTH.AND.&
               (YPOS - DES_RADIUS(NP)) .GT. 0 .AND.&
               (YPOS + DES_RADIUS(NP)) .LT. YLENGTH &
            .OR. &
               DIMN == 3 .AND. &
               (XPOS - DES_RADIUS(NP)) .GT. 0 .AND.&
               (XPOS + DES_RADIUS(NP)) .LT. XLENGTH.AND.&
               (YPOS - DES_RADIUS(NP)) .GT. 0 .AND.&
               (YPOS + DES_RADIUS(NP)) .LT. YLENGTH .AND. &
               (ZPOS - DES_RADIUS(NP)) .GT. 0 .AND.&
               (ZPOS + DES_RADIUS(NP)) .LT. ZLENGTH)THEN

               PEA(NP,2) = .FALSE.
               PC = PC + 1

            ENDIF  ! end if particle has fully entered the system

! check the status of an exiting particle (a particle that has come into
! contact with a mass outlet).
         ELSEIF(PEA(NP,3))THEN
            IF(DIMN == 2 .AND. &
               (XPOS - DES_RADIUS(NP)) .GT. 0 .AND.&
               (XPOS + DES_RADIUS(NP)) .LT. XLENGTH.AND.&
               (YPOS - DES_RADIUS(NP)) .GT. 0 .AND.&
               (YPOS + DES_RADIUS(NP)) .LT. YLENGTH &
            .OR. &
               DIMN == 3 .AND. &
               (XPOS - DES_RADIUS(NP)) .GT. 0 .AND.&
               (XPOS + DES_RADIUS(NP)) .LT. XLENGTH.AND.&
               (YPOS - DES_RADIUS(NP)) .GT. 0 .AND.&
               (YPOS + DES_RADIUS(NP)) .LT. YLENGTH .AND. &
               (ZPOS - DES_RADIUS(NP)) .GT. 0 .AND.&
               (ZPOS + DES_RADIUS(NP)) .LT. ZLENGTH)THEN
! if an exiting particle has altered course so that it is moving back
! into the domain change (no longer on track to exit the domain) and 
! it has fully entered the domain then change its 'exiting' status 
! PEA(NP,3) to false

               PEA(NP,3) = .FALSE.
               PC = PC + 1

            ELSEIF(DIMN == 2 .AND. &
               (XPOS + DES_RADIUS(NP)) .LE. 0 .OR.&
               (XPOS - DES_RADIUS(NP)) .GE. XLENGTH .OR.&
               (YPOS + DES_RADIUS(NP)) .LE. 0 .OR.&
               (YPOS - DES_RADIUS(NP)) .GT. YLENGTH &
            .OR. &
               DIMN == 3 .AND. &
               (XPOS + DES_RADIUS(NP)) .LE. 0 .OR.&
               (XPOS - DES_RADIUS(NP)) .GE. XLENGTH.OR.&
               (YPOS + DES_RADIUS(NP)) .LE. 0 .OR.&
               (YPOS - DES_RADIUS(NP)) .GE. YLENGTH .OR. &
               (ZPOS + DES_RADIUS(NP)) .LE. 0 .OR.&
               (ZPOS - DES_RADIUS(NP)) .GE. ZLENGTH)THEN
! if an exiting particle has fully exited the domain, effectively
! remove the particle (reset PEA(NP,:) to false)       

               PEA(NP,:) = .FALSE.
              
               DES_POS_OLD(NP,:) = ZERO
               DES_POS_NEW(NP,:) = ZERO
               DES_VEL_OLD(NP,:) = ZERO
               DES_VEL_NEW(NP,:) = ZERO
               OMEGA_OLD(NP,:) = ZERO
               OMEGA_NEW(NP,:) = ZERO
               DES_RADIUS(NP) = ZERO
               PMASS(NP) = ZERO
               PVOL(NP) = ZERO
               RO_Sol(NP) = ZERO
               OMOI(NP) = ZERO

               FC(NP,:) = ZERO
               FN(NP,:) = ZERO
               FT(NP,:) = ZERO
               TOW(NP,:) = ZERO

               PN(NP,:) = -1
               PN(NP,1) = 0
               PV(NP,:) = 1
               PFT(NP,:,:) = ZERO
               PPOS(NP,:) = ZERO

! Note that if particle NP has any neighbors then the particle NP will
! still exist in the neighbor's neighbours list.  This information would
! eventually be cleared in the second call to calc_force_des following
! the removal of NP but is forceably removed here to keep the dem
! inlet/outlet code self contained (does not rely on other code)
               NEIGHBOURS(NP,:) = -1
               NEIGHBOURS(NP,1) = 0 
               
! Clear particle NP from any other neighboring particles lists               
               IF (NEIGHBOURS(NP,1) > 0) THEN
                  NLIMNP = NEIGHBOURS(NP,1)+1
                
! Cycle through all neighbours of particle NP
                  DO I = 2, NLIMNP
                     NEIGHNP = NEIGHBOURS(NP,I)

! If any neighbor particle has a lower index than NP then the contact 
! force history will be stored with that particle and needs to be cleared
                     IF (NEIGHNP < NP) THEN
                        IF (PN(NEIGHNP,1) > 0) THEN
                           NLIM = PN(NEIGHNP,1)+1
                           DO J = 2, NLIM
                             NPT = PN(NEIGHNP,J)
                             IF (NPT .NE. NP) CYCLE   ! find particle NP in NEIGHNP list
                             PN(NEIGHNP,J:(MAXNEIGHBORS-1)) = &
                                PN(NEIGHNP,(J+1):MAXNEIGHBORS)
                             PV(NEIGHNP,J:(MAXNEIGHBORS-1)) = &
                                PV(NEIGHNP,(J+1):MAXNEIGHBORS)
                             PFT(NEIGHNP,J:(MAXNEIGHBORS-1),:) = &
                                PFT(NEIGHNP,(J+1):MAXNEIGHBORS,:)
                             PN(NEIGHNP,1) = PN(NEIGHNP,1) -1
                           ENDDO
                        ENDIF
                     ENDIF

                     NLIM_NEIGHNP = NEIGHBOURS(NEIGHNP,1)+1
! Find where particle NP is stored in its neighbours (NEIGHNP) lists
! and remove particle NP from the list
                     DO J = 2, NLIM_NEIGHNP
                        NPT = NEIGHBOURS(NEIGHNP,J)
                        IF (NPT .NE. NP) CYCLE  ! find particle NP in NEIGHNP list
                        NEIGHBOURS(NEIGHNP,1) = NEIGHBOURS(NEIGHNP,1)-1
                        NEIGHBOURS(NEIGHNP,J:(MN-1)) = NEIGHBOURS(NEIGHNP,(J+1):MN)
                     ENDDO
                  ENDDO
               ENDIF  

               PIS = PIS - 1
! In the exiting case do not increment pc

            ENDIF   ! endif particle has reentered/exited of the domain

         ELSE
            PC = PC + 1 
         ENDIF   ! endif PEA(NP,3) or PEA(NP,2)

      ENDDO   ! loop over NP=1,MAX_PIS

      IF (DES_LOC_DEBUG) WRITE(*,1001)

 1000 FORMAT(3X,'---------- START DES_CHECK_PARTICLE ---------->')
 1001 FORMAT(3X,'<---------- END DES_CHECK_PARTICLE ----------') 

! Can be used for debugging purposes      
 9000 FORMAT(F9.4,2X,'NOT NEW  ',I4,1X,6(F9.4,1X))
 9001 FORMAT(F9.4,2X,'RE-ENTER ',I4,1X,6(F9.4,1X))
 9002 FORMAT(F9.4,2X,'EXITED   ',I4)

      RETURN
      END SUBROUTINE DES_CHECK_PARTICLE
