!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C
!     Module name: DES_INIT_NAMELIST                                      C
!     Purpose: DES - initialize the des-namelist                          C                                                                      C
!                                                                         C
!     Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!     Comments: Added some interpolation based inputs                     C
!                                                                         C
!  Keyword Documentation Format:                                          C
!<keyword category="category name" required="true/false"                  C
!                                    legacy="true/false">                 C
!  <description></description>                                            C
!  <arg index="" id="" max="" min=""/>                                    C
!  <dependent keyword="" value="DEFINED"/>                                C
!  <conflict keyword="" value="DEFINED"/>                                 C
!  <valid value="" note="" alias=""/>                                     C
!  <range min="" max="" />                                                C
!  MFIX_KEYWORD=INIT_VALUE                                                C
!</keyword>                                                               C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_INIT_NAMELIST

      USE param1
      USE discretelement
      USE mfix_pic
      USE des_bc
      USE des_ic
      USE des_thermo
      USE des_rxns
      USE pic_bc

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------

!-----------------------------------------------

      INCLUDE 'desnamelist.inc'

!<keyword category="Discrete Element" required="false">
!  <description>Total number of particles to be read in from the 
! user provided particle configuration file. Only valid when  
! GENER_PART_CONFIG is set to False. When GENER_PART_CONFIG is True,
! then Particle count is  automatically calculated by MFIX.
! For an inflow case beginning with no solids inventory, 
! Particles can be specified as 0.  </description>
!  <range min="0" max="+Inf" />
      PARTICLES = UNDEFINED_I
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>maximum number of neighbors per particle.
!  Relevant only for DEM model. </description>
!  <range min="0" max="+Inf" />
      MN = 10
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Only needed if des_continuum_coupled is True
! and DEM model is used.
! Number of times a pure granular simulation is run before
! the coupled DEM simulation is started. It is used to obtain 
! an initial settled configuration. </description>
!  <range min="0" max="+Inf" />
      NFACTOR = 10
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Run one-way coupled simulations. The fluid does not 
! see the particles in terms of drag force. The effect of particle volume
! is still felt by the fluid through non-unity voidage values. 
! </description>
      DES_ONEWAY_COUPLED = .FALSE. 
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Expand the size of the particle arrays by
! an arbitrary factor (multiple of the number of particles). Serves 
! as a knob to allocate more particles than initially specified in the 
! particle configuration file. </description>
!  <dependent keyword="PARTICLES" value="DEFINED"/>
      PARTICLES_FACTOR = 1.2D0
!  <range min="1.0" max="+Inf" />
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Flag to use triangular facet representation for 
! particle/parcel-wall interactions. If PIC model is used, then 
! this flag is forced to true. For DEM model, it will be made default, 
! but currently setting it to false as cohesion model has not been extended 
! to new routines. Expand the size of the particle arrays by
! an arbitrary factor (multiple of the number of particles). Serves 
! as a knob to allocate more particles than initially specified in the 
! particle configuration file. </description>
      USE_STL_DES  = .false.
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description>Use discrete particle model for solids.
! Must be true to invoke DEM, MPPIC, hybrid models.
! Default is DEM model, to turn on MPPIC model, set MPPIC to 
! true as well. </description>
      DISCRETE_ELEMENT = .FALSE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>To switch between pure granular or  coupled simulations of carried and dispersed phase flows. </description>
!  <valid value=".true." note="Performs coupled simulations. "/>
      DES_CONTINUUM_COUPLED = .FALSE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>use an interpolation suite to calculate the drag force on
! each particle based on particle location rather than cell averages.</description>
      DES_INTERP_ON = .FALSE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>use interpolation to compute dispersed phase average fields,
! such as, solids volume fraction, solids velocity fields. if false, the average fields are obtained by simple arithmetic averaging.</description>
!  <valid value=".true." note="if drag is interpolated 
! (i.e., des_interp_on = .true.), then it is forced to 
! true for backward compatibility. 
! Additionally, if MPPIC or Cut-cells are used (DEM or MPPIC), then 
! also the mean field interpolation is forced. "/>
      DES_INTERP_MEAN_FIELDS = .false.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description> Reports mass based on Lagrangian particles 
! and continuum representation. Useful to ensure mass conservation
! between Lagrangian and continuum representations. Recommended 
! use for debugging purposes. </description>
!  <dependent keyword="DES_INTERP_MEAN_FIELDS" value=".TRUE."/>
      DES_REPORT_MASS_INTERP = .false.
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description>Time stepping scheme (relevant to DEM model only). 
! MPPIC is only limited to EULER scheme. </description>
!  <valid value="EULER" note="First-Order Euler Scheme."/>
!  <valid value="ADAMS BASHFORTH" note="Second order ADAMS BASHFORTH scheme"/>
      DES_INTG_METHOD = 'EULER'
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>switch to turn cohesion model (limited to DEM model only) 
! on and off.</description>
!  <conflict keyword="MPPIC" value=".TRUE."/>
      USE_COHESION = .FALSE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Flag to turn on the use hamaker van der waals forces.</description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      VAN_DER_WAALS = .FALSE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Invoke the hybrid scheme.</description>
      DES_CONTINUUM_HYBRID = .FALSE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Flag to set the neighbor search algorithm
! Relevant to DEM model only.</description>
!  <valid value="1" note="N-Square search algorithm (most expensive)"/>
!  <valid value="2-4" note="Grid-Based Neighbor Search (Recommended)"/>
      DES_NEIGHBOR_SEARCH = 4
!</keyword>
!<keyword category="Discrete Element" required="false">
!  <description> Flag to use van der Hoef et al. (2006) 
! model for adjusting the rotation of the 
! contact plane. See the MFIX-DEM documentation. </description>
      USE_VDH_DEM_MODEL = .false. 
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description>Maximum number of steps through a Dem loop
! before a neighbor search will be performed. 
! (search may be called earlier based on other logic).</description>
!  <range min="0.0" max="+Inf" />
      NEIGHBOR_SEARCH_N = 25
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Ratio of the distance (imaginary sphere radius) 
! to particle radius that is allowed before a neighbor search is performed. This works in conjunction with the logic imposed by 
! NEIGHBOR_SEARCH_N in deciding calls to the neighbor search 
! algorithm. </description>
      NEIGHBOR_SEARCH_RAD_RATIO = 1.0D0
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Effectively increase the radius of a particle 
!(multiple of the sum of particle radii) during the  building of  
! particle neighbor list. Relevant to DEM model only. </description>
      FACTOR_RLM = 1.2
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Number of des grid cells in the I- direction. 
! If left undefined, then it is set by MFIX such that its size 
! equals three times the maximum particle diameter with a minimum of 
! 1 cell.   </description>
      DESGRIDSEARCH_IMAX = UNDEFINED_I
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Number of des grid cells in the J- direction. 
! If left undefined, then it is set by MFIX such that its size 
! equals three times the maximum particle diameter with a minimum of 
! 1 cell.   </description>
      DESGRIDSEARCH_JMAX = UNDEFINED_I
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description>Number of des grid cells in the K- direction. 
! If left undefined, then it is set by MFIX such that its size 
! equals three times the maximum particle diameter with a minimum of 
! 1 cell.   </description>
      DESGRIDSEARCH_KMAX = UNDEFINED_I
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description>Collision model for the soft-sphere approach used in 
! DEM model. All models require specifying the following parameters:
! des_en_input, des_en_wall_input, mew, and mew_w. </description>
!  <valid value="LSD" note="The linear spring-dashpot model. Requires 
! additional initialization for kn, kn_w, kt_fac, kt w_fac, des_etat_fac, 
! & des_etat_w_fac. "/>
!  <valid value="HERTZIAN" note="The Hertzian model. Requires 
! additional initialization for  des_et_input, des_et_wall_input, e_young, 
! ew_young, v_poisson, & vw_poisson. "/>
      DES_COLL_MODEL = 'LSD'
!</keyword>

! particle properties
!<keyword category="Discrete Element" required="false">
!  <description>Normal spring constant for inter-particle collisions 
!  needed when using the default (linear spring-dashpot) 
! collision model.</description>
      KN = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Ratio of the tangential spring constant to
! normal spring constant for inter-particle collisions. 
! Use it to specify the tangential spring constant for 
! particle-particle collisions as kt_fac*kn. Required 
! when using the default (linear spring-dashpot) collision 
! model.</description>
!  <dependent keyword="DES_COLL_MODEL" value="LSD"/>
!  <range min="0.0" max="1.0" />
      KT_FAC = 2.d0/7.d0
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Normal spring constant for particle-wall 
! collisions. Needed when using the default (linear spring-dashpot)
! collision model.</description>
      KN_W = UNDEFINED
!</keyword>

   
!<keyword category="Discrete Element" required="false">
! <description>Ratio of the tangential spring constant 
! to normal spring constant for particle-wall collisions. 
! Use it to specify the tangential spring constant for
! particle-wall collisions as kt_w_fac*kn_w. 
! Needed when using the default (linear spring-dashpot) 
! collision model.</description>! 
! <dependent keyword="DES_COLL_MODEL" value="LSD"/>
!  <range min="0.0" max="1.0" />
      KT_W_FAC = 2.d0/7.d0
!</keyword>

!<keyword category="Discrete Element" required="false">
! <description>Inter-particle Columb friction coefficient 
! required for DEM model.</description>
! <range min="0.0" max="1.0" />
      MEW = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
! <description>Particle-wall friction coefficient required for 
! DEM model.</description>
! <range min="0.0" max="1.0" />
      MEW_W = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
! <description>The normal restitution coefficient for
! interparticle collisions that is used to determine the 
! inter-particle normal damping factor. Values are stored as a one
! dimensional array (see MFIX-DEM doc). So
! if MAX=3, then 6 values are needed, which
! are defined as follows: en11 en12 en13 en22
! en23 en33.</description>
! <range min="0.0" max="1.0" />
      DES_EN_INPUT(:) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element" required="false">
! <description>Normal restitution coefficient for particle
! wall collisions that is used to determine the
! particle-wall normal damping factor (see
! cfassign.f for details). Values are stored as
! a one dimensional array. So, if MMAX=3,
! then 3 values are needed, which are
! defined as follows: enw1 enw2 enw3.
! </description>
! <range min="0.0" max="1.0" />
      DES_EN_WALL_INPUT(:) = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
! <description>Tangential restitution coefficient for
! interparticle collisions. Values are stored as
! a one dimensional array. Only needed
! when using the Hertzian collision model.
! </description>
! <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
! <range min="0.0" max="1.0" />
      DES_ET_INPUT(:) = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
! <description>Tangential restitution coefficient for particle
! wall collisions. Values are stored as a one
! dimensional array. Only needed when using
! the Hertzian collision model.
! </description>
! <range min="0.0" max="1.0" />
! <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      DES_ET_WALL_INPUT(:) = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
! <description>Ratio of the tangential damping factor to the 
! normal damping factor for inter-particle collisions. 
! Required for the linear spring-dashpot model for soft-spring collision 
! modelling under DEM. 
! For the Hertzian model, the tangential damping coefficients 
! have to be explicity specified and this variable is not 
! required. </description>
! <dependent keyword="DES_COLL_MODEL" value="LSD"/>
! <range min="0.0" max="1.0" />
! <valid value="UNDEFINED" note="For LSD model, if left undefined, MFIX
! will will revert to default value of 0.5" />
      DES_ETAT_FAC = UNDEFINED 
!</keyword>

!<keyword category="Discrete Element" required="false">
! <description>Ratio of the tangential damping factor to the 
! normal damping factor for particle-wall collisions. 
! Required for the linear spring-dashpot model for soft-spring collision 
! modelling under DEM. 
! For the Hertzian model, the tangential damping coefficients 
! have to be explicity specified and specification of this 
! variable is not required. </description>
! <dependent keyword="DES_COLL_MODEL" value="LSD"/>
! <range min="0.0" max="1.0" />
! <valid value="UNDEFINED" note="For LSD model, if left undefined, MFIX
! will will revert to default value of 0.5" />
      DES_ETAT_W_FAC = UNDEFINED 
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description>Youngs modulus for the wall. Needed when using the Hertzian
! spring-dashpot model for soft-spring collision modelling under DEM.
!</description>
!  <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      EW_YOUNG = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Poisson ratio for the wall. Needed when using the Hertzian
! spring-dashpot model for soft-spring collision modelling under DEM.
!</description>
!  <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      VW_POISSON = UNDEFINED
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description>Youngs modulus for the particle. Needed when using the Hertzian
! spring-dashpot model for soft-spring collision modelling under DEM.
!</description>
!  <arg index="1" id="Discerte solid phase" min="1" max="DES_MMAX"/>
!  <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      E_YOUNG(:DIM_M) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description>Poissons ratio for the particle. Needed when using the Hertzian
! spring-dashpot model for soft-spring collision modelling under DEM.
!</description>
!  <arg index="1" id="Discerte solid phase" min="1" max="DES_MMAX"/>
!  <dependent keyword="DES_COLL_MODEL" value="HERTZIAN"/>
      V_POISSON(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Particle diameter associated with discrete solids 
! phase M</description>
!  <arg index="1" id="Discerte solid phase" min="1" max="DES_MMAX"/>
      DES_D_P0(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Particle density associated with discrete solids phase M</description>
!  <arg index="1" id="Discerte solid phase" min="1" max="DES_MMAX"/>
      DES_RO_S(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Number of solid phases with discrete 
! representation. When DES_CONTINUUM_HYBRID is false, no need to set 
! this variable as it is set equal to MMAX. Only relevant 
! when DES_CONTINUUM_HYBRID  is true so as to differentiate 
! between discrete and continuum representation of solids phase
!</description>
!  <dependent keyword="DES_CONTINUUM_HYBRID" value=".TRUE."/>
      DES_MMAX = UNDEFINED_I
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Lees-Edwards boundary condition to simulate
! homogeneous shear problem with periodic boundary conditions. 
! Not supported in this version. </description>
      DES_LE_BC = .FALSE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Relative velocity needed for Lees-Edwards BC.
! Not supported in this version. </description>
!  <dependent keyword="DES_LE_BC" value=".TRUE."/>
      DES_LE_REL_VEL = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Direction of shear for Lees-Edwards BC. 
! Not supported in this version. </description>
!  <dependent keyword="DES_LE_BC" value=".TRUE."/>
      DES_LE_SHEAR_DIR = UNDEFINED_C
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description>Maximum number of particles that may exist
! within a simulation. This quantity is used for calculating
! the size of arrays for allocation.</description>
      MAX_PIS = UNDEFINED_I
!</keyword>

!<keyword category="Boundary Condition" required="false">
!  <description>logical to force the inlet to operate with an ordered boundary condition.  this may be useful during long simulations or if the inlet appears to be taking a long time to randomly place particles.</description>
      FORCE_ORD_BC = .FALSE.
!</keyword>

! for cohesion: van der waals
!<keyword category="Discrete Element" required="false">
!  <description>hamaker constant used for particle-particle interactions</description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      HAMAKER_CONSTANT = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>hamaker constant used in particle-wall interactions.</description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      WALL_HAMAKER_CONSTANT = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>maximum separation distance above which van der waals forces are not implemented.</description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      VDW_OUTER_CUTOFF = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>minimum separation distance below which van der waals forces are calculated using a surface adhesion model.</description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      VDW_INNER_CUTOFF = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>maximum separation distance above which van der waals forces are not implemented (particle-wall interactions).</description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      WALL_VDW_OUTER_CUTOFF = ZERO 
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>minimum separation distance below which van der waals forces are calculated using a surface adhesion model (particle-wall interactions).</description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      WALL_VDW_INNER_CUTOFF = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Mean radius of surface asperities that influence the cohesive force following a model by rumpf (1990).</description>
!  <dependent keyword="USE_COHESION" value=".TRUE."/>
      Asperities = ZERO
!</keyword>

! J.Musser : des particle initial conditions
      DES_IC_X_w(:) = UNDEFINED
      DES_IC_X_e(:) = UNDEFINED
      DES_IC_Y_s(:) = UNDEFINED
      DES_IC_Y_n(:) = UNDEFINED
      DES_IC_Z_b(:) = UNDEFINED
      DES_IC_Z_t(:) = UNDEFINED
      DES_IC_T_s(:,:) = UNDEFINED
      DES_IC_X_s(:,:,:) = UNDEFINED

!<keyword category="Discrete Element" required="false">
!  <description>Allows writing of discrete particle data to output 
! files. 
! Relevant to both granular and coupled simulations. </description>
      PRINT_DES_DATA = .FALSE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>If print_des_data, this is the frequency at
! which particle data is written out. This only applies to 
! pure granular simulations. For coupled simulation, the 
! output frequency is controlled by spx_dt(1).</description>
!  <dependent keyword="PRINT_DES_DATA" value=".True."/>
!  <dependent keyword="DES_CONTINUUM_COUPLED" value=".False."/>
!  <conflict keyword="DES_CONTINUUM_COUPLED" value=".True."/>
!  <conflict keyword="MPPIC" value=".True."/>
      DES_SPX_DT = LARGE_NUMBER
!</keyword>

!<keyword category="Discrete Element" required="false">
! <description>This is the frequency 
! at which des.res and .res files will be written. 
! This only applies to pure granular simulations. 
! For coupled simulation the restart frequency is controlled
! by Res_dt.</description>
!  <dependent keyword="DES_CONTINUUM_COUPLED" value=".False."/>
!  <conflict keyword="DES_CONTINUUM_COUPLED" value=".True."/>
!  <conflict keyword="MPPIC" value=".True."/>
      DES_RES_DT = LARGE_NUMBER
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>If undefined the default vtp files are written.</description>
!  <valid value="TECPLOT" note="If tecplot is specified then several .dat files are written including DES_DATA.dat & AVG_EPS.dat. See write_des_data.f for details."/>
      DES_OUTPUT_TYPE = UNDEFINED_C
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Print out additional information from dem model.</description>
      DEBUG_DES = .FALSE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Specify particle number for particle level debugging 
! details. </description>
      FOCUS_PARTICLE = 0
!</keyword>


! flag whether to invoke cluster identification algorithm
!<keyword category="Discrete Element" required="false">
!  <description>Flag to turn on on-the-fly cluster data calculations. 
! </description>
!  <conflict keyword="MPPIC" value=".TRUE."/>
      DES_CALC_CLUSTER = .FALSE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>The physical distance relative to a particle for 
! peforming cluster statistics calculations.</description>
!  <conflict keyword="MPPIC" value=".TRUE."/>
!  <dependent keyword="DES_CALC_CLUSTER" value=".TRUE."/>
      CLUSTER_LENGTH_CUTOFF = UNDEFINED
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description>Automatically generate an initial
! particle configuration (position), otherwise use particle input.dat
! Now uses the conventional flag used by continuum solver. Once defined, 
! this feature will determine the total number of particles in the 
! system and their initial placement. 
! The particle_input.dat file (if present in the run directory)
! is ignored.</description>
      GENER_PART_CONFIG = .FALSE.
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description>Turn on snider's version of frictional model. Does not run very stably.</description>
      MPPIC_SOLID_STRESS_SNIDER = .false.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description> First coefficient of restitution for the frictional stress 
! model in the MPPIC model. See the MPPIC documentation for more details.
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      MPPIC_COEFF_EN1 = UNDEFINED
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description> Second coefficient of restitution for the frictional stress 
! model in the MPPIC model. See the MPPIC documentation for more details.
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      MPPIC_COEFF_EN2 = UNDEFINED
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description> Normal coefficient of restitution for parcel-wall collisions
! in the MPPIC model. 
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      MPPIC_COEFF_EN_WALL = UNDEFINED
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description> Tangential coefficient of restitution for
! parcel-wall collisions in the MPPIC model. 
! Currently not implemented in the code. 
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      MPPIC_COEFF_ET_WALL = 1.0
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description> Turn on the implicit treatment for interphase drag force. 
! Valid only for MPPIC model.. 
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      MPPIC_PDRAG_IMPLICIT = .false.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description> Variable to decide if special treatment is needed or not
!in the direction of gravity in the frictional stress tensor
! Valid only for MPPIC model. See the MPPIC documenation. 
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      MPPIC_GRAV_TREATMENT = .true.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description> Flag to print processor level parcel seeding statistics
! for inflow BC with PIC model. 
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      PIC_REPORT_SEEDING_STATS = .false.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description> Flag to print processor level parcel deletion statistics
! for outflow BC with PIC model. Not recommended for production runs. 
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      PIC_REPORT_DELETION_STATS = .false.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description> A run time flag to report minimum value and location
! of gas voidage. This is useful only for debugging and is not recommended for
! production runs. 
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      PIC_REPORT_MIN_EPG = .false. 
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>P_s term in the frictional stress model of Snider</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      PSFAC_FRIC_PIC = 100
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Beta term in the frictional stress model of Snider</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      FRIC_EXP_PIC = 2.5
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>non-singularity term (epsilon) in the frictional stress model of snider</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      FRIC_NON_SING_FAC = 1E-07
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>CFL number used to decide maximum time
! step size for parcels evolution equations.
! Relevant to MPPIC model only. 
!</description>
!  <dependent keyword="MPPIC" value=".TRUE."/>
      CFL_PIC = 0.1
!</keyword>

! J.Musser : des energy equations
!<keyword category="Discrete Element" required="false">
!  <description>solve energy equations for dem.</description>
      DES_ENERGY_EQ = .FALSE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>include particle-gas convection in dem heat transfer model.</description>
      DES_CONV_EQ = .TRUE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>include particle-particle conduction in dem heat transfer model.</description>
      DES_COND_EQ = .TRUE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>include particle--particle radiation in dem heat transfer model.</description>
      DES_RADI_EQ = .FALSE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>include particle-fluid-particle conduction in dem heat transfer model. (requires des_cond_eq set to true.)</description>
      DES_COND_EQ_PFP = .TRUE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>include particle-particle contact conduction in dem heat transfer model. (requires des_cond_eq set to true.)</description>
      DES_COND_EQ_PP  = .TRUE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>specify the nusselt number correlation used for particle-gas convection. (only ranz_1952 is presently included.)</description>
      DES_CONV_CORR = 'RANZ_1952'
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>minimum separation distance between the surfaces of two contacting particles.</description>
      DES_MIN_COND_DIST = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>fluid lens proportion constant used to calculate the radius of the fluid lens that surrounds a particle. used in the particle-fluid-particle conduction model.</description>
      FLPC = 1.0d0/5.0d0
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Emissivity of solids phase M.</description>
!  <arg index="1" id="Discerte solid phase" min="1" max="DES_MMAX"/>
      DES_Em(:DIM_M) = UNDEFINED
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Flag to turn on species equations in the discrete 
!solids phase M.</description>
!  <arg index="1" id="Discerte solid phase" min="1" max="DES_MMAX"/>
      DES_SPECIES_EQ(:DIM_M) = .FALSE.
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Total number of discrete solids species for each phase
!</description>
!  <arg index="1" id="Discerte solid phase" min="1" max="DES_MMAX"/>
      DES_NMAX_s(:DIM_M) = UNDEFINED_I
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Total number of discrete solids species for each phase
!</description>
!  <arg index="1" id="Discerte solid phase" min="1" max="DES_MMAX"/>
!  <arg index="2" id="Number of species" min="1" max="DES_NMAX_sX"/>
      DES_MW_s(:DIM_M, :DIM_N_s) = UNDEFINED
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description>Name of solids phase M, species N as it appears in the
! materials database. 
!</description>
!  <arg index="1" id="Discerte solid phase" min="1" max="DES_MMAX"/>
!  <arg index="2" id="Number of species" min="1" max="DES_NMAX_sX"/>
      DES_SPECIES_s(:DIM_M, :DIM_N_s) = UNDEFINED_C
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Name of solids phase M, species N as it appears in the
! materials database. 
!</description>
!  <arg index="1" id="Discerte solid phase" min="1" max="DES_MMAX"/>
!  <arg index="2" id="Number of species" min="1" max="DES_NMAX_sX"/>
      DES_SPECIES_ALIAS_s(:DIM_M, :DIM_N_s) = UNDEFINED_C
!</keyword>

!<keyword category="Discrete Element" required="false">
!  <description>Reaction model used to calculate the effects of a gas-solids reaction on a particle.</description>
!  <valid value="variable_density" note="constant particle diameter. particle diameter changes to accommodate a loss/gain in mass."/>
!  <valid value="shrinking_particle" note="constant particle density. particles diameter changes to accommodate a loss/gain in mass."/>
      REACTION_MODEL = 'VARIABLE_DENSITY'
!</keyword>


!<keyword category="Discrete Element" required="false">
!  <description> ! Particle shape factor.
!  </description> 
      DES_PSI_s(:DIM_M) = UNDEFINED
!</keyword>



! des wall boundaries: wall velocities. I think they probably 
! defined for the Lees-Edwards BC's
      DES_BC_Uw_s(:,:) = ZERO
      DES_BC_Vw_s(:,:) = ZERO
      DES_BC_Ww_s(:,:) = ZERO


! These need to be inialized to 0, but they are not part of the namelist
      VTP_FINDEX = 0
      TECPLOT_FINDEX = 0

! not a well supported feature and not generic either. So removing 
! from namelists
      DES_CALC_BEDHEIGHT = .FALSE.
      RETURN
      END SUBROUTINE DES_INIT_NAMELIST
