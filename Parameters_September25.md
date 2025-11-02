Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Below are parameter values used for different simulations, as well as initial conditions used in function animation. Some simulations use variants of the published code. The list is going to be updated gradually. 


### 30 September 2025 - Changing the anion/cation proportion in a central gravity field ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.9
#define INITXMAX 1.9	/* x interval for initial condition */
#define INITYMIN -0.9
#define INITYMAX 0.9	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -1.9
#define ADDXMAX 1.9	/* x interval for adding particles */
#define ADDYMIN 1.2
#define ADDYMAX 1.3	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.2
#define BCXMAX 2.2	/* x interval for boundary condition */
#define BCYMIN -2.2
#define BCYMAX 2.2	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 0  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.5 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.014 	    /* parameter controlling radius of particles */
#define MU_B 0.014           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.014         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.014       /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 40           /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 4400      /* number of frames of movie */
#define NVID 100         /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 9        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 1           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 50000.0      /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 2.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 250.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 4.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 0.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 100.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 10000.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 1     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define CHARGE_ADD_B -1.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 1      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 0        /* time at which to add first particle */
#define ADD_PERIOD 6      /* time interval between adding further particles */
#define ADD_TYPE 4         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 800  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 250 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 50.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 23            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 2                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 6           /* number of types shown in graph */
#define RD_INITIAL_COND 51        /* initial condition of particles */
#define REACTION_DIST 2.5         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e10      /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 4              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 1      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.0        /* min proportion of type 1 particles */
#define PROP_MAX 1.0        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 50      /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 29 September 2025 - Can one trap a beam of light in a whispering gallery mode resonator? ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 1      /* set to 1 for a variable index of refraction */
#define IOR 17              /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 98         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.2              /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define WALL_WIDTH_B 0.05   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 1         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 42  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.05       /* defines oscilling beam range */
#define OSCIL_YMID -0.95      /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.009        /* frequency of periodic excitation */
#define AMPLITUDE 2.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.2        /* Courant number in medium B */
#define COURANTB 0.1       /* Courant number */
#define GAMMA 5.0e-6          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 25.0  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 750            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2600         /* number of frames of movie */
#define NVID 12              /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 400      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.5              /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00025  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.015  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 18     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 100.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.7      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.6    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 28 September 2025 - 20% heavy cations of charge 4 and 80% light anions of charge 1 in a central gravity field ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -0.1
#define INITXMAX 0.1	/* x interval for initial condition */
#define INITYMIN -0.1
#define INITYMAX 0.1	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -1.9
#define ADDXMAX 1.9	/* x interval for adding particles */
#define ADDYMIN 1.2
#define ADDYMAX 1.3	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.2
#define BCXMAX 2.2	/* x interval for boundary condition */
#define BCYMIN -2.2
#define BCYMAX 2.2	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 0  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.1 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.01           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.02         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.01       /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 1           /* number of grid point for grid of disks */
#define NGRIDY 1           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 3900      /* number of frames of movie */
#define NVID 100         /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 9        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 50000.0      /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 2.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 250.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 4.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 4.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 0.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 200.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 10000.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 1     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 4.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 4.0   /* charge of added particles */
#define CHARGE_ADD_B -1.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_REGION 1      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 0        /* time at which to add first particle */
#define ADD_PERIOD 6      /* time interval between adding further particles */
#define ADD_TYPE 4         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 800  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 1     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.8    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 250 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 50.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 1      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 23            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 2                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 6           /* number of types shown in graph */
#define RD_INITIAL_COND 51        /* initial condition of particles */
#define REACTION_DIST 2.5         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e10      /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 4              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 50      /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 27 September 2025 - Classics revisited: Noise-reflecting panels ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:**

```
    init_wave_flat(phi, psi, xy_in);
    
    /* add oscillating waves */
    wave_source_x[0] = -1.0;
    wave_source_y[0] = 0.0;
    source_periods[0] = OSCILLATING_SOURCE_PERIOD;
    source_amp[0] = INITIAL_AMP;
    for (source = 0; source < N_SOURCES; source++)
    {
        dperiod = source_periods[source];
        phase = i - (int)(dperiod*(double)((int)((double)i/dperiod)));
        if ((ADD_OSCILLATING_SOURCE)&&(phase == 1)&&(i<MAX_PULSING_TIME))
        {
            printf("Source %i: Adding pulse %i\n", source, add_counter[source]);
            add_counter[source]++;
            if (ALTERNATE_OSCILLATING_SOURCE) sign[source] = -sign[source];
            add_circular_wave(-sign[source]*source_amp[source], wave_source_x[source], wave_source_y[source], phi, psi, xy_in);
        }
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 44         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 0.1	    /* parameter controlling the dimensions of domain */
#define MU 0.2              /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define WALL_WIDTH_B 0.05   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 44  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.1       /* defines oscilling beam range */
#define OSCIL_YMID -0.75      /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.002         /* frequency of periodic excitation */
#define AMPLITUDE 2.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.005      /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number in medium B */
#define COURANTB 0.05      /* Courant number */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 25.0  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 1800         /* number of frames of movie */
#define NVID 8              /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.5              /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00025  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.015  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 15     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 100.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.7      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.6    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 26 September 2025 - 33% heavy cations of charge 2 and 67% light anions of charge 1 in a central gravity field ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -0.1
#define INITXMAX 0.1	/* x interval for initial condition */
#define INITYMIN -0.1
#define INITYMAX 0.1	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -1.9
#define ADDXMAX 1.9	/* x interval for adding particles */
#define ADDYMIN 1.2
#define ADDYMAX 1.3	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.2
#define BCXMAX 2.2	/* x interval for boundary condition */
#define BCYMIN -2.2
#define BCYMAX 2.2	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 0  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.1 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.017 	    /* parameter controlling radius of particles */
#define MU_B 0.012          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.017         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.012         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 1           /* number of grid point for grid of disks */
#define NGRIDY 1           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 3500      /* number of frames of movie */
#define NVID 100         /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 9        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 50000.0      /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 2.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 250.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 0.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 200.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 10000.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 1     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 2.0   /* charge of added particles */
#define CHARGE_ADD_B -1.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_REGION 1      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 0        /* time at which to add first particle */
#define ADD_PERIOD 6      /* time interval between adding further particles */
#define ADD_TYPE 4         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 800  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 1     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.666666667    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 25.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_/*PARTICLE_MASS*/ 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 1      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 23            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 2                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 6           /* number of types shown in graph */
#define RD_INITIAL_COND 51        /* initial condition of particles */
#define REACTION_DIST 2.5         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e10      /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 2              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 50      /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 25 September 2025 - A linear wave crossing a disc with refractive index 2 ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -1.5
#define XMAX 2.5	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 3          /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 0.6	    /* parameter controlling the dimensions of domain */
#define MU 1.15             /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define WALL_WIDTH_B 0.05   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 1         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 1.0       /* defines oscilling beam range */
#define OSCIL_YMID -10.0      /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.015         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.05      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12.5  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2500         /* number of frames of movie */
#define NVID 14             /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 150    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.08       /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 15     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 50.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 0.5      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 24 September 2025 - How does the composition of a U-Cd alloy affect its reaction when bombarded with neutrons? ###

**Program:** `lennardjones.c` 

**Part 1:**

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.8
#define INITXMAX 1.8	/* x interval for initial condition */
#define INITYMIN -1.01
#define INITYMAX 1.05	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.0
#define ADDXMAX -1.9	/* x interval for adding particles */
#define ADDYMIN -0.2
#define ADDYMAX 0.2	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.49 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.015          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.01         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.01       /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 80           /* number of grid point for grid of disks */
#define NGRIDY 40           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2200      /* number of frames of movie */
#define NVID 15          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 5     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 5        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 31         /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 31        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 1.0e8         /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 340.0       /* hue of particles of type 1 */
#define HUE_TYPE2 250.0      /* hue of particles of type 2 */
#define HUE_TYPE3 310.0      /* hue of particles of type 3 */
#define HUE_TYPE4 60.0       /* hue of particles of type 4 */
#define HUE_TYPE5 210.0       /* hue of particles of type 5 */
#define HUE_TYPE6 150.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 6.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 6.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 10.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 235.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 112.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 1.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 500.0        /* initial velocity range */
#define V_INITIAL_ADD 2000.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 0.01   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 30.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 0.0        /* charge of particles of first type */
#define CHARGE_B 0.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 50        /* time at which to add first particle */
#define ADD_PERIOD 100      /* time interval between adding further particles */
#define ADD_TYPE 4         /* type of added particles */
#define N_ADD_PARTICLES 2  /* number of particles to add */
#define FINAL_NOADD_PERIOD 1000  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 250 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 150.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 2              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 1      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 271            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 6           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 3.1         /* maximal distance for reaction to occur */
#define REACTION_PROB 0.1         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.5     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 160.0         /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 15       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 1              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 1      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60      /* size of hashgrid in x direction */
#define HASHY 30      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

**Part 2:**

```
#define TYPE_PROPORTION 0.51 /* proportion of particles of first type */

```

**Part 3:**

```
#define TYPE_PROPORTION 0.55 /* proportion of particles of first type */

```

### 23 September 2025 - A linear wave crossing a disc with refractive index 1.33 ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -1.5
#define XMAX 2.5	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 3          /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 0.6	    /* parameter controlling the dimensions of domain */
#define MU 1.15             /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define WALL_WIDTH_B 0.05   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 1         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 1.0       /* defines oscilling beam range */
#define OSCIL_YMID -10.0      /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.015         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.075     /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12.5  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2600         /* number of frames of movie */
#define NVID 14             /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 150    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.08       /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 11     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 14   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 50.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 0.5      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 22 September 2025 - Melting a Uranium-Cadmium alloy by bombarding it with neutrons ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.8
#define INITXMAX 1.8	/* x interval for initial condition */
#define INITYMIN -1.01
#define INITYMAX 1.05	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.0
#define ADDXMAX -1.9	/* x interval for adding particles */
#define ADDYMIN -0.2
#define ADDYMAX 0.2	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.51 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.015          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.01         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.01       /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 80           /* number of grid point for grid of disks */
#define NGRIDY 40           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2400      /* number of frames of movie */
#define NVID 15          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 5     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 5        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 0         /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 31        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 1.0e8         /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 340.0       /* hue of particles of type 1 */
#define HUE_TYPE2 250.0      /* hue of particles of type 2 */
#define HUE_TYPE3 310.0      /* hue of particles of type 3 */
#define HUE_TYPE4 60.0       /* hue of particles of type 4 */
#define HUE_TYPE5 210.0       /* hue of particles of type 5 */
#define HUE_TYPE6 150.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 6.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 6.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 10.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 235.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 112.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 1.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 500.0        /* initial velocity range */
#define V_INITIAL_ADD 2000.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 0.01   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 30.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 0.0        /* charge of particles of first type */
#define CHARGE_B 0.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 50        /* time at which to add first particle */
#define ADD_PERIOD 100      /* time interval between adding further particles */
#define ADD_TYPE 4         /* type of added particles */
#define N_ADD_PARTICLES 2  /* number of particles to add */
#define FINAL_NOADD_PERIOD 1000  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 250 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 150.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 2              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 1      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 271            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 6           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 3.1         /* maximal distance for reaction to occur */
#define REACTION_PROB 0.1         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.5     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 160.0         /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 15       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 1              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 1      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60      /* size of hashgrid in x direction */
#define HASHY 30      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 21 September 2025 - A single chirp in a parabolic trap ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 192         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 2.3	    /* parameter controlling the dimensions of domain */
#define MU 1.15              /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define WALL_WIDTH_B 0.05   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 44  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.1       /* defines oscilling beam range */
#define OSCIL_YMID -0.75      /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.002         /* frequency of periodic excitation */
#define AMPLITUDE 2.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.005      /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number in medium B */
#define COURANTB 0.05      /* Courant number */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12.5  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2250         /* number of frames of movie */
#define NVID 14             /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 250    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.08       /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 8

#define PLOT_B 0        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 13     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 11   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 100.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 2.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 20 September 2025 - A fictitious, self-regulating nuclear reaction ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.8
#define INITXMAX 1.8	/* x interval for initial condition */
#define INITYMIN -1.01
#define INITYMAX 1.05	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.0
#define ADDXMAX -1.9	/* x interval for adding particles */
#define ADDYMIN -0.2
#define ADDYMAX 0.2	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.8 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.015          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.01         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.01       /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 30           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 1700      /* number of frames of movie */
#define NVID 15          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 5        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 0          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 1.0e8      /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 5000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 340.0       /* hue of particles of type 1 */
#define HUE_TYPE2 250.0      /* hue of particles of type 2 */
#define HUE_TYPE3 310.0      /* hue of particles of type 3 */
#define HUE_TYPE4 60.0       /* hue of particles of type 4 */
#define HUE_TYPE5 210.0       /* hue of particles of type 5 */
#define HUE_TYPE6 150.0      /* hue of particles of type 6 */
#define HUE_TYPE7 110.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 7.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 7.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 10.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 235.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 112.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 1.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 500.0        /* initial velocity range */
#define V_INITIAL_ADD 2000.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 0.01   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 30.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 0.0        /* charge of particles of first type */
#define CHARGE_B 0.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 50        /* time at which to add first particle */
#define ADD_PERIOD 200      /* time interval between adding further particles */
#define ADD_TYPE 4         /* type of added particles */
#define N_ADD_PARTICLES 2  /* number of particles to add */
#define FINAL_NOADD_PERIOD 1200  /* final period where no particles are added */
#define SAFETY_FACTOR 1000.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 250 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 150.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 2              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 1      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 272            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 9                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 7           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 3.1         /* maximal distance for reaction to occur */
#define REACTION_PROB 0.1         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 1.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 125.0         /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 15       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 1              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 1      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60      /* size of hashgrid in x direction */
#define HASHY 30      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 19 September 2025 - Classics revisited: A quadratic resonance diffuser ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:**

```
    init_wave_flat(phi, psi, xy_in);
    
    /* add oscillating waves */
    wave_source_x[0] = -1.0;
    wave_source_y[0] = 0.0;
    source_periods[0] = OSCILLATING_SOURCE_PERIOD;
    source_amp[0] = INITIAL_AMP;
    for (source = 0; source < N_SOURCES; source++)
    {
        dperiod = source_periods[source];
        phase = i - (int)(dperiod*(double)((int)((double)i/dperiod)));
        if ((ADD_OSCILLATING_SOURCE)&&(phase == 1)&&(i<MAX_PULSING_TIME))
        {
            printf("Source %i: Adding pulse %i\n", source, add_counter[source]);
            add_counter[source]++;
            if (ALTERNATE_OSCILLATING_SOURCE) sign[source] = -sign[source];
            add_circular_wave(-sign[source]*source_amp[source], wave_source_x[source], wave_source_y[source], phi, psi, xy_in);
        }
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 46         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 0.1	    /* parameter controlling the dimensions of domain */
#define MU 0.2              /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define WALL_WIDTH_B 0.05   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 44  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.1       /* defines oscilling beam range */
#define OSCIL_YMID -0.75      /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.002         /* frequency of periodic excitation */
#define AMPLITUDE 2.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.005      /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number in medium B */
#define COURANTB 0.05      /* Courant number */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12.5  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2100         /* number of frames of movie */
#define NVID 12             /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.5              /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.002       /* variance of initial condition */
#define INITIAL_WAVELENGTH 0.005     /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 11     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 60.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.7      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.2    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 18 September 2025 - A chain reaction in a denser ²³⁵U sample ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.8
#define INITXMAX 1.8	/* x interval for initial condition */
#define INITYMIN -1.01
#define INITYMAX 1.05	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.0
#define ADDXMAX -1.9	/* x interval for adding particles */
#define ADDYMIN -0.2
#define ADDYMAX 0.2	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.6 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.015          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.01         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.01       /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 30           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2200      /* number of frames of movie */
#define NVID 15          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 9        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 5000000.0      /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 340.0       /* hue of particles of type 1 */
#define HUE_TYPE2 250.0      /* hue of particles of type 2 */
#define HUE_TYPE3 310.0      /* hue of particles of type 3 */
#define HUE_TYPE4 60.0       /* hue of particles of type 4 */
#define HUE_TYPE5 210.0       /* hue of particles of type 5 */
#define HUE_TYPE6 150.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 7.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 7.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 10.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 235.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 112.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 1.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 500.0        /* initial velocity range */
#define V_INITIAL_ADD 2000.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 0.01   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 30.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 0.0        /* charge of particles of first type */
#define CHARGE_B 0.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 50        /* time at which to add first particle */
#define ADD_PERIOD 200      /* time interval between adding further particles */
#define ADD_TYPE 4         /* type of added particles */
#define N_ADD_PARTICLES 2  /* number of particles to add */
#define FINAL_NOADD_PERIOD 1200  /* final period where no particles are added */
#define SAFETY_FACTOR 1000.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 250 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 150.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 2              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 1      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 271            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 6           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 3.1         /* maximal distance for reaction to occur */
#define REACTION_PROB 0.1         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 1.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 200.0         /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 4.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 1              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 1      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60      /* size of hashgrid in x direction */
#define HASHY 30      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 17 September 2025 - Sending chirps into a parabolic trap ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 192         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 2.3	    /* parameter controlling the dimensions of domain */
#define MU 1.15              /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define WALL_WIDTH_B 0.05   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 44  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.1       /* defines oscilling beam range */
#define OSCIL_YMID -0.75      /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.015         /* frequency of periodic excitation */
#define AMPLITUDE 2.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number in medium B */
#define COURANTB 0.05      /* Courant number */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12.5  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 4225         /* number of frames of movie */
#define NVID 14             /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 250    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.08       /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 8

#define PLOT_B 0        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 12     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 18   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 25.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 0.5      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 16 September 2025 - Adding first 33% light cations and then 67% heavy anions in a central gravity field ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -0.1
#define INITXMAX 0.1	/* x interval for initial condition */
#define INITYMIN -0.1
#define INITYMAX 0.1	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -1.9
#define ADDXMAX 1.9	/* x interval for adding particles */
#define ADDYMIN 1.2
#define ADDYMAX 1.3	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.2
#define BCXMAX 2.2	/* x interval for boundary condition */
#define BCYMIN -2.2
#define BCYMAX 2.2	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 0  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.1 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.013 	    /* parameter controlling radius of particles */
#define MU_B 0.013          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.012         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.016         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 1           /* number of grid point for grid of disks */
#define NGRIDY 1           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 5500      /* number of frames of movie */
#define NVID 100         /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 9        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 50000.0      /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 2.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 250.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 1.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 2.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 0.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 200.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 10000.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 1     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_REGION 1      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 0        /* time at which to add first particle */
#define ADD_PERIOD 6      /* time interval between adding further particles */
#define ADD_TYPE 4         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 800  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 1     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.666666667    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 25.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 1      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 23            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 2                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 6           /* number of types shown in graph */
#define RD_INITIAL_COND 51        /* initial condition of particles */
#define REACTION_DIST 2.5         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e10      /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 1              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 50      /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 15 September 2025 - A linear wave crossing a pentagonal prism with refractive index 2 ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -1.5
#define XMAX 2.5	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 8          /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 0.5	    /* parameter controlling the dimensions of domain */
#define MU 1.15             /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define WALL_WIDTH_B 0.05   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 1         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.1       /* defines oscilling beam range */
#define OSCIL_YMID -0.75      /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.015         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.05       /* Courant number in medium B */
#define COURANTB 0.1       /* Courant number */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12.5  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2600         /* number of frames of movie */
#define NVID 14             /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 150    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.08       /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 11     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 14   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 50.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 0.5      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 14 September 2025 - Adding first 25% cations and then 75% anions in a central gravity field ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -0.1
#define INITXMAX 0.1	/* x interval for initial condition */
#define INITYMIN -0.1
#define INITYMAX 0.1	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -1.9
#define ADDXMAX 1.9	/* x interval for adding particles */
#define ADDYMIN 1.2
#define ADDYMAX 1.3	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.2
#define BCXMAX 2.2	/* x interval for boundary condition */
#define BCYMIN -2.2
#define BCYMAX 2.2	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 0  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.1 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.014 	    /* parameter controlling radius of particles */
#define MU_B 0.014          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.014         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.014         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 1           /* number of grid point for grid of disks */
#define NGRIDY 1           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 4600      /* number of frames of movie */
#define NVID 100          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 9        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 50000.0      /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 2.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 250.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 1.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 0.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 200.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 10000.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 1     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_REGION 1      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 0        /* time at which to add first particle */
#define ADD_PERIOD 6      /* time interval between adding further particles */
#define ADD_TYPE 4         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 800  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 1     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.75    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 25.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 1      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 23            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 2                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 6           /* number of types shown in graph */
#define RD_INITIAL_COND 51        /* initial condition of particles */
#define REACTION_DIST 1.85        /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e10      /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 1              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 50      /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 13 September 2025 - A polychromatic beam of light hitting a parabolic trap ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 192         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 2.3	    /* parameter controlling the dimensions of domain */
#define MU 1.15              /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define WALL_WIDTH_B 0.05   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 51  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.1       /* defines oscilling beam range */
#define OSCIL_YMID -0.75      /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.015         /* frequency of periodic excitation */
#define AMPLITUDE 2.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number in medium B */
#define COURANTB 0.05      /* Courant number */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12.5  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2600         /* number of frames of movie */
#define NVID 14             /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 250    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.08       /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 8

#define PLOT_B 0        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 13     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 15   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 25.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 0.5      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 12 September 2025 - Adding first cations and then anions in a central gravity field ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -0.1
#define INITXMAX 0.1	/* x interval for initial condition */
#define INITYMIN -0.1
#define INITYMAX 0.1	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -1.9
#define ADDXMAX 1.9	/* x interval for adding particles */
#define ADDYMIN 1.2
#define ADDYMAX 1.3	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.2
#define BCXMAX 2.2	/* x interval for boundary condition */
#define BCYMIN -2.2
#define BCYMAX 2.2	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 0  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.1 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.014 	    /* parameter controlling radius of particles */
#define MU_B 0.014          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.014         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.014         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 1           /* number of grid point for grid of disks */
#define NGRIDY 1           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 4000      /* number of frames of movie */
#define NVID 100          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 9        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 50000.0      /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 2.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 250.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 1.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 0.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 200.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 10000.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 1     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_REGION 1      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 0        /* time at which to add first particle */
#define ADD_PERIOD 6      /* time interval between adding further particles */
#define ADD_TYPE 4         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 1000  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 1     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 25.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 1      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 23            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 2                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 6           /* number of types shown in graph */
#define RD_INITIAL_COND 51        /* initial condition of particles */
#define REACTION_DIST 1.85        /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e10      /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 1              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 50      /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 11 September 2025 - Using a parabolic trap in reverse ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:**

```
    init_wave_flat(phi, psi, xy_in);
    
    /* add oscillating waves */
    wave_source_x[0] = 0.0;
    wave_source_y[0] = 0.0;
    source_periods[0] = OSCILLATING_SOURCE_PERIOD;
    source_amp[0] = INITIAL_AMP;
    for (source = 0; source < N_SOURCES; source++)
    {
        dperiod = source_periods[source];
        phase = i - (int)(dperiod*(double)((int)((double)i/dperiod)));
        if ((ADD_OSCILLATING_SOURCE)&&(phase == 1)&&(i<MAX_PULSING_TIME))
        {
            printf("Source %i: Adding pulse %i\n", source, add_counter[source]);
            add_counter[source]++;
            if (ALTERNATE_OSCILLATING_SOURCE) sign[source] = -sign[source];
            add_circular_wave(-sign[source]*source_amp[source], wave_source_x[source], wave_source_y[source], phi, psi, xy_in);
        }
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -1.5
#define XMAX 2.5	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 191         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 1.15	    /* parameter controlling the dimensions of domain */
#define MU 2.3              /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define WALL_WIDTH_B 0.05   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 42  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.1       /* defines oscilling beam range */
#define OSCIL_YMID -0.75      /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.015         /* frequency of periodic excitation */
#define AMPLITUDE 2.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number in medium B */
#define COURANTB 0.05      /* Courant number */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12.5  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 1       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 1600         /* number of frames of movie */
#define NVID 10             /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.5             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.001      /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 17     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 40.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 10 September 2025 - Poisoning of a nuclear reactor ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.8
#define INITXMAX 1.8	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.0
#define ADDXMAX -1.9	/* x interval for adding particles */
#define ADDYMIN -0.2
#define ADDYMAX 0.2	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 0  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.53 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.015          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.01         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.01       /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 4400      /* number of frames of movie */
#define NVID 25          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 9        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 400000.0      /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 340.0       /* hue of particles of type 1 */
#define HUE_TYPE2 250.0      /* hue of particles of type 2 */
#define HUE_TYPE3 310.0      /* hue of particles of type 3 */
#define HUE_TYPE4 60.0       /* hue of particles of type 4 */
#define HUE_TYPE5 210.0       /* hue of particles of type 5 */
#define HUE_TYPE6 150.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 10.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 235.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 112.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 1.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 500.0        /* initial velocity range */
#define V_INITIAL_ADD 2000.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 0.01   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 30.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 0.0        /* charge of particles of first type */
#define CHARGE_B 0.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 50        /* time at which to add first particle */
#define ADD_PERIOD 5000      /* time interval between adding further particles */
#define ADD_TYPE 4         /* type of added particles */
#define N_ADD_PARTICLES 2  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 1000.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 250 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 150.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 2              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 1      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 271            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 6           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 3.3         /* maximal distance for reaction to occur */
#define REACTION_PROB 0.1         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 1.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 200.0         /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 4.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 1              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 1      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60      /* size of hashgrid in x direction */
#define HASHY 30      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 09 September 2025 - Bloopers 26: A not very effective invisibility cloak ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -0.75
#define XMAX 3.25	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 20         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 6   /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 6        /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 0.12	    /* parameter controlling the dimensions of domain */
#define MU 0.1	            /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define WALL_WIDTH_B 0.05   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 1         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.1       /* defines oscilling beam range */
#define OSCIL_YMID -0.75      /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.012        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number in medium B */
#define COURANTB 0.05      /* Courant number */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12.5  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 1       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2200         /* number of frames of movie */
#define NVID 10             /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 74     /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.5             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.001      /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 11     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 75.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.5      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 08 September 2025 - Controlling a chain reaction with cadmium atoms ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.8
#define INITXMAX 1.8	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.0
#define ADDXMAX -1.9	/* x interval for adding particles */
#define ADDYMIN -0.1
#define ADDYMAX 0.1	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 0  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.6 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.015          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.01         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.01       /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2000      /* number of frames of movie */
#define NVID 25          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 9        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 400000.0      /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 340.0       /* hue of particles of type 1 */
#define HUE_TYPE2 250.0      /* hue of particles of type 2 */
#define HUE_TYPE3 310.0      /* hue of particles of type 3 */
#define HUE_TYPE4 60.0       /* hue of particles of type 4 */
#define HUE_TYPE5 210.0       /* hue of particles of type 5 */
#define HUE_TYPE6 150.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 10.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 235.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 112.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 1.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 500.0        /* initial velocity range */
#define V_INITIAL_ADD 2000.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 0.01   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 30.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 0.0        /* charge of particles of first type */
#define CHARGE_B 0.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 50        /* time at which to add first particle */
#define ADD_PERIOD 5000      /* time interval between adding further particles */
#define ADD_TYPE 4         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 1000.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 250 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 150.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 2              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 1      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 271            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 6           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 3.3         /* maximal distance for reaction to occur */
#define REACTION_PROB 0.1         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 1.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 50.0         /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 4.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 1              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 1      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60      /* size of hashgrid in x direction */
#define HASHY 30      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 07 September 2025 - Trapping beams between parabolic reflectors is harder for longer wavelength ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 192         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 2.3	    /* parameter controlling the dimensions of domain */
#define MU 1.15              /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define WALL_WIDTH_B 0.05   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 42  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.1       /* defines oscilling beam range */
#define OSCIL_YMID -0.75      /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.015         /* frequency of periodic excitation */
#define AMPLITUDE 2.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number in medium B */
#define COURANTB 0.05      /* Courant number */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12.5  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2500         /* number of frames of movie */
#define NVID 14             /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 250    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.08       /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 8

#define PLOT_B 0        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 13     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 15   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 50.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 06 September 2025 - Chain reaction: Shooting a single neutron into a chunk of 235U ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.8
#define INITXMAX 1.8	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.0
#define ADDXMAX -1.9	/* x interval for adding particles */
#define ADDYMIN -0.1
#define ADDYMAX 0.1	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 0  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.02           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.01         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.01       /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2000      /* number of frames of movie */
#define NVID 25          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 9        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 400000.0      /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 340.0      /* hue of particles of type 2 */
#define HUE_TYPE3 240.0      /* hue of particles of type 3 */
#define HUE_TYPE4 150.0      /* hue of particles of type 4 */
#define HUE_TYPE5 230.0       /* hue of particles of type 5 */
#define HUE_TYPE6 60.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 10.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 235.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 1.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 500.0        /* initial velocity range */
#define V_INITIAL_ADD 2000.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 0.01   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 30.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 0.0        /* charge of particles of first type */
#define CHARGE_B 0.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 50        /* time at which to add first particle */
#define ADD_PERIOD 5000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 1000.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 250 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 150.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 1              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 1      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 27            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 6                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 4           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 3.3         /* maximal distance for reaction to occur */
#define REACTION_PROB 0.1         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 1.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 50.0         /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 4.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 1              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 1      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60      /* size of hashgrid in x direction */
#define HASHY 30      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 05 September 2025 - Classics revisited: Two facing Fresnel lenses ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 45         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 1.35	    /* parameter controlling the dimensions of domain */
#define MU 0.07  	    /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.013    /* width of wall separating lenses */
#define WALL_WIDTH_B 0.01   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 1         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define AMPLITUDE 1.25     /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number in medium B */
#define COURANTB 0.066667  /* Courant number */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12.5  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3

/* Parameters for length and speed of simulation */

#define NSTEPS 1800         /* number of frames of movie */
#define NVID 10              /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 100    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 11     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 100.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.001  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.5      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.12   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 04 September 2025 - Chain reaction: A simple model for nuclear fission ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.8
#define INITXMAX 1.8	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -1.9
#define ADDXMAX 1.9	/* x interval for adding particles */
#define ADDYMIN 1.2
#define ADDYMAX 1.3	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 0  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.02 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.02           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.01         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.01       /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2500      /* number of frames of movie */
#define NVID 25          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 9        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 400000.0      /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 340.0      /* hue of particles of type 2 */
#define HUE_TYPE3 240.0      /* hue of particles of type 3 */
#define HUE_TYPE4 150.0      /* hue of particles of type 4 */
#define HUE_TYPE5 230.0       /* hue of particles of type 5 */
#define HUE_TYPE6 60.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 10.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 235.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 1.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 500.0        /* initial velocity range */
#define V_INITIAL_ADD 0.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 0.01   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 30.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 0.0        /* charge of particles of first type */
#define CHARGE_B 0.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 0        /* time at which to add first particle */
#define ADD_PERIOD 10      /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 2000  /* final period where no particles are added */
#define SAFETY_FACTOR 1000.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 1   /* set to 1 to randomly select sign of added charge */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 250 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 150.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 1              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 1      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 27            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 6                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 4           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 3.3         /* maximal distance for reaction to occur */
#define REACTION_PROB 0.1         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 1.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 50.0         /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 1              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 1      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60      /* size of hashgrid in x direction */
#define HASHY 30      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 03 September 2025 - Trapping a pulse of light between two parabolic reflectors ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 192         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 2.3	    /* parameter controlling the dimensions of domain */
#define MU 1.15              /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define WALL_WIDTH_B 0.05   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 43  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.1       /* defines oscilling beam range */
#define OSCIL_YMID -0.75      /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.03         /* frequency of periodic excitation */
#define AMPLITUDE 2.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number in medium B */
#define COURANTB 0.05      /* Courant number */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12.5  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2950         /* number of frames of movie */
#define NVID 14             /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 250    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.08       /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 8

#define PLOT_B 0        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 13     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 11   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 50.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 02 September 2025 - “Leftons” and “rightons” under weak gravity ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.8
#define INITXMAX 1.8	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -1.9
#define ADDXMAX 1.9	/* x interval for adding particles */
#define ADDYMIN 1.2
#define ADDYMAX 1.3	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.325	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 0  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 9      /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.2  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 6      /* max number of neighbours to be pinned */
#define OBSTACLE_PINNING_TYPE 0     /* type of obstacle pinning, see OP_ in global_ljones */
#define BDRY_PINNING_STEP 4         /* interval between pinned obstacles on boundary */
#define RECOUPLE_OBSTACLES 0        /* set to 1 to reset obstacle coupling */
#define OBSTACLE_RECOUPLE_TYPE 1    /* algorithm for recoupling, see OR_ in global_ljones */
#define OBSTACLE_RECOUPLE_TIME 200    /* time between obstacle recouplings */
#define UNCOUPLE_MAXLENGTH 2.0      /* length at which bonds decouple */
#define COUPLE_MINLENGTH 0.5        /* length at which bonds decouple */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 15     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.5 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 14       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.015 	    /* parameter controlling radius of particles */
#define MU_B 0.015          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.015         /* parameter controlling radius of added particles */
#define MU_ADD_B 0.015         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 5500      /* number of frames of movie */
#define NVID 100          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 0     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 9        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 0    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT 1.0     /* shift in color hue (for some cyclic palettes) */

#define PRINT_PARAMETERS 1  /* set to 1 to print certain parameters */
#define PRINT_TEMPERATURE 0 /* set to 1 to print current temperature */
#define PRINT_ANGLE 0               /* set to 1 to print obstacle orientation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */
#define PRINT_SEGMENTS_FORCE 0      /* set to 1 to print force on segments */
#define PRINT_NPARTICLES 0          /* print number of active particles */
#define PRINT_TYPE_PROP 0           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 50000.0      /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 1.0e-6   /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 1000.0         /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 1.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 0.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 0.01   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 30.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 10000.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 150000.0      /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 180000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 1.0     /* charge of obstacles */
#define OBSTACLE_MASS 1000.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 1.0e10  /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 1.0e8   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 0         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent electric field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define DRAW_E_ARROW 0           /* set to 1 to draw E field arrow */
#define E_ARROW_YSHIFT 0.05      /* vertical position of E field arrow */
#define PRINT_CURRENT 0          /* set to 1 to print electric current (x component) */
#define DRAW_CURRENT_ARROW 0     /* set to 1 to draw current arrow */
#define MAX_CURRENT 200.0       /* current scale */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 10.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.06    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SMOOTH_CONTAINER_DECREASE 1 /* set to 1 to decrease size smoothly at each simulation step */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 800            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.02  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 9     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 0        /* time at which to add first particle */
#define ADD_PERIOD 10      /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 2000  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 1   /* set to 1 to randomly select sign of added charge */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 1000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 7000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 250 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 150.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 2.0  /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 2        /* width of tracer particle trajectory */

#define TRACK_PARTICLE 0          /* set to 1 to track a given particle */
#define TRACKED_PARTICLE 2        /* number of tracked particle */
#define TRACK_INITIAL_TIME 900    /* time when starting to track */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 150       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 0    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */
#define INACTIVATE_SEGMENTS_UNDER_PRESSURE 0    /* set to 1 to inactivate segment groups when limit pressure is reached */
#define SEGMENT_P_INACTIVATE 6.0e7  /* pressure at which to inactivate group */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0    /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
#define KSPRING_BELT 1.0e4          /* spring constant from belt */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X -0.625         /* threshold value for position-dependent type */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing specaial initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 23            /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 2                /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 51        /* initial condition of particles */
#define REACTION_DIST 2.2         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e10      /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 1              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 0     /* set to 1 for reacting particles to become neutral */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 2      /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 2 /* size limitation on smaller cluster */
#define SMALL_NP_MAXSIZE 2      /* limitation on number of partners of particle in smaller cluster */
#define NOTSELECTED_CLUSTER_MAXSIZE 0   /* limit on size of clusters that can merge with non-selected cluster */
#define REPAIR_CLUSTERS 0       /* set to 1 to repair alignment in clusters */
#define REPAIR_MIN_DIST 0.75    /* relative distance below which overlapping polygons are inactivated */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 1.0 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */
#define PLOT_CURRENTS 0     /* set to 1 to make current vs E field plot */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.025    /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 1         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99     /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 5         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 81     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.035            /* radius of partner particle */
#define PARTICLE_MASS_D 2.0  /* mass or partner particle */
#define CHARGE_D 1.0         /* charge of partner particle */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 40      /* size of hashgrid in x direction */
#define HASHY 26      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 01 September 2025 - Classics revisited: A linear wave crossing a Fresnel lens ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -1.0
#define XMAX 3.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 0.8 /* scaling for Julia sets and some other domains */

/* Choice of the billiard table */

#define B_DOMAIN 43         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 2.3    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA -1.1	    /* parameter controlling the dimensions of domain */
#define MU 0.1   	    /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 14           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.013    /* width of wall separating lenses */
#define WALL_WIDTH_B 0.01   /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9  
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.2
#define ISO_YSHIFT_RIGHT 0.15 
#define ISO_SCALE 0.475           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */
/* Physical parameters of wave equation */

#define TWOSPEEDS 1         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number in medium B */
#define COURANTB 0.066667  /* Courant number */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12.5  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3

/* Parameters for length and speed of simulation */

#define NSTEPS 1500         /* number of frames of movie */
#define NVID 10              /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
#define INITIAL_TIME 200    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 300    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 18     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 100.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.75   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 0       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.01  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.5      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.12   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 18      /* width of maze */
#define NYMAZE 9       /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```


