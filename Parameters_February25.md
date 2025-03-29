Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Below are parameter values used for different simulations, as well as initial conditions used in function animation. Some simulations use variants of the published code. The list is going to be updated gradually. 


### 28 February 2025 - Deep impact 9: Throwing 5 times heavier particles into a sea of light particles ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
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

#define INITXMIN -2.0
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX -0.25	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.05
#define ADDXMAX -2.0	  /* x interval for adding particles */
#define ADDYMIN 0.5
#define ADDYMAX 1.0       /* y interval for adding particles */
#define ADDRMIN 4.75 
#define ADDRMAX 6.0     /* r interval for adding particles */

#define BCXMIN -2.1
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -1.8
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 91     /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.12  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 3      /* max number of neighbours to be pinned */
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
#define TYPE_PROPORTION 0.5 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 4.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.015 	    /* parameter controlling radius of particles */
#define MU_B 0.03           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.02         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.075           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 80           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 40
#define NOBSY 24            /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 3600      /* number of frames of movie */
#define NVID 100          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 2000     /* time after which to start saving frames */
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

#define PLOT 13
#define PLOT_B 12        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 7          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 5        /* type of background coloring, see list in global_ljones.c */
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
#define OBSTACLE_AREA_SHADE_FACTOR 80.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT -0.5     /* shift in color hue (for some cyclic palettes) */

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
#define PARTICLE_EMIN 100.0         /* energy of particle with coolest color */
#define PARTICLE_EMAX 500.0         /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 150.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 1   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.5   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.5    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 2.5    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 2.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 400.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 10.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.5     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.5     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 250.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e9    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 500.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 10.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 250    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 500    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0     /* set to 1 to add an electric field */
#define EFIELD 1000.0       /* value of electric field */
#define EFIELD_Y 0.0     /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 0.0       /* value of magnetic field */
#define CHARGE 0.0       /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 5000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 0      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define OBSTACLE_MASS 100.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 5.0e5   /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 2.0e5   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 1         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 0          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 0          /* set to 1 to draw cross on particles of negative charge */
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
#define OBSTACLE_RADIUS 0.015  /* radius of obstacle for circle boundary conditions */
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
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 250     /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 500  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 100    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 40000    /* length of recorded trajectory */
#define TRACER_LUM_FACTOR 100.0     /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 4.0  /* relative mass of tracer particle */
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

#define REACTION_DIFFUSION 0      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 262           /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 1                /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 0         /* initial condition of particles */
#define REACTION_DIST 1.2         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0       /* change of kinetic energy in reaction */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 6              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 1000     /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 10 /* size limitation on smaller cluster */
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
#define DEACIVATE_CLOSE_PAIRS 1 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e10    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 4         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 9      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.1      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 1  /* set to 1 to pair particle of type 1 */
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

#define HASHX 80      /* size of hashgrid in x direction */
#define HASHY 40      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 27 February 2025 - A lighthouse beam in an eight-pointed star ###

**Program:** `particle_trajectory.c` 

**Initial condition in function `animation()`:**

```
    x = LAMBDA;
    y = 0.0;
```

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define XMIN -1.2
#define XMAX 2.8	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1   /* pattern of circles */
#define POLYLINE_PATTERN 17  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */
#define NABSCIRCLES 16       /* number of absorbing circles */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100004        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDA 1.0 	/* parameter controlling shape of domain */
#define MU 0.005        /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 8             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define LAMBDA_B 0.5411961  /* parameter controlling shape of domain (for P_POLYRING) */
#define NPOLY_B 100000           /* number of sides of second polygon */
#define APOLY_B 1.0         /* angle by which to turn second polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define TRAJ_LENGTH 200 /* length of trajectory */
#define PLOT_NMAX 150    /* max length on length plot */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define HEATMAP 0       /* set to 1 to show heat map of particles */
#define PLOT_HEATMAP_AVERAGE 0      /* set to 1 to plot average number of particles in heat map */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 0   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 1   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_ANGLE 1    /* set to 1 to print angle of trajectory 0 */
#define PRINT_TRAJECTORY_SLOPE 0    /* set to 1 to print slope of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 1         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 0         /* set to 1 to use log scale for plot of lengths */
#define LENGTH_PLOT_POLAR 0         /* set to 1 to plot lengths in polar coordinates */
#define MIN_ANGLE 157.5          /* range of angles of trajectory */
#define MAX_ANGLE 202.5         /* range of angles of trajectory */
#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define SLOW_AT_LONG_TRAJ 0     /* set to 1 to slow down movie for long trajectories */
#define ADD_SUCCESS_GALLERY 0   /* set to 1 to add gallery of successful trajectories at end of movie */
#define SUCCESS_GALLERY_FRAMES 25   /* number of frames per success */
#define EXIT_BOTH_WAYS 1        /* set to 1 to add exits to he left to succesful trajectories */

#define NSTEPS 5500       /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 10    /* Color palette, see list in global_pdes.c  */

#define NCOLORS 64       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define COLOR_HUEMIN 60  /* minimal color hue */
#define COLOR_HUEMAX 360 /* maximal color hue */
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define GRAPH_HUE 90.0      /* color hue of graph */
#define SUCCESS_HUE 180.0   /* color hue of success */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 250       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define END_FRAMES  250   /* number of frames at end of movie */

#define NXMAZE 25      /* width of maze */
#define NYMAZE 25      /* height of maze */
#define MAZE_MAX_NGBH 8     /* max number of neighbours of maze cell */
#define RAND_SHIFT 10       /* seed of random number generator */
#define MAZE_XSHIFT -0.7     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.4     /* radius of tounded corners in maze */
#define CLOSE_MAZE 0        /* set to 1 to close maze exits */

```

### 26 February 2025 - Waves in a lattice of circles connected by channels with random width ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:**

```
    /* add oscillating waves */
    for (k=0; k<N_SOURCES; k++)
    {
        wave_source_x[k] = 0.0;
        wave_source_y[k] = 0.0;
    }
    for (source = 0; source < N_SOURCES; source++)
    {
        source_shift[source] = 0;
        if ((ADD_OSCILLATING_SOURCE)&&(i%(OSCILLATING_SOURCE_PERIOD) == source_shift[source])&&(i<MAX_PULSING_TIME))
        {
            if (ALTERNATE_OSCILLATING_SOURCE) sign[source] = -sign[source];
            add_circular_wave(-sign[source]*INITIAL_AMP, wave_source_x[source], wave_source_y[source], phi, psi, xy_in);
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

#define B_DOMAIN 931        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.5    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.25	    /* parameter controlling the dimensions of domain */
#define MU 0.075            /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 9            /* number of grid point for grid of disks */
#define WALL_WIDTH 0.02     /* width of wall separating lenses */
#define WALL_WIDTH_RND 0.5  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */

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
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 61  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.025         /* frequency of periodic excitation */
#define AMPLITUDE 0.5      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.025       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 40.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 13    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 1500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 2750         /* number of frames of movie */
#define NVID 7              /* number of iterations between images displayed on screen */
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

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
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
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.25   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 0.2    /* additional scaling factor for wave amplitude */
#define E_SCALE 15.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.85   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */

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

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
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

### 25 February 2025 - A lighthouse beam in a six-pointed star ###

**Program:** `particle_trajectory.c` 

**Initial condition in function `animation()`:**

```
    x = LAMBDA;
    y = 0.0;
```

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define XMIN -1.2
#define XMAX 2.8	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1   /* pattern of circles */
#define POLYLINE_PATTERN 17  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */
#define NABSCIRCLES 12       /* number of absorbing circles */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100004        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDA 1.0 	/* parameter controlling shape of domain */
#define MU 0.005        /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define LAMBDA_B 0.5773503  /* parameter controlling shape of domain (for P_POLYRING) */
#define NPOLY_B 100000           /* number of sides of second polygon */
#define APOLY_B 1.0         /* angle by which to turn second polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define TRAJ_LENGTH 200 /* length of trajectory */
#define PLOT_NMAX 150    /* max length on length plot */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define HEATMAP 0       /* set to 1 to show heat map of particles */
#define PLOT_HEATMAP_AVERAGE 0      /* set to 1 to plot average number of particles in heat map */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 0   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 1   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_ANGLE 1    /* set to 1 to print angle of trajectory 0 */
#define PRINT_TRAJECTORY_SLOPE 0    /* set to 1 to print slope of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 1         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 0         /* set to 1 to use log scale for plot of lengths */
#define LENGTH_PLOT_POLAR 0         /* set to 1 to plot lengths in polar coordinates */
#define MIN_ANGLE 150.0          /* range of angles of trajectory */
#define MAX_ANGLE 210.0         /* range of angles of trajectory */
#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define SLOW_AT_LONG_TRAJ 0     /* set to 1 to slow down movie for long trajectories */
#define ADD_SUCCESS_GALLERY 0   /* set to 1 to add gallery of successful trajectories at end of movie */
#define SUCCESS_GALLERY_FRAMES 25   /* number of frames per success */
#define EXIT_BOTH_WAYS 1        /* set to 1 to add exits to he left to succesful trajectories */

#define NSTEPS 4400       /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 16    /* Color palette, see list in global_pdes.c  */

#define NCOLORS 64       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define COLOR_HUEMIN 60  /* minimal color hue */
#define COLOR_HUEMAX 360 /* maximal color hue */
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define GRAPH_HUE 90.0      /* color hue of graph */
#define SUCCESS_HUE 180.0   /* color hue of success */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 250       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define END_FRAMES  250   /* number of frames at end of movie */

#define NXMAZE 25      /* width of maze */
#define NYMAZE 25      /* height of maze */
#define MAZE_MAX_NGBH 8     /* max number of neighbours of maze cell */
#define RAND_SHIFT 10       /* seed of random number generator */
#define MAZE_XSHIFT -0.7     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.4     /* radius of tounded corners in maze */
#define CLOSE_MAZE 0        /* set to 1 to close maze exits */

```

### 24 February 2025 - Electric bowling with mass ratio 1000 ###

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

#define INITXMIN -1.6
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.2
#define ADDXMAX -2.0	  /* x interval for adding particles */
#define ADDYMIN -1.0
#define ADDYMAX 1.0       /* y interval for adding particles */
#define ADDRMIN 4.75 
#define ADDRMAX 6.0     /* r interval for adding particles */

#define BCXMIN -2.5
#define BCXMAX 2.5	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -1.8
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 91     /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.12  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 3      /* max number of neighbours to be pinned */
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
#define TYPE_PROPORTION 0.5 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 4.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.03           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.007         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.075           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 20           /* number of grid point for grid of disks */
#define NGRIDY 10          /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 40
#define NOBSY 24            /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 3800      /* number of frames of movie */
#define NVID 250          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 30     /* time after which to start saving frames */
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

#define BOUNDARY_COND 23

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 12        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 7          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 5        /* type of background coloring, see list in global_ljones.c */
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
#define OBSTACLE_AREA_SHADE_FACTOR 80.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT -0.5     /* shift in color hue (for some cyclic palettes) */

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
#define PARTICLE_EMIN 100.0         /* energy of particle with coolest color */
#define PARTICLE_EMAX 3000.0       /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 150.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 1   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.5   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.5    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 2.5    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 2.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define PARTICLE_MASS 2000.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.5     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.5     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 100.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e9    /* harmonic potential of obstacles */
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

#define ADD_EFIELD 1     /* set to 1 to add an electric field */
#define EFIELD 1000.0       /* value of electric field */
#define EFIELD_Y 0.0     /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 0.0       /* value of magnetic field */
#define CHARGE 0.0       /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 5000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 0      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define OBSTACLE_MASS 100.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 5.0e5   /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 2.0e5   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 1         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 0          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 0          /* set to 1 to draw cross on particles of negative charge */
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
#define OBSTACLE_RADIUS 0.015  /* radius of obstacle for circle boundary conditions */
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
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 40     /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 500  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 100    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 40000    /* length of recorded trajectory */
#define TRACER_LUM_FACTOR 100.0     /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 4.0  /* relative mass of tracer particle */
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

#define REACTION_DIFFUSION 0      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 262           /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 1                /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 0         /* initial condition of particles */
#define REACTION_DIST 1.2         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0       /* change of kinetic energy in reaction */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 6              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 1000     /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 10 /* size limitation on smaller cluster */
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
#define DEACIVATE_CLOSE_PAIRS 1 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e10    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 4         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 9      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.1      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 1  /* set to 1 to pair particle of type 1 */
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
#define HASHY 20      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 23 February 2025 - Possibly a metamaterial: Waves in a 12 by 8 hex lattice of connected discs ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:**

```
    init_wave_flat(phi, psi, xy_in);
    
    /* add oscillating waves */
    for (k=0; k<N_SOURCES; k++)
    {
        wave_source_x[k] = 0.0;
        wave_source_y[k] = 0.0;
    }
    for (source = 0; source < N_SOURCES; source++)
    {
        source_shift[source] = 0;
        if ((ADD_OSCILLATING_SOURCE)&&(i%(OSCILLATING_SOURCE_PERIOD) == source_shift[source])&&(i<MAX_PULSING_TIME))
        {
            if (ALTERNATE_OSCILLATING_SOURCE) sign[source] = -sign[source];
            add_circular_wave(-sign[source]*INITIAL_AMP, wave_source_x[source], wave_source_y[source], phi, psi, xy_in);
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

#define B_DOMAIN 94         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.5    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.25	    /* parameter controlling the dimensions of domain */
#define MU 0.075            /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 12           /* number of grid point for grid of disks */
#define NGRIDY 8            /* number of grid point for grid of disks */
#define WALL_WIDTH 0.022    /* width of wall separating lenses */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */

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
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 61  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.025         /* frequency of periodic excitation */
#define AMPLITUDE 0.5      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.025       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 40.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 13    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 1500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 2850         /* number of frames of movie */
#define NVID 7              /* number of iterations between images displayed on screen */
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

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

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
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.25   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 0.2    /* additional scaling factor for wave amplitude */
#define E_SCALE 15.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.85   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */

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

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
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

### 22 February 2025 - A lighthouse beam in a Sinai billiard ###

**Program:** `particle_trajectory.c` 

**Initial condition in function `animation()`:**

```
    x = 0.999*LAMBDA*cos(0.5*PID);
    y = -0.999*LAMBDA*sin(0.5*PID);
```

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define XMIN -1.2
#define XMAX 2.8	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1   /* pattern of circles */
#define POLYLINE_PATTERN 2  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */
#define NABSCIRCLES 4       /* number of absorbing circles */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100004        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDA 1.4 	/* parameter controlling shape of domain */
#define MU 0.02        /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 4             /* number of sides of polygon */
#define APOLY 0.5           /* angle by which to turn polygon, in units of Pi/2 */ 
#define LAMBDA_B 0.577350269    /* parameter controlling shape of domain (for P_POLYRING) */
#define NPOLY_B 100000           /* number of sides of second polygon */
#define APOLY_B 1.0         /* angle by which to turn second polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define TRAJ_LENGTH 500 /* length of trajectory */
#define PLOT_NMAX 450    /* max length on length plot */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define HEATMAP 0       /* set to 1 to show heat map of particles */
#define PLOT_HEATMAP_AVERAGE 0      /* set to 1 to plot average number of particles in heat map */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 0   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 1   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_ANGLE 1    /* set to 1 to print angle of trajectory 0 */
#define PRINT_TRAJECTORY_SLOPE 0    /* set to 1 to print slope of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 1         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 0         /* set to 1 to use log scale for plot of lengths */
#define LENGTH_PLOT_POLAR 0         /* set to 1 to plot lengths in polar coordinates */
#define MIN_ANGLE 90.0          /* range of angles of trajectory */
#define MAX_ANGLE 135.0         /* range of angles of trajectory */
#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define SLOW_AT_LONG_TRAJ 0     /* set to 1 to slow down movie for long trajectories */
#define ADD_SUCCESS_GALLERY 0   /* set to 1 to add gallery of successful trajectories at end of movie */
#define SUCCESS_GALLERY_FRAMES 25   /* number of frames per success */
#define EXIT_BOTH_WAYS 1        /* set to 1 to add exits to he left to succesful trajectories */

#define NSTEPS 6000       /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 15    /* Color palette, see list in global_pdes.c  */

#define NCOLORS 64       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define COLOR_HUEMIN 60  /* minimal color hue */
#define COLOR_HUEMAX 360 /* maximal color hue */
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define GRAPH_HUE 90.0      /* color hue of graph */
#define SUCCESS_HUE 180.0   /* color hue of success */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 250       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define END_FRAMES  250   /* number of frames at end of movie */

#define NXMAZE 25      /* width of maze */
#define NYMAZE 25      /* height of maze */
#define MAZE_MAX_NGBH 8     /* max number of neighbours of maze cell */
#define RAND_SHIFT 10       /* seed of random number generator */
#define MAZE_XSHIFT -0.7     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.4     /* radius of tounded corners in maze */
#define CLOSE_MAZE 0        /* set to 1 to close maze exits */

```

### 21 February 2025 - Electric bowling with mass ratio 100 ###

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

#define INITXMIN -1.6
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.2
#define ADDXMAX -2.0	  /* x interval for adding particles */
#define ADDYMIN -1.0
#define ADDYMAX 1.0       /* y interval for adding particles */
#define ADDRMIN 4.75 
#define ADDRMAX 6.0     /* r interval for adding particles */

#define BCXMIN -2.5
#define BCXMAX 2.5	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -1.8
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 91     /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.12  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 3      /* max number of neighbours to be pinned */
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
#define TYPE_PROPORTION 0.5 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 4.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.03           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.007         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.075           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 20           /* number of grid point for grid of disks */
#define NGRIDY 10          /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 40
#define NOBSY 24            /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 4400      /* number of frames of movie */
#define NVID 200          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 30     /* time after which to start saving frames */
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

#define BOUNDARY_COND 23

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 12        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 7          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 5        /* type of background coloring, see list in global_ljones.c */
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
#define OBSTACLE_AREA_SHADE_FACTOR 80.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT -0.5     /* shift in color hue (for some cyclic palettes) */

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
#define PARTICLE_EMIN 100.0         /* energy of particle with coolest color */
#define PARTICLE_EMAX 3000.0       /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 150.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 1   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.5   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.5    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 2.5    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 2.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define PARTICLE_MASS 200.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0   /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.5     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.5     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 100.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e9    /* harmonic potential of obstacles */
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

#define ADD_EFIELD 1     /* set to 1 to add an electric field */
#define EFIELD 1000.0       /* value of electric field */
#define EFIELD_Y 0.0     /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 0.0       /* value of magnetic field */
#define CHARGE 0.0       /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 5000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 0      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define OBSTACLE_MASS 100.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 5.0e5   /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 2.0e5   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 1         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 0          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 0          /* set to 1 to draw cross on particles of negative charge */
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
#define OBSTACLE_RADIUS 0.015  /* radius of obstacle for circle boundary conditions */
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
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 40     /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 500  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 100    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 40000    /* length of recorded trajectory */
#define TRACER_LUM_FACTOR 100.0     /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 4.0  /* relative mass of tracer particle */
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

#define REACTION_DIFFUSION 0      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 262           /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 1                /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 0         /* initial condition of particles */
#define REACTION_DIST 1.2         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0       /* change of kinetic energy in reaction */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 6              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 1000     /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 10 /* size limitation on smaller cluster */
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
#define DEACIVATE_CLOSE_PAIRS 1 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e10    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 4         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 9      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.1      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 1  /* set to 1 to pair particle of type 1 */
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
#define HASHY 20      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 20 February 2025 - Possibly a metamaterial: Waves in a 15 by 9 lattice of connected discs ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:**

```
    init_wave_flat(phi, psi, xy_in);
    
    /* add oscillating waves */
    for (k=0; k<N_SOURCES; k++)
    {
        wave_source_x[k] = 0.0;
        wave_source_y[k] = 0.0;
    }
    for (source = 0; source < N_SOURCES; source++)
    {
        source_shift[source] = 0;
        if ((ADD_OSCILLATING_SOURCE)&&(i%(OSCILLATING_SOURCE_PERIOD) == source_shift[source])&&(i<MAX_PULSING_TIME))
        {
            if (ALTERNATE_OSCILLATING_SOURCE) sign[source] = -sign[source];
            add_circular_wave(-sign[source]*INITIAL_AMP, wave_source_x[source], wave_source_y[source], phi, psi, xy_in);
        }
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
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

#define B_DOMAIN 93         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.5    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.25	    /* parameter controlling the dimensions of domain */
#define MU 0.075            /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 9            /* number of grid point for grid of disks */
#define WALL_WIDTH 0.022    /* width of wall separating lenses */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */

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
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 61  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.025         /* frequency of periodic excitation */
#define AMPLITUDE 0.5      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.025       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 40.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 13    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 1500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 4450         /* number of frames of movie */
#define NVID 7              /* number of iterations between images displayed on screen */
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

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 13     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 14   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.25   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 0.2    /* additional scaling factor for wave amplitude */
#define E_SCALE 20.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.85   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */

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

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
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

### 19 February 2025 - A lighthouse beam in a hexagonal ring ###

**Program:** `particle_trajectory.c` 

**Initial condition in function `animation()`:**

```
    x = 0.999*LAMBDA;
    y = 0.0;
```

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define XMIN -1.2
#define XMAX 2.8	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1   /* pattern of circles */
#define POLYLINE_PATTERN 2  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */
#define NABSCIRCLES 12       /* number of absorbing circles */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100000        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDA 1.0 	/* parameter controlling shape of domain */
#define MU 0.008        /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define LAMBDA_B 0.577350269    /* parameter controlling shape of domain (for P_POLYRING) */
#define NPOLY_B 6           /* number of sides of second polygon */
#define APOLY_B 1.0         /* angle by which to turn second polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define TRAJ_LENGTH 120 /* length of trajectory */
#define PLOT_NMAX 80    /* max length on length plot */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define HEATMAP 0       /* set to 1 to show heat map of particles */
#define PLOT_HEATMAP_AVERAGE 0      /* set to 1 to plot average number of particles in heat map */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 0   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 1   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_ANGLE 1    /* set to 1 to print angle of trajectory 0 */
#define PRINT_TRAJECTORY_SLOPE 0    /* set to 1 to print slope of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 1         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 0         /* set to 1 to use log scale for plot of lengths */
#define LENGTH_PLOT_POLAR 0         /* set to 1 to plot lengths in polar coordinates */
#define MIN_ANGLE 120.0          /* range of angles of trajectory */
#define MAX_ANGLE 240.0         /* range of angles of trajectory */
#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define SLOW_AT_LONG_TRAJ 0     /* set to 1 to slow down movie for long trajectories */
#define ADD_SUCCESS_GALLERY 0   /* set to 1 to add gallery of successful trajectories at end of movie */
#define SUCCESS_GALLERY_FRAMES 25   /* number of frames per success */
#define EXIT_BOTH_WAYS 1        /* set to 1 to add exits to he left to succesful trajectories */

#define NSTEPS 4500       /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 11    /* Color palette, see list in global_pdes.c  */

#define NCOLORS 64       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define COLOR_HUEMIN 60  /* minimal color hue */
#define COLOR_HUEMAX 360 /* maximal color hue */
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define GRAPH_HUE 90.0      /* color hue of graph */
#define SUCCESS_HUE 180.0   /* color hue of success */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 250       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define END_FRAMES  250   /* number of frames at end of movie */

#define NXMAZE 25      /* width of maze */
#define NYMAZE 25      /* height of maze */
#define MAZE_MAX_NGBH 8     /* max number of neighbours of maze cell */
#define RAND_SHIFT 10       /* seed of random number generator */
#define MAZE_XSHIFT -0.7     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.4     /* radius of tounded corners in maze */
#define CLOSE_MAZE 0        /* set to 1 to close maze exits */

```

### 18 February 2025 - Electric bowling: Charged particles in an electric field colliding with neutral particles ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
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

#define INITXMIN -1.6
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.2
#define ADDXMAX -2.0	  /* x interval for adding particles */
#define ADDYMIN -1.0
#define ADDYMAX 1.0       /* y interval for adding particles */
#define ADDRMIN 4.75 
#define ADDRMAX 6.0     /* r interval for adding particles */

#define BCXMIN -2.5
#define BCXMAX 2.5	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -1.8
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 91     /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.12  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 3      /* max number of neighbours to be pinned */
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
#define TYPE_PROPORTION 0.5 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 4.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.03           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.007         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.075           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 20           /* number of grid point for grid of disks */
#define NGRIDY 10          /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 40
#define NOBSY 24            /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 4400      /* number of frames of movie */
#define NVID 100          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 30     /* time after which to start saving frames */
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

#define BOUNDARY_COND 23

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 12        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 7          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 5        /* type of background coloring, see list in global_ljones.c */
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
#define OBSTACLE_AREA_SHADE_FACTOR 80.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT -0.5     /* shift in color hue (for some cyclic palettes) */

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
#define PARTICLE_EMIN 100.0         /* energy of particle with coolest color */
#define PARTICLE_EMAX 1500.0       /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 150.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 1   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.5   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.5    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 2.5    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 2.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define PARTICLE_MASS 20.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0   /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.5     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.5     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 100.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e9    /* harmonic potential of obstacles */
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

#define ADD_EFIELD 1     /* set to 1 to add an electric field */
#define EFIELD 1000.0       /* value of electric field */
#define EFIELD_Y 0.0     /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 0.0       /* value of magnetic field */
#define CHARGE 0.0       /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 5000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 0      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define OBSTACLE_MASS 100.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 5.0e5   /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 2.0e5   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 1         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 0          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 0          /* set to 1 to draw cross on particles of negative charge */
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
#define OBSTACLE_RADIUS 0.015  /* radius of obstacle for circle boundary conditions */
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
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 50     /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 500  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 100    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 40000    /* length of recorded trajectory */
#define TRACER_LUM_FACTOR 100.0     /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 4.0  /* relative mass of tracer particle */
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

#define REACTION_DIFFUSION 0      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 262           /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 1                /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 0         /* initial condition of particles */
#define REACTION_DIST 1.2         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0       /* change of kinetic energy in reaction */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 6              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 1000     /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 10 /* size limitation on smaller cluster */
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
#define DEACIVATE_CLOSE_PAIRS 1 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e10    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 4         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 9      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.1      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 1  /* set to 1 to pair particle of type 1 */
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
#define HASHY 20      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 17 February 2025 - 3D render of waves in a lattice of connected circles ###

**Program:** `wave_3d.c` 

**Initial condition in function `animation()`:**

```
    init_wave_flat_mod(phi, psi, xy_in);
    
    /* add oscillating waves */
    for (source = 0; source < N_SOURCES; source++)
    {
        wave_source_x[0] = 0.0;
        wave_source_y[0] = 0.0;
        if ((ADD_OSCILLATING_SOURCE)&&(i%(OSCILLATING_SOURCE_PERIOD) == 1)&&(i<MAX_PULSING_TIME))
        {
            if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
            add_circular_wave_mod(-sign*INITIAL_AMP, wave_source_x[source], wave_source_y[source], phi, psi, xy_in);
        }
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 1920          /* number of grid points on x axis */
#define NY 1150          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1        /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 0.8 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 93        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 2   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 181             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.0 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.25   /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.25	    /* parameter controlling the dimensions of domain */
#define MU 0.15             /* parameter controlling the dimensions of domain */
#define NPOLY 4             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 7            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 2000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 20.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 9            /* number of grid point for grid of disks */
#define NGRIDY 5            /* number of grid point for grid of disks */
#define WALL_WIDTH 0.03     /* width of wall separating lenses */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -2.9
#define ISO_XSHIFT_RIGHT 1.4
#define ISO_YSHIFT_LEFT -0.15
#define ISO_YSHIFT_RIGHT -0.15 
#define ISO_SCALE 0.5           /* coordinates for isospectral billiards */


/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 3  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.025        /* frequency of periodic excitation */
#define AMPLITUDE 0.5      /* amplitude of periodic excitation */ 
#define ACHIRP 0.2         /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.08      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 1.0e-6         /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 25     /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define MAX_PULSING_TIME 1000           /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define WAVE_PACKET_RADIUS 20            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

#define PRECOMPUTE_BC 0     /* set to 1 to compute neighbours for Laplacian in advance */

/* Parameters for length and speed of simulation */

#define NSTEPS 2500        /* number of frames of movie */
#define NVID 3             /* number of iterations between images displayed on screen */
#define NSEG 1000          /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* set to 1 to print speed of moving source */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 3         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */
#define ROTATE_VIEW_WHILE_FADE 1    /* set to 1 to keep rotating viewpoint during fade */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define ZPLOT 103     /* wave height */
#define CPLOT 103     /* color scheme */

#define ZPLOT_B 108 
#define CPLOT_B 108        /* plot type for second movie */

#define CHANGE_LUMINOSITY 1     /* set to 1 to let luminosity depend on energy flux intensity */
#define FLUX_WINDOW 30           /* size of averaging window of flux intensity */
#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define FLOOR_ZCOORD 1          /* set to 1 to draw only facets with z not too negative */
#define DRAW_BILLIARD 0         /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 0   /* set to 1 to draw front of boundary after drawing wave */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw construction lines of certain domains */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental, draw outside of billiard in gray */

#define PLOT_SCALE_ENERGY 0.4          /* vertical scaling in energy plot */
#define PLOT_SCALE_LOG_ENERGY 0.5       /* vertical scaling in log energy plot */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 180.0  /* total angle of rotation during simulation */

/* Color schemes */

#define COLOR_PALETTE 11       /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 14     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 0.5   /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSHIFT_AMPLITUDE 0.0   /* additional shift for wave amplitude */
#define VSCALE_ENERGY 3.0      /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0      /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0     /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 40.0        /* scaling factor for energy representation */
#define LOG_SCALE 1.5       /* scaling factor for energy log representation */
#define LOG_SHIFT 0.25      /* shift of colors on log scale */
#define LOG_ENERGY_FLOOR -10.0    /* floor value for log of (total) energy */
#define LOG_MEAN_ENERGY_SHIFT 1.0   /* additional shift for log of mean energy */
#define FLUX_SCALE 4000.0    /* scaling factor for energy flux representation */
#define FLUX_CSCALE 2.0      /* scaling factor for color in energy flux representation */
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 240.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -200.0    /* amplitude of variation of hue for color scheme C_HUE */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

#define DRAW_COLOR_SCHEME 1       /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.5      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.8    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0     /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

#define ADD_POTENTIAL 0         /* set to 1 to add potential to z coordinate */
#define POTENTIAL 10
#define POT_FACT 200.0
#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define MU_B 1.0           /* parameter controlling the dimensions of domain */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define WAVE_PROFILE_X 2.1      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y -1.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 1  /* set to 1 to draw time-average of wave profile squared*/
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define SHADE_2D 0       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
/* end of constants only used by sub_wave and sub_maze */


/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters controlling 3D projection */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-8.0, -6.0, 6.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 0.05    /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 1.8   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0         /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.1           /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.2           /* overall y shift for REP_PROJ_3D representation */

```

### 16 February 2025 - A lighthouse beam in a square with an octagonal hole ###

**Program:** `particle_trajectory.c` 

**Initial condition in function `animation()`:**

```
    x = 0.999*LAMBDA*cos(0.5*PID);
    y = -0.999*LAMBDA*sin(0.5*PID);
```

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define XMIN -1.2
#define XMAX 2.8	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1   /* pattern of circles */
#define POLYLINE_PATTERN 2  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */
#define NABSCIRCLES 12       /* number of absorbing circles */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100000        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDA 1.4 	/* parameter controlling shape of domain */
#define MU 0.008        /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 4             /* number of sides of polygon */
#define APOLY 0.5           /* angle by which to turn polygon, in units of Pi/2 */ 
#define LAMBDA_B 0.4666666667    /* parameter controlling shape of domain (for P_POLYRING) */
#define NPOLY_B 8           /* number of sides of second polygon */
#define APOLY_B 0.25         /* angle by which to turn second polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define TRAJ_LENGTH 120 /* length of trajectory */
#define PLOT_NMAX 110   /* max length on length plot */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define HEATMAP 0       /* set to 1 to show heat map of particles */
#define PLOT_HEATMAP_AVERAGE 0      /* set to 1 to plot average number of particles in heat map */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 0   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 1   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_ANGLE 1    /* set to 1 to print angle of trajectory 0 */
#define PRINT_TRAJECTORY_SLOPE 0    /* set to 1 to print slope of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 1         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 0         /* set to 1 to use log scale for plot of lengths */
#define LENGTH_PLOT_POLAR 0         /* set to 1 to plot lengths in polar coordinates */
#define MIN_ANGLE 90.0          /* range of angles of trajectory */
#define MAX_ANGLE 180.0         /* range of angles of trajectory */
#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define SLOW_AT_LONG_TRAJ 0     /* set to 1 to slow down movie for long trajectories */
#define ADD_SUCCESS_GALLERY 0   /* set to 1 to add gallery of successful trajectories at end of movie */
#define SUCCESS_GALLERY_FRAMES 25   /* number of frames per success */
#define EXIT_BOTH_WAYS 1        /* set to 1 to add exits to he left to succesful trajectories */

#define NSTEPS 4900       /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 14    /* Color palette, see list in global_pdes.c  */

#define NCOLORS 64       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define COLOR_HUEMIN 60  /* minimal color hue */
#define COLOR_HUEMAX 360 /* maximal color hue */
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define GRAPH_HUE 90.0      /* color hue of graph */
#define SUCCESS_HUE 180.0   /* color hue of success */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 250       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define END_FRAMES  250   /* number of frames at end of movie */

#define NXMAZE 25      /* width of maze */
#define NYMAZE 25      /* height of maze */
#define MAZE_MAX_NGBH 8     /* max number of neighbours of maze cell */
#define RAND_SHIFT 10       /* seed of random number generator */
#define MAZE_XSHIFT -0.7     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.4     /* radius of tounded corners in maze */
#define CLOSE_MAZE 0        /* set to 1 to close maze exits */

```

### 15 February 2025 - Electromagnetic bowling with mass ratio 1000 ###

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

#define INITXMIN -1.6
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.2
#define ADDXMAX -2.0	  /* x interval for adding particles */
#define ADDYMIN -1.0
#define ADDYMAX 1.0       /* y interval for adding particles */
#define ADDRMIN 4.75 
#define ADDRMAX 6.0     /* r interval for adding particles */

#define BCXMIN -2.5
#define BCXMAX 2.5	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -1.8
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 91     /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.12  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 3      /* max number of neighbours to be pinned */
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
#define TYPE_PROPORTION 0.5 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 4.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.03           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.007         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.075           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 20           /* number of grid point for grid of disks */
#define NGRIDY 10          /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 40
#define NOBSY 24            /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2200      /* number of frames of movie */
#define NVID 10000          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 30     /* time after which to start saving frames */
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

#define BOUNDARY_COND 23

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 12        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 7          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 5        /* type of background coloring, see list in global_ljones.c */
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
#define OBSTACLE_AREA_SHADE_FACTOR 80.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT -0.5     /* shift in color hue (for some cyclic palettes) */

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
#define PARTICLE_EMIN 100.0         /* energy of particle with coolest color */
#define PARTICLE_EMAX 1500.0       /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 150.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 1   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.5   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.5    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-7    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 2.5    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 2.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define PARTICLE_MASS 2000.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0   /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.5     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.5     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 50.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e9    /* harmonic potential of obstacles */
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

#define ADD_EFIELD 1     /* set to 1 to add an electric field */
#define EFIELD 0.0       /* value of electric field */
#define EFIELD_Y 1500.0       /* value of electric field */
#define ADD_BFIELD 1     /* set to 1 to add a magnetic field */
#define BFIELD 1200.0      /* value of magnetic field */
#define CHARGE 0.0       /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 5000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 0      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define OBSTACLE_MASS 100.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 5.0e5   /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 2.0e5   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 1         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 0          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 0          /* set to 1 to draw cross on particles of negative charge */
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
#define OBSTACLE_RADIUS 0.015  /* radius of obstacle for circle boundary conditions */
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
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 50     /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 500  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 100    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 40000    /* length of recorded trajectory */
#define TRACER_LUM_FACTOR 100.0     /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 4.0  /* relative mass of tracer particle */
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

#define REACTION_DIFFUSION 0      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 262           /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 1                /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 0         /* initial condition of particles */
#define REACTION_DIST 1.2         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0       /* change of kinetic energy in reaction */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 6              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 1000     /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 10 /* size limitation on smaller cluster */
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
#define DEACIVATE_CLOSE_PAIRS 1 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e10    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 4         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 9      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.1      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 1  /* set to 1 to pair particle of type 1 */
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
#define HASHY 20      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 14 February 2025 - Waves of higher frequency in a hex lattice of connected circles ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:**

```
    init_wave_flat(phi, psi, xy_in);
    
    /* add oscillating waves */
    for (k=0; k<N_SOURCES; k++)
    {
        wave_source_x[k] = 0.0;
        wave_source_y[k] = 0.0;
    }
    for (source = 0; source < N_SOURCES; source++)
    {
        source_shift[source] = 0;
        if ((ADD_OSCILLATING_SOURCE)&&(i%(OSCILLATING_SOURCE_PERIOD) == source_shift[source])&&(i<MAX_PULSING_TIME))
        {
            if (ALTERNATE_OSCILLATING_SOURCE) sign[source] = -sign[source];
            add_circular_wave(-sign[source]*INITIAL_AMP, wave_source_x[source], wave_source_y[source], phi, psi, xy_in);
        }
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
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

#define B_DOMAIN 94         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.5    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.25	    /* parameter controlling the dimensions of domain */
#define MU 0.15             /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 9            /* number of grid point for grid of disks */
#define NGRIDY 6            /* number of grid point for grid of disks */
#define WALL_WIDTH 0.022    /* width of wall separating lenses */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */

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
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 61  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.025         /* frequency of periodic excitation */
#define AMPLITUDE 0.5      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.025       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 40.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 13    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 1500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 4400         /* number of frames of movie */
#define NVID 7              /* number of iterations between images displayed on screen */
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

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 14     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 14   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.25   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 0.2    /* additional scaling factor for wave amplitude */
#define E_SCALE 20.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.85   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */

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

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
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

### 13 February 2025 - A lighthouse beam in a pentagon ###

**Program:** `particle_trajectory.c` 

**Initial condition in function `animation()`:**

```
    x = 0.999*LAMBDA;
    y = 0.0;
```

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define XMIN -1.2
#define XMAX 2.8	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1   /* pattern of circles */
#define POLYLINE_PATTERN 5  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */
#define NABSCIRCLES 3       /* number of absorbing circles */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100000        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDA 1.00 	/* parameter controlling shape of domain */
#define MU 0.012        /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define TRAJ_LENGTH 350 /* length of trajectory */
#define PLOT_NMAX 330   /* max length on length plot */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define HEATMAP 0       /* set to 1 to show heat map of particles */
#define PLOT_HEATMAP_AVERAGE 0      /* set to 1 to plot average number of particles in heat map */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 0   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 1   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_ANGLE 1    /* set to 1 to print angle of trajectory 0 */
#define PRINT_TRAJECTORY_SLOPE 0    /* set to 1 to print slope of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 1         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 0         /* set to 1 to use log scale for plot of lengths */
#define LENGTH_PLOT_POLAR 0         /* set to 1 to plot lengths in polar coordinates */
#define MIN_ANGLE 126.0         /* range of angles of trajectory */
#define MAX_ANGLE 180.0         /* range of angles of trajectory */
#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define SLOW_AT_LONG_TRAJ 0     /* set to 1 to slow down movie for long trajectories */
#define ADD_SUCCESS_GALLERY 0   /* set to 1 to add gallery of successful trajectories at end of movie */
#define SUCCESS_GALLERY_FRAMES 25   /* number of frames per success */
#define EXIT_BOTH_WAYS 1        /* set to 1 to add exits to he left to succesful trajectories */

#define NSTEPS 5000       /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 10    /* Color palette, see list in global_pdes.c  */

#define NCOLORS 64       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define COLOR_HUEMIN 60  /* minimal color hue */
#define COLOR_HUEMAX 360 /* maximal color hue */
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define GRAPH_HUE 90.0      /* color hue of graph */
#define SUCCESS_HUE 180.0   /* color hue of success */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 250       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define END_FRAMES  100   /* number of frames at end of movie */

#define NXMAZE 25      /* width of maze */
#define NYMAZE 25      /* height of maze */
#define MAZE_MAX_NGBH 8     /* max number of neighbours of maze cell */
#define RAND_SHIFT 10       /* seed of random number generator */
#define MAZE_XSHIFT -0.7     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.4     /* radius of tounded corners in maze */
#define CLOSE_MAZE 0        /* set to 1 to close maze exits */

```

### 12 February 2025 - Electromagnetic bowling with mass ratio 250 ###

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

#define INITXMIN -1.6
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.2
#define ADDXMAX -2.0	  /* x interval for adding particles */
#define ADDYMIN -1.0
#define ADDYMAX 1.0       /* y interval for adding particles */
#define ADDRMIN 4.75 
#define ADDRMAX 6.0     /* r interval for adding particles */

#define BCXMIN -2.5
#define BCXMAX 2.5	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -1.8
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 91     /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.12  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 3      /* max number of neighbours to be pinned */
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
#define TYPE_PROPORTION 0.5 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 4.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.03           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.007         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.075           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 20           /* number of grid point for grid of disks */
#define NGRIDY 10          /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 40
#define NOBSY 24            /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2200      /* number of frames of movie */
#define NVID 10000          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 30     /* time after which to start saving frames */
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

#define BOUNDARY_COND 23

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 12        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 7          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 5        /* type of background coloring, see list in global_ljones.c */
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
#define OBSTACLE_AREA_SHADE_FACTOR 80.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT -0.5     /* shift in color hue (for some cyclic palettes) */

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
#define PARTICLE_EMIN 100.0         /* energy of particle with coolest color */
#define PARTICLE_EMAX 1500.0       /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 150.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 1   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.5   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.5    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-7    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 2.5    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 2.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define PARTICLE_MASS 500.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0   /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.5     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.5     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 50.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e9    /* harmonic potential of obstacles */
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

#define ADD_EFIELD 1     /* set to 1 to add an electric field */
#define EFIELD 0.0       /* value of electric field */
#define EFIELD_Y 1500.0       /* value of electric field */
#define ADD_BFIELD 1     /* set to 1 to add a magnetic field */
#define BFIELD 1200.0      /* value of magnetic field */
#define CHARGE 0.0       /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 5000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 0      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define OBSTACLE_MASS 100.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 5.0e5   /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 2.0e5   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 1         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 0          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 0          /* set to 1 to draw cross on particles of negative charge */
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
#define OBSTACLE_RADIUS 0.015  /* radius of obstacle for circle boundary conditions */
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
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 50     /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 500  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 100    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 40000    /* length of recorded trajectory */
#define TRACER_LUM_FACTOR 100.0     /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 4.0  /* relative mass of tracer particle */
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

#define REACTION_DIFFUSION 0      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 262           /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 1                /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 0         /* initial condition of particles */
#define REACTION_DIST 1.2         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0       /* change of kinetic energy in reaction */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 6              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 1000     /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 10 /* size limitation on smaller cluster */
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
#define DEACIVATE_CLOSE_PAIRS 1 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e10    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 4         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 9      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.1      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 1  /* set to 1 to pair particle of type 1 */
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
#define HASHY 20      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 11 February 2025 - Waves in a hex lattice of connected circles ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:**

```
    init_wave_flat(phi, psi, xy_in);
    
    /* add oscillating waves */
    for (k=0; k<N_SOURCES; k++)
    {
        wave_source_x[k] = 0.0;
        wave_source_y[k] = 0.0;
    }
    for (source = 0; source < N_SOURCES; source++)
    {
        source_shift[source] = 0;
        if ((ADD_OSCILLATING_SOURCE)&&(i%(OSCILLATING_SOURCE_PERIOD) == source_shift[source])&&(i<MAX_PULSING_TIME))
        {
            if (ALTERNATE_OSCILLATING_SOURCE) sign[source] = -sign[source];
            add_circular_wave(-sign[source]*INITIAL_AMP, wave_source_x[source], wave_source_y[source], phi, psi, xy_in);
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

#define B_DOMAIN 94         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.5    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.25	    /* parameter controlling the dimensions of domain */
#define MU 0.15             /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 9            /* number of grid point for grid of disks */
#define NGRIDY 6            /* number of grid point for grid of disks */
#define WALL_WIDTH 0.0275   /* width of wall separating lenses */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */

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
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 61  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.025         /* frequency of periodic excitation */
#define AMPLITUDE 0.5      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.025       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 40.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 25    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 1500            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 3200         /* number of frames of movie */
#define NVID 7              /* number of iterations between images displayed on screen */
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

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 14   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.25   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 0.2    /* additional scaling factor for wave amplitude */
#define E_SCALE 20.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.85   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */

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

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
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

### 10 February 2025 - A lighthouse beam in an equilateral triangle ###

**Program:** `particle_trajectory.c` 

**Initial condition in function `animation()`:**

```
    x = 0.999*LAMBDA*cos(PI/6.0);
    y = -0.999*LAMBDA*sin(PI/6.0);
```

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define XMIN -1.2
#define XMAX 2.8	/* x interval */
#define YMIN -0.925
#define YMAX 1.325	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1   /* pattern of circles */
#define POLYLINE_PATTERN 5  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */
#define NABSCIRCLES 3       /* number of absorbing circles */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100000        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDA 1.00 	/* parameter controlling shape of domain */
#define MU 0.01         /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define TRAJ_LENGTH 350 /* length of trajectory */
#define PLOT_NMAX 330   /* max length on length plot */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define HEATMAP 0       /* set to 1 to show heat map of particles */
#define PLOT_HEATMAP_AVERAGE 0      /* set to 1 to plot average number of particles in heat map */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 0   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 1   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_ANGLE 1    /* set to 1 to print angle of trajectory 0 */
#define PRINT_TRAJECTORY_SLOPE 0    /* set to 1 to print slope of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 1         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 0         /* set to 1 to use log scale for plot of lengths */
#define LENGTH_PLOT_POLAR 0         /* set to 1 to plot lengths in polar coordinates */
#define MIN_ANGLE 120.0         /* range of angles of trajectory */
#define MAX_ANGLE 180.0         /* range of angles of trajectory */
#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define SLOW_AT_LONG_TRAJ 0     /* set to 1 to slow down movie for long trajectories */
#define ADD_SUCCESS_GALLERY 0   /* set to 1 to add gallery of successful trajectories at end of movie */
#define SUCCESS_GALLERY_FRAMES 25   /* number of frames per success */
#define EXIT_BOTH_WAYS 1        /* set to 1 to add exits to he left to succesful trajectories */

#define NSTEPS 5500       /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 17    /* Color palette, see list in global_pdes.c  */

#define NCOLORS 64       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define COLOR_HUEMIN 120   /* minimal color hue */
#define COLOR_HUEMAX 240 /* maximal color hue */
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define GRAPH_HUE 180.0      /* color hue of graph */
#define SUCCESS_HUE 180.0   /* color hue of success */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 250       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define END_FRAMES  100   /* number of frames at end of movie */

#define NXMAZE 25      /* width of maze */
#define NYMAZE 25      /* height of maze */
#define MAZE_MAX_NGBH 8     /* max number of neighbours of maze cell */
#define RAND_SHIFT 10       /* seed of random number generator */
#define MAZE_XSHIFT -0.7     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.4     /* radius of tounded corners in maze */
#define CLOSE_MAZE 0        /* set to 1 to close maze exits */

```

### 09 February 2025 - Electromagnetic bowling with mass ratio 100 ###

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

#define INITXMIN -1.6
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.2
#define ADDXMAX -2.0	  /* x interval for adding particles */
#define ADDYMIN -1.0
#define ADDYMAX 1.0       /* y interval for adding particles */
#define ADDRMIN 4.75 
#define ADDRMAX 6.0     /* r interval for adding particles */

#define BCXMIN -2.5
#define BCXMAX 2.5	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -1.8
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 91     /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.12  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 3      /* max number of neighbours to be pinned */
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
#define TYPE_PROPORTION 0.5 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 4.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.03           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.007         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.075           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 20           /* number of grid point for grid of disks */
#define NGRIDY 10          /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 40
#define NOBSY 24            /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2200      /* number of frames of movie */
#define NVID 10000          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 30     /* time after which to start saving frames */
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

#define BOUNDARY_COND 23

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 12        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 7          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 5        /* type of background coloring, see list in global_ljones.c */
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
#define OBSTACLE_AREA_SHADE_FACTOR 80.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT -0.5     /* shift in color hue (for some cyclic palettes) */

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
#define PARTICLE_EMIN 100.0         /* energy of particle with coolest color */
#define PARTICLE_EMAX 1500.0       /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 150.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 1   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.5   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.5    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-7    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 2.5    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 2.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define PARTICLE_MASS 200.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0   /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.5     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.5     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 50.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e9    /* harmonic potential of obstacles */
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

#define ADD_EFIELD 1     /* set to 1 to add an electric field */
#define EFIELD 0.0       /* value of electric field */
#define EFIELD_Y 1500.0       /* value of electric field */
#define ADD_BFIELD 1     /* set to 1 to add a magnetic field */
#define BFIELD 1200.0      /* value of magnetic field */
#define CHARGE 0.0       /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 5000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 0      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define OBSTACLE_MASS 100.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 5.0e5   /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 2.0e5   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 1         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 0          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 0          /* set to 1 to draw cross on particles of negative charge */
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
#define OBSTACLE_RADIUS 0.015  /* radius of obstacle for circle boundary conditions */
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
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 50     /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 500  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 100    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 40000    /* length of recorded trajectory */
#define TRACER_LUM_FACTOR 100.0     /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 4.0  /* relative mass of tracer particle */
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

#define REACTION_DIFFUSION 0      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 262           /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 1                /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 0         /* initial condition of particles */
#define REACTION_DIST 1.2         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0       /* change of kinetic energy in reaction */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 6              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 1000     /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 10 /* size limitation on smaller cluster */
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
#define DEACIVATE_CLOSE_PAIRS 1 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e10    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 4         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 9      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.1      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 1  /* set to 1 to pair particle of type 1 */
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
#define HASHY 20      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 08 February 2025 - Waves of higher frequency in a lattice of connected circles ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:**

```
    init_wave_flat(phi, psi, xy_in);
    
    /* add oscillating waves */
    for (k=0; k<N_SOURCES; k++)
    {
        wave_source_x[k] = 0.0;
        wave_source_y[k] = 0.0;
    }
    for (source = 0; source < N_SOURCES; source++)
    {
        source_shift[source] = 0;
        if ((ADD_OSCILLATING_SOURCE)&&(i%(OSCILLATING_SOURCE_PERIOD) == source_shift[source])&&(i<MAX_PULSING_TIME))
        {
            if (ALTERNATE_OSCILLATING_SOURCE) sign[source] = -sign[source];
            add_circular_wave(-sign[source]*INITIAL_AMP, wave_source_x[source], wave_source_y[source], phi, psi, xy_in);
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

#define B_DOMAIN 93         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.5    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.25	    /* parameter controlling the dimensions of domain */
#define MU 0.15             /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 9            /* number of grid point for grid of disks */
#define NGRIDY 5            /* number of grid point for grid of disks */
#define WALL_WIDTH 0.03     /* width of wall separating lenses */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */

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
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 61  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.025         /* frequency of periodic excitation */
#define AMPLITUDE 0.5      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.025       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 40.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 12    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 1000            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 3000         /* number of frames of movie */
#define NVID 7              /* number of iterations between images displayed on screen */
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

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 11     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 17   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.25   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 0.2    /* additional scaling factor for wave amplitude */
#define E_SCALE 5.0        /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.85   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */

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

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
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

### 07 February 2025 - A lighthouse beam in a 100000-gon ###

**Program:** `particle_trajectory.c` 

**Initial condition in function `animation()`:**

``` x = LAMBDA;
    y = 0.0;
```

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define XMIN -1.2
#define XMAX 2.8	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1   /* pattern of circles */
#define POLYLINE_PATTERN 5  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */
#define NABSCIRCLES 1       /* number of absorbing circles */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100000        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDA 1.00 	/* parameter controlling shape of domain */
#define MU 0.03         /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 100000             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define TRAJ_LENGTH 320 /* length of trajectory */
#define PLOT_NMAX 305   /* max length on length plot */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define HEATMAP 0       /* set to 1 to show heat map of particles */
#define PLOT_HEATMAP_AVERAGE 0      /* set to 1 to plot average number of particles in heat map */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 0   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 1   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_ANGLE 1    /* set to 1 to print angle of trajectory 0 */
#define PRINT_TRAJECTORY_SLOPE 0    /* set to 1 to print slope of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 1         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 0         /* set to 1 to use log scale for plot of lengths */
#define LENGTH_PLOT_POLAR 0         /* set to 1 to plot lengths in polar coordinates */
#define MIN_ANGLE 180.0         /* range of angles of trajectory */
#define MAX_ANGLE 270.0         /* range of angles of trajectory */
#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define SLOW_AT_LONG_TRAJ 0     /* set to 1 to slow down movie for long trajectories */
#define ADD_SUCCESS_GALLERY 0   /* set to 1 to add gallery of successful trajectories at end of movie */
#define SUCCESS_GALLERY_FRAMES 25   /* number of frames per success */
#define EXIT_BOTH_WAYS 1        /* set to 1 to add exits to he left to succesful trajectories */

#define NSTEPS 4900       /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 14    /* Color palette, see list in global_pdes.c  */

#define NCOLORS 64       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define COLOR_HUEMIN 120   /* minimal color hue */
#define COLOR_HUEMAX 240 /* maximal color hue */
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define GRAPH_HUE 90.0      /* color hue of graph */
#define SUCCESS_HUE 180.0   /* color hue of success */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 250       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define END_FRAMES  100   /* number of frames at end of movie */

#define NXMAZE 25      /* width of maze */
#define NYMAZE 25      /* height of maze */
#define MAZE_MAX_NGBH 8     /* max number of neighbours of maze cell */
#define RAND_SHIFT 10       /* seed of random number generator */
#define MAZE_XSHIFT -0.7     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.4     /* radius of tounded corners in maze */
#define CLOSE_MAZE 0        /* set to 1 to close maze exits */

```

### 06 February 2025 - Electromagnetic bowling with mass ratio 20 ###

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

#define INITXMIN -1.6
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.2
#define ADDXMAX -2.0	  /* x interval for adding particles */
#define ADDYMIN -1.0
#define ADDYMAX 1.0       /* y interval for adding particles */
#define ADDRMIN 4.75 
#define ADDRMAX 6.0     /* r interval for adding particles */

#define BCXMIN -2.5
#define BCXMAX 2.5	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -1.8
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 91     /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.12  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 3      /* max number of neighbours to be pinned */
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
#define TYPE_PROPORTION 0.5 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 4.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.03           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.007         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.075           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 20           /* number of grid point for grid of disks */
#define NGRIDY 10          /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 40
#define NOBSY 24            /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2200      /* number of frames of movie */
#define NVID 10000          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 30     /* time after which to start saving frames */
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

#define BOUNDARY_COND 23

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 12        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 7          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 5        /* type of background coloring, see list in global_ljones.c */
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
#define OBSTACLE_AREA_SHADE_FACTOR 80.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT -0.5     /* shift in color hue (for some cyclic palettes) */

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
#define PARTICLE_EMIN 100.0         /* energy of particle with coolest color */
#define PARTICLE_EMAX 1500.0       /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 150.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 1   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.5   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.5    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-7    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 2.5    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 2.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define PARTICLE_MASS 40.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0   /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.5     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.5     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 50.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e9    /* harmonic potential of obstacles */
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

#define ADD_EFIELD 1     /* set to 1 to add an electric field */
#define EFIELD 0.0       /* value of electric field */
#define EFIELD_Y 1500.0       /* value of electric field */
#define ADD_BFIELD 1     /* set to 1 to add a magnetic field */
#define BFIELD 1200.0      /* value of magnetic field */
#define CHARGE 0.0       /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 5000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 0      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define OBSTACLE_MASS 100.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 5.0e5   /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 2.0e5   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 1         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 0          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 0          /* set to 1 to draw cross on particles of negative charge */
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
#define OBSTACLE_RADIUS 0.015  /* radius of obstacle for circle boundary conditions */
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
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 50     /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 500  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 100    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 40000    /* length of recorded trajectory */
#define TRACER_LUM_FACTOR 100.0     /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 4.0  /* relative mass of tracer particle */
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

#define REACTION_DIFFUSION 0      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 262           /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 1                /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 0         /* initial condition of particles */
#define REACTION_DIST 1.2         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0       /* change of kinetic energy in reaction */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 6              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 1000     /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 10 /* size limitation on smaller cluster */
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
#define DEACIVATE_CLOSE_PAIRS 1 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e10    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 4         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 9      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.1      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 1  /* set to 1 to pair particle of type 1 */
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
#define HASHY 20      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 05 February 2025 - Waves in a lattice of connected circles ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:**

```
    init_wave_flat(phi, psi, xy_in);
    
    /* add oscillating waves */
    for (k=0; k<N_SOURCES; k++)
    {
        wave_source_x[k] = 0.0;
        wave_source_y[k] = 0.0;
    }
    for (source = 0; source < N_SOURCES; source++)
    {
        source_shift[source] = 0;
        if ((ADD_OSCILLATING_SOURCE)&&(i%(OSCILLATING_SOURCE_PERIOD) == source_shift[source])&&(i<MAX_PULSING_TIME))
        {
            if (ALTERNATE_OSCILLATING_SOURCE) sign[source] = -sign[source];
            add_circular_wave(-sign[source]*INITIAL_AMP, wave_source_x[source], wave_source_y[source], phi, psi, xy_in);
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

#define B_DOMAIN 93         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202  /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.5    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.25	    /* parameter controlling the dimensions of domain */
#define MU 0.15             /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 9            /* number of grid point for grid of disks */
#define NGRIDY 5            /* number of grid point for grid of disks */
#define WALL_WIDTH 0.03     /* width of wall separating lenses */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */

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
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 61  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.025         /* frequency of periodic excitation */
#define AMPLITUDE 0.5      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.025       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 40.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 25    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 1000            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 2500         /* number of frames of movie */
#define NVID 7              /* number of iterations between images displayed on screen */
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

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
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
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.25   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 0.2    /* additional scaling factor for wave amplitude */
#define E_SCALE 5.0        /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.85   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */

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

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
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

### 04 February 2025 - A lighthouse beam in a octagon ###

**Program:** `particle_trajectory.c` 

**Initial condition in function `animation()`:**

```
    x = LAMBDA;
    y = 0.0;
```

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define XMIN -1.2
#define XMAX 2.8	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1   /* pattern of circles */
#define POLYLINE_PATTERN 5  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 50000        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDA 1.00 	/* parameter controlling shape of domain */
#define MU 0.01         /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 8             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define TRAJ_LENGTH 120 /* length of trajectory */
#define PLOT_NMAX 105   /* max length on length plot */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define HEATMAP 0       /* set to 1 to show heat map of particles */
#define PLOT_HEATMAP_AVERAGE 0      /* set to 1 to plot average number of particles in heat map */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 0   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 1   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_ANGLE 1    /* set to 1 to print angle of trajectory 0 */
#define PRINT_TRAJECTORY_SLOPE 0    /* set to 1 to print slope of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 1         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 0         /* set to 1 to use log scale for plot of lengths */
#define LENGTH_PLOT_POLAR 0         /* set to 1 to plot lengths in polar coordinates */
#define MIN_ANGLE 112.5         /* range of angles of trajectory */
#define MAX_ANGLE 180.0         /* range of angles of trajectory */
#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define SLOW_AT_LONG_TRAJ 0     /* set to 1 to slow down movie for long trajectories */
#define ADD_SUCCESS_GALLERY 0   /* set to 1 to add gallery of successful trajectories at end of movie */
#define SUCCESS_GALLERY_FRAMES 25   /* number of frames per success */
#define EXIT_BOTH_WAYS 1        /* set to 1 to add exits to he left to succesful trajectories */

#define NSTEPS 5100       /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 15    /* Color palette, see list in global_pdes.c  */

#define NCOLORS 64       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define COLOR_HUEMIN 120   /* minimal color hue */
#define COLOR_HUEMAX 240 /* maximal color hue */
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define GRAPH_HUE 90.0      /* color hue of graph */
#define SUCCESS_HUE 180.0   /* color hue of success */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 250       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define END_FRAMES  100   /* number of frames at end of movie */

#define NXMAZE 25      /* width of maze */
#define NYMAZE 25      /* height of maze */
#define MAZE_MAX_NGBH 8     /* max number of neighbours of maze cell */
#define RAND_SHIFT 10       /* seed of random number generator */
#define MAZE_XSHIFT -0.7     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.4     /* radius of tounded corners in maze */
#define CLOSE_MAZE 0        /* set to 1 to close maze exits */

```

### 03 February 2025 - Electromagnetic bowling ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
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

#define INITXMIN -1.6
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN -2.2
#define ADDXMAX -2.0	  /* x interval for adding particles */
#define ADDYMIN -1.0
#define ADDYMAX 1.0       /* y interval for adding particles */
#define ADDRMIN 4.75 
#define ADDRMAX 6.0     /* r interval for adding particles */

#define BCXMIN -2.5
#define BCXMAX 2.5	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -1.8
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 0  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 91     /* pattern of obstacles, see list in global_ljones.c */
#define RATTLE_OBSTACLES 0      /* set to 1 to rattle obstacles (for pattern O_SIEVE_B) */
#define OSCILLATE_OBSTACLES 1   /* set to 1 to make obstacles oscillate */ 
#define COUPLE_OBSTACLES 1      /* set to 1 to couple obstacles to neighbours */
#define OBSTACLE_PISC_DISTANCE 0.08  /* minimal distance in Poisson disc process for obstacles, controls density of obstacles */
#define OBSTACLE_COUPLING_DIST 0.12  /* max distance of coupled obstacles */
#define NMAX_OBSTACLE_NEIGHBOURS 8  /* max number of obstacle neighbours */
#define NMAX_OBSTACLE_PINNED 3      /* max number of neighbours to be pinned */
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
#define TYPE_PROPORTION 0.5 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 4.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 7.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.03           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.007         /* parameter controlling radius of added particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.075           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 20           /* number of grid point for grid of disks */
#define NGRIDY 10          /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 40
#define NOBSY 24            /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 5500      /* number of frames of movie */
#define NVID 4000          /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 30     /* time after which to start saving frames */
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

#define BOUNDARY_COND 23

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 12        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 7          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 5        /* type of background coloring, see list in global_ljones.c */
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
#define OBSTACLE_AREA_SHADE_FACTOR 80.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 0        /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_HUE_CLUSTER_SELECTED 90.0    /* Color hue for selected cluster */
#define COLOR_HUE_CLUSTER_NOT_SELECTED 220.0    /* Color hue for selected cluster */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */
#define COLOR_HUESHIFT -0.5     /* shift in color hue (for some cyclic palettes) */

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
#define PARTICLE_EMIN 100.0         /* energy of particle with coolest color */
#define PARTICLE_EMAX 1500.0       /* energy of particle with hottest color */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 150.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 320.0      /* hue of particles of type 0 */
#define HUE_TYPE1 60.0       /* hue of particles of type 1 */
#define HUE_TYPE2 320.0      /* hue of particles of type 2 */
#define HUE_TYPE3 260.0      /* hue of particles of type 3 */
#define HUE_TYPE4 200.0      /* hue of particles of type 4 */
#define HUE_TYPE5 60.0       /* hue of particles of type 5 */
#define HUE_TYPE6 130.0      /* hue of particles of type 6 */
#define HUE_TYPE7 150.0      /* hue of particles of type 7 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 1   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.5   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.5    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-7    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 2.5    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 2.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 1000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 5.0      /* damping coefficient for rotation of particles */
#define PARTICLE_MASS 10.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0   /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.5     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.5     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define V_INITIAL_ADD 50.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.004          /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e9    /* harmonic potential of obstacles */
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

#define ADD_EFIELD 1     /* set to 1 to add an electric field */
#define EFIELD 0.0       /* value of electric field */
#define EFIELD_Y 1500.0       /* value of electric field */
#define ADD_BFIELD 1     /* set to 1 to add a magnetic field */
#define BFIELD 1200.0      /* value of magnetic field */
#define CHARGE 0.0       /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define CHARGE_ADD 1.0   /* charge of added particles */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 5000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 0      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define OBSTACLE_MASS 100.0       /* mass of obstacles, if oscillating */
#define KSPRING_OBSTACLE_OSC 5.0e5   /* spring constant for oscillating obstacles */
#define KSPRING_OBSTACLE_COUPLE 2.0e5   /* spring constant for coupled obstacles */
#define OBSTACLE_HARDCORE 1         /* set to 1 to add "hard core" repulsion between obstacles */
#define KSPRING_OBSTACLE_HARDCORE 1.0e11     /* spring constant for obstacle hard core repulsion */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */
#define EFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */
#define BFIELD_REGION 0          /* space-dependent magnetic field (0 for constant) */

#define ADD_WIND 0          /* set to 1 to add a "wind" friction force */
#define WIND_FORCE 1.35e6    /* force of wind */
#define WIND_YMIN -0.6      /* min altitude of region with wind */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 1.0e5         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 0          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 0          /* set to 1 to draw cross on particles of negative charge */
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
#define OBSTACLE_RADIUS 0.015  /* radius of obstacle for circle boundary conditions */
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
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 150     /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 500  /* final period where no particles are added */
#define SAFETY_FACTOR 4.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 100    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 40000    /* length of recorded trajectory */
#define TRACER_LUM_FACTOR 100.0     /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 4.0  /* relative mass of tracer particle */
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

#define REACTION_DIFFUSION 0      /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 262           /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 1                /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 0         /* initial condition of particles */
#define REACTION_DIST 1.2         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0       /* change of kinetic energy in reaction */
#define COLLISION_TIME 25       /* time during which collisions are shown */
#define COLLISION_RADIUS 2.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 200.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 6              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define CLUSTER_PARTICLES 0     /* set to 1 for particles to form rigid clusters */
#define CLUSTER_MAXSIZE 1000     /* max size of clusters */
#define SMALL_CLUSTER_MAXSIZE 10 /* size limitation on smaller cluster */
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
#define DEACIVATE_CLOSE_PAIRS 1 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e10    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 4         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 5             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 9      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.1      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.035            /* radius of partner particle */
#define PARTICLE_MASS_C 2.0  /* mass or partner particle */
#define CHARGE_C 1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 1   /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 1  /* set to 1 to pair particle of type 1 */
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
#define HASHY 20      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 02 February 2025 - Exciting standing waves in an octagon ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:**

```
    init_wave_flat(phi, psi, xy_in);
    
    /* add oscillating waves */
    for (k=0; k<N_SOURCES; k++)
    {
        angle = DPI*(double)k/(double)N_SOURCES + APOLY*PID;
        wave_source_x[k] = LAMBDA*cos(angle);
        wave_source_y[k] = LAMBDA*sin(angle);
    }
    for (source = 0; source < N_SOURCES; source++)
    {
        source_shift[source] = 0;
        if ((ADD_OSCILLATING_SOURCE)&&(i%(OSCILLATING_SOURCE_PERIOD) == source_shift[source])&&(i<MAX_PULSING_TIME))
        {
            if (ALTERNATE_OSCILLATING_SOURCE) sign[source] = -sign[source];
            add_circular_wave(-sign[source]*INITIAL_AMP, wave_source_x[source], wave_source_y[source], phi, psi, xy_in);
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

#define B_DOMAIN 8         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202   /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.5    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.2	    /* parameter controlling the dimensions of domain */
#define MU 0.5              /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 8             /* number of sides of polygon */
#define APOLY 0.25          /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 25           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.05     /* width of wall separating lenses */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */

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
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 61  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.025         /* frequency of periodic excitation */
#define AMPLITUDE 0.5      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.025       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 40.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 35    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 8                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 1000            /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2200         /* number of frames of movie */
#define NVID 7              /* number of iterations between images displayed on screen */
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

#define INITIAL_AMP 1.25             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 16   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.25   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 0.2    /* additional scaling factor for wave amplitude */
#define E_SCALE 4.0        /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.85   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
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

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
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


### 01 February 2025 - A lighthouse beam in a bestagon ###

**Program:** `particle_trajectory.c` 

**Initial condition in function `animation()`:**

```
    x = LAMBDA;
    y = 0.0;
```

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define XMIN -1.2
#define XMAX 2.8	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1   /* pattern of circles */
#define POLYLINE_PATTERN 5  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 50000        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDA 1.00 	/* parameter controlling shape of domain */
#define MU 0.01         /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define TRAJ_LENGTH 150 /* length of trajectory */
#define PLOT_NMAX 140   /* max length on length plot */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define HEATMAP 0       /* set to 1 to show heat map of particles */
#define PLOT_HEATMAP_AVERAGE 0      /* set to 1 to plot average number of particles in heat map */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 0   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 1   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_ANGLE 1    /* set to 1 to print angle of trajectory 0 */
#define PRINT_TRAJECTORY_SLOPE 0    /* set to 1 to print slope of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 1         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 0         /* set to 1 to use log scale for plot of lengths */
#define LENGTH_PLOT_POLAR 0         /* set to 1 to plot lengths in polar coordinates */
#define MIN_ANGLE 120.0         /* range of angles of trajectory */
#define MAX_ANGLE 180.0         /* range of angles of trajectory */
#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define SLOW_AT_LONG_TRAJ 0     /* set to 1 to slow down movie for long trajectories */
#define ADD_SUCCESS_GALLERY 0   /* set to 1 to add gallery of successful trajectories at end of movie */
#define SUCCESS_GALLERY_FRAMES 25   /* number of frames per success */
#define EXIT_BOTH_WAYS 1        /* set to 1 to add exits to he left to succesful trajectories */

#define NSTEPS 4900       /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 13    /* Color palette, see list in global_pdes.c  */

#define NCOLORS 64       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define COLOR_HUEMIN 120   /* minimal color hue */
#define COLOR_HUEMAX 240 /* maximal color hue */
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define GRAPH_HUE 90.0      /* color hue of graph */
#define SUCCESS_HUE 180.0   /* color hue of success */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 250       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define END_FRAMES  100   /* number of frames at end of movie */

#define NXMAZE 25      /* width of maze */
#define NYMAZE 25      /* height of maze */
#define MAZE_MAX_NGBH 8     /* max number of neighbours of maze cell */
#define RAND_SHIFT 10       /* seed of random number generator */
#define MAZE_XSHIFT -0.7     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.4     /* radius of tounded corners in maze */
#define CLOSE_MAZE 0        /* set to 1 to close maze exits */

```

