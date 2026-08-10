Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Below are parameter values used for different simulations, as well as initial conditions used in function animation. Some simulations use variants of the published code. The list is going to be updated gradually. 


### 31 July 2026 - Increasing the force constant for charged particles on a sphere: Charge ratio 3:1 ###

**Program:** `lennardjones.c` 

**3D part:**

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 0.0
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.25 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 4.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.014 	    /* parameter controlling radius of particles */
#define MU_B 0.009          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2600      /* number of frames of movie */
#define NVID 80          /* number of iterations between images displayed on screen */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 4.0            /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 6.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 100.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1      /* set to 1 to add an electric field */
#define EFIELD 2500000.0  /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 3.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 2500.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 1   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 1.0        /* min proportion of type 1 particles */
#define PROP_MAX 0.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 50    /* initial time without change */
#define PROP_TFINAL 50      /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 53       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 1800
#define NY_SPHERE 1350   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

**2D part:**

```
#define DRAW_SPHERE 0    /* set to 1 to draw 3D sphere */

```

### 30 July 2026 - Sunflower grid wave protections: comparison of two different wavelengths ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** `init_wave_flat_comp(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 183             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_B 183           /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */
#define YMID 1150        /* mid point of display */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 20    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 11      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 11    /* pattern of circles, see list in global_pdes.c */
#define SYMMETRIC_CIRCLE_PATTERNS 1  /* set to 1 to have symmetric patterns in top and bottom */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 4.0    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.028            /* parameter controlling the dimensions of domain */
#define MU_B 0.028          /* parameter controlling the dimensions of domain */
#define MUB 0.028 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MDEPTH_B 10         /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 14           /* number of grid point for grid of disks */
#define NGRIDY 16           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.012    /* width of channels/wall separating lenses */
#define WALL_WIDTH_B 0.012  /* width of channels/wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -1.65  
#define ISO_XSHIFT_RIGHT 0.4
#define ISO_YSHIFT_LEFT -0.05
#define ISO_YSHIFT_RIGHT -0.05 
#define ISO_SCALE 0.85           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCIL_YMID -0.9        /* defines oscilling beam midpoint */
#define OSCILLATION_SCHEDULE 7  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.005        /* frequency of periodic excitation */
#define OMEGA_B 0.015      /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.1       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 15.625  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 2                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define NSOURCES 48         /* number of sources */
#define MAX_PULSING_TIME 10000           /* max time for adding pulses */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 2600      /* number of frames of movie */
#define NVID 14          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100   /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 50   /* number of still frames between movies */
#define END_FRAMES 300   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3

/* Color schemes */

#define COLOR_PALETTE 18      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 18    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 0.8    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE 0.0    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0    /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 2.0     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.2    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/

/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */

#define MEAN_FLUX (PLOT == P_TOTAL_ENERGY_FLUX)||(PLOT_B == P_TOTAL_ENERGY_FLUX)
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */

```

### 29 July 2026 - Increasing the force constant for charged particles on a sphere: Charge ratio 4:1 ###

**Program:** `lennardjones.c` 

**3D part:**

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 0.0
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.25 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 4.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.014 	    /* parameter controlling radius of particles */
#define MU_B 0.009          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2300      /* number of frames of movie */
#define NVID 80          /* number of iterations between images displayed on screen */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 4.0            /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 8.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 100.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1      /* set to 1 to add an electric field */
#define EFIELD 2500000.0  /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 4.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 2500.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 1   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 1.0        /* min proportion of type 1 particles */
#define PROP_MAX 0.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 50    /* initial time without change */
#define PROP_TFINAL 50      /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 53       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 1800
#define NY_SPHERE 1350   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

**2D part:**

```
#define DRAW_SPHERE 0    /* set to 1 to draw 3D sphere */

```

### 28 July 2026 - A high-pass filter based on two square grids with different spacing ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** `init_wave_flat_comp(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 183             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_B 183           /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */
#define YMID 1150        /* mid point of display */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 20    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 14      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 14    /* pattern of circles, see list in global_pdes.c */
#define SYMMETRIC_CIRCLE_PATTERNS 1  /* set to 1 to have symmetric patterns in top and bottom */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 4.0    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.028            /* parameter controlling the dimensions of domain */
#define MU_B 0.028          /* parameter controlling the dimensions of domain */
#define MUB 0.028 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MDEPTH_B 10         /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 14           /* number of grid point for grid of disks */
#define NGRIDY 16           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.012    /* width of channels/wall separating lenses */
#define WALL_WIDTH_B 0.012  /* width of channels/wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -1.65  
#define ISO_XSHIFT_RIGHT 0.4
#define ISO_YSHIFT_LEFT -0.05
#define ISO_YSHIFT_RIGHT -0.05 
#define ISO_SCALE 0.85           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCIL_YMID -0.9        /* defines oscilling beam midpoint */
#define OSCILLATION_SCHEDULE 7  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.005        /* frequency of periodic excitation */
#define OMEGA_B 0.015      /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.1       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 15.625  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 2                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define NSOURCES 48         /* number of sources */
#define MAX_PULSING_TIME 10000           /* max time for adding pulses */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 3600      /* number of frames of movie */
#define NVID 14          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100   /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 50   /* number of still frames between movies */
#define END_FRAMES 300   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3

/* Color schemes */

#define COLOR_PALETTE 11      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 0.8    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE 0.0    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0    /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 2.0     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.2    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/

/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */

#define MEAN_FLUX (PLOT == P_TOTAL_ENERGY_FLUX)||(PLOT_B == P_TOTAL_ENERGY_FLUX)
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */

```

### 27 July 2026 - Heavy cations and light anions in an increasing electric field on a sphere: Charge ratio 4:1 ###

**Program:** `lennardjones.c` 

**3D part:**

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 0.0
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.25 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 4.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.014 	    /* parameter controlling radius of particles */
#define MU_B 0.009          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2500      /* number of frames of movie */
#define NVID 80          /* number of iterations between images displayed on screen */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 8.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 100.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1      /* set to 1 to add an electric field */
#define EFIELD 10000.0  /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 4.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 1      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 2500.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 1.0        /* min proportion of type 1 particles */
#define PROP_MAX 0.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 50    /* initial time without change */
#define PROP_TFINAL 50      /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 53       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 1800
#define NY_SPHERE 1350   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

**2D part:**

```
#define DRAW_SPHERE 0    /* set to 1 to draw 3D sphere */

```

### 26 July 2026 - Randomized square grid wave protections: Comparison of two different wavelengths ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** `init_wave_flat_comp(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 183             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_B 183           /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */
#define YMID 1150        /* mid point of display */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 20    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 2      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 2    /* pattern of circles, see list in global_pdes.c */
#define SYMMETRIC_CIRCLE_PATTERNS 1  /* set to 1 to have symmetric patterns in top and bottom */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 4.0    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.028            /* parameter controlling the dimensions of domain */
#define MU_B 0.028          /* parameter controlling the dimensions of domain */
#define MUB 0.028 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MDEPTH_B 10         /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 18           /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.012    /* width of channels/wall separating lenses */
#define WALL_WIDTH_B 0.012  /* width of channels/wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -1.65  
#define ISO_XSHIFT_RIGHT 0.4
#define ISO_YSHIFT_LEFT -0.05
#define ISO_YSHIFT_RIGHT -0.05 
#define ISO_SCALE 0.85           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCIL_YMID -0.9        /* defines oscilling beam midpoint */
#define OSCILLATION_SCHEDULE 7  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.005        /* frequency of periodic excitation */
#define OMEGA_B 0.015      /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.1       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 15.625  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 2                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define NSOURCES 48         /* number of sources */
#define MAX_PULSING_TIME 10000           /* max time for adding pulses */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 3400      /* number of frames of movie */
#define NVID 14          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100   /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 50   /* number of still frames between movies */
#define END_FRAMES 300   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3

/* Color schemes */

#define COLOR_PALETTE 18      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 0.8    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE 0.0    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0    /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 2.0     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.2    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/

/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */

#define MEAN_FLUX (PLOT == P_TOTAL_ENERGY_FLUX)||(PLOT_B == P_TOTAL_ENERGY_FLUX)
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */

```

### 25 July 2026 - Heavy cations and light anions in an increasing electric field on a sphere: Charge ratio 3:1 ###

**Program:** `lennardjones.c` 

**3D part:**

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 0.0
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.25 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 4.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.014 	    /* parameter controlling radius of particles */
#define MU_B 0.009          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2700      /* number of frames of movie */
#define NVID 80          /* number of iterations between images displayed on screen */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 6.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 500.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1      /* set to 1 to add an electric field */
#define EFIELD 10000.0  /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 3.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 1      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 2500.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 1.0        /* min proportion of type 1 particles */
#define PROP_MAX 0.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 50    /* initial time without change */
#define PROP_TFINAL 50      /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 53       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 1800
#define NY_SPHERE 1350   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

**2D part:**

```
#define DRAW_SPHERE 0    /* set to 1 to draw 3D sphere */

```

### 24 July 2026 - Poisson disc wave protections: Comparison of two different small wavelengths ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** `init_wave_flat_comp(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 183             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_B 183           /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */
#define YMID 1150        /* mid point of display */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 20    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 8      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 8    /* pattern of circles, see list in global_pdes.c */
#define SYMMETRIC_CIRCLE_PATTERNS 1  /* set to 1 to have symmetric patterns in top and bottom */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 4.0    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.028            /* parameter controlling the dimensions of domain */
#define MU_B 0.028          /* parameter controlling the dimensions of domain */
#define MUB 0.028 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MDEPTH_B 10         /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 14           /* number of grid point for grid of disks */
#define NGRIDY 16           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.012    /* width of channels/wall separating lenses */
#define WALL_WIDTH_B 0.012  /* width of channels/wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -1.65  
#define ISO_XSHIFT_RIGHT 0.4
#define ISO_YSHIFT_LEFT -0.05
#define ISO_YSHIFT_RIGHT -0.05 
#define ISO_SCALE 0.85           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCIL_YMID -0.9        /* defines oscilling beam midpoint */
#define OSCILLATION_SCHEDULE 7  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define OMEGA_B 0.02       /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.1       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 15.625  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 2                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define NSOURCES 48         /* number of sources */
#define MAX_PULSING_TIME 10000           /* max time for adding pulses */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 3200      /* number of frames of movie */
#define NVID 14          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100   /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 50   /* number of still frames between movies */
#define END_FRAMES 300   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3

/* Color schemes */

#define COLOR_PALETTE 10      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.2   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0    /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 2.0     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.2    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/

/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */

#define MEAN_FLUX (PLOT == P_TOTAL_ENERGY_FLUX)||(PLOT_B == P_TOTAL_ENERGY_FLUX)
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */

```

### 23 July 2026 - Heavy cations and light anions in an increasing electric field on a sphere ###

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 0.0
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.35 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.4  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.014 	    /* parameter controlling radius of particles */
#define MU_B 0.01           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 3100      /* number of frames of movie */
#define NVID 80          /* number of iterations between images displayed on screen */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 6.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 2.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 500.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1      /* set to 1 to add an electric field */
#define EFIELD 10000.0  /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 2.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 1      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 2500.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 1.0        /* min proportion of type 1 particles */
#define PROP_MAX 0.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 50    /* initial time without change */
#define PROP_TFINAL 50      /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 53       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 1800
#define NY_SPHERE 1350   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

### 22 July 2026 - Poisson disc wave protections: Comparison of two different wavelengths ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** `init_wave_flat_comp(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 183             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_B 183           /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */
#define YMID 1150        /* mid point of display */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 20    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 8      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 8    /* pattern of circles, see list in global_pdes.c */
#define SYMMETRIC_CIRCLE_PATTERNS 1  /* set to 1 to have symmetric patterns in top and bottom */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 4.0    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.028            /* parameter controlling the dimensions of domain */
#define MU_B 0.028          /* parameter controlling the dimensions of domain */
#define MUB 0.028 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MDEPTH_B 10         /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 14           /* number of grid point for grid of disks */
#define NGRIDY 16           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.012    /* width of channels/wall separating lenses */
#define WALL_WIDTH_B 0.012  /* width of channels/wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -1.65  
#define ISO_XSHIFT_RIGHT 0.4
#define ISO_YSHIFT_LEFT -0.05
#define ISO_YSHIFT_RIGHT -0.05 
#define ISO_SCALE 0.85           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCIL_YMID -0.9        /* defines oscilling beam midpoint */
#define OSCILLATION_SCHEDULE 7  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.005        /* frequency of periodic excitation */
#define OMEGA_B 0.01       /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.1       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 15.625  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 2                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define NSOURCES 48         /* number of sources */
#define MAX_PULSING_TIME 10000           /* max time for adding pulses */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 3200      /* number of frames of movie */
#define NVID 14          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100   /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 50   /* number of still frames between movies */
#define END_FRAMES 300   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3

/* Color schemes */

#define COLOR_PALETTE 10      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.2   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0    /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 2.0     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.2    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/

/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.5  /* lower value increases sensitivity of shading */

#define MEAN_FLUX (PLOT == P_TOTAL_ENERGY_FLUX)||(PLOT_B == P_TOTAL_ENERGY_FLUX)
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */

```

### 21 July 2026 - Charged particles in an increasing electric field on a sphere ###

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 0.0
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.7 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.4  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.014 	    /* parameter controlling radius of particles */
#define MU_B 0.014          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2400      /* number of frames of movie */
#define NVID 80          /* number of iterations between images displayed on screen */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 4.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 500.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1      /* set to 1 to add an electric field */
#define EFIELD 10000.0  /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 1      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 1000.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 1.0        /* min proportion of type 1 particles */
#define PROP_MAX 0.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 50    /* initial time without change */
#define PROP_TFINAL 50      /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 53       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 1800
#define NY_SPHERE 1350   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 540.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

### 20 July 2026 - Triangle grid wave protections: Comparison of two different wavelengths ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** `init_wave_flat_comp(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 183             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_B 183           /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */
#define YMID 1150        /* mid point of display */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 20    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 1      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 1    /* pattern of circles, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.75   /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.028            /* parameter controlling the dimensions of domain */
#define MU_B 0.028          /* parameter controlling the dimensions of domain */
#define MUB 0.028 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MDEPTH_B 10         /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 14           /* number of grid point for grid of disks */
#define NGRIDY 16           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.012    /* width of channels/wall separating lenses */
#define WALL_WIDTH_B 0.012  /* width of channels/wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -1.65  
#define ISO_XSHIFT_RIGHT 0.4
#define ISO_YSHIFT_LEFT -0.05
#define ISO_YSHIFT_RIGHT -0.05 
#define ISO_SCALE 0.85           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCIL_YMID -0.9        /* defines oscilling beam midpoint */
#define OSCILLATION_SCHEDULE 7  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.005        /* frequency of periodic excitation */
#define OMEGA_B 0.01       /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.1       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 15.625  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 2                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define NSOURCES 48         /* number of sources */
#define MAX_PULSING_TIME 10000           /* max time for adding pulses */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 3500      /* number of frames of movie */
#define NVID 14          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100   /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 50   /* number of still frames between movies */
#define END_FRAMES 300   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3

/* Color schemes */

#define COLOR_PALETTE 15      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.2   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0    /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 2.0     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.2    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/

/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.5  /* lower value increases sensitivity of shading */

#define MEAN_FLUX (PLOT == P_TOTAL_ENERGY_FLUX)||(PLOT_B == P_TOTAL_ENERGY_FLUX)
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */

```

### 19 July 2026 - Changing the anion/cation proportion of particles on a sphere in an eastward electric field ###

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 0.0
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.4  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.01           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2600      /* number of frames of movie */
#define NVID 80          /* number of iterations between images displayed on screen */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 4.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 500.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1      /* set to 1 to add an electric field */
#define EFIELD 1000000.0  /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 1000.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 1.0        /* min proportion of type 1 particles */
#define PROP_MAX 0.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 50    /* initial time without change */
#define PROP_TFINAL 50      /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 53       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 1800
#define NY_SPHERE 1350   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 540.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

### 18 July 2026 - Square grid wave protections: Comparison of two different wavelengths ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** `init_wave_flat_comp(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 183             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_B 183           /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */
#define YMID 1150        /* mid point of display */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 20    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 0      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 0    /* pattern of circles, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.75   /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.028            /* parameter controlling the dimensions of domain */
#define MU_B 0.028          /* parameter controlling the dimensions of domain */
#define MUB 0.028 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MDEPTH_B 10         /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 14           /* number of grid point for grid of disks */
#define NGRIDY 16           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.012    /* width of channels/wall separating lenses */
#define WALL_WIDTH_B 0.012  /* width of channels/wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -1.65  
#define ISO_XSHIFT_RIGHT 0.4
#define ISO_YSHIFT_LEFT -0.05
#define ISO_YSHIFT_RIGHT -0.05 
#define ISO_SCALE 0.85           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCIL_YMID -0.9        /* defines oscilling beam midpoint */
#define OSCILLATION_SCHEDULE 7  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.005        /* frequency of periodic excitation */
#define OMEGA_B 0.01       /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.1       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 15.625  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 2                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define NSOURCES 48         /* number of sources */
#define MAX_PULSING_TIME 10000           /* max time for adding pulses */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 3200      /* number of frames of movie */
#define NVID 14          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100   /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 50   /* number of still frames between movies */
#define END_FRAMES 300   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3

/* Color schemes */

#define COLOR_PALETTE 11      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.2   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0    /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 2.0     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.2    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/

/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.5  /* lower value increases sensitivity of shading */

#define MEAN_FLUX (PLOT == P_TOTAL_ENERGY_FLUX)||(PLOT_B == P_TOTAL_ENERGY_FLUX)
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */

```

### 17 July 2026 - Changing the cation/anion proportion of particles starting in the center of a maze on the sphere in an electric field ###

**Program:** `lennardjones.c` 

**3D part:**

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 1.37
#define INITYMAX 1.77	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 4.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.01           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2550      /* number of frames of movie */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 4.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 500.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1      /* set to 1 to add an electric field */
#define EFIELD 0.0    /* value of electric field */
#define EFIELD_Y -1500000.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 1000.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 1.0        /* min proportion of type 1 particles */
#define PROP_MAX 0.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 50    /* initial time without change */
#define PROP_TFINAL 50      /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 52       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 1800
#define NY_SPHERE 1350   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 720.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

**2D part:**

```
#define DRAW_SPHERE 0    /* set to 1 to draw 3D sphere */

```


### 16 July 2026 - Comparison of square grid and staggered square grid wave protections with Neumann boundary conditions ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** `init_wave_flat_comp(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 183             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_B 183           /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */
#define YMID 1150        /* mid point of display */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 20    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 0      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 1    /* pattern of circles, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.75   /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.028            /* parameter controlling the dimensions of domain */
#define MU_B 0.028          /* parameter controlling the dimensions of domain */
#define MUB 0.028 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MDEPTH_B 10         /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 14           /* number of grid point for grid of disks */
#define NGRIDY 16           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.012    /* width of channels/wall separating lenses */
#define WALL_WIDTH_B 0.012  /* width of channels/wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -1.65  
#define ISO_XSHIFT_RIGHT 0.4
#define ISO_YSHIFT_LEFT -0.05
#define ISO_YSHIFT_RIGHT -0.05 
#define ISO_SCALE 0.85           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCIL_YMID -0.9        /* defines oscilling beam midpoint */

#define OMEGA 0.007        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.1       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 15.625  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 2                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define NSOURCES 48         /* number of sources */
#define MAX_PULSING_TIME 10000           /* max time for adding pulses */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 2700      /* number of frames of movie */
#define NVID 14          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100   /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 50   /* number of still frames between movies */
#define END_FRAMES 300   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3

/* Color schemes */

#define COLOR_PALETTE 11      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.2   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0    /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 2.0     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.2    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define OSCILLATION_SCHEDULE 3  /* oscillation schedule, see list in global_pdes.c */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/

/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.5  /* lower value increases sensitivity of shading */

#define MEAN_FLUX (PLOT == P_TOTAL_ENERGY_FLUX)||(PLOT_B == P_TOTAL_ENERGY_FLUX)
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */

```

### 15 July 2026 - Changing the cation/anion proportion in a maze 
on the sphere in an electric field (short) ###

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

// #define WINWIDTH 	1760  /* window width */
#define WINWIDTH 	990  /* window width */
#define WINHEIGHT 	990   /* window height */

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 2.7
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 1.4  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.01           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2000      /* number of frames of movie */
#define NVID 80          /* number of iterations between images displayed on screen */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 4.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 500.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1      /* set to 1 to add an electric field */
#define EFIELD 0.0    /* value of electric field */
#define EFIELD_Y -1500000.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 1000.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 1.0        /* min proportion of type 1 particles */
#define PROP_MAX 0.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 50    /* initial time without change */
#define PROP_TFINAL 50      /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 53       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 1800
#define NY_SPHERE 1350   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 540.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

### 14 July 2026 - Comparison of randomized square grid and Poisson disc wave protections with Neumann boundary conditions ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** `init_wave_flat_comp(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 183             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_B 183           /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */
#define YMID 1150        /* mid point of display */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 20    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 2      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 8    /* pattern of circles, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.75   /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.028            /* parameter controlling the dimensions of domain */
#define MU_B 0.028          /* parameter controlling the dimensions of domain */
#define MUB 0.03 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MDEPTH_B 10         /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 14           /* number of grid point for grid of disks */
#define NGRIDY 16           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.012    /* width of channels/wall separating lenses */
#define WALL_WIDTH_B 0.012  /* width of channels/wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -1.65  
#define ISO_XSHIFT_RIGHT 0.4
#define ISO_YSHIFT_LEFT -0.05
#define ISO_YSHIFT_RIGHT -0.05 
#define ISO_SCALE 0.85           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCIL_YMID -0.9        /* defines oscilling beam midpoint */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.1       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 15.625  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 2                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define NSOURCES 48         /* number of sources */
#define MAX_PULSING_TIME 10000           /* max time for adding pulses */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 2800      /* number of frames of movie */
#define NVID 14          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100   /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 50   /* number of still frames between movies */
#define END_FRAMES 300   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3

/* Color schemes */

#define COLOR_PALETTE 18      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 0.9    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.0   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0    /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 2.0     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.2    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define OSCILLATION_SCHEDULE 3  /* oscillation schedule, see list in global_pdes.c */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/

/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.5  /* lower value increases sensitivity of shading */

#define MEAN_FLUX (PLOT == P_TOTAL_ENERGY_FLUX)||(PLOT_B == P_TOTAL_ENERGY_FLUX)
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */

```

### 13 July 2026 - Changing the cation/anion proportion of in a maze on the sphere in an electric field ###

**Program:** `lennardjones.c` 

**3D part:**

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 2.7
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 1.4  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.01           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2750      /* number of frames of movie */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 4.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 500.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1      /* set to 1 to add an electric field */
#define EFIELD 0.0    /* value of electric field */
#define EFIELD_Y -500000.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 1000.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 1.0        /* min proportion of type 1 particles */
#define PROP_MAX 0.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 50    /* initial time without change */
#define PROP_TFINAL 50      /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 51       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 2072
#define NY_SPHERE 1536   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 720.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

**2D part:**

```
#define DRAW_SPHERE 0    /* set to 1 to draw 3D sphere */

```

### 12 July 2026 - Comparison of hex grid and Poisson disc wave protections with Neumann boundary conditions ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** `init_wave_flat_comp(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 183             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_B 183           /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */
#define YMID 1150        /* mid point of display */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 20    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 1      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 8    /* pattern of circles, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.75   /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.028            /* parameter controlling the dimensions of domain */
#define MU_B 0.028          /* parameter controlling the dimensions of domain */
#define MUB 0.03 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MDEPTH_B 10         /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 14           /* number of grid point for grid of disks */
#define NGRIDY 16           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.012    /* width of channels/wall separating lenses */
#define WALL_WIDTH_B 0.012  /* width of channels/wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -1.65  
#define ISO_XSHIFT_RIGHT 0.4
#define ISO_YSHIFT_LEFT -0.05
#define ISO_YSHIFT_RIGHT -0.05 
#define ISO_SCALE 0.85           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCIL_YMID -0.9        /* defines oscilling beam midpoint */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.1       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 15.625  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 2                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define NSOURCES 48         /* number of sources */
#define MAX_PULSING_TIME 10000           /* max time for adding pulses */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 2900      /* number of frames of movie */
#define NVID 14          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100   /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 50   /* number of still frames between movies */
#define END_FRAMES 300   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3

/* Color schemes */

#define COLOR_PALETTE 15      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.2   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0    /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 2.0     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.2    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define OSCILLATION_SCHEDULE 3  /* oscillation schedule, see list in global_pdes.c */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/

/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.5  /* lower value increases sensitivity of shading */

#define MEAN_FLUX (PLOT == P_TOTAL_ENERGY_FLUX)||(PLOT_B == P_TOTAL_ENERGY_FLUX)
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */

```

### 11 July 2026 - Changing the cation/anion proportion of charged particles in a maze on the sphere ###

**Program:** `lennardjones.c` 

**3D part:**

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 2.7
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 1.4  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.01           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 3200      /* number of frames of movie */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 4.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 500.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1      /* set to 1 to add an electric field */
#define EFIELD 0.0    /* value of electric field */
#define EFIELD_Y 20000.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 1000.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 1.0        /* min proportion of type 1 particles */
#define PROP_MAX 0.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 50    /* initial time without change */
#define PROP_TFINAL 50      /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 50       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 2072
#define NY_SPHERE 1536   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 720.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

**2D part:**

```
#define DRAW_SPHERE 0    /* set to 1 to draw 3D sphere */

```

### 10 July 2026 - Comparison of regular and random square grid wave protections with Neumann boundary conditions ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** `init_wave_flat_comp(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 183             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_B 183           /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */
#define YMID 1150        /* mid point of display */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 20    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 0      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 2    /* pattern of circles, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.25   /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.03             /* parameter controlling the dimensions of domain */
#define MU_B 0.03           /* parameter controlling the dimensions of domain */
#define MUB 0.03 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MDEPTH_B 10         /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 16           /* number of grid point for grid of disks */
#define NGRIDY 16           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.012    /* width of channels/wall separating lenses */
#define WALL_WIDTH_B 0.012  /* width of channels/wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -1.65  
#define ISO_XSHIFT_RIGHT 0.4
#define ISO_YSHIFT_LEFT -0.05
#define ISO_YSHIFT_RIGHT -0.05 
#define ISO_SCALE 0.85           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCIL_YMID -0.9        /* defines oscilling beam midpoint */

#define OMEGA 0.012        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.1       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 15.625  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 2                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define NSOURCES 48         /* number of sources */
#define MAX_PULSING_TIME 10000           /* max time for adding pulses */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 2700      /* number of frames of movie */
#define NVID 14          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100   /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 50   /* number of still frames between movies */
#define END_FRAMES 300   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3

/* Color schemes */

#define COLOR_PALETTE 10      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.2   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0    /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 2.0     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.2    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define OSCILLATION_SCHEDULE 3  /* oscillation schedule, see list in global_pdes.c */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/

/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.5  /* lower value increases sensitivity of shading */

#define MEAN_FLUX (PLOT == P_TOTAL_ENERGY_FLUX)||(PLOT_B == P_TOTAL_ENERGY_FLUX)
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */

```

### 09 July 2026 - A 60%-40% mixture of cations and anions in a mazeon a sphere ###

**Program:** `lennardjones.c` 

**3D part:**

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 2.7
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.6 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 1.4  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.01           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2700      /* number of frames of movie */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 4.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 500.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 20000.0    /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 1000.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 0.0        /* min proportion of type 1 particles */
#define PROP_MAX 1.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 250   /* initial time without change */
#define PROP_TFINAL 250     /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 30       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 2072
#define NY_SPHERE 1536   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 720.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

**2D part:**

```
#define DRAW_SPHERE 0    /* set to 1 to draw 3D sphere */

```

### 08 July 2026 - Comparison without interaction of Dirichlet and Neumann boundary conditions for waves crossing a square grid ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** `init_wave_flat_comp(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 183             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_B 183           /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */
#define YMID 1150        /* mid point of display */

// #define XMIN -2.2
// #define XMAX 3.8	/* x interval */
// #define YMIN -1.796875 
// #define YMAX 1.796875	/* y interval for 9/16 aspect ratio */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 202      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 202    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 0      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 0    /* pattern of circles, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.25   /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.045            /* parameter controlling the dimensions of domain */
#define MU_B 0.045          /* parameter controlling the dimensions of domain */
#define MUB 0.045 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MDEPTH_B 10         /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.012    /* width of channels/wall separating lenses */
#define WALL_WIDTH_B 0.012  /* width of channels/wall separating lenses */
#define WALL_WIDTH_RND 0.0  /* proportion of width of width for random arrangements */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define WALL_WIDTH_ASYM 0.75      /* asymmetry of wall width (D_CIRCLE_LATTICE_NONISO) */
#define WALL_WIDTH_ASYM_B 0.75    /* asymmetry of wall width (D_CIRCLE_LATTICE_HEX_NONISO) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

#define ISO_XSHIFT_LEFT -1.65  
#define ISO_XSHIFT_RIGHT 0.4
#define ISO_YSHIFT_LEFT -0.05
#define ISO_YSHIFT_RIGHT -0.05 
#define ISO_SCALE 0.85           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCIL_YMID -0.9        /* defines oscilling beam midpoint */

#define OMEGA 0.006        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.1       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 15.625  /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 2                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define NSOURCES 48         /* number of sources */
#define MAX_PULSING_TIME 10000           /* max time for adding pulses */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 2700      /* number of frames of movie */
#define NVID 14          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100   /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 50   /* number of still frames between movies */
#define END_FRAMES 300   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75             /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 3

/* Color schemes */

#define COLOR_PALETTE 11      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.1   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 1.0    /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 2.0     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.2    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define OSCILLATION_SCHEDULE 3  /* oscillation schedule, see list in global_pdes.c */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/

/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 6     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.5  /* lower value increases sensitivity of shading */

#define MEAN_FLUX (PLOT == P_TOTAL_ENERGY_FLUX)||(PLOT_B == P_TOTAL_ENERGY_FLUX)
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */

```

### 07 July 2026 - A 75%-25% mixture of cations and anions in a maze on a sphere ###

**Program:** `lennardjones.c` 

**3D part:**

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 2.7
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.75 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 1.4  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.01           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2500      /* number of frames of movie */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 4.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 500.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 20000.0    /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 1000.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 0.0        /* min proportion of type 1 particles */
#define PROP_MAX 1.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 250   /* initial time without change */
#define PROP_TFINAL 250     /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 2        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 2072
#define NY_SPHERE 1536   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 540.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

**2D part:**

```
#define DRAW_SPHERE 0    /* set to 1 to draw 3D sphere */

```

### 06 July 2026 - Comparison of Dirichlet and Neumann boundary conditions for waves crossing a square grid ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
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

#define B_DOMAIN 202         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 0   /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 5.0    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.045            /* parameter controlling the dimensions of domain */
#define MU_B 0.45           /* parameter controlling the dimensions of domain */
#define NPOLY 9             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.075    /* width of wall separating lenses */
#define WALL_WIDTH_B 0.1    /* width of wall separating lenses */
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
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 1.2          /* defines oscilling beam range */
#define OSCIL_YMID 0.0          /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.006        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.09        /* Courant number in medium B */
#define COURANTB 0.25      /* Courant number */
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
#define OSCILLATING_SOURCE_PERIOD 7.0   /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 13           /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 3200          /* number of frames of movie */
#define NVID 14              /* number of iterations between images displayed on screen */
#define NSEG 1000            /* number of segments of boundary */
#define INITIAL_TIME 100      /* time after which to start saving frames */
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

#define INITIAL_AMP 1.0              /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0003    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.01   /* wavelength of initial condition */

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
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.7    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 2.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.5  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.7      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 4.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
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

### 05 July 2026 - Positively charged particles in a maze on a sphere ###

**Program:** `lennardjones.c` 

**3D part:**

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 2.7
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */
#define BELT_SPEED1 10.0     /* speed of first conveyor belt */
#define BELT_SPEED2 15.0   /* speed of second conveyor belt */
#define BELT_SPEED3 6.0   /* speed of second conveyor belt */
#define OBSTACLE_OMEGA 300.0  /* obstacle rotation speed */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 1.4  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.01           /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2900      /* number of frames of movie */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 17
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 4.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 500.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 20000.0    /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 1.0        /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 1000.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 0.0        /* min proportion of type 1 particles */
#define PROP_MAX 1.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 250   /* initial time without change */
#define PROP_TFINAL 250     /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
#define RAND_SHIFT 1        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100     /* size of hashgrid in x direction */
#define HASHY 50      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 2072
#define NY_SPHERE 1536   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 540.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

**2D part:**

```
#define DRAW_SPHERE 0    /* set to 1 to draw 3D sphere */

```

### 04 July 2026 - Plane waves crossing a "sunflower" grid of obstacles with Neumann boundary conditions ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
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

#define B_DOMAIN 20         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 11   /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 5.0    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.035            /* parameter controlling the dimensions of domain */
#define MU_B 0.45           /* parameter controlling the dimensions of domain */
#define NPOLY 9             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.075    /* width of wall separating lenses */
#define WALL_WIDTH_B 0.1    /* width of wall separating lenses */
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
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 1.2          /* defines oscilling beam range */
#define OSCIL_YMID 0.0          /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.006        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.09        /* Courant number in medium B */
#define COURANTB 0.25      /* Courant number */
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
#define OSCILLATING_SOURCE_PERIOD 7.0   /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 13           /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 3500          /* number of frames of movie */
#define NVID 14              /* number of iterations between images displayed on screen */
#define NSEG 1000            /* number of segments of boundary */
#define INITIAL_TIME 100      /* time after which to start saving frames */
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

#define INITIAL_AMP 1.0              /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0003    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.01   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 18     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 17   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.7    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 2.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.1  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.7      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 4.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
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

### 03 July 2026 - Particles in a maze on a rotating sphere ###

**Program:** `lennardjones.c` 

**3D part:**

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 2.7
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
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

#define INTERACTION 1        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 1.6  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.022          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 2700      /* number of frames of movie */
#define NVID 75          /* number of iterations between images displayed on screen */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 3          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 4.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 500.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0002           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 15000.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 20000.0    /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 0.0        /* charge of particles of first type */
#define CHARGE_B 0.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 1000.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 1     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 10.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 1   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 5.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 0.0        /* min proportion of type 1 particles */
#define PROP_MAX 1.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 250   /* initial time without change */
#define PROP_TFINAL 250     /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
#define RAND_SHIFT 21        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 120     /* size of hashgrid in x direction */
#define HASHY 60      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 2072
#define NY_SPHERE 1536   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 540.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

**2D part:**

```
#define DRAW_SPHERE 0    /* set to 1 to draw 3D sphere */

```

### 02 July 2026 - Plane waves crossing a randomized square grid of obstacles with Neumann boundary conditions ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
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

#define B_DOMAIN 20         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 2    /* pattern of circles or polygons, see list in global_pdes.c */
#define IMAGE_FILE 5        /* for option D_IMAGE */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 5.0    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */
#define PDISC_CONNECT_FACTOR 1.5    /* controls which discs are connected for D_CIRCLE_LATTICE_POISSON domain */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.03             /* parameter controlling the dimensions of domain */
#define MU_B 0.45           /* parameter controlling the dimensions of domain */
#define NPOLY 9             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.075    /* width of wall separating lenses */
#define WALL_WIDTH_B 0.1    /* width of wall separating lenses */
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
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 1.2          /* defines oscilling beam range */
#define OSCIL_YMID 0.0          /* defines oscilling beam midpoint */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.006        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.09        /* Courant number in medium B */
#define COURANTB 0.25      /* Courant number */
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
#define OSCILLATING_SOURCE_PERIOD 7.0   /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */
#define MAX_PULSING_TIME 13           /* max time for adding pulses */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3
#define BC_NEUMANN 1        /* set to 1 to use Neumann boundary conditions on domain */

/* Parameters for length and speed of simulation */

#define NSTEPS 3300          /* number of frames of movie */
#define NVID 14              /* number of iterations between images displayed on screen */
#define NSEG 1000            /* number of segments of boundary */
#define INITIAL_TIME 100      /* time after which to start saving frames */
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

#define INITIAL_AMP 1.0              /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0003    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.01   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 15     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 17   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define COLOR_RANGE 1.0    /* max range of color (default: 1.0) */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE -0.7    /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 2.0     /* additional scaling factor for wave amplitude */
#define E_SCALE 300.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 1       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.25  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.7      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 4.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.75     /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.0      /* value of y to sample wave profile */
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

### 01 July 2026 - Fewer particles in an 18 by 8 maze on a sphere ###

**Program:** `lennardjones.c` 

**3D part:**

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

#define XMIN 0.0
#define XMAX 6.283185307	/* x interval */
#define YMIN 0.0
#define YMAX 3.141592654	/* y interval for 9/16 aspect ratio */

#define INITXMIN 0.1
#define INITXMAX 6.18	/* x interval for initial condition */
#define INITYMIN 2.7
#define INITYMAX 3.14	/* y interval for initial condition */

#define THERMOXMIN -1.25
#define THERMOXMAX 1.25	/* x interval for initial condition */
#define THERMOYMIN 0.0
#define THERMOYMAX 0.75	/* y interval for initial condition */

#define ADDXMIN 0.1
#define ADDXMAX 0.2	/* x interval for adding particles */
#define ADDYMIN 1.57
#define ADDYMAX 1.57	/* y interval for adding particles */
#define ADDRMIN 2.0 
#define ADDRMAX 2.1     /* r interval for adding particles */

#define BCXMIN 0.0
#define BCXMAX 6.283185307	/* x interval for boundary condition */
#define BCYMIN 0.3
#define BCYMAX 2.841592654	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */
#define OBSYMIN -1.125
#define OBSYMAX 1.125     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 81  /* pattern of circles, see list in global_ljones.c */

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

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 153    /* pattern of repelling segments, see list in global_ljones.c */
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

#define INTERACTION 1        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1      /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 5.0 /* angular frequency of spin-spin interaction for second particle type */
#define MOL_ANGLE_FACTOR 1.0    /* rotation angle for P_MOL_ANGLE color scheme */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.022          /* parameter controlling radius of particles of second type */
#define MU_ADD 0.022        /* parameter controlling radius of added particles */
#define MU_ADD_B 0.022      /* parameter controlling radius of added particles */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 24
#define NOBSY 14           /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */
#define OFSSET_TREES 0.5    /* vertical offset in S_TREES_B */
#define SLOPE_TREES 0.015   /* slope in S_TREES_B (default: 0.3) */
#define SLOPE_TREES_B 0.015   /* slope in S_TREES_B (default: 0.25) */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */
 
#define NSTEPS 3100      /* number of frames of movie */
#define NVID 75          /* number of iterations between images displayed on screen */
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
#define END_FRAMES 200   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 30

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 3          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 3        /* type of background coloring, see list in global_ljones.c */
#define OBSTACLE_COLOR 0    /* type of obstacle, see OC_ in global_ljones.c */
#define SHADE_BG_COLOR_2D 1 /* set to 1 to shade BG color, for option BG_POTENTIAL */
#define SHADE_SCALE_BG_2D 0.1   /* controls 2D shading */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define DRAW_CLUSTER_LINKS 0    /* set to 1 to draw links between particles in cluster */
#define DRAW_OBSTACLE_LINKS 0   /* set to 1 to draw links between interacting obstacles */
#define FILL_OBSTACLE_TRIANGLES 0   /* set to 1 to fill triangles between interacting obstacles */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 5   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 1     /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.999          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */
#define OBSTACLE_AREA_SHADE_FACTOR 75.0     /* controls sensitivity of triangle shade for option FILL_OBSTACLE_TRIANGLES */
#define SHADE_OBSTACLE_FACETS 1     /* set to 1 to shade facets instead of triangles */

/* Color schemes */

#define COLOR_PALETTE 10             /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10        /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 17       /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 10    /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 10 /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11    /* Color palette for different neighbours representation */
#define COLOR_PALETTE_CHARGE 18      /* Color palette for charge representation */
#define COLOR_PALETTE_CLUSTER 14     /* Color palette for cluster representation */
#define COLOR_PALETTE_CLUSTER_SIZE 13 /* Color palette for cluster size representation */
#define COLOR_PALETTE_CLUSTER_SELECTED 11 /* Color palette for selected cluster representation */
#define COLOR_PALETTE_ANGULAR_MOMENTUM 17   /* Color palette for angular momentum */
#define COLOR_PALETTE_CURRENT 17      /* Color palette for current */
#define COLOR_PALETTE_POTENTIAL 11    /* Color palette for direction representation */
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
#define PRINT_NABSORBED 0           /* print number of absorbed particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 350.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMIN 0.0           /* energy of particle with coolest color */
#define PARTICLE_EMAX 100000.0        /* energy of particle with hottest color */
#define PARTICLE_DMIN 200.0         /* energy of particle with largest local density */
#define PARTICLE_DMAX 500.0         /* energy of particle with largest local density */
#define SEGMENT_HUE_MIN 275.0       /* color of original segment */
#define SEGMENT_HUE_MAX 30.0        /* color of saturated segment */
#define OBSTACLE_EMAX 1000000.0         /* energy of obstacle with hottest color */
#define OBSTACLE_VMAX 4.0           /* speed of obstacle with largest luminosity */
#define HUE_TYPE0 280.0      /* hue of particles of type 0 */
#define HUE_TYPE1 70.0       /* hue of particles of type 1 */
#define HUE_TYPE2 100.0      /* hue of particles of type 2 */
#define HUE_TYPE3 140.0      /* hue of particles of type 3 */
#define HUE_TYPE4 180.0       /* hue of particles of type 4 */
#define HUE_TYPE5 220.0       /* hue of particles of type 5 */
#define HUE_TYPE6 260.0      /* hue of particles of type 6 */
#define HUE_TYPE7 300.0      /* hue of particles of type 7 */
#define HUE_TYPE8 330.0      /* hue of particles of type 7 */
#define BG_LOG_EKIN_SHIFT 1.0    /* constant in BG_LOG_EKIN background color scheme */
#define BG_FORCE_SLOPE 1.0e-6    /* constant in BG_FORCE backgound color scheme */
#define BG_CHARGE_SLOPE 1.0     /* constant in BG_CHARGE backgound color scheme (default: 0.5) */
#define BG_POTENTIAL_SLOPE 50.0  /* constant in BG_POTENTIAL background color scheme */
#define CHARGE_HUE_RANGE 0.5    /* range of charge colors */
#define PARTICLE_LMAX 1.5e4     /* angular momentum particle with brightest color */

#define RANDOM_RADIUS 0          /* set to 1 for random particle radius */
#define RANDOM_RADIUS_MIN 0.4    /* min of random particle radius (default 0.75) */
#define RANDOM_RADIUS_RANGE 1.0  /* range of random particle radius (default 0.5) */
#define ADAPT_MASS_TO_RADIUS 0   /* set to positive value to for mass prop to power of radius */
#define ADAPT_DAMPING_TO_RADIUS 0.0   /* set to positive value to for friction prop to power of radius */
#define ADAPT_DAMPING_FACTOR 0.0    /* factor by which damping is adapted to radius */
#define DT_PARTICLE 2.0e-6    /* time step for particle displacement */
#define KREPEL 40.0           /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define SEGMENT_FORCE_EQR 1.0   /* equilibrium distance factor for force from segments (default 1.5) */
#define REPEL_RADIUS 25.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 0.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 0.0      /* damping coefficient for rotation of particles */
#define DAMPING_PAIRS 0.0    /* damping between paired particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 4.0     /* mass of particle of radius MU_B */
#define PARTICLE_ADD_MASS 2.0   /* mass of added particles */
#define PARTICLE_ADD_MASS_B 1.0   /* mass of added particles */
#define PARTICLE_INERTIA_MOMENT 0.1     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.1     /* moment of inertia of second type of particle */
#define V_INITIAL 50.0        /* initial velocity range */
#define V_INITIAL_ADD 4500.0        /* initial velocity range for added particles */
#define OMEGA_INITIAL 100.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */
#define COULOMB_LJ_FACTOR 1.0   /* relative intensity of LJ interaction in I_COULOMB_LJ interaction (default: 0.01) */
#define KCOULOMB_FACTOR 500.0  /* relative intensity of Coulomb interaction in I_COULOMB_LJ (default: 100.0) */
#define COULOMB_ALWAYS_REPEL 1  /* set to 1 to always repel with I_COULOMB_IMAGINARY */
#define OBSTACLE_DAMPING 0.0   /* damping of oscillating obstacles */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001           /* initial inverse temperature */
#define MU_XI 0.005           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 6.0        /* radius in which to count neighbours */
#define BOND_DIST_FACTOR 6.0       /* radius in which to draw bonds */
#define GRAVITY 15000.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define SPHERE_GRAVITY 0       /* set to 1 to have gravity at constant angle wrt sphere */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 1     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 2000.0     /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 100    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */
#define GRAVITY_INITIAL_ANGLE 0.0   /* initial angle for SPHERE_GRAVITY */
#define GRAVITY_DELTA_ANGLE 1440.0   /* increase of angle for SPHERE_GRAVITY */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0      /* set to 1 to add an electric field */
#define EFIELD 20000.0    /* value of electric field */
#define EFIELD_Y 0.0      /* value of electric field */
#define ADD_BFIELD 0      /* set to 1 to add a magnetic field */
#define BFIELD 20000.0       /* value of magnetic field */
#define CHARGE 0.0        /* charge of particles of first type */
#define CHARGE_B 0.0     /* charge of particles of second type */
#define CHARGE_ADD 0.0   /* charge of added particles */
#define CHARGE_ADD_B 0.0   /* charge of added particles */
#define INCREASE_E 0      /* set to 1 to increase electric field */
#define OSCILLATE_E 0     /* set to 1 for oscillating electric field */
#define E_PERIOD 1000      /* period of oscillating electric field */
#define EFIELD_FACTOR 1000.0    /* factor by which to increase electric field */
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

#define ROTATE_SPHERE 0     /* set to 1 to add Coriolis and centripetal force */
#define OMEGA_SPHERE 50.0    /* angular frequency of rotating sphere */
#define CHANGE_OMEGA_SPHERE 0   /* set to 1 to change sphere rotation frequency */
#define OMEGA_SPHERE_FACTOR 4.0    /* change factor of sphere rotation frequency */

#define ROTATION 0          /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 0.25  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 2.0e3         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e5  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 500.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define DRAW_MINUS 1          /* set to 1 to draw cross on particles of negative charge */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_SCHEDULE 3    /* type of temperature schedule, see TS_* in global_ljones */
#define BETA_FACTOR 0.002    /* factor by which to change BETA during simulation */
#define TS_SLOPE 8.5          /* controls speed of change of BETA for TS_TANH schedule (default 1.0) */
#define N_TOSCILLATIONS 1.0   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define INITIAL_CONSTANT_PHASE 200 /* initial phase in which temperature is constant */
#define MIDDLE_CONSTANT_PHASE 0   /* middle phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE 400     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 2     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.3    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_RIN 0.5      /* initial radius of region without coupling */
#define PARTIAL_THERMO_RFIN 1.3     /* final radius of region without coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 100.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_REGION 0      /* shape of add regions, cf ADD_* in global_ljones */
#define ADD_TIME 20        /* time at which to add first particle */
#define ADD_PERIOD 10000      /* time interval between adding further particles */
#define ADD_TYPE 1         /* type of added particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 100  /* final period where no particles are added */
#define SAFETY_FACTOR 10.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */
#define ADD_ALTERNATE_CHARGE 0   /* set to 1 to randomly select sign of added charge */
#define TIME_DEPENDENT_ADD_CHARGE 0     /* set to 1 to have added charge depend on time */
#define ALTERNATE_CHARGE_PROPORTION 0.5    /* proportion of particles of opposite charge */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 6000    /* number of tracer particles */
#define TRACER_STEPS 5           /* number of tracer steps recorded between images */
#define TRAJECTORY_LENGTH 4000    /* length of recorded trajectory */
#define TRAJECTORY_DRAW_LENGTH 1000 /* length of drawn trajectory */
#define TRACER_LUM_FACTOR 100.0    /* controls luminosity decrease of trajectories with time */
#define TRACER_PARTICLE_MASS 1.0  /* relative mass of tracer particle */
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

#define POSITION_DEPENDENT_TYPE 0   /* set to PDIC_* to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define POSITION_DEP_SIGN -1.0      /* sign in position dependence condition */
#define POSITION_DEP_X 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_Y 1.5         /* threshold value for position-dependent type */
#define POSITION_DEP_MASS_RATIO 5.0    /* position-dependent mass factor */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define SPECIAL_IC 0              /* set to 1 for choosing special initial condition RD_INITIAL_COND */
#define REACTION_DIFFUSION 0     /* set to 1 to simulate a chemical reaction (particles may change type) */
#define REACTION_MAX_TIME 100000     /* time after which no reactions take place */  
#define RD_REACTION 22             /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 8                /* number of types in reaction-diffusion equation */
#define RD_PLOT_TYPES 8           /* number of types shown in graph */
#define RD_INITIAL_COND 2         /* initial condition of particles */
#define REACTION_DIST 2.8         /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0         /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.0     /* probability controlling dissociation reaction */ 
#define KILLING_PROB 0.0015       /* probability of enzymes being killed */
#define DELTAMAX 0.1              /* max orientation difference for pairing polygons */
#define CENTER_COLLIDED_PARTICLES 1  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 0            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN -2.0e3       /* change of kinetic energy in reaction */
#define CORRECT_EQUILIBRIUM_POSITION 1  /* set to 1 to nudge particle dist towards eq dist */
#define NUDGE_FACTOR 0.0005      /* factor by which to correct particle distance */
#define COLLISION_TIME 35       /* time during which collisions are shown */
#define COLLISION_RADIUS 3.0    /* radius of discs showing collisions, in units of MU */
#define DELTAVMAX 500.0         /* maximal deltav allowed for pairing molecules */
#define AGREGMAX 3              /* maximal number of partners for CHEM_AGGREGATION reaction */
#define AGREG_DECOUPLE 12       /* minimal number of partners to decouple from thermostat */
#define NEUTRALIZE_REACTING_PARTICLES 1     /* set to 1 for reacting particles to become neutral */
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
#define PROP_MIN 0.0        /* min proportion of type 1 particles */
#define PROP_MAX 1.0        /* max proportion of type 1 particles */
#define PROP_TINITIAL 250   /* initial time without change */
#define PROP_TFINAL 250     /* final time without change */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define RANDOMIZE_ANGLE 0   /* set to 1 for random orientation */
#define DEACIVATE_CLOSE_PAIRS 0 /* set to 1 to test for closeness to other particles */
#define PAIR_SAFETY_FACTOR 1.2  /* distance to deactivate divided by sum of radii */
#define THIRD_TYPE_PROPORTION 1.0   /* proportion of third type pairings, for certain pairing types */

#define KSPRING_PAIRS 5.0e9    /* spring constant for pair interaction */
#define KTORQUE_PAIRS 1.0e10   /* constant for angular coupling in pair interaction */
#define KTORQUE_PAIR_ANGLE 0.0    /* constant for coupling between orientation in pairs */
#define NPARTNERS 2         /* number of partners of particles - for DNA, set NPARTNERS_DNA */
#define NPARTNERS_DNA 8     /* number of partners of particles, case of DNA, should be at least 8 */
#define NARMS 4             /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE 99      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between ions for POLY_WATER case */
#define PAIR_DRATIO 1.0      /* ratio between equilibrium distance and radius (default: 1.0) */
#define MU_C 0.022            /* radius of partner particle */
#define PARTICLE_MASS_C 1.0  /* mass or partner particle */
#define CHARGE_C -1.0         /* charge of partner particle */
#define CLUSTER_COLOR_FACTOR 40  /* factor for initialization of cluster colors */
#define ALTERNATE_POLY_CHARGE 0  /* set to 1 for alternating charges in molecule */
#define SECONDARY_PAIRING 0     /* set to 1 to pair with secondary partners, experimental */
#define DNA_RIGIDITY 0.5     /* controls rigidity for POLY_DNA_DOUBLE pairs, default = 1 */

#define PAIR_TYPEB_PARTICLES 0  /* set to 1 to pair particle of type 1 */
#define NPARTNERS_B 18         /* number of partners of particles */
#define NARMS_B 1               /* number of "arms" for certain paring types */ 
#define PAIRING_TYPE_B 5     /* type of pairing, see POLY_ in global_ljones.c */
#define MU_D 0.022            /* radius of partner particle */
#define PARTICLE_MASS_D 1.0  /* mass or partner particle */
#define CHARGE_D -1.0         /* charge of partner particle */

#define ADD_ABSORBERS 0     /* set to 1 to add absorbing discs */
#define ABSORBER_PATTERN 3  /* pattern of absorbers, see AP_* in global_ljones */
#define ABSORBER_X 0.0
#define ABSORBER_Y 0.0      /* coordinates of first absorber */
#define ABSORBER_R 0.015     /* radius of absorber */
#define ABSORBER_PDIST 0.4  /* parameter of Poisson disc process */

#define ADD_POTENTIAL_SPHERE 0  /* add potential for gradient force field on sphere */
#define DRAW_POTENTIAL_SPHERE 1 /* draw sphere radius depending on potential */
#define SPHERE_POTENTIAL 2      /* type of sphere potential */
#define SPHERE_POT_PATTERN 3    /* pattern of local minma of SPP_WELLS sphere potential */
#define PLANET_DEM 4            /* planet DEM used for SPP_PLANET */
#define POT_SPHERE_AMP 1.0      /* amplitude in definition of potential on sphere */
#define POT_SPHERE_RADIUS 0.1   /* radius in definition of potential on sphere */
#define POT_SPHERE_SMOOTH 0.5   /* smoothing of potential well */
#define POT_SPHERE_STRENGTH 2.5e4    /* coefficient of gradient force */

#define NXMAZE 18     /* width of maze */
#define NYMAZE 8      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
#define RAND_SHIFT 7        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.015    /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 120     /* size of hashgrid in x direction */
#define HASHY 60      /* size of hashgrid in y direction */
#define HASHMAX 100   /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

/* constants related to evolution on a sphere */
#define SPHERE 1        /* set to 1 to compute evolution in spherical geometry */
#define SIN_THETA_REG 0.01   /* regularization of sin(theta) for motion on sphere */
#define POLAR_PADDING 0.01   /* region around poles that belong to the same hashcell */
#define DRAW_SPHERE 1    /* set to 1 to draw 3D sphere */
#define DRAW_ELLIPSES_ON_SPHERE 1   /* set to 1 to draw ellipses instead of circles on sphere in 2D */
#define NX_SPHERE 2072
#define NY_SPHERE 1536   /* number of points on sphere */
#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define FLIPX -1.0             /* set to -1 to flip left/right */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D -0.0          /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.35       /* limit on cosine of normal to shown facets */
#define RSCALE_POTENTIAL 0.15   /* radial scaling of potential */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 540.0   /* total angle of rotation during simulation */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */
#define DRAW_POLAR_AXIS 1   /* set to 1 to draw polar axis */

double light[3] = {-0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-2.0, 3.0, 2.0};    /* location of observer for REP_PROJ_3D representation */ 

```

**2D part:**

```
#define DRAW_SPHERE 0    /* set to 1 to draw 3D sphere */

```

