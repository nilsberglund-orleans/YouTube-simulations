 ### Parameter values for YouTube simulations ###

Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Below are parameter values used for different simulations, as well as initial conditions used in 
function animation. Some simulations use variants of the published code. The list is going to be 
updated gradually. 


### 31 December 23 - Compressing molecular caterpillars ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.95
#define INITXMAX 1.95	/* x interval for initial condition */
#define INITYMIN -1.1
#define INITYMAX 1.1	/* y interval for initial condition */

#define ADDXMIN 1.9
#define ADDXMAX 2.0	/* x interval for adding particles */
#define ADDYMIN -0.9
#define ADDYMAX 0.9	/* y interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 10.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.012          /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
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
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 4200      /* number of frames of movie */
#define NVID 30         /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 1
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 0        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11   /* Color palette for different neighbours representation */

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
#define PARTICLE_EMAX 20000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 1000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0     /* set to 1 to add an electric field */
#define EFIELD 200000.0   /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 10000.0       /* value of magnetic field */
#define CHARGE 1.5      /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 1000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 10000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 10.0  /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 1   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.25      /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1600            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 3     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define KSPRING_PAIRS 1.0e10    /* spring constant for pair interaction */
#define NPARTNERS 2         /* number of partners of particles */
#define PAIRING_TYPE 1      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between anions for POLY_WATER case */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100   /* size of hashgrid in x direction */
#define HASHY 50    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 30 December 23 - Trying to grow snowflakes with periodic boundary conditions ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.9
#define INITXMAX 1.9	/* x interval for initial condition */
#define INITYMIN -1.9
#define INITYMAX 13.9	/* y interval for initial condition */

#define ADDXMIN -1.9
#define ADDXMAX 1.9	/* x interval for adding particles */
#define ADDYMIN 1.2
#define ADDYMAX 3.0	/* y interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -2.0
#define BCYMAX 14.0	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 10.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.012 	    /* parameter controlling radius of particles */
#define MU_B 0.007          /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 20           /* number of grid point for grid of disks */
#define NGRIDY 30           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */
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
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 5        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 0          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 1        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11   /* Color palette for different neighbours representation */

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
#define PARTICLE_EMAX 20000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 30.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 1000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 16.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU_B */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0002         /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 2000.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0     /* set to 1 to add an electric field */
#define EFIELD 200000.0   /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 10000.0       /* value of magnetic field */
#define CHARGE -2.0      /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 1000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 10000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 10.0  /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1200            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_COUPLING 1   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 3     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define KSPRING_PAIRS 1.0e10    /* spring constant for pair interaction */
#define NPARTNERS 2         /* number of partners of particles */
#define PAIRING_TYPE 2      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between anions for POLY_WATER case */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100   /* size of hashgrid in x direction */
#define HASHY 400    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 29 December 23 - Evanescent coupling of two optical fibers ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** 

```
    init_wave_flat(phi, psi, xy_in);
    
    if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == 1))
    {
        if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
        add_circular_wave(sign, -1.7, 0.92, phi, psi, xy_in);
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 1      /* set to 1 for a variable index of refraction */
#define IOR 15               /* choice of index of refraction, see list in global_pdes.c */
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

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 522        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 103   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.25	    /* parameter controlling the dimensions of domain */
#define MU 0.03             /* parameter controlling the dimensions of domain */
#define MU_B 0.42           /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY -0.666666666666          /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.04       /* Courant number */
#define COURANTB 0.12      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0005           /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 15     /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 2       /* type of wave packet sources */
#define N_WAVE_PACKETS 2                /* number of wave packets */
#define WAVE_PACKET_RADIUS 25            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2500       /* number of frames of movie */
#define NVID 18           /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
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

#define INITIAL_AMP 0.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.000012   /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.05   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 4        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 17     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 125.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.25    /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 0.5     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.5   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

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
**Part 2:**

```
    #define MU 0.1
    
    add_circular_wave(sign, -1.7, 0.92, phi, psi, xy_in);
```

### 28 December 23 - Terminal Frost: Trying to grow some snowflakes  ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.9
#define INITXMAX 1.9	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 11.125	/* y interval for initial condition */

#define ADDXMIN -1.9
#define ADDXMAX 1.9	/* x interval for adding particles */
#define ADDYMIN 1.2
#define ADDYMAX 3.0	/* y interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 11.5	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 10.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.012 	    /* parameter controlling radius of particles */
#define MU_B 0.007          /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 20           /* number of grid point for grid of disks */
#define NGRIDY 30           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 5000      /* number of frames of movie */
#define NVID 100         /* number of iterations between images displayed on screen */
#define NSEG 25          /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 5        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 0          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 1        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11   /* Color palette for different neighbours representation */

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
#define PARTICLE_EMAX 20000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 30.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 1500.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 16.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU_B */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.001         /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 2000.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0     /* set to 1 to add an electric field */
#define EFIELD 200000.0   /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 10000.0       /* value of magnetic field */
#define CHARGE -2.0      /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 1000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 10000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 1  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 20.0  /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1200            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_COUPLING 1   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 3     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 1.0    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.0   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 1  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define KSPRING_PAIRS 1.0e10    /* spring constant for pair interaction */
#define NPARTNERS 2         /* number of partners of particles */
#define PAIRING_TYPE 2      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between anions for POLY_WATER case */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100   /* size of hashgrid in x direction */
#define HASHY 300    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 27 December 23 - What if 14 active volcanoes erupted at the same time? ###

**Program:** `wave_sphere.c` 

**Initial condition in function `animation()`:** 

```
    /* Kilauea, Hawaii, USA */
    init_circular_wave_sphere(155.2867*PI/180.0, 19.42109*PI/180.0, phi, psi, xy_in, wsphere);
    
    /* Mount Etna, Sicily, Italy */
    add_circular_wave_sphere(1.0, -14.9950*PI/180.0, 37.7550*PI/180.0, phi, psi, xy_in, wsphere);
    
    /* Eyjafjallajökull, Iceland */
    add_circular_wave_sphere(1.0, 19.61333*PI/180.0, 63.62000*PI/180.0, phi, psi, xy_in, wsphere);

    /* Mount Vesuvius, Italy */
    add_circular_wave_sphere(1.0, -14.433*PI/180.0, 40.817*PI/180.0, phi, psi, xy_in, wsphere);

    /* Sakurajima, Japan */
    add_circular_wave_sphere(1.0, -130.65806*PI/180.0, 31.58056*PI/180.0, phi, psi, xy_in, wsphere);

    /* Mount Merapi, Indonesia */
    add_circular_wave_sphere(1.0, -110.44611*PI/180.0, -7.54139*PI/180.0, phi, psi, xy_in, wsphere);
    
    /* Ulawun, Papua New Guinea */
    add_circular_wave_sphere(1.0, -151.33333*PI/180.0, -5.05000*PI/180.0, phi, psi, xy_in, wsphere);
    
    /* Santa María, Guatemala */
    add_circular_wave_sphere(1.0, 91.55167*PI/180.0, 14.75556*PI/180.0, phi, psi, xy_in, wsphere);
    
    /*  Mount Erebus, Antarctica */
    add_circular_wave_sphere(1.0, -167.15333*PI/180.0, -77.52972*PI/180.0, phi, psi, xy_in, wsphere);
    
    /*  Piton de la Fournaise, Reunion island, Indian Ocean */
    add_circular_wave_sphere(0.75, -55.70889*PI/180.0, -21.24250*PI/180.0, phi, psi, xy_in, wsphere);

    /*  Soufrière Hills, Montserrat */
    add_circular_wave_sphere(0.75, 62.1773*PI/180.0, 16.7103*PI/180.0, phi, psi, xy_in, wsphere);
    
    /*  Axial Seamount, Pacific Ocean */
    add_circular_wave_sphere(0.5, 130.00*PI/180.0, 45.95*PI/180.0, phi, psi, xy_in, wsphere);
    
    /*  Havre Seamount, New Zealand */
    add_circular_wave_sphere(0.5, 179.968611*PI/180.0, -31.120278*PI/180.0, phi, psi, xy_in, wsphere);

    /*  Boomerang Seamount, Indian Ocean */
    add_circular_wave_sphere(0.5, -77.825*PI/180.0, -37.721*PI/180.0, phi, psi, xy_in, wsphere);

```

**2D version:**

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 2560          /* number of grid points on x axis */
#define NY 1280           /* number of grid points on y axis */

#define DPOLE 20         /* safety distance to poles */
#define SMOOTHPOLE 0.1     /* smoothing coefficient at poles */
#define ZERO_MERIDIAN 180.0     /* choice of zero meridian (will be at left/right boundary of 2d plot) */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 0        /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 0.6  /* scaling for Julia sets */
#define JULIA_ROT -20.0       /* rotation of Julia set, in degrees */
#define JULIA_RE 0.5    
#define JULIA_IM 0.462    /* parameters for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 84         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 33   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 10              /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.0 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.1              /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 7            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 2000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 20.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30            /* number of grid point for grid of disks */
#define NGRIDY 18            /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 3  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.001       /* frequency of periodic excitation */
#define AMPLITUDE 0.8     /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0       /* damping of periodic excitation */
#define COURANT 0.05       /* Courant number */
#define COURANTB 0.002    /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 1.0e-2         /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 25    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 0  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define WAVE_PACKET_RADIUS 20            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

#define PRECOMPUTE_BC 0     /* set to 1 to compute neighbours for Laplacian in advance */

/* Parameters for length and speed of simulation */

#define NSTEPS 2400       /* number of frames of movie */
#define NVID 4            /* number of iterations between images displayed on screen */
#define NSEG 1000          /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* set to 1 to print speed of moving source */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 3         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 500   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 1.0            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.002  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define ZPLOT 103     /* wave height */
#define CPLOT 103     /* color scheme */

#define ZPLOT_B 108 
#define CPLOT_B 108        /* plot type for second movie */

#define CHANGE_LUMINOSITY 1     /* set to 1 to let luminosity depend on energy flux intensity */
#define FLUX_WINDOW 30          /* size of averaging window of flux intensity */
#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of plot */
#define SHADE_3D 0              /* set to 1 to change luminosity according to normal vector */
#define SHADE_2D 1              /* set to 1 to change luminosity according to normal vector to plane */
#define SHADE_WAVE 1            /* set to 1 to have luminosity depend on wave height */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define FLOOR_ZCOORD 1          /* set to 1 to draw only facets with z not too negative */
#define DRAW_BILLIARD 0         /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 0   /* set to 1 to draw front of boundary after drawing wave */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw construction lines of certain domains */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental, draw outside of billiard in gray */
#define SHADE_SCALE_2D 10.0     /* controls "depth" of 2D shading */
#define COS_LIGHT_MIN 0.0       /* controls angle-dependence of 2D shading */
#define COS_LIGHT_MAX 0.8       /* controls angle-dependence of 2D shading */

#define PLOT_SCALE_ENERGY 0.4          /* vertical scaling in energy plot */
#define PLOT_SCALE_LOG_ENERGY 0.5       /* vertical scaling in log energy plot */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 
#define PLOT_2D 1               /* switch to 2D representation, equirectangular projection */
#define PHISHIFT 0.0            /* shift of phi in 2D plot (in degrees) */
#define FLOODING 1              /* set to 1 to draw waves above altitude (for Earth representations) */

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE -50.0   /* total angle of rotation during simulation */
#define ROTATE_VIEW_WHILE_FADE 1    /* set to 1 to keep rotating viewpoint during fade */

#define VIEWPOINT_TRAJ 2    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */

/* Color schemes */

#define COLOR_PALETTE 11      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 16     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define COLOR_OUT_R 1.0    /* color outside domain */
#define COLOR_OUT_G 1.0    
#define COLOR_OUT_B 1.0    

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.5   /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSCALE_ENERGY 5.0     /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0      /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0    /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 200.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.25      /* scaling factor for energy log representation */
#define LOG_SHIFT 0.5      /* shift of colors on log scale */
#define LOG_ENERGY_FLOOR -10.0    /* floor value for log of (total) energy */
#define LOG_MEAN_ENERGY_SHIFT 1.0   /* additional shift for log of mean energy */
#define FLUX_SCALE 1200.0    /* scaling factor for energy flux representation */
#define FLUX_CSCALE 2.0      /* scaling factor for color in energy flux representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 240.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -200.0    /* amplitude of variation of hue for color scheme C_HUE */

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
#define RAND_SHIFT 1        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

#define DRAW_COLOR_SCHEME 0       /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 6.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 6.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0     /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

#define ADD_POTENTIAL 0         /* set to 1 to add potential to z coordinate */
#define POTENTIAL 10
#define POT_FACT 20.0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters controlling 3D projection */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {-0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-6.0, -3.0, 4.5};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define ADD_DEM 1               /* add DEM (digital elevation model) */
#define ADD_NEGATIVE_DEM 0      /* add DEM with bathymetric data */
#define RSCALE_DEM 0.075          /* scaling factor of radial component for DEM */
#define SMOOTH_DEM 0            /* set to 1 to smoothen DEM (to make altitude less constant) */
#define DEM_SMOOTH_STEPS 5      /* number of smoothening steps */
#define DEM_SMOOTH_HEIGHT 0.5   /* relative height below which to smoothen */
#define DEM_MAXHEIGHT 9000.0     /* max height of DEM (estimated from Everest/Olympus Mons) */
#define PLANET_SEALEVEL 0.0      /* sea level for flooded planet */
#define VENUS_NODATA_FACTOR 0.5     /* altitude to assign to DEM points without data (fraction of mean altitude) */

#define RSCALE 0.025             /* scaling factor of radial component */
#define RMAX 1.5               /* max value of radial component */
#define RMIN 0.5               /* min value of radial component */
#define Z_SCALING_FACTOR 0.8   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.2   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0         /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D 0.0           /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.0           /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.1        /* limit on cosine of normal to shown facets */

```

**3D version:**

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 2560          /* number of grid points on x axis */
#define NY 1280           /* number of grid points on y axis */

#define DPOLE 20         /* safety distance to poles */
#define SMOOTHPOLE 0.1     /* smoothing coefficient at poles */
#define ZERO_MERIDIAN 180.0     /* choice of zero meridian (will be at left/right boundary of 2d plot) */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 0        /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 0.6  /* scaling for Julia sets */
#define JULIA_ROT -20.0       /* rotation of Julia set, in degrees */
#define JULIA_RE 0.5    
#define JULIA_IM 0.462    /* parameters for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 84         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 33   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 10              /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.0 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.1              /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 7            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 2000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 20.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30            /* number of grid point for grid of disks */
#define NGRIDY 18            /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 3  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.001       /* frequency of periodic excitation */
#define AMPLITUDE 0.8     /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0       /* damping of periodic excitation */
#define COURANT 0.05       /* Courant number */
#define COURANTB 0.002    /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 1.0e-2         /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 25    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 0  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define WAVE_PACKET_RADIUS 20            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

#define PRECOMPUTE_BC 0     /* set to 1 to compute neighbours for Laplacian in advance */

/* Parameters for length and speed of simulation */

#define NSTEPS 2400       /* number of frames of movie */
#define NVID 4            /* number of iterations between images displayed on screen */
#define NSEG 1000          /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* set to 1 to print speed of moving source */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 3         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 500   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 1.0            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.002  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define ZPLOT 103     /* wave height */
#define CPLOT 103     /* color scheme */

#define ZPLOT_B 108 
#define CPLOT_B 108        /* plot type for second movie */

#define CHANGE_LUMINOSITY 1     /* set to 1 to let luminosity depend on energy flux intensity */
#define FLUX_WINDOW 30          /* size of averaging window of flux intensity */
#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define SHADE_2D 0              /* set to 1 to change luminosity according to normal vector to plane */
#define SHADE_WAVE 1            /* set to 1 to have luminosity depend on wave height */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define FLOOR_ZCOORD 1          /* set to 1 to draw only facets with z not too negative */
#define DRAW_BILLIARD 0         /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 0   /* set to 1 to draw front of boundary after drawing wave */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw construction lines of certain domains */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental, draw outside of billiard in gray */
#define SHADE_SCALE_2D 10.0     /* controls "depth" of 2D shading */
#define COS_LIGHT_MIN 0.0       /* controls angle-dependence of 2D shading */
#define COS_LIGHT_MAX 0.8       /* controls angle-dependence of 2D shading */

#define PLOT_SCALE_ENERGY 0.4          /* vertical scaling in energy plot */
#define PLOT_SCALE_LOG_ENERGY 0.5       /* vertical scaling in log energy plot */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 
#define PLOT_2D 0               /* switch to 2D representation, equirectangular projection */
#define PHISHIFT 0.0            /* shift of phi in 2D plot (in degrees) */
#define FLOODING 1              /* set to 1 to draw waves above altitude (for Earth representations) */

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 720.0   /* total angle of rotation during simulation */
#define ROTATE_VIEW_WHILE_FADE 1    /* set to 1 to keep rotating viewpoint during fade */

#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */

/* Color schemes */

#define COLOR_PALETTE 11      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 16     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define COLOR_OUT_R 1.0    /* color outside domain */
#define COLOR_OUT_G 1.0    
#define COLOR_OUT_B 1.0    

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.5   /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSCALE_ENERGY 5.0     /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0      /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0    /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 200.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.25      /* scaling factor for energy log representation */
#define LOG_SHIFT 0.5      /* shift of colors on log scale */
#define LOG_ENERGY_FLOOR -10.0    /* floor value for log of (total) energy */
#define LOG_MEAN_ENERGY_SHIFT 1.0   /* additional shift for log of mean energy */
#define FLUX_SCALE 1200.0    /* scaling factor for energy flux representation */
#define FLUX_CSCALE 2.0      /* scaling factor for color in energy flux representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 240.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -200.0    /* amplitude of variation of hue for color scheme C_HUE */

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
#define RAND_SHIFT 1        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

#define DRAW_COLOR_SCHEME 1       /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 6.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 6.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0     /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

#define ADD_POTENTIAL 0         /* set to 1 to add potential to z coordinate */
#define POTENTIAL 10
#define POT_FACT 20.0
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters controlling 3D projection */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {-0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-6.0, -3.0, 4.5};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define ADD_DEM 1               /* add DEM (digital elevation model) */
#define ADD_NEGATIVE_DEM 0      /* add DEM with bathymetric data */
#define RSCALE_DEM 0.075          /* scaling factor of radial component for DEM */
#define SMOOTH_DEM 0            /* set to 1 to smoothen DEM (to make altitude less constant) */
#define DEM_SMOOTH_STEPS 5      /* number of smoothening steps */
#define DEM_SMOOTH_HEIGHT 0.5   /* relative height below which to smoothen */
#define DEM_MAXHEIGHT 9000.0     /* max height of DEM (estimated from Everest/Olympus Mons) */
#define PLANET_SEALEVEL 0.0      /* sea level for flooded planet */
#define VENUS_NODATA_FACTOR 0.5     /* altitude to assign to DEM points without data (fraction of mean altitude) */

#define RSCALE 0.025             /* scaling factor of radial component */
#define RMAX 1.5               /* max value of radial component */
#define RMIN 0.5               /* min value of radial component */
#define Z_SCALING_FACTOR 0.8   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.2   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0         /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D 0.0           /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.0           /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.7        /* limit on cosine of normal to shown facets */

```

### 26 December 23 - Return of the water molecules ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.125
#define INITYMAX 1.125	/* y interval for initial condition */

#define ADDXMIN 1.9
#define ADDXMAX 2.0	/* x interval for adding particles */
#define ADDYMIN -0.9
#define ADDYMAX 0.9	/* y interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.016 	    /* parameter controlling radius of particles */
#define MU_B 0.009          /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 55           /* number of grid point for grid of disks */
#define NGRIDY 27           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 3600      /* number of frames of movie */
#define NVID 30         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 0        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */

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
#define PARTICLE_EMAX 20000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 60.0     /* hue of particles of type 0 */
#define HUE_TYPE1 280.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 100.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 3000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 5.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU_B */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.00002          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0     /* set to 1 to add an electric field */
#define EFIELD 200000.0   /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 10000.0       /* value of magnetic field */
#define CHARGE 2.0      /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 1000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 10000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 1  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 5.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1200            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 3  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define KSPRING_PAIRS 1.0e10    /* spring constant for pair interaction */
#define NPARTNERS 2         /* number of partners of particles */
#define PAIRING_TYPE 2      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 104.45    /* angle (in degrees) between anions for POLY_WATER case */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100   /* size of hashgrid in x direction */
#define HASHY 50    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 25 December 23 - A Christmas wave guide (short version) ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** 

```
    init_wave_flat(phi, psi, xy_in);
    
    if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == 1))
    {
        if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
        add_circular_wave(sign, 0.0, -3.1*LAMBDA, phi, psi, xy_in);
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 1      /* set to 1 for a variable index of refraction */
#define IOR 14               /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1150  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 2300          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -1.197916667
#define XMAX 1.197916667	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 73        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 103   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.32	    /* parameter controlling the dimensions of domain */
#define MU 0.48             /* parameter controlling the dimensions of domain */
#define MU_B 0.42           /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY -0.666666666666          /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.04       /* Courant number */
#define COURANTB 0.12      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0001           /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 15     /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 2       /* type of wave packet sources */
#define N_WAVE_PACKETS 2                /* number of wave packets */
#define WAVE_PACKET_RADIUS 25            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 1100       /* number of frames of movie */
#define NVID 15           /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
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

#define INITIAL_AMP 0.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.000012   /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.05   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 5        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 17     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 50.0      /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.25     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 0.5     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.5   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

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

### 25 December 23 - A Christmas wave guide ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** 

```
    init_wave_flat(phi, psi, xy_in);
    
    if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == 1))
    {
        if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
        add_circular_wave(sign, 0.0, -3.1*LAMBDA, phi, psi, xy_in);
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 1      /* set to 1 for a variable index of refraction */
#define IOR 14               /* choice of index of refraction, see list in global_pdes.c */
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

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 73        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 103   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.32	    /* parameter controlling the dimensions of domain */
#define MU 0.48             /* parameter controlling the dimensions of domain */
#define MU_B 0.42           /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY -0.666666666666          /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.06       /* Courant number */
#define COURANTB 0.18      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0001           /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 8     /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 2       /* type of wave packet sources */
#define N_WAVE_PACKETS 2                /* number of wave packets */
#define WAVE_PACKET_RADIUS 25            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 3750       /* number of frames of movie */
#define NVID 25           /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
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

#define INITIAL_AMP 0.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.000012   /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.05   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 5

#define PLOT_B 0        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 13     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 17     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 50.0      /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.25     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 0.5     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.5   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

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

### 24 December 23 - Let it snow! (short verion) ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	900  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -1.125
#define XMAX 1.125	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.125
#define INITXMAX 1.125	/* x interval for initial condition */
#define INITYMIN 1.125
#define INITYMAX 2.5	/* y interval for initial condition */

#define ADDXMIN -1.125
#define ADDXMAX 1.125	/* x interval for adding particles */
#define ADDYMIN 1.2
#define ADDYMAX 4.0	/* y interval for adding particles */

#define BCXMIN -1.5
#define BCXMAX 1.5	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 4.5	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 28      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.07	    /* parameter controlling the dimensions of domain */
#define MU 0.006 	    /* parameter controlling radius of particles */
#define MU_B 0.009          /* parameter controlling radius of particles of second type */
#define NPOLY 3            /* number of sides of polygon */
#define APOLY -1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 1           /* number of grid point for grid of disks */
#define NGRIDY 2           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */
#define NTREES 15           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 1200      /* number of frames of movie */
#define NVID 60         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 3   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 0        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11   /* Color palette for different neighbours representation */

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
#define PARTICLE_HUE_MIN 300.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 20000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 1     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 30.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 100.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 16.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU_B */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.00002          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 5000.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0     /* set to 1 to add an electric field */
#define EFIELD 200000.0   /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 10000.0       /* value of magnetic field */
#define CHARGE -2.0      /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 1000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 10000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 5.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1200            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 3  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 5.0e8        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define KSPRING_PAIRS 1.0e10    /* spring constant for pair interaction */
#define NPARTNERS 2         /* number of partners of particles */
#define PAIRING_TYPE 2      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 90.0    /* angle (in degrees) between anions for POLY_WATER case */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60   /* size of hashgrid in x direction */
#define HASHY 30    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 24 December 23 - Let it snow! ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN 1.125
#define INITYMAX 2.5	/* y interval for initial condition */

#define ADDXMIN -2.0
#define ADDXMAX 2.0	/* x interval for adding particles */
#define ADDYMIN 1.2
#define ADDYMAX 4.0	/* y interval for adding particles */

#define BCXMIN -2.5
#define BCXMAX 2.5	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 4.5	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 28      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1        /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.075	    /* parameter controlling the dimensions of domain */
#define MU 0.006 	    /* parameter controlling radius of particles */
#define MU_B 0.009          /* parameter controlling radius of particles of second type */
#define NPOLY 3            /* number of sides of polygon */
#define APOLY -1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 1           /* number of grid point for grid of disks */
#define NGRIDY 2           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */
#define NTREES 20           /* number of trees in S_TREES */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 4000      /* number of frames of movie */
#define NVID 30         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 3   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 0  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 0        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */
#define COLOR_PALETTE_PRESSURE 11   /* Color palette for different neighbours representation */

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
#define PARTICLE_HUE_MIN 300.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 20000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 1     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 30.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 100.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 16.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU_B */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.00002          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 5000.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0     /* set to 1 to add an electric field */
#define EFIELD 200000.0   /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 10000.0       /* value of magnetic field */
#define CHARGE -2.0      /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 1000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 10000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 5.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1200            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 3  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 5.0e8        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 0    /* set to 1 to form particles pairs */
#define KSPRING_PAIRS 1.0e10    /* spring constant for pair interaction */
#define NPARTNERS 2         /* number of partners of particles */
#define PAIRING_TYPE 2      /* type of pairing, see POLY_ in global_ljones.c */
#define PARTNER_ANGLE 90.0    /* angle (in degrees) between anions for POLY_WATER case */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60   /* size of hashgrid in x direction */
#define HASHY 30    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 23 December 23 - Side by side comparison of two optical fibers ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** 

```
    init_wave_flat(phi, psi, xy_in);
    
    if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == 1))
    {
        if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
        add_circular_wave(sign, -1.5, 0.8, phi, psi, xy_in);
        add_circular_wave(sign, -1.5, -0.35, phi, psi, xy_in);
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 13               /* choice of index of refraction, see list in global_pdes.c */
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

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 521        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 103   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.48             /* parameter controlling the dimensions of domain */
#define MU_B 0.42           /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY -0.666666666666          /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.06       /* Courant number */
#define COURANTB 0.18       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.007            /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 8     /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 2       /* type of wave packet sources */
#define N_WAVE_PACKETS 2                /* number of wave packets */
#define WAVE_PACKET_RADIUS 25            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2700       /* number of frames of movie */
#define NVID 20           /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
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

#define INITIAL_AMP 0.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.000012   /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.05   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 5        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 17     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 40.0      /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.25     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.5     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.4   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

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
### 22 December 23 - Molecular caterpillars ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.125
#define INITYMAX 1.125	/* y interval for initial condition */

#define ADDXMIN 1.9
#define ADDXMAX 2.0	/* x interval for adding particles */
#define ADDYMIN -0.9
#define ADDYMAX 0.9	/* y interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.012          /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
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
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 2400      /* number of frames of movie */
#define NVID 30         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 0        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */

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
#define PARTICLE_EMAX 10000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 60.0     /* hue of particles of type 0 */
#define HUE_TYPE1 280.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 3000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.00004          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0     /* set to 1 to add an electric field */
#define EFIELD 200000.0   /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 10000.0       /* value of magnetic field */
#define CHARGE 1.5      /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 1000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 10000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 1  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 5.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1200            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 3  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define KSPRING_PAIRS 1.0e10    /* spring constant for pair interaction */
#define NPARTNERS 2         /* number of partners of particles */
#define PAIRING_TYPE 1      /* type of pairing, see POLY_ in global_ljones.c */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100   /* size of hashgrid in x direction */
#define HASHY 50    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 21 December 23 - Optical fibers, single-mode vs multimode ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** 

```
    init_wave_flat(phi, psi, xy_in);
    
    if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == 1))
    {
        if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
        add_circular_wave(sign, -1.5, 0.6, phi, psi, xy_in);
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 13               /* choice of index of refraction, see list in global_pdes.c */
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

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 52        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 103   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.48             /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY -0.666666666666          /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.06       /* Courant number */
#define COURANTB 0.18       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.007            /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 8     /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 2       /* type of wave packet sources */
#define N_WAVE_PACKETS 2                /* number of wave packets */
#define WAVE_PACKET_RADIUS 25            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2400       /* number of frames of movie */
#define NVID 20           /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
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

#define INITIAL_AMP 0.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.000012   /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.05   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 5        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 17     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 60.0      /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.25     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.5     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.6   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

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

**Part 2:**

```
#define MU 0.42             /* parameter controlling the dimensions of domain */

```

### 20 December 23 - Triatomic molecules under temperature cycling ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.125
#define INITYMAX 1.125	/* y interval for initial condition */

#define ADDXMIN 1.9
#define ADDXMAX 2.0	/* x interval for adding particles */
#define ADDYMIN -0.9
#define ADDYMAX 0.9	/* y interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.012 	    /* parameter controlling radius of particles */
#define MU_B 0.01           /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 70           /* number of grid point for grid of disks */
#define NGRIDY 35           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 4100      /* number of frames of movie */
#define NVID 30         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 0        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */

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
#define PARTICLE_EMAX 10000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 60.0     /* hue of particles of type 0 */
#define HUE_TYPE1 280.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 3000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 2.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.00004          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0     /* set to 1 to add an electric field */
#define EFIELD 200000.0   /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 10000.0       /* value of magnetic field */
#define CHARGE 2.0      /* charge of particles of first type */
#define CHARGE_B -1.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 1000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 10000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 1  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 5.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1200            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 3  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define KSPRING_PAIRS 1.0e10    /* spring constant for pair interaction */
#define NPARTNERS 2         /* number of partners of particles */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 100   /* size of hashgrid in x direction */
#define HASHY 50    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 19 December 23 - Flooding the Moon (short version) ###

**Program:** `earth_sealevel.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

/* General geometrical parameters */

#define WINWIDTH 	1150  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 2048           /* number of grid points on x axis */
#define NY 2048           /* number of grid points on y axis */

#define DPOLE 20         /* safety distance to poles */
#define SMOOTHPOLE 0.1     /* smoothing coefficient at poles */
#define ZERO_MERIDIAN 180.0     /* choice of zero meridian (will be at left/right boundary of 2d plot) */

#define XMIN -1.041666667
#define XMAX 1.041666667	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 0        /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 0.25     /* scaling for Julia sets */
#define JULIA_ROT 180.0       /* rotation of Julia set, in degrees */
#define JULIA_RE -0.77145    
#define JULIA_IM -0.10295    /* parameters for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 87         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 33   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 10              /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.0 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.1             /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 7            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 2000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 20.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30            /* number of grid point for grid of disks */
#define NGRIDY 18            /* number of grid point for grid of disks */

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

#define OMEGA 0.001       /* frequency of periodic excitation */
#define AMPLITUDE 0.8     /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0       /* damping of periodic excitation */
#define COURANT 0.05       /* Courant number */
#define COURANTB 0.002    /* Courant number in medium B */
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

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 30    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define WAVE_PACKET_RADIUS 20            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

#define PRECOMPUTE_BC 0     /* set to 1 to compute neighbours for Laplacian in advance */

/* Parameters for length and speed of simulation */

#define NSTEPS 1200       /* number of frames of movie */
#define NVID 4            /* number of iterations between images displayed on screen */
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

/* Parameters of initial condition */

#define INITIAL_AMP 0.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.001  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.05  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define ZPLOT 103     /* wave height */
#define CPLOT 103     /* color scheme */

#define ZPLOT_B 108 
#define CPLOT_B 108        /* plot type for second movie */

#define CHANGE_LUMINOSITY 1     /* set to 1 to let luminosity depend on energy flux intensity */
#define FLUX_WINDOW 30          /* size of averaging window of flux intensity */
#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define SHADE_2D 0              /* set to 1 to change luminosity according to normal vector to plane */
#define SHADE_WAVE 1            /* set to 1 to have luminosity depend on wave height */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define FLOOR_ZCOORD 1          /* set to 1 to draw only facets with z not too negative */
#define DRAW_BILLIARD 0         /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 0   /* set to 1 to draw front of boundary after drawing wave */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw construction lines of certain domains */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental, draw outside of billiard in gray */
#define SHADE_SCALE_2D 10.0     /* controls "depth" of 2D shading */
#define COS_LIGHT_MIN 0.0       /* controls angle-dependence of 2D shading */
#define COS_LIGHT_MAX 0.8       /* controls angle-dependence of 2D shading */

#define PLOT_SCALE_ENERGY 0.4          /* vertical scaling in energy plot */
#define PLOT_SCALE_LOG_ENERGY 0.5       /* vertical scaling in log energy plot */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 
#define PLOT_2D 0               /* switch to 2D representation, equirectangular projection */
#define PHISHIFT 0.0            /* shift of phi in 2D plot (in degrees) */
#define FLOODING 1              /* set to 1 to draw waves above altitude (for Earth representations) */

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0   /* total angle of rotation during simulation */
#define ROTATE_VIEW_WHILE_FADE 0    /* set to 1 to keep rotating viewpoint during fade */

#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */

/* Color schemes */

#define COLOR_PALETTE 11      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 16     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define COLOR_OUT_R 1.0    /* color outside domain */
#define COLOR_OUT_G 1.0    
#define COLOR_OUT_B 1.0    

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.0   /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSCALE_ENERGY 4.0     /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0      /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0    /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 150.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.5      /* shift of colors on log scale */
#define LOG_ENERGY_FLOOR -10.0    /* floor value for log of (total) energy */
#define LOG_MEAN_ENERGY_SHIFT 1.0   /* additional shift for log of mean energy */
#define FLUX_SCALE 1200.0    /* scaling factor for energy flux representation */
#define FLUX_CSCALE 2.0      /* scaling factor for color in energy flux representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 240.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -200.0    /* amplitude of variation of hue for color scheme C_HUE */

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

#define DRAW_COLOR_SCHEME 0       /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 3.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0     /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

#define ADD_POTENTIAL 0         /* set to 1 to add potential to z coordinate */
#define POTENTIAL 10
#define POT_FACT 20.0
/* end of constants only used by sub_wave and sub_maze */


/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters controlling 3D projection */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {-0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-7.0, -2.0, 2.5};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */


#define ADD_DEM 1               /* add DEM (digital elevation model) */
#define ADD_NEGATIVE_DEM 0      /* add DEM with bathymetric data */
#define RSCALE_DEM 0.1          /* scaling factor of radial component for DEM */
#define SMOOTH_DEM 1            /* set to 1 to smoothen DEM (to make altitude less constant) */
#define DEM_SMOOTH_STEPS 5      /* number of smoothening steps */
#define DEM_SMOOTH_HEIGHT 1.0   /* relative height below which to smoothen */
#define DEM_MAXHEIGHT 15000     /* max height of DEM */
#define DEM_MAXDEPTH -10000     /* max depth of DEM */
#define PLANET_SEALEVEL 0.0      /* sea level for flooded planet */
#define VENUS_NODATA_FACTOR 0.5     /* altitude to assign to DEM points without data (fraction of mean altitude) */

#define MAXHEIGHT 5000.0        /* max sea level height */
#define MAXDEPTH 4000.0         /* max sea level depth */

#define RSCALE 0.01             /* scaling factor of radial component */
#define RMAX 10.0               /* max value of radial component */
#define RMIN 0.5               /* min value of radial component */
#define Z_SCALING_FACTOR 0.8    /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.2   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define X_SCALING_FACTOR 1.8   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define Y_SCALING_FACTOR 2.2   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0         /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D 0.0           /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.0           /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.5        /* limit on cosine of normal to shown facets */

```

### 18 December 23 - Interacting polar diatomic molecules ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.95
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.075
#define INITYMAX 1.125	/* y interval for initial condition */

#define ADDXMIN 1.9
#define ADDXMAX 2.0	/* x interval for adding particles */
#define ADDYMIN -0.9
#define ADDYMAX 0.9	/* y interval for adding particles */


#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 1.0 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.014          /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
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
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 4200      /* number of frames of movie */
#define NVID 30         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 5        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 0          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 4        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */

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
#define PARTICLE_EMAX 10000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 3000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.5   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.00004          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1     /* set to 1 to add an electric field */
#define EFIELD 200000.0   /* value of electric field */
#define ADD_BFIELD 1     /* set to 1 to add a magnetic field */
#define BFIELD 10000.0       /* value of magnetic field */
#define CHARGE -1.0      /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 1000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 10000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 1  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 10.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1200            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 3  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define PAIR_PARTICLES 1    /* set to 1 to form particles pairs */
#define KSPRING_PAIRS 1.0e10    /* spring constant for pair interaction */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 90   /* size of hashgrid in x direction */
#define HASHY 45    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 17 December 23 - Combining a biconcave and a convex lens ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 1      /* set to 1 for a variable index of refraction */
#define IOR 13               /* choice of index of refraction, see list in global_pdes.c */
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

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 69        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 103   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.8	    /* parameter controlling the dimensions of domain */
#define MU 0.915            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY -0.666666666666          /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.04      /* Courant number in medium B */
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
#define OSCILLATING_SOURCE_PERIOD 20    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 2       /* type of wave packet sources */
#define N_WAVE_PACKETS 2                /* number of wave packets */
#define WAVE_PACKET_RADIUS 25            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3

/* Parameters for length and speed of simulation */

#define NSTEPS 2600       /* number of frames of movie */
#define NVID 7            /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 300      /* time after which to start saving frames */
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

#define INITIAL_AMP 1.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00001    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 5        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 17     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.25     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.4   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 1   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

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

### 16 December 23 - Crystals or proteins? ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.01	/* x interval for initial condition */
#define INITYMIN -1.1
#define INITYMAX 1.17	/* y interval for initial condition */

#define ADDXMIN 1.9
#define ADDXMAX 2.0	/* x interval for adding particles */
#define ADDYMIN -0.9
#define ADDYMAX 0.9	/* y interval for adding particles */


#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

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

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.01           /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 90           /* number of grid point for grid of disks */
#define NGRIDY 45           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 4500      /* number of frames of movie */
#define NVID 25         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 0

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 0        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */

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
#define PARTICLE_EMAX 30000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 3000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0     /* set to 1 to add an electric field */
#define EFIELD 1000000.0      /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 5333.33333333333       /* value of magnetic field */
#define CHARGE -1.0      /* charge of particles of first type */
#define CHARGE_B 1.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 1000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 200.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1200            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 3  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 0      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 80   /* size of hashgrid in x direction */
#define HASHY 40    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 15 December 23 - A biconcave lens ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 1      /* set to 1 for a variable index of refraction */
#define IOR 12               /* choice of index of refraction, see list in global_pdes.c */
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

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 68        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 103   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.8	    /* parameter controlling the dimensions of domain */
#define MU 0.9              /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY -0.666666666666          /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.01         /* frequency of periodic excitation */
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
#define OSCILLATING_SOURCE_PERIOD 20    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 2       /* type of wave packet sources */
#define N_WAVE_PACKETS 2                /* number of wave packets */
#define WAVE_PACKET_RADIUS 25            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3

/* Parameters for length and speed of simulation */

#define NSTEPS 2400       /* number of frames of movie */
#define NVID 7            /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 300      /* time after which to start saving frames */
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

#define INITIAL_AMP 1.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00001    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 5        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 18     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.25     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.4   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 1   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

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

### 14 December 23 - Changing the anion/cation proportion in a mixture ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.01	/* x interval for initial condition */
#define INITYMIN -1.1
#define INITYMAX 1.17	/* y interval for initial condition */

#define ADDXMIN 1.9
#define ADDXMAX 2.0	/* x interval for adding particles */
#define ADDYMIN -0.9
#define ADDYMAX 0.9	/* y interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.1 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.012          /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
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
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 2900      /* number of frames of movie */
#define NVID 25         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 13        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 0        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */

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
#define PRINT_TYPE_PROP 1           /* print type proportion */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 30000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 3000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.5   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1     /* set to 1 to add an electric field */
#define EFIELD 1000000.0      /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 5.33333333333333       /* value of magnetic field */
#define CHARGE -1.0      /* charge of particles of first type */
#define CHARGE_B 2.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 1000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 200.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1200            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 3  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define CHANGE_TYPES 1      /* set to 1 to change type proportion in course of simulation */
#define PROP_MIN 0.1        /* min proportion of type 1 particles */
#define PROP_MAX 0.9        /* max proportion of type 1 particles */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60   /* size of hashgrid in x direction */
#define HASHY 30    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 13 December 23 - Linear waves crossing a lens angled  at 45 and 60 degrees ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 10               /* choice of index of refraction, see list in global_pdes.c */
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

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 67        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 103   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.6	    /* parameter controlling the dimensions of domain */
#define MU 1.11             /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY -0.5          /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.058      /* Courant number in medium B */
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
#define OSCILLATING_SOURCE_PERIOD 20    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 2       /* type of wave packet sources */
#define N_WAVE_PACKETS 2                /* number of wave packets */
#define WAVE_PACKET_RADIUS 25            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3

/* Parameters for length and speed of simulation */

#define NSTEPS 2500       /* number of frames of movie */
#define NVID 7            /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 300      /* time after which to start saving frames */
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

#define INITIAL_AMP 1.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00001    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 5        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 17     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.25     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.4   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 1   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

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

**Part 2:** 

```
#define APOLY -0.666666666666          /* angle by which to turn polygon, in units of Pi/2 */ 
```

### 12 December 23 - Ions in an increasing electric field ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.01	/* x interval for initial condition */
#define INITYMIN -1.1
#define INITYMAX 1.17	/* y interval for initial condition */

#define ADDXMIN 1.9
#define ADDXMAX 2.0	/* x interval for adding particles */
#define ADDYMIN -0.9
#define ADDYMAX 0.9	/* y interval for adding particles */


#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.35 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.01 	    /* parameter controlling radius of particles */
#define MU_B 0.015          /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
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
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 2700      /* number of frames of movie */
#define NVID 25         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 5        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 0          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 2        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */

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
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 30000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 3000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.5   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1     /* set to 1 to add an electric field */
#define EFIELD 5.33333333333333      /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 5.33333333333333       /* value of magnetic field */
#define CHARGE -1.0      /* charge of particles of first type */
#define CHARGE_B 2.0     /* charge of particles of second type */
#define INCREASE_E 1     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 1000000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 200.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1200            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 3  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60   /* size of hashgrid in x direction */
#define HASHY 30    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 11 December 23 - Two lenses with a high index of refraction ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** 

```
    init_wave_flat(phi, psi, xy_in);
    
    if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == 1))
    {
        if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
        add_circular_wave(sign, -1.75, 0.16, phi, psi, xy_in);
        add_circular_wave(sign, -1.75, -0.16, phi, psi, xy_in);
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 1      /* set to 1 for a variable index of refraction */
#define IOR 10               /* choice of index of refraction, see list in global_pdes.c */
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

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 63        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 103   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define MU 0.7              /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY -0.666666666666          /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
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
#define OSCILLATING_SOURCE_PERIOD 20    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 2       /* type of wave packet sources */
#define N_WAVE_PACKETS 2                /* number of wave packets */
#define WAVE_PACKET_RADIUS 25            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 3600       /* number of frames of movie */
#define NVID 7            /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 300      /* time after which to start saving frames */
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

#define INITIAL_AMP 1.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00001    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 5        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 17     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.25     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.4   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 1   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

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

### 10 December 23 - Ions in an electric field under temperature cycling ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.01	/* x interval for initial condition */
#define INITYMIN -1.1
#define INITYMAX 1.17	/* y interval for initial condition */

#define ADDXMIN 1.9
#define ADDXMAX 2.0	/* x interval for adding particles */
#define ADDYMIN -0.9
#define ADDYMAX 0.9	/* y interval for adding particles */


#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.33 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.012 	    /* parameter controlling radius of particles */
#define MU_B 0.018          /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 85           /* number of grid point for grid of disks */
#define NGRIDY 43           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 2400      /* number of frames of movie */
#define NVID 25         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 5        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 0          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 2        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */

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
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 30000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 3000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.5   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.00002          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1     /* set to 1 to add an electric field */
#define EFIELD 500000.0      /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 5.33333333333333       /* value of magnetic field */
#define CHARGE -2.0      /* charge of particles of first type */
#define CHARGE_B 3.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 100000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 1  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 200.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1200            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 3  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60   /* size of hashgrid in x direction */
#define HASHY 30    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 09 December 23 - Linear waves crossing an angled lens ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 10               /* choice of index of refraction, see list in global_pdes.c */
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

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 67        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 103   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.6	    /* parameter controlling the dimensions of domain */
#define MU 1.11             /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY -0.1          /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.01         /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.058      /* Courant number in medium B */
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
#define OSCILLATING_SOURCE_PERIOD 20    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 2       /* type of wave packet sources */
#define N_WAVE_PACKETS 2                /* number of wave packets */
#define WAVE_PACKET_RADIUS 25            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3

/* Parameters for length and speed of simulation */

#define NSTEPS 2500       /* number of frames of movie */
#define NVID 7            /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 300      /* time after which to start saving frames */
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

#define INITIAL_AMP 1.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00001    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 5        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 17     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.25     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.4   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 1   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 1     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

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

**Part 2:**

```
    #define APOLY -0.2          /* angle by which to turn polygon, in units of Pi/2 */ 
```


### 08 December 23 - Bloopers 12: Compressing ions and electrons with mass ratio 100 ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.01	/* x interval for initial condition */
#define INITYMIN -1.1
#define INITYMAX 1.17	/* y interval for initial condition */

#define ADDXMIN 1.9
#define ADDXMAX 2.0	/* x interval for adding particles */
#define ADDYMIN -0.9
#define ADDYMAX 0.9	/* y interval for adding particles */


#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

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

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.008 	    /* parameter controlling radius of particles */
#define MU_B 0.018          /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 85           /* number of grid point for grid of disks */
#define NGRIDY 43           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 4500      /* number of frames of movie */
#define NVID 25         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 13
#define PLOT_B 5        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 0          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 2        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */

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
#define PRINT_NPARTICLES 1          /* print number of active particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 30000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 3000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 100.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.000005          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e11    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 2.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0     /* set to 1 to add an electric field */
#define EFIELD 100000.0      /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 5.33333333333333       /* value of magnetic field */
#define CHARGE -1.0      /* charge of particles of first type */
#define CHARGE_B 2.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 100000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 200.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 1   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1200            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 3  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60   /* size of hashgrid in x direction */
#define HASHY 30    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 07 December 23 - A binary zone plate ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** 

```
    init_wave_flat(phi, psi, xy_in);
    
    if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == 1))
    {
        if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
        add_circular_wave(sign, -1.75, 0.25, phi, psi, xy_in);
        add_circular_wave(sign, -1.75, -0.25, phi, psi, xy_in);
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 10               /* choice of index of refraction, see list in global_pdes.c */
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

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 65        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 103   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.0	    /* parameter controlling the dimensions of domain */
#define MU 0.075             /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 10           /* number of grid point for grid of disks */

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
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.024        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.063      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT -400.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 20    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 2       /* type of wave packet sources */
#define N_WAVE_PACKETS 2                /* number of wave packets */
#define WAVE_PACKET_RADIUS 25            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2500       /* number of frames of movie */
#define NVID 7            /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 300      /* time after which to start saving frames */
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

#define INITIAL_AMP 1.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00001    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 5        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 17     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.25     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.4   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 1   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

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

### 06 December 23 - Injecting electrons in a crystal of ions ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.01	/* x interval for initial condition */
#define INITYMIN -1.1
#define INITYMAX 1.17	/* y interval for initial condition */

#define ADDXMIN 1.9
#define ADDXMAX 2.0	/* x interval for adding particles */
#define ADDYMIN -0.9
#define ADDYMAX 0.9	/* y interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 1   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.55 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.008 	    /* parameter controlling radius of particles */
#define MU_B 0.018          /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 10         /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 30
#define NOBSY 20            /* obstacles for O_HEX obstacle pattern */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 4750      /* number of frames of movie */
#define NVID 25         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 5        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 4        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */

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
#define PRINT_NPARTICLES 1          /* print number of active particles */
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 30000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */
#define BG_FORCE_SLOPE 7.5e-8   /* contant in BG_FORCE backgound color scheme*/

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 3000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 25.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.00002          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e10    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0     /* set to 1 to add an electric field */
#define EFIELD 100000.0      /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 5.33333333333333       /* value of magnetic field */
#define CHARGE -1.0      /* charge of particles of first type */
#define CHARGE_B 2.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 100000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 1000.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 200.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1000            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 4       /* time interval between adding further particles */
#define N_ADD_PARTICLES 3  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 30   /* size of hashgrid in x direction */
#define HASHY 20    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 05 December 23 - Flooding Mercury ###

**Program:** `earth_sealevel.c` 

**Part 1 (3D)** 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 2304           /* number of grid points on x axis */
#define NY 1152           /* number of grid points on y axis */

#define DPOLE 20         /* safety distance to poles */
#define SMOOTHPOLE 0.1     /* smoothing coefficient at poles */
#define ZERO_MERIDIAN 0.0     /* choice of zero meridian (will be at left/right boundary of 2d plot) */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 0        /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 0.25     /* scaling for Julia sets */
#define JULIA_ROT 180.0       /* rotation of Julia set, in degrees */
#define JULIA_RE -0.77145    
#define JULIA_IM -0.10295    /* parameters for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 89         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 33   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 10              /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.0 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.1             /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 7            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 2000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 20.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30            /* number of grid point for grid of disks */
#define NGRIDY 18            /* number of grid point for grid of disks */

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

#define OMEGA 0.001       /* frequency of periodic excitation */
#define AMPLITUDE 0.8     /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0       /* damping of periodic excitation */
#define COURANT 0.05       /* Courant number */
#define COURANTB 0.002    /* Courant number in medium B */
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

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 30    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define WAVE_PACKET_RADIUS 20            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

#define PRECOMPUTE_BC 0     /* set to 1 to compute neighbours for Laplacian in advance */

/* Parameters for length and speed of simulation */

#define NSTEPS 3500       /* number of frames of movie */
#define NVID 4            /* number of iterations between images displayed on screen */
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

/* Parameters of initial condition */

#define INITIAL_AMP 0.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.001  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.05  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define ZPLOT 103     /* wave height */
#define CPLOT 103     /* color scheme */

#define ZPLOT_B 108 
#define CPLOT_B 108        /* plot type for second movie */

#define CHANGE_LUMINOSITY 1     /* set to 1 to let luminosity depend on energy flux intensity */
#define FLUX_WINDOW 30          /* size of averaging window of flux intensity */
#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define SHADE_2D 0              /* set to 1 to change luminosity according to normal vector to plane */
#define SHADE_WAVE 1            /* set to 1 to have luminosity depend on wave height */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define FLOOR_ZCOORD 1          /* set to 1 to draw only facets with z not too negative */
#define DRAW_BILLIARD 0         /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 0   /* set to 1 to draw front of boundary after drawing wave */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw construction lines of certain domains */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental, draw outside of billiard in gray */
#define SHADE_SCALE_2D 10.0     /* controls "depth" of 2D shading */
#define COS_LIGHT_MIN 0.0       /* controls angle-dependence of 2D shading */
#define COS_LIGHT_MAX 0.8       /* controls angle-dependence of 2D shading */

#define PLOT_SCALE_ENERGY 0.4          /* vertical scaling in energy plot */
#define PLOT_SCALE_LOG_ENERGY 0.5       /* vertical scaling in log energy plot */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 
#define PLOT_2D 0               /* switch to 2D representation, equirectangular projection */
#define PHISHIFT 0.0            /* shift of phi in 2D plot (in degrees) */
#define FLOODING 1              /* set to 1 to draw waves above altitude (for Earth representations) */

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 720.0   /* total angle of rotation during simulation */
#define ROTATE_VIEW_WHILE_FADE 0    /* set to 1 to keep rotating viewpoint during fade */

#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */

/* Color schemes */

#define COLOR_PALETTE 11      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 16     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define COLOR_OUT_R 1.0    /* color outside domain */
#define COLOR_OUT_G 1.0    
#define COLOR_OUT_B 1.0    

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.0   /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSCALE_ENERGY 4.0     /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0      /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0    /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 150.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.5      /* shift of colors on log scale */
#define LOG_ENERGY_FLOOR -10.0    /* floor value for log of (total) energy */
#define LOG_MEAN_ENERGY_SHIFT 1.0   /* additional shift for log of mean energy */
#define FLUX_SCALE 1200.0    /* scaling factor for energy flux representation */
#define FLUX_CSCALE 2.0      /* scaling factor for color in energy flux representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 240.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -200.0    /* amplitude of variation of hue for color scheme C_HUE */

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

#define DRAW_COLOR_SCHEME 0       /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 3.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0     /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

#define ADD_POTENTIAL 0         /* set to 1 to add potential to z coordinate */
#define POTENTIAL 10
#define POT_FACT 20.0
/* end of constants only used by sub_wave and sub_maze */


/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters controlling 3D projection */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, 0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {0.0, 8.0, 2.5};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define ADD_DEM 1               /* add DEM (digital elevation model) */
#define ADD_NEGATIVE_DEM 0      /* add DEM with bathymetric data */
#define RSCALE_DEM 0.15          /* scaling factor of radial component for DEM */
#define SMOOTH_DEM 1            /* set to 1 to smoothen DEM (to make altitude less constant) */
#define DEM_SMOOTH_STEPS 5      /* number of smoothening steps */
#define DEM_SMOOTH_HEIGHT 1.0   /* relative height below which to smoothen */
#define DEM_MAXHEIGHT 9860      /* max height of DEM */
#define DEM_MAXDEPTH -10000     /* max depth of DEM */
#define PLANET_SEALEVEL 0.0      /* sea level for flooded planet */
#define VENUS_NODATA_FACTOR 0.5     /* altitude to assign to DEM points without data (fraction of mean altitude) */

#define MAXHEIGHT 5000.0        /* max sea level height */
#define MAXDEPTH 4000.0         /* max sea level depth */

#define RSCALE 0.01             /* scaling factor of radial component */
#define RMAX 10.0               /* max value of radial component */
#define RMIN 0.5               /* min value of radial component */
#define Z_SCALING_FACTOR 0.8    /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.2   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define X_SCALING_FACTOR 2.2   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define Y_SCALING_FACTOR 2.2   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0         /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D 0.0           /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.0           /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.5        /* limit on cosine of normal to shown facets */

```

**Part 2 (2D)** 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 2304           /* number of grid points on x axis */
#define NY 1152           /* number of grid points on y axis */

#define DPOLE 20         /* safety distance to poles */
#define SMOOTHPOLE 0.1     /* smoothing coefficient at poles */
#define ZERO_MERIDIAN 0.0     /* choice of zero meridian (will be at left/right boundary of 2d plot) */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 0        /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 0.25     /* scaling for Julia sets */
#define JULIA_ROT 180.0       /* rotation of Julia set, in degrees */
#define JULIA_RE -0.77145    
#define JULIA_IM -0.10295    /* parameters for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 89         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 33   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 10              /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.0 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.1             /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 7            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 2000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 20.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30            /* number of grid point for grid of disks */
#define NGRIDY 18            /* number of grid point for grid of disks */

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

#define OMEGA 0.001       /* frequency of periodic excitation */
#define AMPLITUDE 0.8     /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0       /* damping of periodic excitation */
#define COURANT 0.05       /* Courant number */
#define COURANTB 0.002    /* Courant number in medium B */
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

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 30    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define WAVE_PACKET_RADIUS 20            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

#define PRECOMPUTE_BC 0     /* set to 1 to compute neighbours for Laplacian in advance */

/* Parameters for length and speed of simulation */

#define NSTEPS 3500       /* number of frames of movie */
#define NVID 4            /* number of iterations between images displayed on screen */
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

/* Parameters of initial condition */

#define INITIAL_AMP 0.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.001  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.05  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define ZPLOT 103     /* wave height */
#define CPLOT 103     /* color scheme */

#define ZPLOT_B 108 
#define CPLOT_B 108        /* plot type for second movie */

#define CHANGE_LUMINOSITY 1     /* set to 1 to let luminosity depend on energy flux intensity */
#define FLUX_WINDOW 30          /* size of averaging window of flux intensity */
#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of plot */
#define SHADE_3D 0              /* set to 1 to change luminosity according to normal vector */
#define SHADE_2D 1              /* set to 1 to change luminosity according to normal vector to plane */
#define SHADE_WAVE 1            /* set to 1 to have luminosity depend on wave height */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define FLOOR_ZCOORD 1          /* set to 1 to draw only facets with z not too negative */
#define DRAW_BILLIARD 0         /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 0   /* set to 1 to draw front of boundary after drawing wave */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw construction lines of certain domains */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental, draw outside of billiard in gray */
#define SHADE_SCALE_2D 10.0     /* controls "depth" of 2D shading */
#define COS_LIGHT_MIN 0.0       /* controls angle-dependence of 2D shading */
#define COS_LIGHT_MAX 0.8       /* controls angle-dependence of 2D shading */

#define PLOT_SCALE_ENERGY 0.4          /* vertical scaling in energy plot */
#define PLOT_SCALE_LOG_ENERGY 0.5       /* vertical scaling in log energy plot */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 
#define PLOT_2D 1               /* switch to 2D representation, equirectangular projection */
#define PHISHIFT 0.0            /* shift of phi in 2D plot (in degrees) */
#define FLOODING 1              /* set to 1 to draw waves above altitude (for Earth representations) */

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 720.0   /* total angle of rotation during simulation */
#define ROTATE_VIEW_WHILE_FADE 0    /* set to 1 to keep rotating viewpoint during fade */

#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */

/* Color schemes */

#define COLOR_PALETTE 11      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 16     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define COLOR_OUT_R 1.0    /* color outside domain */
#define COLOR_OUT_G 1.0    
#define COLOR_OUT_B 1.0    

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.0   /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSCALE_ENERGY 4.0     /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0      /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0    /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 150.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.5      /* shift of colors on log scale */
#define LOG_ENERGY_FLOOR -10.0    /* floor value for log of (total) energy */
#define LOG_MEAN_ENERGY_SHIFT 1.0   /* additional shift for log of mean energy */
#define FLUX_SCALE 1200.0    /* scaling factor for energy flux representation */
#define FLUX_CSCALE 2.0      /* scaling factor for color in energy flux representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 240.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -200.0    /* amplitude of variation of hue for color scheme C_HUE */

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

#define DRAW_COLOR_SCHEME 0       /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 3.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0     /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

#define ADD_POTENTIAL 0         /* set to 1 to add potential to z coordinate */
#define POTENTIAL 10
#define POT_FACT 20.0
/* end of constants only used by sub_wave and sub_maze */


/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters controlling 3D projection */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, 0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {0.0, 8.0, 2.5};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define ADD_DEM 1               /* add DEM (digital elevation model) */
#define ADD_NEGATIVE_DEM 0      /* add DEM with bathymetric data */
#define RSCALE_DEM 0.15          /* scaling factor of radial component for DEM */
#define SMOOTH_DEM 1            /* set to 1 to smoothen DEM (to make altitude less constant) */
#define DEM_SMOOTH_STEPS 5      /* number of smoothening steps */
#define DEM_SMOOTH_HEIGHT 1.0   /* relative height below which to smoothen */
#define DEM_MAXHEIGHT 9860      /* max height of DEM */
#define DEM_MAXDEPTH -10000     /* max depth of DEM */
#define PLANET_SEALEVEL 0.0      /* sea level for flooded planet */
#define VENUS_NODATA_FACTOR 0.5     /* altitude to assign to DEM points without data (fraction of mean altitude) */

#define MAXHEIGHT 5000.0        /* max sea level height */
#define MAXDEPTH 4000.0         /* max sea level depth */

#define RSCALE 0.01             /* scaling factor of radial component */
#define RMAX 10.0               /* max value of radial component */
#define RMIN 0.5               /* min value of radial component */
#define Z_SCALING_FACTOR 0.8    /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.2   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define X_SCALING_FACTOR 2.2   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define Y_SCALING_FACTOR 2.2   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0         /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D 0.0           /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.0           /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.5        /* limit on cosine of normal to shown facets */

```

### 04 December 23 - The Hall effect in a crystal with an increasing magnetic field ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.01	/* x interval for initial condition */
#define INITYMIN -1.1
#define INITYMAX 1.17	/* y interval for initial condition */

#define ADDXMIN 1.9
#define ADDXMAX 2.0	/* x interval for adding particles */
#define ADDYMIN -0.9
#define ADDYMAX 0.9	/* y interval for adding particles */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.25
#define BCYMAX 1.25	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 1   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 6  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.55 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.008 	    /* parameter controlling radius of particles */
#define MU_B 0.018          /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 10         /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */
#define NOBSX 30
#define NOBSY 18            /* obstacles for O_HEX obstacle pattern */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 3800      /* number of frames of movie */
#define NVID 25         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 2   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 5        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 0          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 2        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */

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
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 30000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 3000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 25.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL -50.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.00002          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e10    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 1     /* set to 1 to add an electric field */
#define EFIELD 100000.0      /* value of electric field */
#define ADD_BFIELD 1     /* set to 1 to add a magnetic field */
#define BFIELD 5.33333333333333       /* value of magnetic field */
#define CHARGE -1.0      /* charge of particles of first type */
#define CHARGE_B 2.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 100000.0    /* factor by which to increase electric field */
#define INCREASE_B 1     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */
#define CHARGE_OBSTACLES 1      /* set to 1 for obstacles to be charged */
#define OBSTACLE_CHARGE 3.0     /* charge of obstacles */
#define KCOULOMB_OBSTACLE 100.0   /* Coulomb force constant for charged obstacles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 200.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1000            /* time before end of sim at which to restore size */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 6       /* time interval between adding further particles */
#define N_ADD_PARTICLES 2  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 30   /* size of hashgrid in x direction */
#define HASHY 20    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 03 December 23 - The effect of an obstruction on a two-lenses set-up ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** 

```
    init_wave_flat_comp(phi, psi, xy_in);
    
    if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == 1))
    {
        if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
        add_circular_wave(sign, -1.75, 0.16, phi, psi, xy_in);
        add_circular_wave(sign, -1.75, -0.16, phi, psi, xy_in);
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 1      /* set to 1 for a variable index of refraction */
#define IOR 11               /* choice of index of refraction, see list in global_pdes.c */
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

#define B_DOMAIN 64      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 64    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 13      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 13     /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define MU 0.0              /* parameter controlling the dimensions of domain */
#define MUB 0.4 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 20            /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */

#define OMEGA 0.024        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.063      /* Courant number in medium B */
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

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 20    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define NSOURCES 48         /* number of sources */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 4

/* Parameters for length and speed of simulation */

#define NSTEPS 2400      /* number of frames of movie */
#define NVID 7            /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100    /* time after which to start saving frames */
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

#define INITIAL_AMP 1.0            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00001    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 5

/* Color schemes */

#define COLOR_PALETTE 17     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0     /* scaling factor for energy representation */
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
#define COLORBAR_RANGE 1.5    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define OSCILLATION_SCHEDULE 3  /* oscillation schedule, see list in global_pdes.c */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
/* end of constants only used by sub_wave and sub_maze */

```

### 02 December 23 - Compressing charged molecules: charge density and force field ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1 /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 0  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1600  /* window width */
#define WINHEIGHT 	900   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.01	/* x interval for initial condition */
#define INITYMIN -1.1
#define INITYMAX 1.17	/* y interval for initial condition */

#define ADDXMIN -1.97
#define ADDXMAX -0.8	/* x interval for adding particles */
#define ADDYMIN -1.125
#define ADDYMAX 1.125	/* y interval for adding particles */


#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1  /* pattern of circles, see list in global_ljones.c */

#define ADD_INITIAL_PARTICLES 0 /* set to 1 to add a second type of particles */
#define CIRCLE_PATTERN_B 1  /* pattern of circles for additional particles */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 4  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 0    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 27      /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 3        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 3      /* shape of second rocket */
#define NOZZLE_SHAPE 6        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 6      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 1         /* set to 1 to have two types of particles */
#define TYPE_PROPORTION 0.55 /* proportion of particles of first type */
#define TWOTYPE_CONFIG 0    /* choice of types, see TTC_ list in global_ljones.c */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 12       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 12     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 4.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.8  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.008 	    /* parameter controlling radius of particles */
#define MU_B 0.018          /* parameter controlling radius of particles of second type */
#define NPOLY 40            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define AWEDGE 0.5          /* opening angle of wedge, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 85           /* number of grid point for grid of disks */
#define NGRIDY 43           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 4000      /* number of frames of movie */
#define NVID 25         /* number of iterations between images displayed on screen */
#define NSEG 25         /* number of segments of boundary of circles */
#define INITIAL_TIME 0     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 150     /* time after which to start moving obstacle */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 4   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 5        /* plot type for second movie */

/* Background color depending on particle properties */

#define COLOR_BACKGROUND 1  /* set to 1 to color background */
#define BG_COLOR 2          /* type of background coloring, see list in global_ljones.c */
#define BG_COLOR_B 4        /* type of background coloring, see list in global_ljones.c */

#define DRAW_BONDS 1    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */
#define N_PARTICLE_COLORS 200   /* number of colors for P_NUMBER color scheme */
#define INITIAL_POS_TYPE 0      /* type of initial position dependence */
#define ERATIO 0.995          /* ratio for time-averaging in P_EMEAN color scheme */
#define DRATIO 0.995          /* ratio for time-averaging in P_DIRECT_EMEAN color scheme */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */
#define COLOR_PALETTE_EKIN 10   /* Color palette for kinetic energy */
#define COLOR_PALETTE_ANGLE 10  /* Color palette for angle representation */
#define COLOR_PALETTE_DIRECTION 17   /* Color palette for direction representation */
#define COLOR_PALETTE_INITIAL_POS 0  /* Color palette for initial position representation */
#define COLOR_PALETTE_DIFFNEIGH 10   /* Color palette for different neighbours representation */

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
#define FORCE_FACTOR 0.1            /* factor controlling length of force vector */

/* particle properties */

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 3000.0        /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 60.0      /* hue of particles of type 1 */
#define HUE_TYPE2 140.0      /* hue of particles of type 2 */
#define HUE_TYPE3 200.0     /* hue of particles of type 3 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 50.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 3000.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 5000.0  /* damping coefficient of particles during initial phase */
#define DAMPING_ROT 100.0      /* dampint coefficient for rotation of particles */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 25.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */
#define VICSEK_VMIN 1.0    /* minimal speed of particles in Vicsek model */
#define VICSEK_VMAX 40.0    /* minimal speed of particles in Vicsek model */

#define V_INITIAL_TYPE 0    /* type of initial speed distribution (see VI_ in global_ljones.c) */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.00004          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e10    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.2       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0          /* horizontal gravity acting on all particles */
#define CIRCULAR_GRAVITY 0     /* set to 1 to have gravity directed to center */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */
#define KSPRING_VICSEK 0.2   /* spring constant for I_VICSEK_SPEED interaction */
#define VICSEK_REPULSION 10.0    /* repulsion between particles in Vicsek model */

#define ADD_EFIELD 0     /* set to 1 to add an electric field */
#define EFIELD 50000.0      /* value of electric field */
#define ADD_BFIELD 0     /* set to 1 to add a magnetic field */
#define BFIELD 5.33333333333333       /* value of magnetic field */
#define CHARGE -1.0      /* charge of particles of first type */
#define CHARGE_B 3.0     /* charge of particles of second type */
#define INCREASE_E 0     /* set to 1 to increase electric field */
#define EFIELD_FACTOR 100000.0    /* factor by which to increase electric field */
#define INCREASE_B 0     /* set to 1 to increase magnetic field */
#define BFIELD_FACTOR 1000.0    /* factor by which to increase magnetic field */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0         /* force constant in angular dynamics */
#define KTORQUE_BOUNDARY 1.0e6  /* constant in torque from the boundary */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 10.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 200.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 0   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 0    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 1   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 1000            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.15  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 1        /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only some particles to thermostat */
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.25   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0   /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 10       /* time interval between adding further particles */
#define N_ADD_PARTICLES 2  /* number of particles to add */
#define FINAL_NOADD_PERIOD 0  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define ROTATION_SCHEDULE 0         /* time-dependence of rotation angle, see ROT_* in global_ljones.c */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 300       /* initial time without rotation */
#define ROTATE_FINAL_TIME 300       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.5     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX -2.0*PI              /* maximal rotation speed */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 20   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 0    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */
#define SHOW_SEGMENTS_PRESSURE 1    /* set to 1 to show (averaged) pressure acting on segments */
#define SEGMENT_PMAX 7.5e7        /* pressure of segment with hottest color */
#define P_AVRG_FACTOR 0.02      /* factor in computation of mean pressure */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 500.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 0           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 5.0e11       /* harmonic potential between segment groups */
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

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_REACTION 16          /* type of reaction, see list in global_ljones.c */
#define RD_TYPES 5              /* number of types in reaction-diffusion equation */
#define RD_INITIAL_COND 9       /* initial condition of particles */
#define REACTION_DIST 4.0       /* maximal distance for reaction to occur */
#define REACTION_PROB 1.0        /* probability controlling reaction term */ 
#define DISSOCIATION_PROB 0.002  /* probability controlling dissociation reaction */ 
#define CENTER_COLLIDED_PARTICLES 0  /* set to 1 to recenter particles upon reaction (may interfere with thermostat) */
#define EXOTHERMIC 1            /* set to 1 to make reaction exo/endothermic */
#define DELTA_EKIN 2000.0        /* change of kinetic energy in reaction */
#define COLLISION_TIME 15       /* time during which collisions are shown */

#define CHANGE_RADIUS 0         /* set to 1 to change particle radius during simulation */
#define MU_RATIO 0.666666667    /* ratio by which to increase radius */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PLOT_PARTICLE_NUMBER 0      /* set to 1 to make of plot of particle number over time */
#define PARTICLE_NB_PLOT_FACTOR 0.5 /* expected final number of particles over initial number */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.25    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 12      /* width of maze */
#define NYMAZE 12      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 4        /* seed of random number generator */
#define MAZE_XSHIFT 0.5     /* horizontal shift of maze */
#define MAZE_WIDTH 0.01     /* width of maze walls */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 60   /* size of hashgrid in x direction */
#define HASHY 30    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define LIMIT_ENERGY 0      /* set to 1 to limit energy, when there is no thermostat */

```

### 01 December 23 - Comparison of two-lenses set-ups with different numerical apertures ###

**Program:** `wave_comparison.c` 

**Initial condition in function `animation()`:** 
```
    init_wave_flat_comp(phi, psi, xy_in);
    
    if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == 1))
    {
        if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
        add_circular_wave(sign, -1.75, 0.16, phi, psi, xy_in);
        add_circular_wave(sign, -1.75, -0.16, phi, psi, xy_in);
    }
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 1      /* set to 1 for a variable index of refraction */
#define IOR 10               /* choice of index of refraction, see list in global_pdes.c */
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

#define B_DOMAIN 63      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 63    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 13      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 13     /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */
#define RANDOM_POLY_ANGLE_B 0 /* set to 1 to randomize angle of polygons */

#define XDEP_POLY_ANGLE 0   /* set to 1 to rotate polygons depending on x coordinate */
#define XDEP_POLY_ANGLE_B 0   /* set to 1 to rotate polygons depending on x coordinate */
#define POLY_ROTATION_ANGLE -0.645 /* rotation angle for |x|=1 in units of Pi/2 */
#define HEX_NONUNIF_COMPRESSSION 0.15 /* compression factor for HEX_NONUNIF pattern */
#define HEX_NONUNIF_COMPRESSSION_B -0.15 /* compression factor for HEX_NONUNIF pattern */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define MU 0.7             /* parameter controlling the dimensions of domain */
#define MUB 0.2 	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define APOLY_B 2.0         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 20            /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */

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

#define TWOSPEEDS 1         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */

#define OMEGA 0.024        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.063      /* Courant number in medium B */
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

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 20    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define NSOURCES 48         /* number of sources */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 4

/* Parameters for length and speed of simulation */

#define NSTEPS 2600      /* number of frames of movie */
#define NVID 7            /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 0    /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20   /* number of still frames between movies */
#define END_FRAMES 100   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 1.0            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00001    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 5

/* Color schemes */

#define COLOR_PALETTE 17     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 2.5     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.5    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define OSCILLATION_SCHEDULE 3  /* oscillation schedule, see list in global_pdes.c */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
/* end of constants only used by sub_wave and sub_maze */

```



