/*********************************************************************************/
/*                                                                               */
/*  Animation of interacting particles in a planar domain                        */
/*                                                                               */
/*  N. Berglund, november 2021                                                   */
/*                                                                               */
/*  UPDATE 24/04: distinction between damping and "elasticity" parameters        */
/*  UPDATE 27/04: new billiard shapes, bug in color scheme fixed                 */
/*  UPDATE 28/04: code made more efficient, with help of Marco Mancini           */
/*                                                                               */
/*  Feel free to reuse, but if doing so it would be nice to drop a               */
/*  line to nils.berglund@univ-orleans.fr - Thanks!                              */
/*                                                                               */
/*  compile with                                                                 */
/*  gcc -o lennardjones lennardjones.c                                           */
/* -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp        */
/*                                                                               */
/*  OMP acceleration may be more effective after executing                       */
/*  export OMP_NUM_THREADS=2 in the shell before running the program             */
/*                                                                               */
/*  To make a video, set MOVIE to 1 and create subfolder tif_ljones              */
/*  It may be possible to increase parameter PAUSE                               */
/*                                                                               */
/*  create movie using                                                           */
/*  ffmpeg -i lj.%05d.tif -vcodec libx264 lj.mp4                                 */
/*                                                                               */
/*********************************************************************************/

#include <math.h>
#include <string.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <unistd.h>
#include <sys/types.h>
#include <tiffio.h>     /* Sam Leffler's libtiff library. */
#include <omp.h>
#include <time.h>

#define MOVIE 0        /* set to 1 to generate movie */
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
#define NTREES 15           /* number of trees in S_TREES */

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

#define NO_WRAP_BC ((BOUNDARY_COND != BC_PERIODIC)&&(BOUNDARY_COND != BC_PERIODIC_CIRCLE)&&(BOUNDARY_COND != BC_PERIODIC_TRIANGLE)&&(BOUNDARY_COND != BC_KLEIN)&&(BOUNDARY_COND != BC_PERIODIC_FUNNEL)&&(BOUNDARY_COND != BC_BOY)&&(BOUNDARY_COND != BC_GENUS_TWO))
#define PERIODIC_BC ((BOUNDARY_COND == BC_PERIODIC)||(BOUNDARY_COND == BC_PERIODIC_CIRCLE)||(BOUNDARY_COND == BC_PERIODIC_FUNNEL)||(BOUNDARY_COND == BC_PERIODIC_TRIANGLE))
#define TWO_OBSTACLES ((SEGMENT_PATTERN == S_TWO_CIRCLES_EXT)||(SEGMENT_PATTERN == S_TWO_ROCKETS))
#define COMPUTE_EMEAN ((PLOT == P_EMEAN)||(PLOT_B == P_EMEAN)||(PLOT == P_DIRECT_EMEAN)||(PLOT_B == P_DIRECT_EMEAN))
#define COMPUTE_DIRMEAN ((PLOT == P_DIRECT_EMEAN)||(PLOT_B == P_DIRECT_EMEAN))

double xshift = 0.0;      /* x shift of shown window */
double xspeed = 0.0;      /* x speed of obstacle */
double ylid = 0.9;        /* y coordinate of lid (for BC_RECTANGLE_LID b.c.) */
double vylid = 0.0;       /* y speed coordinate of lid (for BC_RECTANGLE_LID b.c.) */
double xwall = 0.0;       /* x coordinate of wall (for BC_RECTANGLE_WALL b.c.) */
double vxwall = 0.0;      /* x speed of wall (for BC_RECTANGLE_WALL b.c.) */
double angular_speed = 0.0;    /* angular speed of rotating segments */
double xtrack = 0.0;      /* traking coordinate */
double ytrack = 0.0;      /* traking coordinate */
double xsegments[2] = {SEGMENTS_X0, -SEGMENTS_X0};  /* x coordinate of segments (for option MOVE_BOUNDARY) */
double ysegments[2] = {SEGMENTS_Y0, SEGMENTS_Y0};  /* y coordinate of segments (for option MOVE_BOUNDARY) */
double vxsegments[2] = {SEGMENTS_VX0, SEGMENTS_VX0};       /* vx coordinate of segments (for option MOVE_BOUNDARY) */
double vysegments[2] = {SEGMENTS_VY0, SEGMENTS_VY0};       /* vy coordinate of segments (for option MOVE_BOUNDARY) */
int thermostat_on = 1;    /* thermostat switch used when VARY_THERMOSTAT is on */
double cosangle[NSEG+1];      
double sinangle[NSEG+1];      /* precomputed trig functions of angles to draw circles faster */

#define THERMOSTAT_ON ((THERMOSTAT)&&((!VARY_THERMOSTAT)||(thermostat_on)))

#include "global_ljones.c"
#include "sub_maze.c"
#include "sub_hashgrid.c"
#include "sub_lj.c"

FILE *lj_time_series, *lj_final_position;

/*********************/
/* animation part    */
/*********************/


double repel_schedule(int i)
{
    static double kexponent;
    static int first = 1;
    double krepel;
    
    if (first) 
    {
        kexponent = log(KREPEL_FACTOR)/(double)(INITIAL_TIME + NSTEPS);
        first = 0;
    }
    krepel = KREPEL*exp(kexponent*(double)i);
    printf("krepel = %.3lg\n", krepel);
    return(krepel);
}


double efield_schedule(int i)
{
    static double efactor;
    static int first = 1;
    double efield;
    
    if (first) 
    {
        efactor = EFIELD_FACTOR/(double)(INITIAL_TIME + NSTEPS);
        first = 0;
    }
    efield = EFIELD*(double)i*efactor;
    printf("E = %.3lg\n", efield);
    return(efield);
}


double bfield_schedule(int i)
{
    static double bfactor;
    static int first = 1;
    double bfield;
    
    if (first) 
    {
        bfactor = BFIELD_FACTOR/(double)(INITIAL_TIME + NSTEPS);
        first = 0;
    }
    bfield = BFIELD*(double)i*bfactor;
    printf("B = %.3lg\n", bfield);
    return(bfield);
}


double temperature_schedule(int i)
{
    static double bexponent, omega, bexp2;
    static int first = 1, t1, t2, t3;
    double beta;
    
    if (first)
    {
        t1 = NSTEPS - MIDDLE_CONSTANT_PHASE - FINAL_DECREASE_PHASE - FINAL_CONSTANT_PHASE;
        t2 = NSTEPS - FINAL_DECREASE_PHASE - FINAL_CONSTANT_PHASE;
        t3 = NSTEPS - FINAL_CONSTANT_PHASE;
        bexponent = log(BETA_FACTOR)/(double)(t1);
        omega = N_TOSCILLATIONS*DPI/(double)(t1);
        bexp2 = -log(BETA_FACTOR)/(double)(FINAL_DECREASE_PHASE);
        first = 0;
    }
    if (i < INITIAL_TIME) beta = BETA;
    else if (i < INITIAL_TIME + t1)
    {
        beta = BETA*exp(bexponent*(double)(i - INITIAL_TIME));
        if (!NO_OSCILLATION) beta = beta*2.0/(1.0 + cos(omega*(double)(i - INITIAL_TIME)));
    }
    else if (i < INITIAL_TIME + t2) beta = BETA*BETA_FACTOR;
    else if (i < INITIAL_TIME + t3)
    {
        beta = BETA*exp(bexp2*(double)(i - INITIAL_TIME - t3));
    }
    else beta = BETA;
    printf("beta = %.3lg\n", beta);
    return(beta);
}

double container_size_schedule(int i)
{
    if ((i < INITIAL_TIME)||(i > INITIAL_TIME + NSTEPS - RESTORE_TIME)) return(INITXMIN);
    else 
        return(INITXMIN + (1.0-COMPRESSION_RATIO)*(INITXMAX-INITXMIN)*(double)(i-INITIAL_TIME)/(double)(NSTEPS-RESTORE_TIME));
}

double obstacle_schedule_old(int i)
{
    double time;
    static double t1 = 0.5, t2 = 0.75, t3 = 0.875;
    
    if (i < INITIAL_TIME) return(OBSTACLE_XMIN);
    else 
    {
        time = (double)(i-INITIAL_TIME)/(double)(NSTEPS);
        
        if (time < t1) return(OBSTACLE_XMIN + (OBSTACLE_XMAX - OBSTACLE_XMIN)*time/t1);
        else if (time < t2) return(OBSTACLE_XMIN + (OBSTACLE_XMAX - OBSTACLE_XMIN)*(time - t1)/(t2 - t1));
        else if (time < t3) return(OBSTACLE_XMIN + (OBSTACLE_XMAX - OBSTACLE_XMIN)*(time - t2)/(t3 - t2));
        else return(OBSTACLE_XMAX);
    }
}

double obstacle_schedule(int i)
{
    double time, acceleration = 40.0;
    double x;
//     static double t1 = 0.5, t2 = 0.75, t3 = 0.875;
    
    if (i < INITIAL_TIME) return(OBSTACLE_XMIN);
    else 
    {
        time = (double)(i-INITIAL_TIME)/(double)(NSTEPS);
        
        x = OBSTACLE_XMIN + 0.5*acceleration*time*time;
        xspeed = acceleration*time;
//         while (x > OBSXMAX) x += OBSXMIN - OBSXMAX;
        return(x);
    }
}

double obstacle_schedule_smooth(int i, int j)
{
    double time, acceleration = 50.0;
    double x;
//     static double t1 = 0.5, t2 = 0.75, t3 = 0.875;
    
    if (i < INITIAL_TIME) return(OBSTACLE_XMIN);
    else 
    {
        time = ((double)(i-INITIAL_TIME) + (double)j/(double)NVID)/(double)(NSTEPS);
        
        x = OBSTACLE_XMIN + 0.5*acceleration*time*time;
        xspeed = acceleration*time;
//         while (x > OBSXMAX) x += OBSXMIN - OBSXMAX;
        return(x);
    }
}

double gravity_schedule(int i, int j)
{
    double time, gravity, x, y;
    
    switch (GRAVITY_SCHEDULE){
    
        case (G_INCREASE_RELEASE):
        {
            if ((i < INITIAL_TIME + GRAVITY_INITIAL_TIME)||(i > NSTEPS + INITIAL_TIME - GRAVITY_RESTORE_TIME)) return(GRAVITY);
            else 
            {
                time = ((double)(i - INITIAL_TIME - GRAVITY_INITIAL_TIME) 
                + (double)j/(double)NVID)/(double)(NSTEPS - GRAVITY_RESTORE_TIME - GRAVITY_INITIAL_TIME);
                gravity = GRAVITY*(1.0 + time*(GRAVITY_FACTOR - 1.0));
                return(gravity);
            }
            break;
        }
        
        case (G_INCREASE_DECREASE): 
        {
            if ((i < INITIAL_TIME + GRAVITY_INITIAL_TIME)||(i > NSTEPS + INITIAL_TIME - GRAVITY_RESTORE_TIME)) return(GRAVITY);
            else 
            {
                time = ((double)(i - INITIAL_TIME - GRAVITY_INITIAL_TIME) 
                + (double)j/(double)NVID)/(double)(NSTEPS - GRAVITY_RESTORE_TIME - GRAVITY_INITIAL_TIME);
                x = 2.0 - cos(DPI*time);
                y = 0.5*((GRAVITY_FACTOR - 1.0)*x + 3.0 - GRAVITY_FACTOR);
                gravity = GRAVITY*y;
                return(gravity);
            }
            break;
        }
    }
}

double rotation_angle(double phase)
{
    /* case of rotating hourglass */
//     while (phase > DPI) phase -= DPI;
//     return(phase - 0.5*sin(2.0*phase));
    
    /* case of centrifuge */
//     while (phase > 1.0) phase -= 1.0;
//     phase *= DPI;
//     angular_speed = 0.5*OMEGAMAX*(1.0 - cos(phase));
//     return(0.5*OMEGAMAX*(phase - sin(phase)));
    
    /* case of centrifuge remaining at constant speed for a while */
    switch (ROTATION_SCHEDULE) {
        case (ROT_SPEEDUP_SLOWDOWN):
        {
            if (phase < ROTATE_CHANGE_TIME)
            {
                return(0.5*OMEGAMAX*(phase - (ROTATE_CHANGE_TIME/PI)*sin(phase*PI/ROTATE_CHANGE_TIME)));
            }
            else if (phase < 1.0 - ROTATE_CHANGE_TIME)
            {
                return(0.5*OMEGAMAX*(2.0*phase - ROTATE_CHANGE_TIME));
            }
            else
            {
                return(0.5*OMEGAMAX*(2.0 - 2.0*ROTATE_CHANGE_TIME + phase - 1.0 + (ROTATE_CHANGE_TIME/PI)*sin((1.0-phase)*PI/ROTATE_CHANGE_TIME)));
            }
        }
        case (ROT_BACK_FORTH):
        {
            return(OMEGAMAX*(1.0 - cos(DPI*phase))/DPI);
        }
    }
}

double rotation_schedule(int i)
{
    double phase;
    static int imin = INITIAL_TIME + ROTATE_INITIAL_TIME, imax = INITIAL_TIME + NSTEPS - ROTATE_FINAL_TIME;
    
    if (i < imin) 
    {
        angular_speed = 0.0;
        return(0.0);
    }
    else
    {
        if (i > imax) i = imax;
        phase = (DPI/(double)PERIOD_ROTATE_BOUNDARY)*(double)(i - imin);
        return(rotation_angle(phase));
    }
}

double rotation_schedule_smooth(int i, int j)
{
    double phase, angle, phase1, angle1;
    static int imin = INITIAL_TIME + ROTATE_INITIAL_TIME, imax = INITIAL_TIME + NSTEPS - ROTATE_FINAL_TIME;
    
    if (i < imin) 
    {
        angular_speed = 0.0;
        return(0.0);
    }
    else
    {
        if (i > imax) 
        {
            angle = rotation_angle(1.0);
            angular_speed = 0.0;
        }
        else 
        {
            phase = (1.0/(double)(imax - imin))*((double)(i - imin) + (double)j/(double)NVID);
            angle = rotation_angle(phase);
            
            phase1 = (1.0/(double)(imax - imin))*((double)(i + 1 - imin) + (double)j/(double)NVID);
            angle1 = rotation_angle(phase1);
            angular_speed = 25.0*(angle1 - angle);
        }
        
        return(angle);
    }
}

int thermostat_schedule(int i)
{
    if (i < INITIAL_TIME) return(1);
    else return(0);
}

double radius_schedule(int i)
{
    return(1.0 + (MU_RATIO - 1.0)*(double)i/(double)NSTEPS);
}

double evolve_particles(t_particle particle[NMAXCIRCLES], t_hashgrid hashgrid[HASHX*HASHY], 
                        double qx[NMAXCIRCLES], double qy[NMAXCIRCLES], double qangle[NMAXCIRCLES],
                        double px[NMAXCIRCLES], double py[NMAXCIRCLES], double pangle[NMAXCIRCLES], 
                        double beta, int *nactive, int *nsuccess, int *nmove, int *ncoupled, int initial_phase)
{
    double a, totalenergy = 0.0, damping, direction, dmean;  
    static double b = 0.25*SIGMA*SIGMA*DT_PARTICLE/MU_XI, xi = 0.0;
    int j, move, ncoup;
    
    if (initial_phase) damping = INITIAL_DAMPING;
    else damping = DAMPING;
    
    #pragma omp parallel for private(j,xi,totalenergy,a,move)
    for (j=0; j<ncircles; j++) if (particle[j].active)
    {
        particle[j].vx = px[j] + 0.5*DT_PARTICLE*particle[j].fx;
        particle[j].vy = py[j] + 0.5*DT_PARTICLE*particle[j].fy;
        particle[j].omega = pangle[j] + 0.5*DT_PARTICLE*particle[j].torque;
        
        px[j] = particle[j].vx + 0.5*DT_PARTICLE*particle[j].fx;
        py[j] = particle[j].vy + 0.5*DT_PARTICLE*particle[j].fy;
        pangle[j] = particle[j].omega + 0.5*DT_PARTICLE*particle[j].torque;
        
        particle[j].energy = (px[j]*px[j] + py[j]*py[j])*particle[j].mass_inv;
        
        if (COMPUTE_EMEAN)
            particle[j].emean = ERATIO*particle[j].emean + (1.0-ERATIO)*particle[j].energy;
        
        if (COMPUTE_DIRMEAN)
        {
            direction = argument(particle[j].vx, particle[j].vy);
            dmean = particle[j].dirmean;
            if (dmean < direction - PI) dmean += DPI;
            else if (dmean > direction + PI) dmean -= DPI;
            particle[j].dirmean = DRATIO*dmean + (1.0-DRATIO)*direction;
            if (particle[j].dirmean < 0.0) particle[j].dirmean += DPI;
            else if (particle[j].dirmean > DPI) particle[j].dirmean -= DPI;
        }
        
        if ((COUPLE_ANGLE_TO_THERMOSTAT)&&(particle[j].thermostat))
            particle[j].energy += pangle[j]*pangle[j]*particle[j].inertia_moment_inv;
        
        qx[j] = particle[j].xc + 0.5*DT_PARTICLE*px[j]*particle[j].mass_inv;
        qy[j] = particle[j].yc + 0.5*DT_PARTICLE*py[j]*particle[j].mass_inv;
        qangle[j] = particle[j].angle + 0.5*DT_PARTICLE*pangle[j]*particle[j].inertia_moment_inv;
        
        if ((THERMOSTAT_ON)&&(particle[j].thermostat))
        {
            px[j] *= exp(- 0.5*DT_PARTICLE*xi);
            py[j] *= exp(- 0.5*DT_PARTICLE*xi);
        }
        if ((COUPLE_ANGLE_TO_THERMOSTAT)&&(particle[j].thermostat)) 
            pangle[j] *= exp(- 0.5*DT_PARTICLE*xi);
    }
                
    /* compute kinetic energy */
//     *nactive = 0;
    ncoup = 1;
    for (j=0; j<ncircles; j++)
        if ((particle[j].active)&&(particle[j].thermostat))
        {
            totalenergy += particle[j].energy;
            ncoup++;
//             *nactive++;
        }
    totalenergy *= DIMENSION_FACTOR;    /* normalize energy to take number of degrees of freedom into account */
    if (THERMOSTAT_ON)
    {
        /* TODO - fix nactive vs ncoupled */
//         a = DT_PARTICLE*(totalenergy - (double)*nactive/beta)/MU_XI;
        a = DT_PARTICLE*(totalenergy - (double)ncoup/beta)/MU_XI;
        a += SIGMA*sqrt(DT_PARTICLE)*gaussian();
        xi = (xi + a - b*xi)/(1.0 + b);
    }
    
    move = 0;
    for (j=0; j<ncircles; j++) if (particle[j].active) 
    {
        if ((THERMOSTAT_ON)&&(particle[j].thermostat))
        {
            px[j] *= exp(- 0.5*DT_PARTICLE*xi);
            py[j] *= exp(- 0.5*DT_PARTICLE*xi);
            if (!COUPLE_ANGLE_TO_THERMOSTAT) pangle[j] *= exp(- DT_PARTICLE*DAMPING_ROT);
        }
        else 
        {
            px[j] *= exp(- DT_PARTICLE*damping);
            py[j] *= exp(- DT_PARTICLE*damping);
            pangle[j] *= exp(- DT_PARTICLE*DAMPING_ROT);
        }
        if ((THERMOSTAT_ON)&&(COUPLE_ANGLE_TO_THERMOSTAT)&&(particle[j].thermostat))
            pangle[j] *= exp(- 0.5*DT_PARTICLE*xi);
        
        particle[j].xc = qx[j] + 0.5*DT_PARTICLE*px[j]*particle[j].mass_inv;
        particle[j].yc = qy[j] + 0.5*DT_PARTICLE*py[j]*particle[j].mass_inv;
        particle[j].angle = qangle[j] + 0.5*DT_PARTICLE*pangle[j]*particle[j].inertia_moment_inv;
        
//         particle[j].vx = px[j] + 0.5*DT_PARTICLE*particle[j].fx;
//         particle[j].vy = py[j] + 0.5*DT_PARTICLE*particle[j].fy;
//         particle[j].omega = pangle[j] + 0.5*DT_PARTICLE*particle[j].torque;

        /* TO DO: move this to function wrap_particle */
        if ((BOUNDARY_COND == BC_PERIODIC_CIRCLE)||(BOUNDARY_COND == BC_PERIODIC_FUNNEL)||(BOUNDARY_COND == BC_PERIODIC_TRIANGLE))
        {
//                       if (particle[j].xc < BCXMIN) 
            if (particle[j].xc < xshift + BCXMIN) 
            {
                particle[j].xc += BCXMAX - BCXMIN;
                if (RESAMPLE_Y) 
                {
                    *nmove++;
                    if (resample_particle(j, NTRIALS, particle) == 1) 
                    {
                        px[j] = particle[j].vx;
                        py[j] = particle[j].vy;
                        update_hashgrid(particle, hashgrid, 0);
                        *nsuccess++;
                    }
                }
            }
//                       else if (particle[j].xc > BCXMAX) 
            else if (particle[j].xc > xshift + BCXMAX) 
            {
                particle[j].xc += BCXMIN - BCXMAX;
            }
            if (particle[j].yc > BCYMAX) particle[j].yc += BCYMIN - BCYMAX;
            else if (particle[j].yc < BCYMIN) particle[j].yc += BCYMAX - BCYMIN;
        }
        else if (!NO_WRAP_BC)
        {
            move += wrap_particle(&particle[j], &px[j], &py[j]);
        }
//                 if (move > 0) 
//                 {
//                     compute_relative_positions(particle, hashgrid);
//                     update_hashgrid(particle, hashgrid, 0);   /* REDUNDANT ? */
//                 }
    }
    
    *ncoupled = ncoup;
    return(totalenergy);
}


void evolve_lid(double fboundary)
{
    double force;
    
    force = fboundary - GRAVITY*LID_MASS;
    if (ylid > BCYMAX + LID_WIDTH) force -= KSPRING_BOUNDARY*(ylid - BCYMAX - LID_WIDTH);
    vylid += force*DT_PARTICLE/LID_MASS;
    ylid += vylid*DT_PARTICLE;
}

void evolve_wall(double fboundary)
{
    double force;
    
    force = fboundary;
    if (xwall > BCYMAX - WALL_WIDTH) force -= KSPRING_BOUNDARY*(xwall - BCYMAX + WALL_WIDTH);
    else if (xwall < BCYMIN + WALL_WIDTH) force += KSPRING_BOUNDARY*(BCYMIN + WALL_WIDTH - xwall);
    
    force -= vxwall*WALL_FRICTION;
    
    vxwall += fboundary*DT_PARTICLE/WALL_MASS;
    
    if (vxwall > WALL_VMAX) vxwall = WALL_VMAX;
    else if (vxwall < -WALL_VMAX) vxwall = -WALL_VMAX;
    
    xwall += vxwall*DT_PARTICLE;
//     printf("fboundary = %.3lg, xwall = %.3lg, vxwall = %.3lg\n", fboundary, xwall, vxwall);
}


void evolve_segments(t_segment segment[NMAXSEGMENTS], int time)
{
    int i, nactive = 0, group;
    double fx[2] = {0.0, 0.0}, fy[2] = {0.0, 0.0}, x, y, padding = 3.0*MU, mass2 = SEGMENTS_MASS;
    
    if (SEGMENT_PATTERN == S_TWO_CIRCLES_EXT) mass2 = SEGMENTS_MASS*TWO_CIRCLES_RADIUS_RATIO;
    
    for (group=0; group<2; group++)
    {
        fx[group] = 0.0;
        fy[group] = 0.0;
    }
    for (i=0; i<nsegments; i++) if (segment[i].active)
    {
        group = segment[i].group;
        fx[group] += segment[i].fx;
        fy[group] += segment[i].fy;
        nactive++;
        if (BOUNDARY_COND == BC_SCREEN) /* add force from simulation boundary */
        {
            x = 0.5*(segment[i].x1 + segment[i].x2);
            y = 0.5*(segment[i].y1 + segment[i].y2);
            if (x < XMIN + padding) fx[group] += KSPRING_BOUNDARY*(XMIN + padding - x);
            else if (x > XMAX - padding) fx[group] -= KSPRING_BOUNDARY*(x - XMAX + padding);
            if (y < YMIN + padding) fy[group] += KSPRING_BOUNDARY*(YMIN + padding - y);
            else if (y > YMAX - padding) fy[group] -= KSPRING_BOUNDARY*(y - YMAX + padding);
        }
        else if ((BOUNDARY_COND == BC_REFLECT_ABS)||(BOUNDARY_COND == BC_REFLECT_ABS_BOTTOM)) 
        /* add force from simulation boundary */
        {
             y = 0.5*(segment[i].y1 + segment[i].y2);
             if (y < YMIN) fy[group] += KSPRING_BOUNDARY*(YMIN - y);
        }
        if (group == 0) fy[group] -= GRAVITY*SEGMENTS_MASS;
        else fy[group] -= GRAVITY*mass2;
    }
    if (nactive > 0) for (group=0; group<2; group++)
    {
        fx[group] = fx[group]/(double)nactive;
        fy[group] = fy[group]/(double)nactive;
    }
    if (FLOOR_FORCE) 
    {   
        if (fx[0] > FMAX) fx[0] = FMAX;
        else if (fx[0] < -FMAX) fx[0] = -FMAX;
        if (fy[0] > FMAX) fy[0] = FMAX;
        else if (fy[0] < -FMAX) fy[0] = -FMAX;
    }
    vxsegments[0] += fx[0]*DT_PARTICLE/SEGMENTS_MASS;
    vysegments[0] += fy[0]*DT_PARTICLE/SEGMENTS_MASS;
    xsegments[0] += vxsegments[0]*DT_PARTICLE;
    ysegments[0] += vysegments[0]*DT_PARTICLE;
    if (TWO_OBSTACLES)
    {
        if (FLOOR_FORCE)
        {
            if (fx[1] > FMAX) fx[1] = FMAX;
            else if (fx[1] < -FMAX) fx[1] = -FMAX;
            if (fy[1] > FMAX) fy[1] = FMAX;
            else if (fy[1] < -FMAX) fy[1] = -FMAX;
        }
        vxsegments[1] += fx[1]*DT_PARTICLE/mass2;
        vysegments[1] += fy[1]*DT_PARTICLE/mass2;
        xsegments[1] += vxsegments[1]*DT_PARTICLE;
        ysegments[1] += vysegments[1]*DT_PARTICLE;
    }
    
    /* add some damping if y coordinate is small (for lunar landing) */
    if (DAMP_SEGS_AT_NEGATIVE_Y)
        for (group=0; group<2; group++) 
            if (ysegments[group] < 0.1) 
            {
                vysegments[group] *= exp(-DAMPING*DT_PARTICLE);
                vxsegments[group] *= exp(-DAMPING*DT_PARTICLE);
            }

    /* to avoid numerical instabilities */
    for (group=0; group<2; group++) 
    {
        if (xsegments[group] + 1.0 > BCXMAX) 
        {
            xsegments[group] = BCXMAX - 1.0;
            vxsegments[group] = 0.0;
        }
        if ((RELEASE_ROCKET_AT_DEACTIVATION)&&((BOUNDARY_COND == BC_REFLECT_ABS)||(BOUNDARY_COND == BC_ABSORBING))) 
        {
//             ysegments[group] = SEGMENTS_Y0;
            if (time < SEGMENT_DEACTIVATION_TIME) vysegments[group] = 0.0;
            else if ((ysegments[group] < SEGMENTS_Y0)&&(vysegments[group] < 0.0)) 
                vysegments[group] = -0.5*vysegments[group];
        }
    }
        
    translate_segments(segment, xsegments, ysegments);
}


void evolve_segment_groups(t_segment segment[NMAXSEGMENTS], int time, t_group_segments segment_group[NMAXGROUPS])
/* new version of evolve_segments that takes the group structure into account */
{
    double fx[NMAXGROUPS], fy[NMAXGROUPS], torque[NMAXGROUPS], dx[NMAXGROUPS], dy[NMAXGROUPS], dalpha[NMAXGROUPS];
    double x, y, dx0, dy0, padding, proj, distance, f, xx[2], yy[2], xmean = 0.0, ymean = 0.0, ymax = 0.0;
    int i, j, k, group = 0;
    static double maxdepth, saturation_depth, xmax;
    static int first = 1;
    
    if (first)
    {
        xmax = XMAX - TRACK_X_PADDING;
        if ((PLOT_SPEEDS)||(PLOT_TRAJECTORIES)) xmax -= 1.8;
        first = 0;
    }
    
    maxdepth = 0.5*GROUP_WIDTH;
    saturation_depth = 0.1*GROUP_WIDTH;
    padding = 0.1;
    
    for (group=0; group<ngroups; group++)
    {
        fx[group] = 0.0;
        fy[group] = 0.0;
        torque[group] = 0.0;
    }
    
    /* only groups of segments of index 1 or larger are mobile */
    for (i=0; i<nsegments; i++) if ((segment[i].active)&&(segment[i].group > 0))
    {
        group = segment[i].group;
        
        fx[group] += segment[i].fx;
        fy[group] += segment[i].fy;
        torque[group] += segment[i].torque;
        
        dx0 = segment[i].xc - segment_group[group].xc;
        dy0 = segment[i].yc - segment_group[group].yc;
        torque[group] += dx0*segment[i].fy - dy0*segment[i].fx;
                
        if (BOUNDARY_COND == BC_SCREEN) /* add force from simulation boundary */
        {
            x = 0.5*(segment[i].x1 + segment[i].x2);
            y = 0.5*(segment[i].y1 + segment[i].y2);
            if (x < XMIN + padding) fx[group] += KSPRING_BOUNDARY*(XMIN + padding - x);
            else if (x > XMAX - padding) fx[group] -= KSPRING_BOUNDARY*(x - XMAX + padding);
            if (y < YMIN + padding) fy[group] += KSPRING_BOUNDARY*(YMIN + padding - y);
            else if (y > YMAX - padding) fy[group] -= KSPRING_BOUNDARY*(y - YMAX + padding);
        }
        else if ((BOUNDARY_COND == BC_REFLECT_ABS)||(BOUNDARY_COND == BC_REFLECT_ABS_BOTTOM)) 
        /* add force from simulation boundary */
        {
             y = 0.5*(segment[i].y1 + segment[i].y2);
             if (y < YMIN) fy[group] += KSPRING_BOUNDARY*(YMIN - y);
        }
        
        /* repulsion between different groups */
        if (GROUP_REPULSION) for (j=0; j<nsegments; j++) if ((segment[j].active)&&(segment[j].group != group))
        {
            xx[0] = segment[j].x1;
            yy[0] = segment[j].y1;
            xx[1] = segment[j].x2;
            yy[1] = segment[j].y2;
            for (k=0; k<2; k++)
            {
                x = xx[k];
                y = yy[k];
                proj = (segment[i].ny*(x - segment[i].x1) - segment[i].nx*(y - segment[i].y1))/segment[i].length;
                if ((proj > 0.0)&&(proj < 1.0))
                {
                    distance = segment[i].nx*x + segment[i].ny*y - segment[i].c;
                    if ((distance > -maxdepth)&&(distance < 0.0))
                    {
                        if (distance < -saturation_depth) distance = -saturation_depth;
                        f = KSPRING_GROUPS*(-distance);
                        segment[j].fx += f*segment[i].nx;
                        segment[j].fy += f*segment[i].ny;
                        segment[j].torque += (x - segment[i].xc)*f*segment[i].ny - (y - segment[i].yc)*f*segment[i].nx;
                    
                        fx[group] -= f*segment[i].nx;
                        fy[group] -= f*segment[i].ny;
                        torque[group] -= (x - segment[i].xc)*f*segment[i].ny - (y - segment[i].yc)*f*segment[i].nx;
                    }
                }
            }
        }
    }
    
    if (GROUP_G_REPEL) for (i=0; i<ngroups; i++) for (j=i+1; j<ngroups; j++) 
    {
        x = segment_group[j].xc - segment_group[i].xc;
        y = segment_group[j].yc - segment_group[i].yc;
        distance = module2(x, y);
        
        if (distance < GROUP_G_REPEL_RADIUS)
        {
            if (distance < 0.1*GROUP_G_REPEL_RADIUS) distance = 0.1*GROUP_G_REPEL_RADIUS;
            f = KSPRING_GROUPS*(GROUP_G_REPEL_RADIUS - distance);
            fx[j] += f*x/distance;
            fy[j] += f*y/distance;
            fx[i] -= f*x/distance;
            fy[i] -= f*y/distance;
        }
    }
        
    if (FLOOR_FORCE) for (group=1; group<ngroups; group++)
    {   
        if (fx[group] > FMAX) fx[group] = FMAX;
        else if (fx[group] < -FMAX) fx[group] = -FMAX;
        if (fy[group] > FMAX) fy[group] = FMAX;
        else if (fy[group] < -FMAX) fy[group] = -FMAX;
    }
    
    for (group=1; group<ngroups; group++)
    {
        fy[group] -= GRAVITY*segment_group[group].mass;
        fx[group] += GRAVITY_X*segment_group[group].mass;
        
        segment_group[group].vx += fx[group]*DT_PARTICLE/segment_group[group].mass;
        segment_group[group].vy += fy[group]*DT_PARTICLE/segment_group[group].mass;
        segment_group[group].omega += torque[group]*DT_PARTICLE/segment_group[group].moment_inertia;

        segment_group[group].vx *= exp(- DT_PARTICLE*SEGMENT_GROUP_DAMPING);
        segment_group[group].vy *= exp(- DT_PARTICLE*SEGMENT_GROUP_DAMPING);
        segment_group[group].omega *= exp(- DT_PARTICLE*SEGMENT_GROUP_DAMPING);

        dx[group] = segment_group[group].vx*DT_PARTICLE;
        dy[group] = segment_group[group].vy*DT_PARTICLE;
        dalpha[group] = segment_group[group].omega*DT_PARTICLE;
        
        segment_group[group].xc += dx[group];
        segment_group[group].yc += dy[group];
        segment_group[group].angle += dalpha[group];
        
//         printf("group %i: (dx, dy) = (%.3lg, %.3lg)\n", group, dx[group], dy[group]);
    }
    
        
    for (i=0; i<nsegments; i++) if ((segment[i].active)&&(segment[i].group > 0))
    {
        group = segment[i].group;
        
        translate_one_segment(segment, i, dx[group], dy[group]);
        rotate_one_segment(segment, i, dalpha[group],  segment_group[group].xc,  segment_group[group].yc);
    }
    
    if (TRACK_SEGMENT_GROUPS) 
    {
        /* compute mean position */
        for (group=1; group<ngroups; group++)
        {
            xmean += segment_group[group].xc;
            ymean += segment_group[group].yc;
            if (segment_group[group].yc > ymax) ymax = segment_group[group].yc;
        }
        xmean = xmean/((double)(ngroups-1));
        ymean = ymean/((double)(ngroups-1));
        
        /* bias towards ymax */
        ymean = 0.75*ymax + 0.25*ymean;
        
        if (ymean > ytrack) ytrack = ymean;
        if (xmean > xmax) 
            xtrack = xmean - xmax;
        else if (xmean < XMIN + TRACK_X_PADDING) 
            xtrack = xmean - XMIN - TRACK_X_PADDING;
    }
}


void animation()
{
    double time, scale, diss, rgb[3], dissip, gradient[2], x, y, dx, dy, dt, xleft, xright, 
            a, b, length, fx, fy, force[2], totalenergy = 0.0, pos[2], prop, vx, xi = 0.0, torque, torque_ij, pleft = 0.0, pright = 0.0, entropy[2], speed_ratio, ymin, ymax, delta_energy, speed, ratio = 1.0, ratioc, cum_etot = 0.0, emean = 0.0, radius_ratio, t;
    double *qx, *qy, *px, *py, *qangle, *pangle, *pressure, *obstacle_speeds;
    int i, j, k, n, m, s, ij[2], i0, iplus, iminus, j0, jplus, jminus, p, q, p1, q1, p2, q2, total_neighbours = 0, 
        min_nb, max_nb, close, wrapx = 0, wrapy = 0, nadd_particle = 0, nmove = 0, nsuccess = 0, 
        tracer_n[N_TRACER_PARTICLES], traj_position = 0, traj_length = 0, move = 0, old, m0, floor, nthermo, wall = 0,
        group, gshift, n_total_active = 0, ncollisions = 0, ncoupled = 1;
    int *particle_numbers;
    static int imin, imax;
    static short int first = 1;
    t_particle *particle;
    t_obstacle *obstacle;
    t_segment *segment;
    t_group_segments *segment_group;     
    t_tracer *trajectory;
    t_group_data *group_speeds;
    t_collision *collisions;
    t_hashgrid *hashgrid;
    t_lj_parameters params;
    char message[100];

    ratioc = 1.0 - ratio;
    
    /* parameter values, grouped in a structure to simplify parameter printing */
    params.beta = BETA;
    params.krepel = KREPEL; 
    params.xmincontainer = BCXMIN;
    params.xmaxcontainer = BCXMAX; 
    params.fboundary = 0.0;
    params.gravity = GRAVITY;
    params.radius = MU;
    
    particle = (t_particle *)malloc(NMAXCIRCLES*sizeof(t_particle));    /* particles */  
    if (ADD_FIXED_OBSTACLES) obstacle = (t_obstacle *)malloc(NMAXOBSTACLES*sizeof(t_obstacle));    /* circular obstacles */  
    if (ADD_FIXED_SEGMENTS) 
    {
        segment = (t_segment *)malloc(NMAXSEGMENTS*sizeof(t_segment));         /* linear obstacles */  
        segment_group = (t_group_segments *)malloc(NMAXGROUPS*sizeof(t_group_segments));
    }

    if (TRACER_PARTICLE) 
        trajectory = (t_tracer *)malloc(TRAJECTORY_LENGTH*N_TRACER_PARTICLES*sizeof(t_tracer));
    
    hashgrid = (t_hashgrid *)malloc(HASHX*HASHY*sizeof(t_hashgrid));    /* hashgrid */      
            
    qx = (double *)malloc(NMAXCIRCLES*sizeof(double));
    qy = (double *)malloc(NMAXCIRCLES*sizeof(double));
    px = (double *)malloc(NMAXCIRCLES*sizeof(double));
    py = (double *)malloc(NMAXCIRCLES*sizeof(double));
    qangle = (double *)malloc(NMAXCIRCLES*sizeof(double));
    pangle = (double *)malloc(NMAXCIRCLES*sizeof(double));
    pressure = (double *)malloc(N_PRESSURES*sizeof(double));
    
    if (REACTION_DIFFUSION) 
    {
        collisions = (t_collision *)malloc(2*NMAXCOLLISIONS*sizeof(t_collision));
        for (i=0; i<2*NMAXCOLLISIONS; i++) collisions[i].time = 0;
    }
    
    if (SAVE_TIME_SERIES) 
    {
        lj_time_series = fopen("lj_time_series.dat", "w");
        lj_final_position = fopen("lj_final_position.dat", "w");
    }
    
    if (ADD_FIXED_OBSTACLES) init_obstacle_config(obstacle);
    if (ADD_FIXED_SEGMENTS) init_segment_config(segment);
    
    if ((MOVE_SEGMENT_GROUPS)&&(ADD_FIXED_SEGMENTS)) 
    {
        for (i=0; i<ngroups; i++) init_segment_group(segment, i, segment_group);
        group_speeds = (t_group_data *)malloc(ngroups*(INITIAL_TIME + NSTEPS)*sizeof(t_group_data));
    }
    
    /* initialise array of trig functions to speed up drawing particles */
    init_angles();
    
    /* initialise positions and radii of circles */
    init_particle_config(particle);
    
    /* add some particles, beta */
    if (ADD_INITIAL_PARTICLES) add_particle_config(particle, -0.6, 1.6, -1.0, 1.0, MU_B);
    
    init_hashgrid(hashgrid);
    
    xshift = OBSTACLE_XMIN;
    speed_ratio = (double)(25*NVID)*DT_PARTICLE;
    
    if (RECORD_PRESSURES) for (i=0; i<N_PRESSURES; i++) pressure[i] = 0.0;
    if (PLOT_SPEEDS) obstacle_speeds = (double *)malloc(2*ngroups*(INITIAL_TIME + NSTEPS)*sizeof(double));
    if (PLOT_PARTICLE_NUMBER) 
        particle_numbers = (int *)malloc((NSTEPS+1)*(RD_TYPES+1)*sizeof(int));
    
    
    
    printf("Initializing configuration\n");

    params.nactive = initialize_configuration(particle, hashgrid, obstacle, px, py, pangle, tracer_n, segment);
    
    printf("%i active particles\n", params.nactive);
     
//     xi = 0.0;
    
//     for (i=0; i<ncircles; i++)
//     {
//         printf("Particle %i at (%.3f, %.3f) of energy %.3f\n", i, particle[i].xc, particle[i].yc, particle[i].energy);
//     }
    sleep(1);
        
    update_hashgrid(particle, hashgrid, 1);
    printf("Updated hashgrid\n");
    compute_relative_positions(particle, hashgrid);
    printf("Computed relative positions\n");
    blank();
//     glColor3f(0.0, 0.0, 0.0);

    glutSwapBuffers();
    
    sleep(SLEEP1);

    for (i=0; i<=INITIAL_TIME + NSTEPS; i++)
    {
        printf("Computing frame %d\n",i);
        
        if (INCREASE_KREPEL) params.krepel = repel_schedule(i);
        if (INCREASE_BETA) params.beta = temperature_schedule(i);  
        if (INCREASE_E) params.efield = efield_schedule(i);
        if (INCREASE_B) params.bfield = bfield_schedule(i);
        if (DECREASE_CONTAINER_SIZE) 
        {
            params.xmincontainer = container_size_schedule(i);
            if (SYMMETRIC_DECREASE) params.xmaxcontainer = -container_size_schedule(i);
        }
        if ((ROTATE_BOUNDARY)&&(!SMOOTH_ROTATION)) rotate_segments(segment, rotation_schedule(i));
        if (VARY_THERMOSTAT) 
        {
            thermostat_on = thermostat_schedule(i);
            printf("Termostat: %i\n", thermostat_on);
        }
        
         
        /* deactivate some segments */
        if ((ADD_FIXED_SEGMENTS)&&(DEACTIVATE_SEGMENT)&&(i == INITIAL_TIME + SEGMENT_DEACTIVATION_TIME + 1))
            for (j=0; j<nsegments; j++) if (segment[j].inactivate) segment[j].active = 0;
            
        /* recolor particles in case if P_INITIAL_POS color code */
        if ((i <= INITIAL_TIME-1)&&(i%10 == 0)&&((PLOT == P_INITIAL_POS)||(PLOT_B == P_INITIAL_POS)))
        {
            printf("Recoloring particles\n"); 
            ymin = particle[0].yc;
            ymax = particle[0].yc;
            for (j=1; j<ncircles; j++) if (particle[j].active)
            {
                if (particle[j].yc < ymin) ymin = particle[j].yc;
                if (particle[j].yc > ymax) ymax = particle[j].yc;
            }
            for (j=0; j<ncircles; j++) if (particle[j].active)
                particle[j].color_hue = 360.0*(particle[j].yc - ymin)/(ymax - ymin);
        }
	
	blank();
        
        params.fboundary = 0.0;
        pleft = 0.0;
        pright = 0.0;
        if (RECORD_PRESSURES) for (j=0; j<N_PRESSURES; j++) pressure[j] = 0.0;
        
//         printf("evolving particles\n");
        for(n=0; n<NVID; n++)
        {
            if (MOVE_OBSTACLE) 
            {
                params.xmincontainer = obstacle_schedule_smooth(i, n);
                xshift = params.xmincontainer;
            }
            if ((ROTATE_BOUNDARY)&&(SMOOTH_ROTATION)) 
                rotate_segments(segment, rotation_schedule_smooth(i,n));
            if (ROTATE_BOUNDARY) 
            {
                params.omega = angular_speed;
                params.angle = rotation_schedule_smooth(i,n);
            }
            
            if (INCREASE_GRAVITY) params.gravity = gravity_schedule(i,n); 
            if ((BOUNDARY_COND == BC_RECTANGLE_WALL)&&(i < INITIAL_TIME + WALL_TIME)) wall = 1;
            else wall = 0;
            
            if ((MOVE_BOUNDARY)||(MOVE_SEGMENT_GROUPS)||(PRINT_SEGMENTS_FORCE)) for (j=0; j<nsegments; j++)
            {
                segment[j].fx = 0.0;
                segment[j].fy = 0.0;
                segment[j].torque = 0.0;
            }
            
            compute_relative_positions(particle, hashgrid);
            update_hashgrid(particle, hashgrid, 0);
            
            /* compute forces on particles */
            for (j=0; j<ncircles; j++) if (particle[j].active)
            {
                particle[j].fx = 0.0;
                particle[j].fy = 0.0;
                particle[j].torque = 0.0;
                
                /* compute force from other particles */
                compute_particle_force(j, params.krepel, particle, hashgrid);
                
                /* take care of boundary conditions */
                params.fboundary += compute_boundary_force(j, particle, obstacle, segment, params.xmincontainer, params.xmaxcontainer, &pleft, &pright, pressure, wall, params.krepel);

                /* align velocities in case of Vicsek models */
//                 if (VICSEK_INT) 
                if ((VICSEK_INT)&&(!particle[j].close_to_boundary))
                {
                    speed = module2(particle[j].vx,particle[j].vy);
                    if ((VICSEK_VMIN > 0.0)&&(speed < VICSEK_VMIN)) speed = VICSEK_VMIN;
                    if ((INTERACTION == I_VICSEK_SHARK)&&(particle[j].type == 2)) speed *= 1.75;
                    if (speed > VICSEK_VMAX) speed = 0.5*(speed + VICSEK_VMAX);
                    particle[j].vx = speed*cos(particle[j].angle);
                    particle[j].vy = speed*sin(particle[j].angle);
                    
                    speed = module2(px[j],py[j]);
                    if ((VICSEK_VMIN > 0.0)&&(speed < VICSEK_VMIN)) speed = VICSEK_VMIN;
                    if ((INTERACTION == I_VICSEK_SHARK)&&(particle[j].type == 2)) speed *= 1.75;
                    if (speed > VICSEK_VMAX) speed = 0.5*(speed + VICSEK_VMAX);
                    px[j] = speed*cos(particle[j].angle);
                    py[j] = speed*sin(particle[j].angle);
                    
                    
                }

                /* add gravity */
                if (INCREASE_GRAVITY) particle[j].fy -= params.gravity/particle[j].mass_inv;
                else if (CIRCULAR_GRAVITY)
                {
                    particle[j].fx -= GRAVITY*particle[j].xc/particle[j].mass_inv;
                    particle[j].fy -= GRAVITY*particle[j].yc/particle[j].mass_inv;
                }
                else 
                {
                    particle[j].fy -= GRAVITY/particle[j].mass_inv;
                    particle[j].fx += GRAVITY_X/particle[j].mass_inv;
                }
                
                /* add electric force */
                if (ADD_EFIELD) 
                {
                    if (INCREASE_E) particle[j].fx += params.efield*particle[j].charge;
                    else particle[j].fx += EFIELD*particle[j].charge;
                }
                
                /* add magnetic force */
                if (ADD_BFIELD)
                {
                    particle[j].fx += params.bfield*particle[j].charge*particle[j].vy*particle[j].mass_inv;
                    particle[j].fy -= params.bfield*particle[j].charge*particle[j].vx*particle[j].mass_inv;
                }
                
                if (FLOOR_FORCE)
                {
                    if (particle[j].fx > FMAX) particle[j].fx = FMAX;
                    if (particle[j].fx < -FMAX) particle[j].fx = -FMAX;
                    if (particle[j].fy > FMAX) particle[j].fy = FMAX;
                    if (particle[j].fy < -FMAX) particle[j].fy = -FMAX;
                    if (particle[j].torque > FMAX) particle[j].torque = FMAX;
                    if (particle[j].torque < -FMAX) particle[j].torque = -FMAX;
                }
            }
            
            /* timestep of thermostat algorithm */
            totalenergy = evolve_particles(particle, hashgrid, qx, qy, qangle, px, py, pangle, params.beta, &params.nactive, &nsuccess, &nmove, &ncoupled, i < INITIAL_TIME);
            
            
            /* evolution of lid coordinate */
            if (BOUNDARY_COND == BC_RECTANGLE_LID) evolve_lid(params.fboundary);
            if (BOUNDARY_COND == BC_RECTANGLE_WALL) 
            {
                if (i < INITIAL_TIME + WALL_TIME) evolve_wall(params.fboundary);
                else xwall = 0.0;
            }
            if ((MOVE_BOUNDARY)&&(i > OBSTACLE_INITIAL_TIME)) evolve_segments(segment, i);
            
            if ((MOVE_SEGMENT_GROUPS)&&(i > INITIAL_TIME + SEGMENT_DEACTIVATION_TIME)) evolve_segment_groups(segment, i, segment_group);
            
//             if ((MOVE_SEGMENT_GROUPS)&&(i > OBSTACLE_INITIAL_TIME)) evolve_segment_groups(segment, i, segment_group);
        } /* end of for (n=0; n<NVID; n++) */
        
        if ((i>INITIAL_TIME)&&(SAVE_TIME_SERIES))
        {
            n_total_active = 0;
            for (j=0; j<ncircles; j++) if (particle[j].active) n_total_active++;
            fprintf(lj_time_series, "%i\n", n_total_active);
            for (j=0; j<ncircles; j++) if (particle[j].active)
            {
                fprintf(lj_time_series, "%.4f\n", particle[j].xc);
                fprintf(lj_time_series, "%.4f\n", particle[j].yc);
            }
        }

//         printf("evolved particles\n");

        if (PLOT_SPEEDS) /* record speeds of segments */
        {
            gshift = INITIAL_TIME + NSTEPS;
            if (MOVE_SEGMENT_GROUPS) for (group = 1; group < ngroups; group++)
            {
                group_speeds[(group-1)*gshift + i].xc = segment_group[group].xc;
                group_speeds[(group-1)*gshift + i].yc = segment_group[group].yc;
                group_speeds[(group-1)*gshift + i].vx = segment_group[group].vx*speed_ratio;
                group_speeds[(group-1)*gshift + i].vy = segment_group[group].vy*speed_ratio;
                group_speeds[(group-1)*gshift + i].omega = segment_group[group].omega*speed_ratio;
            }
            else
            {
                obstacle_speeds[i] = vysegments[0];
                obstacle_speeds[INITIAL_TIME + NSTEPS + i] = vysegments[1];
            }
        }
        
        if (MOVE_BOUNDARY) 
            printf("segment[%i]: (fx, fy) = (%.3lg, %.3lg), torque = %.3lg)\n", i, fx, fy, torque);
        
        if (MOVE_SEGMENT_GROUPS) for (group=1; group<ngroups; group++)
            printf("segments position [%i] (%.3lg, %.3lg) angle %.3lg\n speed (%.3lg, %.3lg) omega %.3lg\n", 
                   group, segment_group[group].xc, segment_group[group].yc, segment_group[group].angle, segment_group[group].vx, segment_group[group].vy, segment_group[group].omega);
        
//         if ((PARTIAL_THERMO_COUPLING)) 
        if ((PARTIAL_THERMO_COUPLING)&&(i>N_T_AVERAGE)) 
        {
            nthermo = partial_thermostat_coupling(particle, xshift + PARTIAL_THERMO_SHIFT, segment);
            printf("%i particles coupled to thermostat out of %i active\n", nthermo, params.nactive);
            params.mean_energy = compute_mean_energy(particle);
        }
        else params.mean_energy = totalenergy/(double)ncircles;
        
        if (CENTER_PX) center_momentum(px);
        if (CENTER_PY) center_momentum(py);
        if (CENTER_PANGLE) center_momentum(pangle);
        if (FLOOR_OMEGA) floor = floor_momentum(pangle);
            
//         printf("pressure left %.5lg, pressure right %.5lg\n", pleft, pright);
//         for (j=0; j<N_PRESSURES; j++) printf("pressure[%i] = %.5lg\n", j, pressure[j]);
        
        /* reset angular values to [0, 2 Pi) */
        if (ROTATION) for (j=0; j<ncircles; j++) 
        {
            while (particle[j].angle > DPI) particle[j].angle -= DPI;
            while (particle[j].angle < 0.0) particle[j].angle += DPI;
        }
        
        /* update tracer particle trajectory */
        if ((TRACER_PARTICLE)&&(i > INITIAL_TIME)) 
        {
            for (j=0; j<N_TRACER_PARTICLES; j++)
            {
                trajectory[j*TRAJECTORY_LENGTH + traj_position].xc = particle[tracer_n[j]].xc;
                trajectory[j*TRAJECTORY_LENGTH + traj_position].yc = particle[tracer_n[j]].yc;
            }
            traj_position++;
            if (traj_position >= TRAJECTORY_LENGTH) traj_position = 0;
            traj_length++;
            if (traj_length >= TRAJECTORY_LENGTH) traj_length = TRAJECTORY_LENGTH - 1;
//             for (j=0; j<traj_length; j++)
//                 printf("Trajectory[%i] = (%.3lg, %.3lg)\n", j, trajectory[j].xc, trajectory[j].yc);
        }
                
        printf("Mean kinetic energy: %.3f\n", totalenergy/(double)ncircles); 
        printf("Kinetic energy by coupled particle: %.3f\n",  totalenergy/(double)ncoupled);
        
        if ((!THERMOSTAT)&&(LIMIT_ENERGY))
        {
            if (cum_etot > 0.0)
            {
                emean = cum_etot/(double)(i+1);
                if (totalenergy > 10.0*emean)
                {
                    reset_energy(particle, px, py, totalenergy, emean);
                    totalenergy = 0.0;
                    for (j=0; j<ncircles; j++) if (particle[j].active) 
                        totalenergy += particle[i].energy;
                    printf("Reset mean kinetic energy: %.3f\n", totalenergy/(double)ncircles);
                }
            }
            printf("Emean: %.3f\n", emean/(double)ncircles);
            cum_etot += totalenergy;
        }
        
        printf("Boundary force: %.3f\n", params.fboundary/(double)(ncircles*NVID)); 
        if (RESAMPLE_Y) printf("%i succesful moves out of %i trials\n", nsuccess, nmove);
        if (INCREASE_GRAVITY) printf("Gravity: %.3f\n", params.gravity);
                
        total_neighbours = 0;
        min_nb = 100;
        max_nb = 0;
        for (j=0; j<ncircles; j++) if (particle[j].active)
        {
            total_neighbours += particle[j].neighb;
            if (particle[j].neighb > max_nb) max_nb = particle[j].neighb; 
            if (particle[j].neighb < min_nb) min_nb = particle[j].neighb; 
        }
//         printf("Mean number of neighbours: %.3f\n", (double)total_neighbours/(double)ncircles); 
//         printf("Min number of neighbours: %i\n", min_nb); 
//         printf("Max number of neighbours: %i\n", max_nb); 
        
        blank();
        /* case of reaction-diffusion equation */
        if ((i > INITIAL_TIME)&&(REACTION_DIFFUSION)) 
        {
            ncollisions = update_types(particle, collisions, ncollisions, particle_numbers, i - INITIAL_TIME - 1,  &delta_energy);
            if (EXOTHERMIC) params.beta *= 1.0/(1.0 + delta_energy/totalenergy);
            params.nactive = 0;
            for (j=0; j<ncircles; j++) if (particle[j].active)
            {
                params.nactive++;
                qx[j] = particle[j].xc;
                qy[j] = particle[j].yc;
                px[j] = particle[j].vx;
                py[j] = particle[j].vy;
            }
//             draw_collisions(collisions, ncollisions);
        }
        
        /* case of varying type proportion */
        if (CHANGE_TYPES) 
        {
            t = (double)i/(double)NSTEPS;
            params.prop = PROP_MIN*(1.0-t) + PROP_MAX*t;
            change_type_proportion(particle, params.prop);
        }
 
        if (TRACER_PARTICLE) draw_trajectory(trajectory, traj_position, traj_length);
        draw_particles(particle, PLOT, params.beta, collisions, ncollisions, BG_COLOR, hashgrid);
        draw_container(params.xmincontainer, params.xmaxcontainer, obstacle, segment, wall);

        /* add a particle */
        if ((ADD_PARTICLES)&&(i > ADD_TIME)&&((i - INITIAL_TIME - ADD_TIME)%ADD_PERIOD == 1)&&(i < NSTEPS - FINAL_NOADD_PERIOD))
        {
            for (k=0; k<N_ADD_PARTICLES; k++)
            nadd_particle = add_particles(particle, px, py, nadd_particle);
//             params.nactive = nadd_particle;
            params.nactive = 0;
            for (j=0; j<ncircles; j++) if (particle[j].active) params.nactive++;
        }
        
        /* change particle radius */
        if (CHANGE_RADIUS)
        {
            radius_ratio = radius_schedule(i+1)/radius_schedule(i);
            printf("Particle radius factor %.5lg\t", radius_schedule(i+1));
            for (j=0; j<ncircles; j++) particle[j].radius *= radius_ratio;
            printf("Particle 0 radius %.5lg\n", particle[0].radius);
            params.radius *= radius_ratio; 
        }
        
        /* compute force on segments */
        if (PRINT_SEGMENTS_FORCE) compute_segments_force(&params, segment);
        
        update_hashgrid(particle, hashgrid, 1);
        
        if (PRINT_PARAMETERS) print_parameters(params, PRINT_LEFT, pressure, 1);
        if ((BOUNDARY_COND == BC_EHRENFEST)||(BOUNDARY_COND == BC_RECTANGLE_WALL)) 
            print_ehrenfest_parameters(particle, pleft, pright);
        else if (PRINT_PARTICLE_NUMBER) 
        {
            if (REACTION_DIFFUSION) 
                print_particle_types_number(particle, RD_TYPES);
            else print_particle_number(ncircles);
        }
        
        if ((i > INITIAL_TIME + WALL_TIME)&&(PRINT_ENTROPY)) 
        {
            compute_entropy(particle, entropy);
            printf("Entropy 1 = %.5lg, Entropy 2 = %.5lg\n", entropy[0], entropy[1]);
            print_entropy(entropy);
        }
        
        if (PLOT_SPEEDS) draw_speed_plot(group_speeds, i);
        if (PLOT_TRAJECTORIES) draw_trajectory_plot(group_speeds, i);
        
        else if (PRINT_PARTICLE_SPEEDS) print_particles_speeds(particle);
        else if (PRINT_SEGMENTS_SPEEDS)
        {
            if (MOVE_BOUNDARY) print_segments_speeds(vxsegments, vysegments);
            else print_segment_group_speeds(segment_group);
        }
        if ((i > INITIAL_TIME)&&(PLOT_PARTICLE_NUMBER))
        {
            count_particle_number(particle, particle_numbers, i - INITIAL_TIME);
            draw_particle_nb_plot(particle_numbers, i - INITIAL_TIME);
        }
        
	if (!((NO_EXTRA_BUFFER_SWAP)&&(MOVIE))) glutSwapBuffers();

	if (MOVIE)
        {
            if (i >= INITIAL_TIME) 
            {
                if ((TIME_LAPSE)&&(TIME_LAPSE_FIRST))
                {
                    if ((TIME_LAPSE)&&((i - INITIAL_TIME)%TIME_LAPSE_FACTOR == 0)&&(!DOUBLE_MOVIE))
                    {
                        save_frame_lj();
                    }
                    save_frame_lj_counter(NSTEPS/TIME_LAPSE_FACTOR + MID_FRAMES + i - INITIAL_TIME);
                }
                else
                {
                    save_frame_lj();
                    if ((TIME_LAPSE)&&((i - INITIAL_TIME)%TIME_LAPSE_FACTOR == 0)&&(!DOUBLE_MOVIE))
                    {
                        save_frame_lj_counter(NSTEPS + END_FRAMES + (i - INITIAL_TIME)/TIME_LAPSE_FACTOR);
                    }
                }
            }
            else printf("Initial phase time %i of %i\n", i, INITIAL_TIME);
            
            
            if ((i >= INITIAL_TIME)&&(DOUBLE_MOVIE))
            {
                if (TRACER_PARTICLE) draw_trajectory(trajectory, traj_position, traj_length);
                draw_particles(particle, PLOT_B, params.beta, collisions, ncollisions, BG_COLOR_B, hashgrid);
                draw_container(params.xmincontainer, params.xmaxcontainer, obstacle, segment, wall);
                if (PRINT_PARAMETERS) print_parameters(params, PRINT_LEFT, pressure, 0);
                if (PLOT_SPEEDS) draw_speed_plot(group_speeds, i);
                if (PLOT_TRAJECTORIES) draw_trajectory_plot(group_speeds, i);
                if (BOUNDARY_COND == BC_EHRENFEST) print_ehrenfest_parameters(particle, pleft, pright);
                else if (PRINT_PARTICLE_NUMBER) 
                {
                    if (REACTION_DIFFUSION) 
                        print_particle_types_number(particle, RD_TYPES);
                    else print_particle_number(ncircles);
                }
                else if (PRINT_PARTICLE_SPEEDS) print_particles_speeds(particle);
                else if (PRINT_SEGMENTS_SPEEDS) print_segment_group_speeds(segment_group);
//                     print_segments_speeds(vxsegments, vysegments);
                glutSwapBuffers();
                save_frame_lj_counter(NSTEPS + MID_FRAMES + 1 + counter);
                counter++;
            }
            else if (NO_EXTRA_BUFFER_SWAP) glutSwapBuffers();

            /* it seems that saving too many files too fast can cause trouble with the file system */
            /* so this is to make a pause from time to time - parameter PAUSE may need adjusting   */
            if (i % PAUSE == PAUSE - 1)
            {
                printf("Making a short pause\n");
                sleep(PSLEEP);
                s = system("mv lj*.tif tif_ljones/");
            }
        }

    }
    
    if (SAVE_TIME_SERIES)
    {
        n_total_active = 0;
        for (j=0; j<ncircles; j++) if (particle[j].active) n_total_active++;
        fprintf(lj_final_position, "%i\n", n_total_active);
        for (j=0; j<ncircles; j++) if (particle[j].active)
        {
            fprintf(lj_final_position, "%i\n", j);
            fprintf(lj_final_position, "%.4f\n", particle[j].xc);
            fprintf(lj_final_position, "%.4f\n", particle[j].yc);
        }
    }

    if (MOVIE) 
    {
        if (DOUBLE_MOVIE) 
        {
            blank();
            if (TRACER_PARTICLE) draw_trajectory(trajectory, traj_position, traj_length);
            draw_particles(particle, PLOT, params.beta, collisions, ncollisions, BG_COLOR, hashgrid); 
            draw_container(params.xmincontainer, params.xmaxcontainer, obstacle, segment, wall);
            if (PRINT_PARAMETERS) print_parameters(params, PRINT_LEFT, pressure, 1);
            if (PLOT_SPEEDS) draw_speed_plot(group_speeds, i);
            if (PLOT_TRAJECTORIES) draw_trajectory_plot(group_speeds, i);
            if (BOUNDARY_COND == BC_EHRENFEST) print_ehrenfest_parameters(particle, pleft, pright);
            else if (PRINT_PARTICLE_NUMBER) 
            {
                if (REACTION_DIFFUSION) 
                    print_particle_types_number(particle, RD_TYPES);
                else print_particle_number(ncircles);
            }
            else if (PRINT_PARTICLE_SPEEDS) print_particles_speeds(particle);
            else if (PRINT_SEGMENTS_SPEEDS) print_segment_group_speeds(segment_group);
//                 print_segments_speeds(vxsegments, vysegments);
//             glutSwapBuffers();
        }
        for (i=0; i<MID_FRAMES; i++) 
        {
            save_frame_lj();
            if (!NO_EXTRA_BUFFER_SWAP) glutSwapBuffers();
        }
        glutSwapBuffers();
        if (DOUBLE_MOVIE) 
        {
            if (TRACER_PARTICLE) draw_trajectory(trajectory, traj_position, traj_length);
            draw_particles(particle, PLOT_B, params.beta, collisions, ncollisions, BG_COLOR_B, hashgrid);
            draw_container(params.xmincontainer, params.xmaxcontainer, obstacle, segment, wall);
            if (PRINT_PARAMETERS) print_parameters(params, PRINT_LEFT, pressure, 0);
            if (PLOT_SPEEDS) draw_speed_plot(group_speeds, i);
            if (PLOT_TRAJECTORIES) draw_trajectory_plot(group_speeds, i);
            if (BOUNDARY_COND == BC_EHRENFEST) print_ehrenfest_parameters(particle, pleft, pright);
            else if (PRINT_PARTICLE_NUMBER) 
            {
                if (REACTION_DIFFUSION) 
                    print_particle_types_number(particle, RD_TYPES);
                else print_particle_number(ncircles);
            }
            else if (PRINT_PARTICLE_SPEEDS) print_particles_speeds(particle);
            else if (PRINT_SEGMENTS_SPEEDS) print_segment_group_speeds(segment_group);
//                 print_segments_speeds(vxsegments, vysegments);
            if (!((NO_EXTRA_BUFFER_SWAP)&&(MOVIE))) glutSwapBuffers();
        }
        if ((TIME_LAPSE)&&(!DOUBLE_MOVIE)) 
        {
            for (i=0; i<END_FRAMES; i++) 
                save_frame_lj_counter(NSTEPS + MID_FRAMES + NSTEPS/TIME_LAPSE_FACTOR + i);
        }
        else for (i=0; i<END_FRAMES; i++) save_frame_lj_counter(NSTEPS + MID_FRAMES + 1 + counter + i);
        
        s = system("mv lj*.tif tif_ljones/");
    }
    
    params.nactive = 0;
    for (j=0; j<ncircles; j++) if (particle[j].active) params.nactive++;
    printf("%i active particles\n", params.nactive);  
    
//     printf("1\n");
    free(particle);
//     printf("2\n");
    if (ADD_FIXED_OBSTACLES) free(obstacle);
//     printf("3\n");
    if (ADD_FIXED_SEGMENTS) 
    {
        free(segment);
        free(segment_group);
    }
//     printf("4\n");
    if (MOVE_SEGMENT_GROUPS) free(group_speeds);
//     printf("5\n");
    if (TRACER_PARTICLE) free(trajectory);
//     printf("6\n");
    if (PLOT_SPEEDS) free(obstacle_speeds);
//     printf("7\n");
    if (PLOT_PARTICLE_NUMBER) free(particle_numbers);
//     printf("8\n");
    free(hashgrid);
//     printf("9\n");
    free(qx);
//     printf("10\n");
    free(qy);
//     printf("11\n");
    free(px);
//     printf("12\n");
    free(py);
//     printf("13\n");
    free(qangle);
//     printf("14\n");
    free(pangle);
//     printf("15\n");
    free(pressure);
//     printf("16\n");
    if (REACTION_DIFFUSION) free(collisions);
    
    if (SAVE_TIME_SERIES) 
    {
        fclose(lj_time_series);
        fclose(lj_final_position);
    }
}


void display(void)
{
    time_t rawtime;
    struct tm * timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    
    glPushMatrix();

    blank();
    glutSwapBuffers();
    blank();
    glutSwapBuffers();

    animation();
    sleep(SLEEP2);

    glPopMatrix();

    glutDestroyWindow(glutGetWindow());

    printf("Start local time and date: %s", asctime(timeinfo));
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    printf("Current local time and date: %s", asctime(timeinfo));
}


int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(WINWIDTH,WINHEIGHT);
    glutCreateWindow("Particles with Lennard-Jones interaction in a planar domain");

    init();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}

