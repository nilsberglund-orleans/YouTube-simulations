### Parameter values for YouTube simulations ###

Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Below are parameter values used for different simulations, as well as initial conditions used in 
function animation. Some simulations use variants of the published code. The list is going to be 
updated gradually. 



### 13 March 22 - Colliding vortices in the Rock-Paper-Scissors reaction-diffusion equation ###

**Program:** `bz.c` 

**Initial condition in function `animation()`:** `init_random(0.5, 0.4, phi, psi, chi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 1920          /* number of grid points on x axis */
#define NY 1000          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 99    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define MU 0.8	            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MRATIO 5            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15            /* number of grid point for grid of disks */
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
/* xy_in_billiard and draw_billiard in sub_wave.c */

/* Physical patameters of wave equation */

#define DT 0.0000005

#define VISCOSITY 0.05

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */

#define EPSILON 0.8     /* time scale separation */
#define DELTA 0.1       /* time scale separation */
#define BZQ 0.0008       /* parameter in BZ equation */
#define BZF 1.2          /* parameter in BZ equation */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.1     /* noise intensity */

#define CHANGE_VISCOSITY 1      /* set to 1 to change the viscosity in the course of the simulation */
#define ADJUST_INTSTEP 1        /* set to 1 to decrease integration step when viscosity increases */
#define VISCOSITY_INITIAL_TIME 100  /* initial time during which viscosity remains constant */
#define VISCOSITY_FACTOR 200.0   /* factor by which to change viscosity */
#define VISCOSITY_MAX 2.0        /* max value of viscosity beyond which NVID is increased and integration step is decrase, 
                                    for numerical stability */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 1500      /* number of frames of movie */
#define NVID 10           /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */

#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define INITIAL_TIME 5  /* initial still time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 50   /* number of still frames at end of movie */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Plot type */

#define PLOT 21

#define PLOT_B 25        /* plot type for second movie */

#define P_RGB 20        /* RGB plot */
#define P_POLAR 21      /* polar angle associated with RBG plot */
#define P_GRADIENT 22   /* gradient of polar angle */
#define P_GRADIENTX 23  /* direction of gradient of u */
#define P_GRADIENT_INTENSITY 24  /* gradient and intensity of polar angle */
#define P_CURL 25       /* curl of polar angle */

#define PRINT_TIME 0    /* set to 1 to print running time */
#define PRINT_VISCOSITY 1    /* set to 1 to print v */

#define DRAW_FIELD_LINES 0  /* set to 1 to draw field lines */
#define FIELD_LINE_WIDTH 1  /* width of field lines */
#define N_FIELD_LINES 120   /* number of field lines */
#define FIELD_LINE_FACTOR 120 /* factor controlling precision when computing origin of field lines */

/* Color schemes, see list in global_pdes.c  */

// #define COLOR_PALETTE 11     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE 14     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 4   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 1.0   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 0.00001   /* scaling factor for curl representation */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0      /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.0 

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 0.55    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 12 March 22 - Interacting pentagons on the projective plane ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 2    /* factor of time-lapse movie */


/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.85
#define INITXMAX 1.85	/* x interval for initial condition */
#define INITYMIN -0.9
#define INITYMAX 0.9	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 1   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 0  /* pattern of obstacles, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.8 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 0  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */

#define INTERACTION 3       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 1.8   /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.045 	    /* parameter controlling radius of particles */
#define MU_B 0.02427051     /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 3           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 3750      /* number of frames of movie */
#define NVID 100          /* number of iterations between images displayed on screen */
#define NSEG 150         /* number of segments of boundary */
#define INITIAL_TIME 0    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 4   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 13

/* Plot type, see list in global_ljones.c  */

#define PLOT 4
#define PLOT_B 3        /* plot type for second movie */

#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */

/* Color schemes */

#define COLOR_PALETTE 0     /* Color palette, see list in global_ljones.c  */

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

/* particle properties */

#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 5.0e2           /* energy of particle with hottest color */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 10.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 6.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 0.1   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.02     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.002            /* initial inverse temperature */
#define MU_XI 0.05           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e8    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e8    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.5       /* radius in which to count neighbours */
#define GRAVITY 0.0         /* gravity acting on all particles */

#define ROTATION 1           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 600.0         /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 7.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 1  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 500.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define FINAL_CONSTANT_PHASE 0  /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 700            /* time before end of sim at which to restore size */

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
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 25        /* time at which to add first particle */
#define ADD_PERIOD 20       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 250  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define TRAJECTORY_LENGTH 6000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 0.1    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define FLOOR_FORCE 0      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */

#define HASHX 30   /* size of hashgrid in x direction */
#define HASHY 20   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 11 March 22 - Interacting dipoles on the projective plane ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 2    /* factor of time-lapse movie */


/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.85
#define INITXMAX 1.85	/* x interval for initial condition */
#define INITYMIN -0.9
#define INITYMAX 0.9	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 1   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 0  /* pattern of obstacles, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.8 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 0  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */

#define INTERACTION 5       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 2.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 1.8   /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.045 	    /* parameter controlling radius of particles */
#define MU_B 0.02427051     /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 3           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 2700      /* number of frames of movie */
#define NVID 200          /* number of iterations between images displayed on screen */
#define NSEG 150         /* number of segments of boundary */
#define INITIAL_TIME 0    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 4   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 13

/* Plot type, see list in global_ljones.c  */

#define PLOT 4
#define PLOT_B 3        /* plot type for second movie */

#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */

/* Color schemes */

#define COLOR_PALETTE 0     /* Color palette, see list in global_ljones.c  */

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

/* particle properties */

#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 5.0e2           /* energy of particle with hottest color */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 5.0e-7    /* time step for particle displacement */
#define KREPEL 10.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 6.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 0.1   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.02     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.002            /* initial inverse temperature */
#define MU_XI 0.05           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e8    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e8    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.5       /* radius in which to count neighbours */
#define GRAVITY 0.0         /* gravity acting on all particles */

#define ROTATION 1           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 600.0         /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 7.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 1  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 500.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define FINAL_CONSTANT_PHASE 0  /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 700            /* time before end of sim at which to restore size */

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
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 25        /* time at which to add first particle */
#define ADD_PERIOD 20       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 250  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define TRAJECTORY_LENGTH 6000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 0.1    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define FLOOR_FORCE 0      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */

#define HASHX 30   /* size of hashgrid in x direction */
#define HASHY 20   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 10 March 22 - Direction of rotation in the Rock-Paper-Scissors reaction-diffusion equation ###

**Program:** `bz.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);`

```
#define MOVIE 1        /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 1920          /* number of grid points on x axis */
#define NY 1000          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 99    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define MU 0.8	            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MRATIO 5            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15            /* number of grid point for grid of disks */
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
/* xy_in_billiard and draw_billiard in sub_wave.c */

/* Physical patameters of wave equation */

#define DT 0.0000005

#define VISCOSITY 1.0
// #define GAMMA_INVERSE 1.0

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */

#define EPSILON 0.8     /* time scale separation */
#define DELTA 0.1       /* time scale separation */
#define BZQ 0.0008       /* parameter in BZ equation */
#define BZF 1.2          /* parameter in BZ equation */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.1     /* noise intensity */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 2500      /* number of frames of movie */
#define NVID 10           /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */

#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define INITIAL_TIME 5  /* initial still time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 50   /* number of still frames at end of movie */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Plot type */

#define PLOT 23

#define PLOT_B 22        /* plot type for second movie */

#define P_RGB 20        /* RGB plot */
#define P_POLAR 21      /* polar angle associated with RBG plot */
#define P_GRADIENT 22   /* gradient of polar angle */
#define P_GRADIENTX 23  /* gradient of polar angle */

#define PRINT_TIME 0    /* set to 1 to print running time */

#define DRAW_FIELD_LINES 0  /* set to 1 to draw field lines */
#define FIELD_LINE_WIDTH 1  /* width of field lines */
#define N_FIELD_LINES 120   /* number of field lines */
#define FIELD_LINE_FACTOR 120 /* factor controlling precision when computing origin of field lines */

/* Color schemes, see list in global_pdes.c  */

#define COLOR_PALETTE 18     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 4   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 1.0   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 2.0       /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0      /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.0 

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 0.55    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 9 March 22 - Crystal formation on the projective plane, with smaller particles ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 2    /* factor of time-lapse movie */


/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.85
#define INITXMAX 1.85	/* x interval for initial condition */
#define INITYMIN -0.9
#define INITYMAX 0.9	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 1   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 0  /* pattern of obstacles, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.8 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 0  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.75  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.02 	    /* parameter controlling radius of particles */
#define MU_B 0.02427051     /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 3           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 2400      /* number of frames of movie */
#define NVID 200          /* number of iterations between images displayed on screen */
#define NSEG 150         /* number of segments of boundary */
#define INITIAL_TIME 0    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 4   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 13

/* Plot type, see list in global_ljones.c  */

#define PLOT 3
#define PLOT_B 0        /* plot type for second movie */

#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */

/* particle properties */

#define PARTICLE_HUE_MIN 330.0      /* color of original particle */
#define PARTICLE_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_EMAX 5.0e2           /* energy of particle with hottest color */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 5.0e-7    /* time step for particle displacement */
#define KREPEL 10.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 6.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 0.1   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.25     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.002            /* initial inverse temperature */
#define MU_XI 0.05           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e8    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e8    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 5.0       /* radius in which to count neighbours */
#define GRAVITY 0.0         /* gravity acting on all particles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 100.0         /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 1  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 500.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define FINAL_CONSTANT_PHASE 0  /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 700            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.1  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 25        /* time at which to add first particle */
#define ADD_PERIOD 20       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 250  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define TRAJECTORY_LENGTH 6000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 0.1    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define FLOOR_FORCE 0      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */

#define HASHX 35   /* size of hashgrid in x direction */
#define HASHY 25   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 8 March 22 - Crystal formation on the projective plane ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 2    /* factor of time-lapse movie */


/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.85
#define INITXMAX 1.85	/* x interval for initial condition */
#define INITYMIN -0.9
#define INITYMAX 0.9	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 1   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 0  /* pattern of obstacles, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.8 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 0  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.75  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.025 	    /* parameter controlling radius of particles */
#define MU_B 0.02427051     /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 3           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 2400      /* number of frames of movie */
#define NVID 200          /* number of iterations between images displayed on screen */
#define NSEG 150         /* number of segments of boundary */
#define INITIAL_TIME 100    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 4   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 13

/* Plot type, see list in global_ljones.c  */

#define PLOT 0
#define PLOT_B 3        /* plot type for second movie */

#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */

/* particle properties */

#define PARTICLE_HUE_MIN 330.0      /* color of original particle */
#define PARTICLE_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_EMAX 5.0e2           /* energy of particle with hottest color */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 5.0e-7    /* time step for particle displacement */
#define KREPEL 10.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 0.1   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.25     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.002            /* initial inverse temperature */
#define MU_XI 0.05           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e8    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e8    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 5.0       /* radius in which to count neighbours */
#define GRAVITY 0.0         /* gravity acting on all particles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 100.0         /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 1  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 500.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 0        /* set to 1 to have exponential BETA change only */
#define FINAL_CONSTANT_PHASE 0  /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 700            /* time before end of sim at which to restore size */

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
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 25        /* time at which to add first particle */
#define ADD_PERIOD 20       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 250  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define TRAJECTORY_LENGTH 6000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 0.1    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define FLOOR_FORCE 0      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */

#define HASHX 30   /* size of hashgrid in x direction */
#define HASHY 20   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 7 March 22 - Comparison of a Fresnel lens and a circular lens ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat_comp(phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 1920          /* number of grid points on x axis */
#define NY 1000          /* number of grid points on y axis */
#define YMID 500        /* mid point of display */
#define XMIN -0.5
#define XMAX 3.5	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 43      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 47    /* choice of domain shape, see list in global_pdes.c */

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

#define LAMBDA -1.1	    /* parameter controlling the dimensions of domain */
#define MU 0.1   	    /* parameter controlling the dimensions of domain */
#define MUB 0.1 	    /* parameter controlling the dimensions of domain */
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
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */

#define OMEGA 0.004        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define COURANT 0.02       /* Courant number */
#define COURANTB 0.01154     /* Courant number in medium B */
#define GAMMA 0.0      /* damping factor in wave equation */
#define GAMMAB 0.0          /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-6      /* damping factor on boundary */
#define KAPPA 0.0       /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4       /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0       /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 0

/* Parameters for length and speed of simulation */

#define NSTEPS 2500      /* number of frames of movie */
#define NVID 25          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 20    /* time after which to start saving frames */
#define COMPUTE_ENERGIES 0  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20   /* number of still frames between movies */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0003  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.02  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 1

/* Color schemes */

#define COLOR_PALETTE 18     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define BLACK_TEXT 1     /* set to 1 to write text in black instead of white */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 200.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.5     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0     /* shift of colors on log scale */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.5    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 2.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

```

### 6 March 22 - Mixing of two Lennard-Jones fluids ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 2    /* factor of time-lapse movie */


/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -1.125
#define INITYMAX 1.125	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 0  /* pattern of obstacles, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.8 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 0  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 4.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.015 	    /* parameter controlling radius of particles */
#define MU_B 0.02427051     /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 3           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 3900      /* number of frames of movie */
#define NVID 300          /* number of iterations between images displayed on screen */
#define NSEG 150         /* number of segments of boundary */
#define INITIAL_TIME 0    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 4   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 0

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 6        /* plot type for second movie */

#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */

/* particle properties */

#define PARTICLE_HUE_MIN 330.0      /* color of original particle */
#define PARTICLE_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_EMAX 1.5e3           /* energy of particle with hottest color */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 0.1   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.25     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.001            /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e8    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e8    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 5.0       /* radius in which to count neighbours */
#define GRAVITY 0.0         /* gravity acting on all particles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 100.0         /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 0.1   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define FINAL_CONSTANT_PHASE 0  /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 700            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 1   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 2.0  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 25        /* time at which to add first particle */
#define ADD_PERIOD 20       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 250  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define TRAJECTORY_LENGTH 6000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 0.1    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define POSITION_DEPENDENT_TYPE 1   /* set to 1 to make particle type depend on initial position */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */

#define HASHX 25   /* size of hashgrid in x direction */
#define HASHY 15   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 5 March 22 - A circular lens ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define HIGHRES 1        /* set to 1 if resolution of grid is double that of displayed image */

#define WINWIDTH 	720   /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 1440          /* number of grid points on x axis */
#define NY 1440          /* number of grid points on y axis */

#define XMIN -0.425
#define XMAX 1.825	/* x interval  */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 47       /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 201   /* pattern of circles or polygons, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA -1.1	    /* parameter controlling the dimensions of domain */
#define MU 0.2              /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.3333333333333333           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 36           /* number of grid point for grid of disks */
#define NGRIDY 6           /* number of grid point for grid of disks */

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

#define OMEGA 0.005        /* frequency of periodic excitation */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define COURANT 0.03       /* Courant number */
#define COURANTB 0.015      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0          /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 100    /* period of oscillating source */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3

/* Parameters for length and speed of simulation */

#define NSTEPS 1300        /* number of frames of movie */
#define NVID 30          /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 1    /* width of billiard boundary */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100    /* number of still frames at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00025  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.015  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 0        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 18     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 150.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0     /* shift of colors on log scale */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 4 March 22 - A Brownian particle in a Lennard-Jones fluid ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.9
#define INITXMAX 1.9	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 0  /* pattern of obstacles, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.8 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 0  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.015 	    /* parameter controlling radius of particles */
#define MU_B 0.02427051     /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 3           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 6000      /* number of frames of movie */
#define NVID 375          /* number of iterations between images displayed on screen */
#define NSEG 150         /* number of segments of boundary */
#define INITIAL_TIME 50    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 4   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 3

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 6        /* plot type for second movie */

#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */

/* particle properties */

#define PARTICLE_HUE_MIN 330.0      /* color of original particle */
#define PARTICLE_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_EMAX 1.5e3           /* energy of particle with hottest color */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 20.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 0.1   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.25     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.001            /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e8    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e8    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 5.0       /* radius in which to count neighbours */
#define GRAVITY 0.0         /* gravity acting on all particles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 100.0         /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 0.1   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define FINAL_CONSTANT_PHASE 0  /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 700            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 1   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 2.0  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 25        /* time at which to add first particle */
#define ADD_PERIOD 20       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 250  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 1   /* set to 1 to have a tracer particle */
#define TRAJECTORY_LENGTH 6000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 0.1    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */

#define HASHX 25   /* size of hashgrid in x direction */
#define HASHY 15   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 3 March 22 - Onset of spirals in the stochastic FitzHugh-Nagumo equations ###

**Program:** `fhn.c` 

**Initial condition in function `animation()`:** `init_random(0.0, 1.0, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 1920          /* number of grid points on x axis */
#define NY 1000          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 99    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define MU 0.8	            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MRATIO 5            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15            /* number of grid point for grid of disks */
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
/* xy_in_billiard and draw_billiard in sub_wave.c */

/* Physical patameters of wave equation */

#define DT 0.00000015

#define VISCOSITY 5.0
#define GAMMA_INVERSE 1.0

#define EPSILON 0.1 /* time scale separation */
#define FHNA 1.0    /* parameter in FHN equation */
#define FHNC -0.01    /* parameter in FHN equation */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 1     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.04     /* noise intensity */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 5000      /* number of frames of movie */
#define NVID 2          /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */

#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define INITIAL_TIME 25  /* initial still time */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Plot type */

#define PLOT 11

#define PRINT_TIME 0    /* set to 1 to print running time */

#define DRAW_FIELD_LINES 0  /* set to 1 to draw field lines */
#define FIELD_LINE_WIDTH 1  /* width of field lines */
#define N_FIELD_LINES 120   /* number of field lines */
#define FIELD_LINE_FACTOR 120 /* factor controlling precision when computing origin of field lines */

/* Color schemes, see list in global_pdes.c  */

#define COLOR_PALETTE 18     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 4   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 1.0   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0      /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.0 

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.0    /* scale of color scheme bar */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 2 March 22 - Bernoulli’s principle, illustrated by a Lennard-Jones gas moving through a funnel ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.9
#define INITXMAX 2.9	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define BCXMIN -3.0
#define BCXMAX 3.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 0  /* pattern of obstacles, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.8 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 0  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.012 	    /* parameter controlling radius of particles */
#define MU_B 0.02427051     /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 3           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 3200      /* number of frames of movie */
#define NVID 250          /* number of iterations between images displayed on screen */
#define NSEG 150         /* number of segments of boundary */
#define INITIAL_TIME 200    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 4   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 6

/* Plot type, see list in global_ljones.c  */

#define PLOT 0
#define PLOT_B 6        /* plot type for second movie */

#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */

/* particle properties */

#define PARTICLE_HUE_MIN 330.0      /* color of original particle */
#define PARTICLE_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_EMAX 1.5e3           /* energy of particle with hottest color */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-7     /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 5.0e2   /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 0.1   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.25     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.001            /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e8    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e8    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 5.0       /* radius in which to count neighbours */
#define GRAVITY 0.0         /* gravity acting on all particles */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 100.0         /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 1  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 0.1   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define FINAL_CONSTANT_PHASE 0  /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 700            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 1     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 1   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 1         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 2.0  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 1   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define MAX_PRESSURE 1.0e8   /* pressure shown in "hottest" color */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 25        /* time at which to add first particle */
#define ADD_PERIOD 20       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 250  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */

#define HASHX 25   /* size of hashgrid in x direction */
#define HASHY 15   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 1 March 22 - Larger spirals in the Rock-Paper-Scissors reaction-diffusion equation ###

**Program:** `bz.c` 

**Initial condition in function `animation()`:** `init_random(0.5, 0.4, phi, psi, chi, xy_in);`

```
#define MOVIE 1        /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 1920          /* number of grid points on x axis */
#define NY 1000          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 99    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define MU 0.8	            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MRATIO 5            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15            /* number of grid point for grid of disks */
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
/* xy_in_billiard and draw_billiard in sub_wave.c */

/* Physical patameters of wave equation */

#define DT 0.0000001

#define VISCOSITY 10.0

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */

#define EPSILON 0.8     /* time scale separation */
#define DELTA 0.1       /* time scale separation */
#define BZQ 0.0008       /* parameter in BZ equation */
#define BZF 1.2          /* parameter in BZ equation */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.1     /* noise intensity */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 5000      /* number of frames of movie */
#define NVID 20          /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */

#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define INITIAL_TIME 5  /* initial still time */
#define END_FRAMES 50   /* number of still frames at end of movie */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Plot type */

#define PLOT 0

#define PRINT_TIME 0    /* set to 1 to print running time */

#define DRAW_FIELD_LINES 0  /* set to 1 to draw field lines */
#define FIELD_LINE_WIDTH 1  /* width of field lines */
#define N_FIELD_LINES 120   /* number of field lines */
#define FIELD_LINE_FACTOR 120 /* factor controlling precision when computing origin of field lines */

/* Color schemes, see list in global_pdes.c  */

#define COLOR_PALETTE 18     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 4   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 1.0   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 4.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0      /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.0   

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 0.55    /* scale of color scheme bar */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```


